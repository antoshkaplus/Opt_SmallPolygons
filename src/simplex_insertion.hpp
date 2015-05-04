
#ifndef SIMPLEX_INSERTION
#define SIMPLEX_INSERTION


#include "util.hpp"


struct SimplexInsertionError : logic_error {
    
    using logic_error::logic_error;
    
    vector<Tour> tours;
    vector<i::Point> points;
    
};


struct SimplexInsertion {
protected:
    struct Item {
        Edge edge;
        City city;
        double profit;
        
        Item(Edge e, City c, double p) 
        : edge(e), city(c), profit(p) {}
        
        bool operator<(const Item t) const {
            return profit < t.profit;
        }
    };
    
    const std::vector<i::Point> *points;
    const grid::Grid<City> *closest_cities; 
    const grid::Grid<double> *edge_distance;
    
    Count city_count;
    // get matrix for each point get 10-20 closest, sort them
    // row - each city, col - closest cities
    grid::Grid<char> edge_exists;
    std::vector<bool> visited;
    std::priority_queue<Item> insertion_queue;
    double param_exclude_edge = 0;
    double param_area_ = 0;
    
    double Distance(City c_0, City c_1) {
        return (*edge_distance)(c_0, c_1);
    }
    
    bool Exists(Edge e) {
        return edge_exists(e[0], e[1]);
    }
    
    double Profit(Edge e, City v) {
        double b = Distance(e[0], v);
        double c = Distance(e[1], v);
        double a = Distance(e[0], e[1]);
        return b + c + param_exclude_edge * a + param_area_ * sqrt(HeronFormula(a, b, c));
    }
    
    void Remove(Edge e) {
        edge_exists(e[0], e[1]) = edge_exists(e[1], e[0]) = false;
    }
    
    
    void Add(Edge e) {
        edge_exists(e[0], e[1]) = edge_exists(e[1], e[0]) = true; 
    }
    
    
    void UpdateQueue(City new_v, City e_0, City e_1) {
        auto& cc = *closest_cities; 
        for (int j = 0; j < cc.col_count(); ++j) {
            int k = cc(new_v, j);
            if (visited[k]) continue;
            // queue is max
            insertion_queue.emplace(Edge{e_0, new_v}, k,
                                    -Profit({e_0, new_v}, k));
            insertion_queue.emplace(Edge{e_1, new_v}, k, 
                                    -Profit({e_1, new_v}, k));
        }
    }
    
    vector<Tour> BuildTours() {
        for (int i = 0; i < city_count; ++i) {
            if (!visited[i]) {
                throw logic_error("city not visited");
            }
        }
        fill(visited.begin(), visited.end(), false);
        vector<Tour> tours;
        for (int k = 0; k < city_count; ++k) {
            if (visited[k]) continue;
            City first = k;
            City cur = k;
            City prev = -1;
            tours.resize(tours.size()+1);
            do {
                tours.back().push_back(cur);
                for (int j = 0; j < city_count; ++j) {
                    if (Exists(Edge{cur, j}) && j != prev) {
                        prev = cur;
                        cur = j;
                        break;
                    }
                }
            } while (cur != first);
            for (auto t : tours.back()) {
                visited[t] = true;
            }
        }
        if (find(visited.begin(), visited.end(), false) != visited.end()) {
            throw logic_error("tours are not correct");
        }
        return tours; 
    }
    
    void InitializeWithStartingTours(const vector<Tour>& starting_tours) {
        for (const Tour& t : starting_tours) {
            for (Index i = 0; i < t.size(); ++i) {
                Count sz = t.size();
                City c = t[i];
                City prev = t[(i-1+sz) % sz];
                City next = t[(i+1) % sz];
                Add({c, next});
                visited[c] = true;
                UpdateQueue(c, prev, next);
            }
        }
        
    }

public:
    
    void set_param_exclude_edge(double p) {
        param_exclude_edge = p;
    }
    
    void set_param_area(double p) {
        param_area_ = p;
    }
    
    double param_area() const {
        return param_area_;
    }
    
    virtual std::vector<Tour> Solve(const std::vector<i::Point>& points,
                            const grid::Grid<City>& closest_cities,
                            const grid::Grid<double>& edge_distance, 
                            const vector<Tour>& starting_tours) {
        this->edge_distance = &edge_distance;
        this->points = &points;
        this->closest_cities = &closest_cities;
        city_count = points.size();
        // can be reused
        visited.resize(city_count);
        fill(visited.begin(), visited.end(), false);
        edge_exists.resize(city_count, city_count);
        edge_exists.fill(false);
        InitializeWithStartingTours(starting_tours);
        while (!insertion_queue.empty()) {
            Item t = insertion_queue.top();
            insertion_queue.pop();
            if (visited[t.city] || !Exists(t.edge)) continue;
            Remove(t.edge);
            Add(Edge{t.edge[0], t.city});
            Add(Edge{t.edge[1], t.city});
            visited[t.city] = true;
            UpdateQueue(t.city, t.edge[0], t.edge[1]);
        }
        return BuildTours();
    }
    
};


#endif
