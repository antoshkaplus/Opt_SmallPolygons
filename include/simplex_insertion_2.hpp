
#ifndef SIMPLEX_INSERTION_2
#define SIMPLEX_INSERTION_2


#include "simplex_insertion.hpp"


struct SimplexInsertion_2 : SimplexInsertion {
    
    double Area(City c_0, City c_1, City c_2) {
        return HeronFormula(Distance(c_0, c_1), Distance(c_1, c_2), Distance(c_0, c_2));
    }
    
public:
    
    void set_param_profit(double p) {
        param_profit = p;
    }
    
    
    std::vector<Tour> Solve(const std::vector<i::Point>& points,
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
        
        vector<double> insertion_loss(city_count, numeric_limits<double>::max());
        vector<bool> improved(city_count, false); 
        
        InitializeWithStartingTours(starting_tours);
        while (!insertion_queue.empty()) {
            Item t = insertion_queue.top();
            insertion_queue.pop();
            if (!Exists(t.edge)) {
                continue;
            }
            auto area = Area(t.city, t.edge[0], t.edge[1]);
            if (visited[t.city]) {
                if (area < insertion_loss[t.city]) {
                    //improved[t.city]
                    insertion_loss[t.city] = area;
                }
                continue;
            }
            insertion_loss[t.city] = area;
            Remove(t.edge);
            Add(Edge{t.edge[0], t.city});
            Add(Edge{t.edge[1], t.city});
            visited[t.city] = true;
            UpdateQueue(t.city, t.edge[0], t.edge[1]);
        }
        
        while (true) {
            visited.resize(city_count);
            fill(visited.begin(), visited.end(), false);
            edge_exists.resize(city_count, city_count);
            edge_exists.fill(false);
            InitializeWithStartingTours(starting_tours);
            
            // now we start again
            while (!insertion_queue.empty()) {
                Item t = insertion_queue.top();
                insertion_queue.pop();
                if (!Exists(t.edge)) {
                    continue;
                }
                auto area = Area(t.city, t.edge[0], t.edge[1]);
                if (visited[t.city] || area > insertion_loss[t.city]) {
                    continue;
                }
                Remove(t.edge);
                Add(Edge{t.edge[0], t.city});
                Add(Edge{t.edge[1], t.city});
                visited[t.city] = true;
                UpdateQueue(t.city, t.edge[0], t.edge[1]);
            }
            
            try {
                return BuildTours();
            } catch (logic_error& e) {
            
            }
            
            for (int i = 0; i < city_count; ++i) {
                if (!visited[i]) {
                    insertion_loss[i] = numeric_limits<double>::max();
                }
            }
        }
        
        return BuildTours();
    }
    
};

#endif