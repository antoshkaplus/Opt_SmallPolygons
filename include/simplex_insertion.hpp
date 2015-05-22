
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
    AdjacentEdges* tours;
    
    Count city_count;
    // get matrix for each point get 10-20 closest, sort them
    // row - each city, col - closest cities
    std::priority_queue<Item> insertion_queue;
    double param_exclude_edge = 0;
    double param_area_ = 0;
    
    double Distance(City c_0, City c_1) {
        return (*edge_distance)(c_0, c_1);
    }
    
    double Profit(const Edge& e, City v) {
        double b = Distance(e[0], v);
        double c = Distance(e[1], v);
        double a = Distance(e[0], e[1]);
        return b + c + param_exclude_edge * a + param_area_ * sqrt(HeronFormula(a, b, c));
    }
    
    void UpdateQueue(City new_v) {
        const auto& ts = *tours; 
        auto& cc = *closest_cities; 
        for (int j = 0; j < cc.col_count(); ++j) {
            int k = cc(new_v, j);
            if (ts.Visited(k)) continue;
            // queue is max
            for (const Edge& e : {Edge{{ts.Prev(new_v), new_v}}, Edge{{new_v, ts.Next(new_v)}}}) {
                insertion_queue.emplace(e, k, -Profit(e, k));
            }
        }
    }
        
    void InitializeWithStartingTours() {
        auto& ts = *tours; 
        for (City c = 0; c < city_count; ++c) {
            if (!ts.Visited(c)) continue;
            UpdateQueue(c);
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
    
    // first we are going to create IncompleteTour => Tour
    virtual void Solve(const std::vector<i::Point>& points,
                       const grid::Grid<City>& closest_cities,
                       const grid::Grid<double>& edge_distance, 
                       AdjacentEdges& tours
                       ) {
        this->edge_distance = &edge_distance;
        this->points = &points;
        this->closest_cities = &closest_cities;
        this->tours = &tours;
        auto& ts = tours; 
        city_count = points.size();
        // can be reused
        InitializeWithStartingTours();
        while (!insertion_queue.empty()) {
            Item t = insertion_queue.top();
            insertion_queue.pop();
            if (ts.Visited(t.city) || !ts.EdgeExists(t.edge)) continue;
            ts.InsertCity(t.city, t.edge);
            UpdateQueue(t.city);
        }
    }
    
};


#endif
