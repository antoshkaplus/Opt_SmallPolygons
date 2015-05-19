
#ifndef SMALL_POLYGONS_SOLVER
#define SMALL_POLYGONS_SOLVER


#include "ant/core/core.hpp"
#include "simplex_insertion.hpp"


using namespace chrono;


class SmallPolygonsSolver {
private:
    
    SimplexInsertion* simplex_insertion;
    unsigned time_millis;
    
    grid::Grid<City> closest_cities;
    Count max_tour_count;
    Count city_count;
    vector<Point> points;
    vector<i::Point> points_int;
    // index of tour to which edge corresponds
    grid::Grid<double> edge_distance;
    grid::Grid<char> edge_valid;
    unordered_set<uint64_t> excluded_starting_tours;
    
    double best_area;
    bool should_improve = false;
    
    vector<Tour> GenerateStartingTours() {
        Seed seed = std::chrono::system_clock::now().time_since_epoch().count();
        auto starting_tours = ::GenerateStartingTours(closest_cities, excluded_starting_tours, max_tour_count, seed);
        auto adj_edges = ConstructAdjacentEdges(starting_tours, city_count);
        auto city_groups = ConstructCityGroups(starting_tours, city_count);
        auto inters = FindIntersectionsBetweenTours(points_int, closest_cities, 
                                                    city_groups, starting_tours, adj_edges);
        if (!inters.empty()) {
            throw logic_error("starting tours have intersections");
        }
        return starting_tours;
    }
    
    
    void InitializePoints(const vector<int>& point_coordinates) {
        points.resize(city_count);
        points_int.resize(city_count);
        for (int i = 0; i < city_count; ++i) {
            points[i].x = points_int[i].x = point_coordinates[2*i];
            points[i].y = points_int[i].y = point_coordinates[2*i+1];
        }
    }
    
public:
    // return zero index indices of vertices for each polygon
    // they space-separated
    vector<string> choosePolygons(vector<int> point_coordinates, int max_tour_count) {
        time_t t = ant::GetMillisCount();
        
        excluded_starting_tours.clear();
        this->max_tour_count = max_tour_count;
        city_count = point_coordinates.size()/2;
        InitializePoints(point_coordinates);
        edge_distance = ComputeEdgeDistance(points);
        closest_cities = ComputeClosestCities(edge_distance, 30);
        edge_valid = ConstructValidEdges(closest_cities);
        // should assign only those that are closest
        vector<Tour> best_sol_tours;
        best_area = numeric_limits<double>::max();
        
        Index iter = 0;
        while (true) {
            unsigned time_past = GetMillisCount() - t;
            if (time_past >= time_millis) {
                goto out;
            }
            if (time_past > time_millis/10 && best_area == numeric_limits<double>::max()) {
                simplex_insertion->set_param_area(simplex_insertion->param_area()/2);
                should_improve = false;
            } 
            if (time_past >= time_millis/2) {
                should_improve = false;
            }
            
            auto starting_tours = GenerateStartingTours();
            for (auto& s : starting_tours) {
                sort(s.begin(), s.end());
            }
            if (starting_tours.empty()) throw logic_error("starting tours are empty");
            // created initial solution probably infeasible
            vector<Tour> sol_tours = simplex_insertion->Solve(
                                                points_int, closest_cities, 
                                                edge_distance, starting_tours);
            ++iter;
            auto adj_edges = ConstructAdjacentEdges(sol_tours, city_count);
            
            if (should_improve) {
                vector<Index> city_groups = ConstructCityGroups(sol_tours, city_count);
                Improve(closest_cities, edge_distance, adj_edges, sol_tours, city_groups, 5);
            }
            
            bool success = MakeSolutionValid(points_int, closest_cities, edge_valid, adj_edges, excluded_starting_tours, sol_tours);
            
            if (!success) continue;
            
            double area = ::ComputeArea(points_int, sol_tours);
            if (area < best_area) {
                best_sol_tours = sol_tours;
                best_area = area;
            }
        }
    out:
        if (best_area == numeric_limits<double>::max()) {
            throw logic_error("solution not found");
        }
        cout << "itercount: " << iter << endl;
        cout << "area: " << best_area << endl << endl;
        vector<string> sol;
        for (auto& t : best_sol_tours) { 
            string result = "";
            for (int k : t) {
                result += std::to_string(k) + " ";
            }
            result.erase(result.end()-1);
            sol.push_back(result);
        }
        
        return sol;
    }
    
    double area() const {
        return best_area;
    }
    
    void set_simplex_insertion(SimplexInsertion& simplex_insertion) {
        this->simplex_insertion = &simplex_insertion;
    }
    
    void set_time_millis(unsigned time) {
        time_millis = time;
    }
    
    void set_should_improve(bool improve) {
        should_improve = improve;
    }
};



#endif

