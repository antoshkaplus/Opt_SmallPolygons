

#ifndef UTIL_SMALL_POLY
#define UTIL_SMALL_POLY

#include <fstream>
#include <iostream>
#include <random>
#include <map>
#include <unordered_set>
#include <queue>

#include "ant/core/core.hpp"
#include "ant/geometry/d2.h"
#include "ant/optimization/tsp/tsp.hpp"
#include "ant/optimization/tsp/tsp_ant_colony.hpp"
#include "ant/optimization/tsp/tsp_random_insertion.hpp"
#include "ant/optimization/tsp/tsp_farthest_insertion.hpp"
#include "ant/optimization/tsp/tsp_simplex_insertion.hpp"


using namespace std;
using namespace ant;
template<typename T> using Grid = grid::Grid<T>;
using namespace ant::geometry::d2::f;
using namespace ant::opt::tsp;
using Segment = ant::geometry::d2::i::Segment;
using namespace ant::geometry::d2;
using namespace ant::linalg;


using Polygons = vector<vector<City>>;
using Intersections = vector<pair<TSP::Edge, TSP::Edge>>;
using Seed = std::default_random_engine::result_type;



void MakeEachTourFeasible(AdjacentEdges& edges, vector<City> tour_begin);






struct MakeTourFeasibleError : runtime_error {
    
    Intersections inters;
    Tour tour;
    vector<i::Point> points;
    
    using runtime_error::runtime_error;
    
};

struct IntersectionsBetweenToursError : runtime_error {
    
    Intersections inters;
    vector<Tour> tours;
    vector<i::Point> points;
    
    using runtime_error::runtime_error;
    
};


void InsertCityBetween(Tour& target, City what, Edge where);

bool ValidityCheck(const vector<int>& s);
void ValidityCheck(const AdjacentEdges& adj_edgs, 
                   const vector<Tour>& tours, 
                   const vector<City>& city_groups);


AdjacentEdges ConstructAdjacentEdges(const Polygons& ss, Count city_count);
vector<Index> ConstructCityGroups(const Polygons& ss, Count city_count);

vector<Edge> ToursToEdges(const vector<Tour>& tours);

Intersections FindIntersectionsForTour(const vector<i::Point>& ps,
                                       const grid::Grid<City>& near,
                                       const Tour& tour,
                                       const AdjacentEdges& adj,
                                       const vector<bool>& visited);

Intersections FindIntersectionsBetweenTours(
                                            const vector<i::Point>& ps, 
                                            const grid::Grid<int>& near,
                                            const vector<Index>& city_corresp, 
                                            const vector<Tour>& tours,
                                            const AdjacentEdges& adj);


grid::Grid<City> ComputeClosestCities(const grid::Grid<double>& ds, 
                                      Count closest_count);
                                      
vector<Tour> GenerateStartingTours(const grid::Grid<City>& closest_cities, 
                                   const unordered_set<uint64_t>& excluded_starting_tours,
                                   Count max_tour_count, 
                                   Seed seed);
                                   

// this one should be usuful when not that many points (not that many possible polygon)
template<class Distance>
vector<Tour> GenerateStartingTours(const vector<i::Point>& ps,
                                   const Distance& ds,
                                   const grid::Grid<City>& closest_cities, 
                                   const unordered_set<uint64_t>& excluded_starting_tours,
                                   Count max_tour_count, 
                                   Seed seed) {
    std::default_random_engine rng(seed);
    Count city_count = closest_cities.row_count();
    vector<City> cities(city_count);
    iota(cities.begin(), cities.end(), 0);
    vector<Tour> result;
    vector<bool> visited(city_count, false);
    vector<City> buf(3);
    while (!cities.empty() && result.size() != max_tour_count) {
        // size of cities is changing
        std::uniform_int_distribution<City> distr(0, cities.size()-1);
        Index i = distr(rng);
        City k_0 = cities[i];
        swap(cities.back(), cities[i]);
        cities.pop_back();
        if (visited[k_0]) continue;
        City k_1 = closest_cities(k_0, 0);
        City k_2 = closest_cities(k_0, 1);
        bool b_1 = visited[k_1];
        bool b_2 = visited[k_2];
        
        if (b_1 && b_2) continue;
        
        if (b_1 || b_2) {
            //double area = HeronFormula(ds(k_0, k_1), ds(k_1, k_2), ds(k_2, k_0));
            City k_3 = closest_cities(k_0, 2);
            bool b = visited[k_3];
            if (b) continue;
            if (b_1) {
                PointInsideTriangle<i::Point> tr(ps[k_0], ps[k_3], ps[k_2]);
                if (tr.IsInside(ps[k_1])) { //|| area < HeronFormula(ds(k_0, k_3), ds(k_3, k_2), ds(k_2, k_0))) 
                    continue;
                }
                k_1 = k_3;
            } else { // b_2 visited
                PointInsideTriangle<i::Point> tr(ps[k_0], ps[k_3], ps[k_1]);
                if (tr.IsInside(ps[k_2])) { //|| area < HeronFormula(ds(k_0, k_1), ds(k_1, k_3), ds(k_3, k_0))) 
                    continue;
                }
                k_2 = k_3;
            }
        }
        // not a fact that k_i is inside cities
        // could be removed before
        buf[0] = k_0;
        buf[1] = k_1;
        buf[2] = k_2;
        sort(buf.begin(), buf.end());
        if (excluded_starting_tours.count(Hash(buf[0], buf[1], buf[2], 0)) == 1) {
            continue;
        }
        visited[k_0] = visited[k_1] = visited[k_2] = true;
        result.push_back(buf);
    }
    return result;
}



                                   
vector<bool> ComputeVisitedCities(const vector<Tour>& tours, Count city_count);

void SwitchEdgesInTour(Tour& tour, Edge& e_0, Edge& e_1, AdjacentEdges& adj_edges);

bool MakeTourFeasible(const vector<i::Point>& points,
                      const grid::Grid<City>& closest_cities,
                      const grid::Grid<char>& valid_edge,
                      AdjacentEdges& adj_edges,
                      unordered_set<uint64_t>& excluded_tours,
                      Tour& tour,
                      Count max_iteration_count);
                      

grid::Grid<double> ComputeEdgeDistance(const vector<Point>& ps);                      
int ComputeArea(const vector<i::Point>& ps, const vector<Tour>& ts);

grid::Grid<char> ConstructValidEdges(const grid::Grid<City>& closest_cities);

string TourToString(const Tour& t);

bool ResolveIntersectionsBetweenTwoTours(const vector<i::Point>& ps, 
                                         const grid::Grid<City>& near, 
                                         vector<Tour>& tours, 
                                         AdjacentEdges& adj_edges,
                                         unordered_set<uint64_t>& excluded_tours,
                                         Count max_iteration_count);
                                        
bool EdgeEqual(const Edge& e_0, const Edge& e_1);
Index EdgeCommonCity(const Edge& e_0, const Edge& e_1);

void ReplaceCity_2(AdjacentEdges& adj_edges, Tour& a, City c_a, Tour& b, City c_b);
// need to update those
void ReplaceCity(AdjacentEdges& adj_edges, Tour& source, Tour& target, City what, Edge where);

bool MakeSolutionValid(const vector<i::Point>& points_int, 
                       const grid::Grid<City>& closest_cities, 
                       const grid::Grid<char>& edge_valid, 
                       AdjacentEdges& adj_edges, 
                       unordered_set<uint64_t>& excluded_tours,
                       vector<Tour>& tours);
                    
void Improve(const grid::Grid<City>& closest_cities, 
             const grid::Grid<double>& ds, 
             AdjacentEdges& adj_edgs, 
             vector<Tour>& tours, 
             vector<City>& city_groups, 
             Count closest_count);








#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Sweep_line_2_algorithms.h>





#endif



