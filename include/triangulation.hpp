//
//  triangulation.hpp
//  SmallPolygon
//
//  Created by Anton Logunov on 7/13/15.
//
//

#ifndef SmallPolygon_triangulation_hpp
#define SmallPolygon_triangulation_hpp

#include "ant/geometry/triangle/delaunay_triangulation_i.hpp"
#include "ant/geometry/d2.hpp"

#include "util.hpp"

using namespace std;
using namespace ant;
using namespace ant::geometry::d2::i;
using namespace ant::geometry::d2;


class Triangulation {
    
    void Solve(std::vector<i::Point>& points, Count tour_count) {
        auto r = CircumRectangle(points);
        auto t = CircumTriangle(r);
        
        DelaunayTriangulation tr;
        auto trs = tr.Compute(points, t);
        
        vector<Count> adjacent_triangles(points.size(), 0);
        for (auto& t : trs) {
            for (Index i = 0; i < 3; ++i) {
                ++adjacent_triangles[t[i]]; 
            }
        }
        
        // this should be like a data structure... need to think about it.
        
        unordered_map<Edge, std::array<trianlge::Triangle>> adj_trs;
        vector<bool> visited(points.size(), false);
        
        auto area = [](const triangle::Triangle& t_0, const triangle::Triangle& t_1) {
            return 0.;
        } 
        
        // first choose starting triangles
        // just need to sort indices ????
        vector<Index> trs_sorted;
        iota(trs_sorted.begin(), trs_sorted.end(), 0);
        sort(trs_sorted.begin(), trs_sorted.end(), [](Index i_0, Index i_1) {
            return area(trs[i_0]) < area(trs[i_1]);
        });
        
        // first just take tour_count triangles... but play by the rules
        // all points should be unvisited
        auto can_use = [&](const triangle::Triangle& t) {
            for (auto i = 0; i < 3; ++i) {
                if (visited[t[i]]) {
                    return false;
                }
            }
            return true;
        };
 
        auto visit = [&](const triangle::Triangle& t) {
            for (auto i = 0; i < 3; ++i) {
                visited[t[i]] = true;
            }
        };
 
        // go inside with triangle
        std::priority_queue<Index> q;
        
        vector<Index> used_triangles;
                                                  
        Count use = 0;
        Index current = 0;
        while (use != tour_count) {
            if (can_use(trs[current])) {
                visit(trs[current]);
                used_triangles.push_back(current);
                ++use;
            }
            ++current;
        }
        
        while (true) {
            auto i = q.top();
            q.pop();
            if (!can_use(i)) {
                continue;
            }
            used_triangles.push_back(i);
            // put adjacent triangles into the queue
        }
        
        
        // after that I need somehow construct polygons out of this shit 
        // will  do it by edges
        
        // need to build tours first from edges
        
        AdjacentEdges adjacent_edges;
        
    }

};


#endif
