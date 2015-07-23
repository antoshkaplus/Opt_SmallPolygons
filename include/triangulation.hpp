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
#include "ant/geometry/triangle/adjacent_triangles.hpp"
#include "ant/geometry/d2.hpp"

#include "util.hpp"

using namespace std;
using namespace ant;
using namespace ant::geometry::d2::i;
using namespace ant::geometry::d2;
using namespace ant::geometry;


class Triangulation {
public:
    vector<Tour> Solve(std::vector<i::Point>& points, Count tour_count) {
        auto r = CircumRectangle(points);
        r.origin.x -= 1;
        r.origin.y -= 1;
        r.size.width += 2;
        r.size.height += 2;
        auto t = CircumTriangle(r);
        
        DelaunayTriangulation tr;
        auto trs = tr.Compute(points, t);
        
        vector<Count> adjacent_triangles(points.size(), 0);
        for (auto& t : trs) {
            for (Index i = 0; i < 3; ++i) {
                ++adjacent_triangles[t[i]]; 
            }
        }
        
        ant::geometry::triangle::AdjacentTrianglesIndex adj_trs;
        for (Index i = 0; i < trs.size(); ++i) {
            auto& t = trs[i];
            for (auto e : t.Edges()) {
                adj_trs.Insert(e, i);
            } 
        }
        
        
        // this should be like a data structure... need to think about it.
        
        // visit and visited related to points
        
        vector<bool> visited(points.size(), false);
        
        auto area = [&](const triangle::Triangle& t) {
            return std::abs(CrossProduct(points[t[0]], points[t[1]], points[t[2]]));
        };
        
        // first choose starting triangles
        // just need to sort indices ????
        vector<Index> trs_sorted(trs.size());
        iota(trs_sorted.begin(), trs_sorted.end(), 0);
        sort(trs_sorted.begin(), trs_sorted.end(), [&](Index i_0, Index i_1) {
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
        
        auto count_visited = [&](const triangle::Triangle& t) {
            Count n = 0;
            for (auto i = 0; i < 3; ++i) {
                if (visited[t[i]]) n++;
            }
            return n;
        };
        
        auto comp_by_area = [&](Index i_0, Index i_1) {
            return area(trs[i_0]) > area(trs[i_1]);
        };
 
        // go inside with triangle
        std::priority_queue<Index, std::vector<Index>, decltype(comp_by_area)> q(comp_by_area);
        
        vector<Index> used_triangles;
                                                  
        Count use = 0;
        Index current = 0;
        tour_count = 1;
        while (use != tour_count) {
            auto& t = trs[trs_sorted[current]];
            if (can_use(t)) {
                visit(t);
                used_triangles.push_back(trs_sorted[current]);
                // need to find out adjacent triangles and push them into queue
                for (auto e : t.Edges()) {
                    Index a = adj_trs.another(e, trs_sorted[current]);
                    if (a < 0) continue;
                    q.push(a);
                } 
                ++use;
            }
            ++current;
        }
        
        auto can_add = [&](Index i) {
            auto& t = trs[i];
            bool free = false;
            for (Index k = 0; k < 3; ++k) {
                if (!visited[t[k]]) {
                    free = true;
                    break;
                }
            }
            return free;
        };
        
        
        while (!q.empty()) {
            auto i = q.top();
            q.pop();
            if (!can_add(i)) {
                continue;
            }
            visit(trs[i]);
            used_triangles.push_back(i);
            // put adjacent triangles into the queue
            for (auto e : trs[i].Edges()) {
                // index of triangle 
                Index a = adj_trs.another(e, i);
                if (a < 0 || visited[trs[a].Third(e)]) continue;
                q.push(a);
            } 
        }
        
        for (Index i = 0; i < visited.size(); ++i) {
            if (!visited[i]) cout << i << endl;
            
        }
        
        for (auto& t : trs) {
            if (count_visited(t) == 2) {
                cout << t << endl;
            }
        }
        
        Count c = count(visited.begin(), visited.end(), false); 
        if (c > 0) {
            cout << "we failed" << endl;
        }
        
        
        
        for (auto& t : trs) {
            for (Index k = 0; k < 3; ++k) {
                if (!visited[t[k]]) {
                    cout << t << endl;
                }
            }
        }
        
        
        ant::geometry::triangle::AdjacentTrianglesIndex adj_trs_i;
        for (auto u : used_triangles) {
            for (auto e : trs[u].Edges()) {
                adj_trs_i.Insert(e, u);
            }
        }
        using Edge = std::array<Index, 2>;
        vector<Edge> edges;
        for (auto p : adj_trs_i) {
            if (p.second[0] == -1) {
                edges.push_back({{p.first[0], p.first[1]}});
            } else if (p.second[1] == -1) {
                edges.push_back({{p.first[0], p.first[1]}});
            } 
        }
        
        
        AdjacentItems<Index, Index, -1> adj_ns;
        for (auto e : edges) {
            adj_ns.Insert(e[0], e[1]);
            adj_ns.Insert(e[1], e[0]);
        }
        
        
        vector<Tour> tours;
        while (!adj_ns.Empty()) {
            auto b = adj_ns.begin();
            Tour t;
            Index B = b->first;
            t.push_back(b->second[0]);
            t.push_back(b->first);
            Index stop = b->second[0];
            Index prev = b->first;
            Index next = b->second[1];
            while (next != stop) {
                t.push_back(next);
                Index cur = adj_ns.another(next, prev);
                adj_ns.Remove(next);
                Index prev = next;
                Index next = cur;
            }
            adj_ns.Remove(B);
            tours.push_back(t);
        }
        
        return tours;
        // now need to construct tours from edges
    }

};


#endif
