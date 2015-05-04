
#include "util.hpp"



bool ValidityCheck(const vector<int>& s) {
    // feasibility check
    vector<bool> visited(s.size(), false);
    for (City c : s) {
        if (c < 0 || c >= s.size()) {
            cout << "city is out of range";
            return false;
        }
        if (visited[c]) {
            cout << "should not visit one city twice";
            return false;
        }
        visited[c] = true;
    }
    return true;
}

void ValidityCheck(const AdjacentEdges& adj_edgs, 
                   const vector<Tour>& tours, 
                   const vector<City>& city_groups) {
    Count city_count = city_groups.size();
    vector<bool> visited(city_count, false);
    for (int i = 0; i < tours.size(); ++i) {
        auto& t = tours[i];
        assert(t.size() >= 3);
        for (auto k = 0; k < t.size(); ++k) {
            int c = t[k];
            assert(city_groups[c] == i);
            assert(!visited[c]);
            visited[c] = true;
            int c_prev = t[(k+t.size()-1)%t.size()];
            int c_next = t[(k+1)%t.size()];
            auto& a = adj_edgs[c];
            assert((a.first == c_prev && a.second == c_next) || 
                   (a.first == c_next && a.second == c_prev));
        }
    }
    for (int i = 0; i < city_count; ++i) {
        assert(visited[i]);
    }              
}


AdjacentEdges ConstructAdjacentEdges(const Polygons& ss, Count city_count) {
    AdjacentEdges res(city_count);
    for (auto& s : ss) {
        for (int i = 0; i < s.size(); ++i) {
            int prev = (i + s.size() -1) % s.size();
            int next = (i+1) % s.size();
            res[s[i]] = {s[prev], s[next]};
        }
    }
    return res;
}


vector<Index> ConstructCityGroups(const Polygons& ss, Count city_count) {
    vector<Index> r(city_count);
    for (int i = 0; i < ss.size(); ++i) {
        auto& s = ss[i];
        for (auto c : s) {
            r[c] = i;
        }
    }
    return r;
}


vector<Edge> ToursToEdges(const vector<Tour>& tours) {
    vector<Edge> edges;
    auto sz = 0;
    for (auto& t : tours) {
        sz += t.size();
    }
    edges.reserve(sz);
    for (auto& t : tours) {
        auto s = TSP::tourToEdges(t);
        edges.insert(edges.end(), s.begin(), s.end());
    }
    return edges;
}


Intersections FindIntersectionsForTour(const vector<i::Point>& ps,
                                       const grid::Grid<City>& near,
                                       const Tour& tour,
                                       const AdjacentEdges& adj,
                                       const vector<bool>& visited) {
    unordered_set<uint64_t> inters_hashes;
    Intersections inters;
    array<City, 4> ee;
    for (int i = 0; i < tour.size(); ++i) {
        Edge e{tour[i], tour[(i+1) % tour.size()]};
        if (e[0] > e[1]) swap(e[0], e[1]);
        for (auto e_c : e) {
            for (int j = 0; j < near.col_count(); ++j) {
                int f = near(e_c, j);
                if (!visited[f]) continue;
                for (auto r : {adj[f].first, adj[f].second}) {
                    Edge e_2{f, r};
                    if (e_2[0] > e_2[1]) swap(e_2[0], e_2[1]);
                    
                    if (e[0] == e_2[0] || e[0] == e_2[1] ||
                        e[1] == e_2[0] || e[1] == e_2[1] ||
                        !Segment{ps[e[0]], ps[e[1]]}.IntersectOrLie(
                            Segment{ps[e_2[0]], ps[e_2[1]]})) {
                            
                            continue;
                        }
                    merge(e.begin(), e.end(), e_2.begin(), e_2.end(), ee.begin());
                    auto h = Hash(ee[0], ee[1], ee[2], ee[3]);
                    if (inters_hashes.count(h) == 1) {
                        continue;
                    }
                    inters_hashes.insert(h);
                    inters.emplace_back(e, e_2); 
                }
            }
        }
    }
    return inters;
}


// believe that tours have all the vertices inside 
Intersections FindIntersectionsBetweenTours(
                    const vector<i::Point>& ps, 
                    const grid::Grid<int>& near,
                    const vector<Index>& city_corresp, 
                    const vector<Tour>& tours,
                    const AdjacentEdges& adj) {
    unordered_set<uint64_t> inters_hashes;
    Intersections inters;
    array<City, 4> ee;
    // should be construct
    auto visited = ComputeVisitedCities(tours, ps.size());
    for (int t_ind = 0; t_ind < tours.size(); ++t_ind) {
        auto& t = tours[t_ind];
        for (int c_ind = 0; c_ind < t.size(); ++c_ind) {
            Edge e{t[c_ind], t[(c_ind+1) % t.size()]};
            if (e[0] > e[1]) swap(e[0], e[1]);
            
            for (auto e_c : e) {
                for (int j = 0; j < near.col_count(); ++j) {
                    int f = near(e_c, j);
                    if (!visited[f] || city_corresp[e_c] == city_corresp[f]) {
                        continue;
                    }
                    for (auto r : {adj[f].first, adj[f].second}) {
                        
                        Edge e_2{f, r};
                        if (e_2[0] > e_2[1]) swap(e_2[0], e_2[1]);
                        
                        if (!Segment{ps[e[0]], ps[e[1]]}.IntersectOrLie(
                                    Segment{ps[e_2[0]], ps[e_2[1]]})) {
                            continue;
                        }
                        merge(e.begin(), e.end(), e_2.begin(), e_2.end(), ee.begin());
                        auto h = Hash(ee[0], ee[1], ee[2], ee[3]);
                        if (inters_hashes.count(h) == 1) {
                            continue;
                        }
                        inters_hashes.insert(h);
                        inters.emplace_back(e, e_2); 
                    }
                }
            }    
        }
    }
    return inters;
}


grid::Grid<City> ComputeClosestCities(const grid::Grid<double>& ds, 
                                      Count closest_count) {
    // one element is going to be excluded everytime
    Count k = min(closest_count, ds.row_count()-1);
    Count city_count = ds.row_count();
    grid::Grid<City> closest_cities(city_count, closest_count);
    std::vector<City> cities;
    for (City c = 0; c < city_count; ++c) {
        cities.resize(city_count);
        iota(cities.begin(), cities.end(), 0);
        std::swap(cities[c], cities.back());
        cities.pop_back();
        // could do it quicker too
        partial_sort(cities.begin(), cities.begin() + k, cities.end(), [=, &ds](City c_0, City c_1) {
            return ds(c, c_0) < ds(c, c_1);
        });
        cities.resize(k);
        for (int i = 0; i < k; ++i) {
            closest_cities(c, i) = cities[i];
        }
    }
    return closest_cities;
}


vector<Tour> GenerateStartingTours(const grid::Grid<City>& closest_cities, 
                                   const unordered_set<uint64_t>& excluded_starting_tours,
                                   Count max_tour_count, Seed seed) {
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
        if (visited[k_1] || visited[k_2]) {
            // here could probably try 2 and 3 ???
            continue;
        }
        // not a fact that k_i is inside cities
        // could be removed before
        buf[0] = k_0;
        buf[1] = k_1;
        buf[2] = k_2;
        sort(buf.begin(), buf.end());
        auto h = Hash(buf[0], buf[1], buf[2], 0);
        if (excluded_starting_tours.count(h) == 1) {
            continue;
        }
        visited[k_0] = visited[k_1] = visited[k_2] = true;
        result.push_back(buf);
    }
    return result;
}


vector<bool> ComputeVisitedCities(const vector<Tour>& tours, Count city_count) {
    vector<bool> visited(city_count, false);
    for (auto& t : tours) {
        for (auto c : t) {
            visited[c] = true;
        }
    }
    return visited;
}


// visited should be size of number of cities
void InitVisitedCities(const vector<Tour>& tours, vector<bool>& visited) {
    fill(visited.begin(), visited.end(), false);
    for (auto& t : tours) {
        for (auto c : t) {
            visited[c] = true;
        }
    }
}



// tour may be less than full amount
// O(n) anyway
// need to return new edges for some possible checks
void SwitchEdgesInTour(Tour& tour, Edge& e_0, Edge& e_1, AdjacentEdges& adj_edges) {
    Count n = tour.size();
    array<Index, 2> oe_0;
    array<Index, 2> oe_1;
    for (auto i : {0, 1}) {
        oe_0[i] = find(tour.begin(), tour.end(), e_0[i]) - tour.begin();
        oe_1[i] = find(tour.begin(), tour.end(), e_1[i]) - tour.begin(); 
    }
    if (oe_0[1] != (oe_0[0] + 1) % n) {
        swap(oe_0[0], oe_0[1]);
        swap(e_0[0], e_0[1]);
    }
    if (oe_1[1] != (oe_1[0] + 1) % n) {
        swap(oe_1[0], oe_1[1]);
        swap(e_1[0], e_1[1]);
    }
    
    // updating adj_edges
    auto& p_0_0 = adj_edges[e_0[0]];
    if (p_0_0.first == e_0[1]) {
        p_0_0.first = e_1[0];
    } else {
        p_0_0.second = e_1[0];
    }
    auto& p_0_1 = adj_edges[e_0[1]]; 
    if (p_0_1.first == e_0[0]) {
        p_0_1.first = e_1[1];
    } else {
        p_0_1.second = e_1[1];
    }
    
    auto& p_1_0 = adj_edges[e_1[0]];
    if (p_1_0.first == e_1[1]) {
        p_1_0.first = e_0[0];
    } else {
        p_1_0.second = e_0[0];
    }
    auto& p_1_1 = adj_edges[e_1[1]]; 
    if (p_1_1.first == e_1[0]) {
        p_1_1.first = e_0[1];
    } else {
        p_1_1.second = e_0[1];
    }
    // end updating adj_edges
    
    if (oe_0[0] > oe_1[1]) {
        swap(e_0[0], e_1[1]);
        reverse(tour.begin() + oe_1[1], tour.begin() + oe_0[0] + 1);    
    }
    else { // oe_1[0] > oe_0[1]
        swap(e_1[0], e_0[1]);
        reverse(tour.begin() + oe_0[1], tour.begin() + oe_1[0] + 1);
    }
} 


                      
bool MakeTourFeasible(const vector<i::Point>& points,
                      const grid::Grid<City>& closest_cities,
                      const grid::Grid<char>& valid_edge,
                      AdjacentEdges& adj_edges,
                      unordered_set<uint64_t>& excluded_tours,
                      Tour& tour,
                      Count max_iteration_count) {
    if (tour.size() == 3 && i::Collinear(points[tour[0]], points[tour[1]], points[tour[2]])) {
        excluded_tours.insert(Hash(tour[0], tour[1], tour[2], 0));
        return false;
    }
    
    default_random_engine rng(std::chrono::system_clock::now().time_since_epoch().count());
    vector<bool> vs(points.size());
    InitVisitedCities({tour}, vs);
    Index iteration = 0;
    while (iteration++ < max_iteration_count) {
        // need to do it everytime
        Intersections inters = FindIntersectionsForTour(points, 
            closest_cities, tour, adj_edges, vs);
        if (inters.empty()) return true;
        
        shuffle(inters.begin(), inters.end(), rng);
        Tour new_tour;
        for (auto i : inters) {
            auto& e_0 = i.first;
            if (adj_edges[e_0[0]].first != e_0[1] && adj_edges[e_0[0]].second != e_0[1]) {
                continue;
            }
            auto& e_1 = i.second;
            if (adj_edges[e_1[0]].first != e_1[1] && adj_edges[e_1[0]].second != e_1[1]) {
                continue;
            }
            
            new_tour = tour;
            // need to keep adj_edges valid to check if exists on top
            SwitchEdgesInTour(new_tour, e_0, e_1, adj_edges);
            auto c_0 = e_0[0];
            auto k_0 = e_0[1];
            auto c_1 = e_1[0];
            auto k_1 = e_1[1];
            if (!valid_edge(c_0, k_0) || !valid_edge(c_1, k_1)) {
                return false;
            }
            // hopefully reduced number of intersections
            tour = new_tour;
        }
    }
    return false;
}



grid::Grid<double> ComputeEdgeDistance(const vector<Point>& ps) {
    Count city_count = ps.size();
    grid::Grid<double> h(city_count, city_count);
    for (int r = 0; r < city_count; ++r) {
        for (int c = 0; c < city_count; ++c) {
            h(r, c) = h(c, r) = ps[r].distance(ps[c]);
        }
    }
    return h;
}


int ComputeArea(const vector<i::Point>& ps, const vector<Tour>& ts) {
    int res = 0;
    for (auto t : ts) {
        res += i::ShoelaceFormula(ps, t);
    }
    return res;
}


grid::Grid<char> ConstructValidEdges(const grid::Grid<City>& closest_cities) {
    auto c_count = closest_cities.row_count();
    grid::Grid<char> res(c_count, c_count);
    res.fill(false);
    for (auto i = 0; i < c_count; ++i) {
        for (auto j = 0; j < closest_cities.col_count(); ++j) {
            auto k = closest_cities(i, j);
            res(i, k) = res(k, i) = true;
        }
    }
    return res;
}


string TourToString(const Tour& t) {
    string result = "";
    for (int k : t) {
        result += std::to_string(k) + " ";
    }
    result.erase(result.end()-1);
    return result;
}

bool ResolveIntersectionsBetweenTwoTours(
                        const vector<i::Point>& ps, 
                        const grid::Grid<City>& near, 
                        vector<Tour>& tours, 
                        AdjacentEdges& adj_edges,
                        unordered_set<uint64_t>& excluded_tours,
                        Count max_interation_count) {
    
    Count city_count = ps.size();
    auto city_groups = ConstructCityGroups(tours, city_count);
    Intersections inters;
    vector<bool> visited(tours.size());
    
    using P = pair<Edge, Edge>;
    const auto& g = city_groups; 
    
    // city group equality
    auto eq = [&](const P& p_0, const P& p_1) {
        return g[p_0.first[0]] == g[p_1.first[0]] && g[p_0.second[0]] == g[p_1.second[0]];
    };
    auto edge_exists = [&](const Edge& e) {
        return adj_edges[e[0]].first == e[1] || adj_edges[e[0]].second == e[1];
    };
    
    Index iteration = 0;
    while (iteration++ < max_interation_count) {
        inters = FindIntersectionsBetweenTours(ps, near, city_groups, tours, adj_edges);
        
        if (inters.empty()) return true;
        
        fill(visited.begin(), visited.end(), false);
        
        for (auto i : inters) {
            if (city_groups[i.first[0]] > city_groups[i.second[0]]) {
                swap(i.first, i.second);
            }
        }
        sort(inters.begin(), inters.end(), [&](const P& i_0, const P& i_1) {
            return g[i_0.first[0]] < g[i_1.first[0]] ||
            (g[i_0.first[0]] == g[i_1.first[0]] && g[i_0.second[0]] < g[i_1.second[0]]);
        });
 
        for (int i = 0; i < inters.size(); i += 2) {
            if (i+2 < inters.size() && eq(inters[i], inters[i+2])) {
                // some times angle of one figure intersect angle of another figure
                // don't really want to mess with it
                return false;
            }
            auto& a_0 = inters[i].first;
            auto& a_1 = inters[i+1].first;
            Index a_i = g[a_0[0]];
            if (visited[a_i]) {
                if (!(edge_exists(a_0) && edge_exists(a_1) && city_groups[a_1[0]] == a_i)) {
                    continue;
                }
            }
            auto& b_0 = inters[i].second;
            auto& b_1 = inters[i+1].second;
            Index b_i = g[b_0[0]];
            if (visited[b_i]) {
                if (!(edge_exists(b_0) && edge_exists(b_1) && city_groups[b_1[0]] == b_i)) {
                    continue;
                }
            }
            if (EdgeEqual(a_0, a_1)) {
                int c = EdgeCommonCity(b_0, b_1);
                if (c != -1) {
                    if (tours[b_i].size() == 3) {
                        auto& s = tours[b_i]; 
                        excluded_tours.insert(Hash(s[0], s[1], s[2], 0));
                        return false;
                    }
                    city_groups[c] = a_i;
                    ReplaceCity(adj_edges, tours[b_i], tours[a_i], c, a_0);
                    visited[a_i] = visited[b_i] = true;
                    continue;
                }
            }
            if (EdgeEqual(b_0, b_1)) {
                int c = EdgeCommonCity(a_0, a_1);
                if (c != -1) {
                    if (tours[a_i].size() == 3) {
                        auto& s = tours[a_i]; 
                        excluded_tours.insert(Hash(s[0], s[1], s[2], 0));
                        return false;
                    }
                    city_groups[c] = b_i;
                    ReplaceCity(adj_edges, tours[a_i], tours[b_i], c, b_0);
                    visited[a_i] = visited[b_i] = true;
                    continue;
                }
            }
            City c_a = EdgeCommonCity(a_0, a_1);
            City c_b = EdgeCommonCity(b_0, b_1);
            if (c_a != -1 && c_b != -1) {
                swap(city_groups[c_a], city_groups[c_b]);
                visited[a_i] = visited[b_i] = true;
                ReplaceCity_2(adj_edges, tours[a_i], c_a, tours[b_i], c_b);
                continue;
            }
            
            // a lot of shit between // need to look deeper
            return false;
        }
    }
    return false;
}

bool EdgeEqual(const Edge& e_0, const Edge& e_1) {
    return (e_0[0] == e_1[0] && e_0[1] == e_1[1]) ||
            (e_0[0] == e_1[1] && e_0[1] == e_1[0]);
}


// returns -1 if not found
Index EdgeCommonCity(const Edge& e_0, const Edge& e_1) {
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            if (e_0[i] == e_1[j]) return e_0[i]; 
        }
    }
    return -1;
}

void ReplaceCity(AdjacentEdges& adj_edges, Tour& source, Tour& target, City what, Edge where) {
    // have some problems when tours in the same tour
    
    // start update adj edges
    auto old = adj_edges[what];
    adj_edges[what] = {where[0], where[1]};
    auto& od_0 = adj_edges[old.first];
    auto& od_1 = adj_edges[old.second]; 
    if (od_0.first == what) {
        od_0.first = old.second;
    } else {
        od_0.second = old.second;   
    }
    if (od_1.first == what) {
        od_1.first = old.first;
    } else {
        od_1.second = old.first;   
    }
    // end update adj edges
    
    auto& w_0 = adj_edges[where[0]];
    if (w_0.first == where[1]) w_0.first = what;
    else {
        w_0.second = what;
    }
    auto& w_1 = adj_edges[where[1]];
    if (w_1.first == where[0]) w_1.first = what;
    else {
        w_1.second = what;
    }

    remove(source.begin(), source.end(), what);
    source.pop_back();
    
    InsertCityBetween(target, what, where);
}


void ReplaceCity_2(AdjacentEdges& adj_edges, Tour& a, City c_a, Tour& b, City c_b) {
    // start update adj edges
    auto& old_a = adj_edges[c_a];
    auto& old_b = adj_edges[c_b];
    {
        auto& od_0 = adj_edges[old_a.first];
        auto& od_1 = adj_edges[old_a.second]; 
        if (od_0.first == c_a) {
            od_0.first = c_b;
        } else {
            od_0.second = c_b;   
        }
        if (od_1.first == c_a) {
            od_1.first = c_b;
        } else {
            od_1.second = c_b;   
        }
    }
    {
        auto& od_0 = adj_edges[old_b.first];
        auto& od_1 = adj_edges[old_b.second]; 
        if (od_0.first == c_b) {
            od_0.first = c_a;
        } else {
            od_0.second = c_a;   
        }
        if (od_1.first == c_b) {
            od_1.first = c_a;
        } else {
            od_1.second = c_a;   
        }
    }
    swap(old_a, old_b);
    
    Index i_a = find(a.begin(), a.end(), c_a) - a.begin();
    Index i_b = find(b.begin(), b.end(), c_b) - b.begin();
    
    swap(a[i_a], b[i_b]);
}


void InsertCityBetween(Tour& target, City what, Edge where) {
    Index i = find(target.begin(), target.end(), where[0]) - target.begin();
    assert(target[i] == where[0]);
    Index next = (i+1)%target.size();
    if (target[next] != where[1]) {
        next = i;
    }
    target.insert(target.begin() + next, what);
}

bool MakeSolutionValid(const vector<i::Point>& points_int, 
                       const grid::Grid<City>& closest_cities, 
                       const grid::Grid<char>& edge_valid, 
                       AdjacentEdges& adj_edges, 
                       unordered_set<uint64_t>& excluded_tours,
                       vector<Tour>& tours) {
                       vector<int> groups =ConstructCityGroups(tours, points_int.size());
    
    for (auto& t : tours) {
        // adj_edges aren't changing for the particular tour
        if (!MakeTourFeasible(points_int, closest_cities, edge_valid, adj_edges, excluded_tours, t, 10)) {
            return false;
        } 
    }
    // need to update those dudes
    return ResolveIntersectionsBetweenTwoTours(
                            points_int, closest_cities, 
                            tours, adj_edges, 
                            excluded_tours, 10);
                      
}




// pass three cities, return area of triangle
void Improve(const grid::Grid<City>& closest_cities, 
             const grid::Grid<double>& ds, 
             AdjacentEdges& adj_edgs, 
             vector<Tour>& tours, 
             vector<City>& city_groups, 
             Count closest_count) {
    
    auto area = [&] (City c_0, City c_1, City c_2) {
        return HeronFormula(ds(c_0, c_1), ds(c_1, c_2), ds(c_0, c_2));
    };
    
    Count city_count = city_groups.size();
    for (City c = 0; c < city_count; ++c) {
        // can't move vertex
        if (tours[city_groups[c]].size() == 3) continue;
        
        auto& c_p = adj_edgs[c];
        double c_area = area(c, c_p.first, c_p.second);
        if (isnan(c_area)) continue; 
        double min_area = c_area;
        Edge min_e;
        for (Index i = 0; i < closest_count; ++i) {
            City q = closest_cities(c, i);
            if (city_groups[q] == city_groups[c]) continue;
            auto& c_q = adj_edgs[q];
            for (auto qq : {c_q.first, c_q.second}) {
                //if (qq == c) continue;
                double q_area = area(c, q, qq);
                if (isnan(q_area)) continue;
                if (q_area < min_area) {
                    min_area = q_area;
                    min_e = {q, qq};
                }
            }
        }
        if (min_area == c_area) continue;
        // now need to exchange
        Index c_t = city_groups[c];
        Index min_t = city_groups[min_e[0]];
        city_groups[c] = min_t;
        ReplaceCity(adj_edgs, tours[c_t], tours[min_t], c, min_e);
    }
}
















