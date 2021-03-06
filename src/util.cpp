
#include "util.hpp"





uniform_int_distribution<> POINTS_DISTRIBUTION(0, 699);
default_random_engine RNG;



// get factor of 2 for going through everything multiple times sometimes
void MakeTourFeasible(const vector<i::Point>& points,
                      const Grid<City>& closest_cities, 
                      const vector<City>& city_group,
                      AdjacentEdges& tours, 
                      City tour_start
                      ) {
    auto& ts = tours;
    auto& ps = points;
    City cur = tour_start;
    while (cur != tour_start) { 
        start:
        City next = ts.Next(cur);
        Edge e = {{cur, next}};
        Index cur_group = city_group[cur];
        for (auto e_c : e) {
            for (auto j = 0; j < closest_cities.col_count(); ++j) {
                auto k = closest_cities(e_c, j);
                if (city_group[k] != cur_group) {
                    continue;
                }
                for (const Edge& e_2 : {Edge{{ts.Prev(k), k}}, Edge{{k, ts.Next(k)}}}) {
                    if (!i::Segment{ps[e[0]], ps[e[1]]}.IntersectOrLie(
                            i::Segment{ps[e_2[0]], ps[e_2[1]]})
                        || e[0] == e_2[0] || e[0] == e_2[1] ||
                           e[1] == e_2[0] || e[1] == e_2[1]) {
                    
                        continue;
                    }
                    // got intersection
                    ts.FlipInTour(e, e_2);
                    // don't know will it help 
                    tour_start = cur;
                    goto start;
                }
            }
        }
        cur = next;
    }
}



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
            assert((a[0] == c_prev && a[1] == c_next) || 
                   (a[0] == c_next && a[1] == c_prev));
        }
    }
    for (int i = 0; i < city_count; ++i) {
        assert(visited[i]);
    }              
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
                for (auto r : {adj[f][0], adj[f][1]}) {
                    Edge e_2{f, r};
                    if (e_2[0] > e_2[1]) swap(e_2[0], e_2[1]);
                    
                    if (e[0] == e_2[0] || e[0] == e_2[1] ||
                        e[1] == e_2[0] || e[1] == e_2[1] ||
                        !i::Segment{ps[e[0]], ps[e[1]]}.IntersectOrLie(
                            i::Segment{ps[e_2[0]], ps[e_2[1]]})) {
                            
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
                    for (auto r : {adj[f][0], adj[f][1]}) {
                        
                        Edge e_2{f, r};
                        if (e_2[0] > e_2[1]) swap(e_2[0], e_2[1]);
                        
                        if (!i::Segment{ps[e[0]], ps[e[1]]}.IntersectOrLie(
                                    i::Segment{ps[e_2[0]], ps[e_2[1]]})) {
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


Intersections FindIntersections(const vector<i::Point>& ps, 
                       const grid::Grid<int>& near,
                       const AdjacentEdges& adj) {
    Intersections inters;
    Count city_count = adj.CityCount();
    for (int i = 0; i < city_count; ++i) {
        Edge e{{i, adj.Next(i)}};
        if (e[0] > e[1]) swap(e[0], e[1]);
        for (auto e_c : e) {
            for (int j = 0; j < near.col_count(); ++j) {
                int f = near(e_c, j);
                for (auto r : {adj[f][0], adj[f][1]}) {
                    Edge e_2{{f, r}};
                    if (e_2[0] > e_2[1]) swap(e_2[0], e_2[1]);
                    if (e[0] == e_2[0] || e[0] == e_2[1] ||
                            e[1] == e_2[0] || e[1] == e_2[1] ||
                            !i::Segment{ps[e[0]], ps[e[1]]}.IntersectOrLie(
                                i::Segment{ps[e_2[0]], ps[e_2[1]]})) {
                        continue;
                    }
                    inters.emplace_back(e, e_2); 
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


grid::Grid<double> ComputeEdgeDistance(const vector<Point>& ps) {
    Count city_count = ps.size();
    grid::Grid<double> h(city_count, city_count);
    for (int r = 0; r < city_count; ++r) {
        for (int c = 0; c < city_count; ++c) {
            h(r, c) = h(c, r) = ps[r].Distance(ps[c]);
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

//bool ResolveIntersectionsBetweenTwoTours(
//                        const vector<i::Point>& ps, 
//                        const grid::Grid<City>& near, 
//                        vector<Tour>& tours, 
//                        AdjacentEdges& adj_edges,
//                        unordered_set<uint64_t>& excluded_tours,
//                        Count max_interation_count) {
//    
//    Count city_count = ps.size();
//    auto city_groups = ConstructCityGroups(tours, city_count);
//    Intersections inters;
//    vector<bool> visited(tours.size());
//    
//    using P = pair<Edge, Edge>;
//    const auto& g = city_groups; 
//    
//    // city group equality
//    auto eq = [&](const P& p_0, const P& p_1) {
//        return g[p_0.first[0]] == g[p_1.first[0]] && g[p_0.second[0]] == g[p_1.second[0]];
//    };
//    auto edge_exists = [&](const Edge& e) {
//        return adj_edges[e[0]][0] == e[1] || adj_edges[e[0]][1] == e[1];
//    };
//    
//    Index iteration = 0;
//    while (iteration++ < max_interation_count) {
//        inters = FindIntersectionsBetweenTours(ps, near, city_groups, tours, adj_edges);
//        
//        if (inters.empty()) return true;
//        
//        fill(visited.begin(), visited.end(), false);
//        
//        for (auto i : inters) {
//            if (city_groups[i.first[0]] > city_groups[i.second[0]]) {
//                swap(i.first, i.second);
//            }
//        }
//        sort(inters.begin(), inters.end(), [&](const P& i_0, const P& i_1) {
//            return g[i_0.first[0]] < g[i_1.first[0]] ||
//            (g[i_0.first[0]] == g[i_1.first[0]] && g[i_0.second[0]] < g[i_1.second[0]]);
//        });
// 
//        for (int i = 0; i < inters.size(); i += 2) {
//            if (i+2 < inters.size() && eq(inters[i], inters[i+2])) {
//                // some times angle of one figure intersect angle of another figure
//                // don't really want to mess with it
//                return false;
//            }
//            auto& a_0 = inters[i].first;
//            auto& a_1 = inters[i+1].first;
//            Index a_i = g[a_0[0]];
//            if (visited[a_i]) {
//                if (!(edge_exists(a_0) && edge_exists(a_1) && city_groups[a_1[0]] == a_i)) {
//                    continue;
//                }
//            }
//            auto& b_0 = inters[i].second;
//            auto& b_1 = inters[i+1].second;
//            Index b_i = g[b_0[0]];
//            if (visited[b_i]) {
//                if (!(edge_exists(b_0) && edge_exists(b_1) && city_groups[b_1[0]] == b_i)) {
//                    continue;
//                }
//            }
//            if (EdgeEqual(a_0, a_1)) {
//                int c = EdgeCommonCity(b_0, b_1);
//                if (c != -1) {
//                    if (tours[b_i].size() == 3) {
//                        auto& s = tours[b_i]; 
//                        excluded_tours.insert(Hash(s[0], s[1], s[2], 0));
//                        return false;
//                    }
//                    city_groups[c] = a_i;
//                    ReplaceCity(adj_edges, tours[b_i], tours[a_i], c, a_0);
//                    visited[a_i] = visited[b_i] = true;
//                    continue;
//                }
//            }
//            if (EdgeEqual(b_0, b_1)) {
//                int c = EdgeCommonCity(a_0, a_1);
//                if (c != -1) {
//                    if (tours[a_i].size() == 3) {
//                        auto& s = tours[a_i]; 
//                        excluded_tours.insert(Hash(s[0], s[1], s[2], 0));
//                        return false;
//                    }
//                    city_groups[c] = b_i;
//                    ReplaceCity(adj_edges, tours[a_i], tours[b_i], c, b_0);
//                    visited[a_i] = visited[b_i] = true;
//                    continue;
//                }
//            }
//            City c_a = EdgeCommonCity(a_0, a_1);
//            City c_b = EdgeCommonCity(b_0, b_1);
//            if (c_a != -1 && c_b != -1) {
//                swap(city_groups[c_a], city_groups[c_b]);
//                visited[a_i] = visited[b_i] = true;
//                ReplaceCity_2(adj_edges, tours[a_i], c_a, tours[b_i], c_b);
//                continue;
//            }
//            
//            // a lot of shit between // need to look deeper
//            return false;
//        }
//    }
//    return false;
//}

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


void InsertCityBetween(Tour& target, City what, Edge where) {
    Index i = find(target.begin(), target.end(), where[0]) - target.begin();
    assert(target[i] == where[0]);
    Index next = (i+1)%target.size();
    if (target[next] != where[1]) {
        next = i;
    }
    target.insert(target.begin() + next, what);
}

//bool MakeSolutionValid(const vector<i::Point>& points_int, 
//                       const grid::Grid<City>& closest_cities, 
//                       const grid::Grid<char>& edge_valid, 
//                       AdjacentEdges& adj_edges, 
//                       unordered_set<uint64_t>& excluded_tours,
//                       vector<Tour>& tours) {
//                       vector<int> groups =ConstructCityGroups(tours, points_int.size());
//    
//    for (auto& t : tours) {
//        // adj_edges aren't changing for the particular tour
//        if (!MakeTourFeasible(points_int, closest_cities, edge_valid, adj_edges, excluded_tours, t, 10)) {
//            return false;
//        } 
//    }
//    // need to update those dudes
//    return ResolveIntersectionsBetweenTwoTours(
//                            points_int, closest_cities, 
//                            tours, adj_edges, 
//                            excluded_tours, 10);
//                      
//}




// pass three cities, return area of triangle
//void Improve(const grid::Grid<City>& closest_cities, 
//             const grid::Grid<double>& ds, 
//             AdjacentEdges& adj_edgs, 
//             vector<Tour>& tours, 
//             vector<City>& city_groups, 
//             Count closest_count) {
//    
//    auto area = [&] (City c_0, City c_1, City c_2) {
//        return HeronFormula(ds(c_0, c_1), ds(c_1, c_2), ds(c_0, c_2));
//    };
//    
//    Count city_count = city_groups.size();
//    for (City c = 0; c < city_count; ++c) {
//        // can't move vertex
//        if (tours[city_groups[c]].size() == 3) continue;
//        
//        auto& c_p = adj_edgs[c];
//        double c_area = area(c, c_p[0], c_p[1]);
//        if (isnan(c_area)) continue; 
//        double min_area = c_area;
//        Edge min_e;
//        for (Index i = 0; i < closest_count; ++i) {
//            City q = closest_cities(c, i);
//            if (city_groups[q] == city_groups[c]) continue;
//            auto& c_q = adj_edgs[q];
//            for (auto qq : {c_q[0], c_q[1]}) {
//                //if (qq == c) continue;
//                double q_area = area(c, q, qq);
//                if (isnan(q_area)) continue;
//                if (q_area < min_area) {
//                    min_area = q_area;
//                    min_e = {q, qq};
//                }
//            }
//        }
//        if (min_area == c_area) continue;
//        // now need to exchange
//        Index c_t = city_groups[c];
//        Index min_t = city_groups[min_e[0]];
//        city_groups[c] = min_t;
//        ReplaceCity(adj_edgs, tours[c_t], tours[min_t], c, min_e);
//    }
//}

vector<i::Point> GenerateSample(int point_count) {
    vector<i::Point> res(point_count);
    unordered_set<int> ps;
    for (auto& p : res) {
        while (true) {
            p.set(POINTS_DISTRIBUTION(RNG), POINTS_DISTRIBUTION(RNG));
            int key = p.x * POINTS_DISTRIBUTION.max() + p.y;
            if (ps.find(key) == ps.end()) {
                ps.insert(key);
                break;
            }
        }
    }
    return res;
}















