//
//  main.cpp
//  FragileMirrors
//
//  Created by Anton Logunov on 5/5/13.
//  Copyright (c) 2013 Anton Logunov. All rights reserved.
//


#include "small_polygons_solver.hpp"
#include "small_polygons.hpp"
#include "simplex_insertion.hpp"
#include "triangulation.hpp"
#include "ant/optimization/optimization.h"


//
//struct Test {
//    vector<int> ps;
//    int n;
//    
//    Test() {}
//    Test(vector<int> ps, int n) : ps(ps), n(n) {}
//};
//
//uniform_int_distribution<> POINTS_DISTRIBUTION(0, 699);
//default_random_engine RNG;
//
//vector<i::Point> GenerateTest(int point_count) {
//    vector<i::Point> res(point_count);
//    unordered_set<int> ps;
//    for (auto& p : res) {
//        while (true) {
//            p.set(POINTS_DISTRIBUTION(RNG), POINTS_DISTRIBUTION(RNG));
//            int key = p.x * POINTS_DISTRIBUTION.max() + p.y;
//            if (ps.find(key) == ps.end()) {
//                ps.insert(key);
//                break;
//            }
//        }
//    }
//    return res;
//}
//
Sample GenerateSample() {
    uniform_int_distribution<> distr(20, 1500);
    default_random_engine rng;
    int Np = 3000;// 2*distr(rng);
    vector<int> ps(Np);
    for (int& p : ps) {
        p = POINTS_DISTRIBUTION(rng);
    }
    return {ps, 10}; 
}

string to_string(vector<int> inds) {
    string res;
    for (int i : inds) {
        res += to_string(i) + " ";
    }
    return res;
}

vector<int> to_vector_int(const vector<i::Point>& ps) {
    vector<int> res;
    res.reserve(2*ps.size());
    for (auto p : ps) {
        res.push_back(p.x);
        res.push_back(p.y);
    }
    return res;
}


struct Objective {
    vector<double> areas;
    
    Objective() {}
    
    Objective(const vector<double>& areas) : areas(areas) {}  
    
    bool operator<(Objective& obj) {
        Count c = 0;
        for (int i = 0; i < areas.size(); ++i) {
            if (areas[i] < obj.areas[i]) {
                ++c;
            }
        }
        return c > areas.size()/2;
    }
};

void TestRange(const pair<Count, Count>& range, Count count, double param_profit, unsigned test_time) {
    default_random_engine rng;
    uniform_int_distribution<> distr(range.first, range.second);
    uniform_int_distribution<> distr_poly(2, 20);
    SimplexInsertion si;
    si.set_param_exclude_edge(param_profit);
    SmallPolygonsSolver sol;
    sol.set_simplex_insertion(si);
    sol.set_time_millis(test_time);
    for (int i = 0; i < count; ++i) {
        Sample t{
            to_vector_int(GenerateSample(distr(rng))), 
            distr_poly(rng)
        };
        sol.choosePolygons(t.ps, t.n);
    }
}


void FindParamArea(unsigned test_time) {
    ofstream outprof("param_area.txt");
    Count test_count = 30;
    vector<pair<Count, Count>> groups = {/*{20, 99}, {100, 499},*/ {500, 1500}};
    uniform_int_distribution<> distr_poly(2, 20);
    
    SmallPolygons sm;
    sm.set_time_millis(test_time);
    for (auto g : groups) {
        uniform_int_distribution<> distr(g.first, g.second);
        default_random_engine rng;
        
        vector<Sample> ts;
        for (int i = 0; i < test_count; ++i) {
            ts.emplace_back(
                            to_vector_int(GenerateSample(distr(rng))),
                            distr_poly(rng)); 
        }
        
        auto f = [&, test_count](double param) {
            Objective obj;
            // try to create multiple threads here
            sm.insertion()->set_param_area(param);
            for (int i = 0; i < test_count; ++i) {
                sm.choosePolygons(ts[i].ps, ts[i].n);
                obj.areas.push_back(sm.area());
            }
            outprof << param << endl;
            return obj;
        };
        outprof << "started group " << g.first << " " << g.second << endl;
        double param  = opt::GoldenSectionSearch(0, 1.5, f, 0.1);
        outprof << "param: " << param;
    }
}


void FindParamExcludeEdge(unsigned test_time) {
    ofstream outprof("param_exclude_edge.txt");
    Count test_count = 30;
    vector<pair<Count, Count>> groups = {{20, 99}, {100, 499}, {500, 1500}};
    uniform_int_distribution<> distr_poly(2, 20);
    
    SimplexInsertion ins;
    SmallPolygonsSolver sol;
    sol.set_simplex_insertion(ins);
    sol.set_time_millis(test_time);
    for (auto g : groups) {
        uniform_int_distribution<> distr(g.first, g.second);
        default_random_engine rng;
        
        vector<Sample> ts;
        for (int i = 0; i < test_count; ++i) {
            ts.emplace_back(
                to_vector_int(GenerateSample(distr(rng))),
                distr_poly(rng)); 
        }
        
        auto f = [&, test_count](double param) {
            Objective obj;
            // try to create multiple threads here
            ins.set_param_exclude_edge(param);
            for (int i = 0; i < test_count; ++i) {
                sol.choosePolygons(ts[i].ps, ts[i].n);
                obj.areas.push_back(sol.area());
            }
            outprof << param << endl;
            return obj;
        };
        outprof << "started group " << g.first << " " << g.second << endl;
        double param  = opt::GoldenSectionSearch(0, 1, f, 0.1);
        outprof << "param: " << param;
    }

}

void TestShouldImprove() {
    default_random_engine rng;
    uniform_int_distribution<> distr(20, 99);
    uniform_int_distribution<> distr_poly(2, 20);
    SmallPolygons sm;
    sm.set_time_millis(2000);
    Count improve_better = 0;
    Count equal = 0;
    for (int i = 0; i < 100; ++i) {
        auto v = to_vector_int(GenerateSample(distr(rng)));
        auto N = distr_poly(rng);
        sm.set_should_improve(true);
        sm.choosePolygons(v, N);
        auto a_t = sm.area();
        sm.set_should_improve(false);
        sm.choosePolygons(v, N);
        auto a_f = sm.area();
        if (a_t < a_f) {
            ++improve_better;
            cout << "win" << endl;
        } else if (a_t == a_f) {
            ++equal;
            cout << "equal" << endl; 
        }
    }
    cout << "improve: " << improve_better << "; equal: " << equal << "; from: " << 100 << endl;
}




void ReadTest(ifstream& in, vector<int>& ps, int& N) {
    int Np;
    in >> Np;
    ps.resize(Np);
    for (int i = 0; i < Np; ++i) {
        in >> ps[i];
    }
    in >> N;
}

void WriteTest(ofstream& out, const vector<string>& ts) {
    out << ts.size() << endl;
    for (auto s : ts) {
        out << s << endl;
    }
}



int main(int argc, const char * argv[])
{
    try {
        if (argc == 1) {
            TestShouldImprove();
            return 0;
        }
        if (argc == 2) {
            // first arg where test cases stored
            string root(argv[1]);
            for (int t = 1; t <= 10; ++t) {
                string filename = root + "/output_" + to_string(t) + ".txt";
                ifstream in(filename);
                vector<int> points;
                int N;
                ReadTest(in, points, N);
                in.close();
                SmallPolygons poly;
                vector<string> ss = poly.choosePolygons(points, N);
                string out_name = root + "/out_" + to_string(t) + ".txt";
                ofstream out(out_name);
                WriteTest(out, ss);
                out.close();
            }  
            return 0;  
        } 
        if (argc == 3) {
            // first arg - input.txt
            // second arg - output.txt
            vector<int> ps;
            int N;
            ifstream in(argv[1]);
            ReadTest(in, ps, N);
            in.close();
            SmallPolygons poly;
            poly.set_time_millis(2000);
            vector<string> ss = poly.choosePolygons(ps, N);
            ofstream out(argv[2]);
            WriteTest(out, ss);
            out.close();
            return 0;
        }
        
    } catch (MakeTourFeasibleError& error) {
        ofstream out_points("error_points.txt");
        auto v = to_vector_int(error.points);
        out_points << v.size() << endl;
        for (auto i : v) {
            out_points << i << endl;
        }
        ofstream out_tour("error_tour.txt");
        out_tour << 1 << endl; // one tour
        out_tour << TourToString(error.tour) << endl;
    } catch (IntersectionsBetweenToursError& error) {
        ofstream out_points("error_points.txt");
        auto v = to_vector_int(error.points);
        out_points << v.size() << endl;
        for (auto i : v) {
            out_points << i << endl;
        }
        ofstream out_tour("error_tour.txt");
        out_tour << error.tours.size() << endl; // one tour
        for (auto& t : error.tours) {
            out_tour << TourToString(t) << endl;
        }
    }
    
//    int Np;
//    cin >> Np;
//    vector<int> ps(Np);
//    for (int i = 0; i < Np; ++i) {
//        cin >> ps[i];
//    } 
//    int N;
//    cin >> N;
//    t.n = N;
//    t.ps = ps;
//    t = GenerateTest();
//    
//    SmallPolygons polys;
//    auto ss = polys.choosePolygons(t.ps, t.n);
//    ofstream out("log.txt");
//    cout << ss.size() << endl;
//    out << ss.size() << endl;
//    for (auto& s : ss) {
//        cout << s << endl;
//        out << s << endl;
//    } 
//    cout.flush();
//    out.flush();
    return 0;
}

