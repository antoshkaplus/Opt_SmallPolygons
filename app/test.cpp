

// try to create really big field
// take 30 points like it was before 

// and compare find intersections with our new algorithm.
// know the time

#include "ant/geometry/bentley_ottmann.hpp"
#include "ant/core/core.hpp"
#include "simplex_insertion.hpp"


using namespace std;
using namespace ant;
using namespace geometry;
using namespace d2;


int main(int argc, const char * argv[]) {
    int N = 1500;
    auto sample_i = GenerateSample(N);
    vector<f::Point> sample_f(N);
    transform(sample_i.begin(), sample_i.end(), sample_f.begin(), [](const i::Point& p) {
        return f::Point(p.x, p.y);
    });
    SimplexInsertion si;
    auto edge_distances = ComputeEdgeDistance(sample_f);
    auto closest_cities = ComputeClosestCities(edge_distances, 30);
    AdjacentEdges tours(N);
    auto ts = ::GenerateStartingTours(closest_cities, unordered_set<uint64_t>(), 20, 0);
    tours.AddTourRange(ts);
    si.Solve(sample_i, closest_cities, edge_distances, tours);
    
    auto& ps = sample_f;
    vector<f::Point> shortened_ps(2*N);
    vector<array<Index, 2>> segs(N);
    double a = 1.e-7;
    for (auto i = 0; i < N; ++i) {
        shortened_ps[2*i] = (1-a)*ps[i] + a*ps[tours.Next(i)];
        shortened_ps[2*i+1] = a*ps[i] + (1-a)*ps[tours.Next(i)];
        segs[i] = {{2*i, 2*i+1}};
    }
    auto func = [&](Index s_0, Index s_1) {
        return f::Intersection(
                   f::Segment{shortened_ps[segs[s_0][0]], shortened_ps[segs[s_0][1]]}, 
                   f::Segment{shortened_ps[segs[s_1][0]], shortened_ps[segs[s_1][1]]});
    };
    
    
    // finding intersecctions with bentley ottmann
    auto t = GetMillisCount();
    BentleyOttmann<f::Point, decltype(func)> finder_1;
    finder_1.FindIntersections(shortened_ps, segs, func);
    auto d = GetMillisCount() - t;
    cout << "BentleyOttmann: " << d << endl;
    
    // finding intersections near brute force
    t = GetMillisCount();
    FindIntersections(sample_i, closest_cities, tours);
    d = GetMillisCount() - t;
    cout << "NearBruteForce: " << d << endl;
    
    
    
    
    
}
