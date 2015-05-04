
#include <fstream>
#include <iostream>
#include <random>
#include <map>
#include <unordered_set>
#include <array>
#include <stack>
#include <assert.h>
#include <algorithm>
#include <chrono>


namespace ant {
    
    using Int = int;
    // index and count should be of type int
    // because 
    // 1) unsigned types increase probability of making a bug
    // 2) lesser type will create problem of casting or being too little
    // 3) bigger type impossible to iterate through
    // the only thing is unsigned integers is good for bitwise operations
    using Count = int; 
    using Index = int;
    
    using Long = int64_t;
    using Float = float;
    using Double = double;
    
}



namespace ant {
    
    namespace geometry {
        
        namespace d2 {
            
            
            
            template<class T>
            struct Size {
                Size() : Size(0, 0) {}
                Size(T width, T height)
                : width(width), height(height) {}
                
                void set(T width, T height) {
                    this->width = width;
                    this->height = height;
                }
                T area() const {
                    return width*height;
                }
                T perimeter() const {
                    return 2*(height + width); 
                }
                bool isCovering(const Size<T>& s) const {
                    return width >= s.width && height >= s.height;
                }
                void swap() {
                    std::swap(height, width);
                }
                
                Size swapped() const {
                    return Size(height, width);
                }
                
                T width, height;
            };
            
            template<class T>
            bool operator==(const Size<T>& s_0, const Size<T>& s_1) {
                return s_0.width == s_1.width && s_0.height == s_1.height;
            }
            template<class T>
            bool operator!=(const Size<T>& s_0, const Size<T>& s_1) {
                return s_0.width != s_1.width || s_0.height != s_1.height;
            }
            
            
            
            namespace i {
                
                using Size = d2::Size<size_t>;
                
                struct Point {
                    Point() {}
                    Point(Int x, Int y) : x(x), y(y) {}
                    void set(Int x, Int y) {
                        this->x = x;
                        this->y = y;
                    }
                    void swap() {
                        std::swap(x, y);
                    }
                    Int x, y;
                };
                
                struct Segment {
                    Segment() {}
                    Segment(const Point& fst, const Point& snd)
                    : fst(fst), snd(snd) {}
                    Point fst, snd;
                    
                    void Swap() {
                        std::swap(fst, snd);
                    }
                    
                    Segment Swapped() const {
                        return Segment(snd, fst);
                    }
                    
                    bool Lie(Point q) const
                    {
                        return (q.x <= std::max(fst.x, snd.x) && q.x >= std::min(fst.x, snd.x) &&
                                q.y <= std::max(fst.y, snd.y) && q.y >= std::min(fst.y, snd.y));
                    }
                    
                    bool Intersect(const Segment s) const {
                        int o1 = Orientation(fst, snd, s.fst);
                        int o2 = Orientation(fst, snd, s.snd);
                        int o3 = Orientation(s.fst, s.snd, fst);
                        int o4 = Orientation(s.fst, s.snd, snd);
                        
                        // General case
                        if (o1 != o2 && o3 != o4)
                            return true;
                        
                        // Special Cases
                        // p1, q1 and p2 are colinear and p2 lies on segment p1q1
                        if (o1 == 0 && Lie(s.fst)) return true;
                        
                        // p1, q1 and p2 are colinear and q2 lies on segment p1q1
                        if (o2 == 0 && Lie(s.snd)) return true;
                        
                        // p2, q2 and p1 are colinear and p1 lies on segment p2q2
                        if (o3 == 0 && s.Lie(fst)) return true;
                        
                        // p2, q2 and q1 are colinear and q1 lies on segment p2q2
                        if (o4 == 0 && s.Lie(snd)) return true;
                        
                        return false; // Doesn't fall in any of the above cases
                    }
                    
                private:
                    // To find orientation of ordered triplet (p, q, r).
                    // The function returns following values
                    // 0 --> p, q and r are colinear
                    // 1 --> Clockwise
                    // 2 --> Counterclockwise
                    int Orientation(Point p, Point q, Point r) const
                    {
                        // See 10th slides from following link for derivation of the formula
                        // http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
                        int val = (q.y - p.y) * (r.x - q.x) -
                        (q.x - p.x) * (r.y - q.y);
                        if (val == 0) return 0;  // colinear
                        return (val > 0) ? 1 : 2; // clock or counterclock wise
                    }
                    
                };
                
                bool operator==(const Point& p_0, const Point& p_1);
                bool operator!=(const Point& p_0, const Point& p_1);
                
                struct Rectangle {
                    Rectangle() {}
                    Rectangle(const Point& origin, const Size& size) 
                    : origin(origin), size(size) {}
                    Rectangle(Int x, Int y, Int width, Int height) 
                    : origin(x, y), size(width, height) {}
                    
                    void set(Int x, Int y, size_t width, size_t height) {
                        origin.set(x, y);
                        size.set(width, height);
                    }
                    
                    void set(const Point& origin, const Point& diag) {
                        this->origin = origin;
                        size.set(diag.x - origin.x, diag.y - origin.y);
                    }
                    
                    void swap() {
                        origin.swap();
                        size.swap();
                    }
                    
                    size_t area() const {
                        return size.area();
                    }
                    
                    size_t perimeter() const {
                        return size.perimeter(); 
                    }
                    
                    bool isIntersect(const Rectangle& r) const {
                        return origin.x < r.origin.x + r.size.width  && origin.x + size.width  > r.origin.x &&
                        origin.y < r.origin.y + r.size.height && origin.y + size.height > r.origin.y;
                    }
                    
                    Point origin;
                    Size size;
                };
                
                bool operator==(const Rectangle& r_0, const Rectangle& r_1);
                bool operator!=(const Rectangle& r_0, const Rectangle& r_1);    
                
                // works only for simple polygons    
                int ShoelaceFormula(const std::vector<Point>& ps);
                int ShoelaceFormula(const std::vector<Point>& ps, const std::vector<Index>& order);
                
                
                struct Polygon {
                    std::vector<Point> points;
                    
                    double ComputeArea() const {
                        
                        
                        
                        return 0;
                    }
                };
                
            } // namespace i
            
            
            namespace f {
                
                
                struct Point {
                    Point() : Point(0, 0) {}
                    Point(Float x, Float y) : x(x), y(y) {}
                    
                    Float distance(const Point& p) const {
                        Float 
                        dx = p.x - x,
                        dy = p.y - y;
                        return sqrt(dx*dx + dy*dy);
                    } 
                    
                    static Float distance(const Point& p_0, const Point& p_1) {
                        return p_0.distance(p_1);
                    }
                    
                    void set(Float x, Float y) {
                        this->x = x;
                        this->y = y;
                    }
                    
                    Float x, y;
                };
                
                Point& operator+=(Point& p_0, const Point& p_1);    
                Point& operator/=(Point& p_0, Float f);
                Point operator+(Point p_0, const Point& p_1);
                Point operator/(Point p_0, Float f);
                
                
                struct Indent {
                    Indent() : Indent(0, 0) {}
                    Indent(Float dx, Float dy) : dx(dx), dy(dy) {}
                    
                    Indent& operator+=(const Indent& d) {
                        dx += d.dx;
                        dy += d.dy;
                        return *this;
                    }
                    
                    Float distance() const {
                        return sqrt(dx*dx + dy*dy);
                    }
                    
                    Indent normed() const {
                        auto d = distance();
                        return {dx/d, dy/d};
                    }
                    
                    Float dx, dy;
                };
                
                Indent& operator/=(Indent& i, Float f);
                Indent& operator*=(Indent& i, Float f); 
                Indent operator/(Indent i, Float f);
                Indent operator*(Indent i, Float f);
                Indent operator+(Indent i_0, Indent i_1);
                
                Point operator+(Indent i, Point p);
                Point operator-(Indent i, Point p); 
                Point operator+(Point p, Indent i);
                Point operator-(Point p, Indent i); 
                
                Point& operator+=(Point& p, Indent i);
                Indent operator-(const Point& p_0, const Point& p_1); 
                Indent operator*(Float d, Indent i);
                
                
                using Size = d2::Size<double>;
                
                struct Line {
                    Line() : Line(0, 0, 0) {}
                    Line(double a, double b, double c) : a(a), b(b), c(c) {}
                    Line(const Point& p_0, const Point& p_1) {
                        a = p_1.y - p_0.y;
                        b = p_0.x - p_1.x;
                        c = p_0.x*(p_0.y - p_1.y) + p_0.y*(p_1.x - p_0.x);
                    }
                    
                    double a, b, c;
                };
                
                struct Circle {
                    Circle() : radius(0) {}
                    Circle(Point center, double radius) : center(center), radius(radius) {}
                    Point center;
                    double radius;
                };
                
                struct Rectangle {
                    Rectangle() : origin(0, 0), size(0, 0) {}
                    Rectangle(Float x, Float y, Float width, Float height) 
                    : origin(x, y), size(width, height) {}
                    Rectangle(Point origin, Size size) 
                    : origin(origin), size(size) {}
                    
                    bool isInside(const Point& p) const {
                        return p.x >= origin.x && p.y >= origin.y && 
                        p.x <= origin.x+size.width && p.y <= origin.y+size.height;
                    }
                    
                    Point origin;
                    Size size; 
                };
                
                bool operator==(const Point& p_0, const Point& p_1);
                std::ostream& operator<<(std::ostream& output, const Point& p);
                std::pair<Point, Point> circleLineIntersection(const Circle& circle, const Line& line);
                
            } // namespace f
            
            
            
            // oab
            // Returns a positive value, if p_0, p_1, p_2 makes a counter-clockwise turn,
            // negative for clockwise turn, and zero if the points are collinear.
            template<class P>
            double CrossProduct(const P& p_0, const P& p_1, const P& p_2) {
                return (p_1.x - p_0.x)*(p_2.y - p_0.y) - (p_1.y - p_0.y)*(p_2.x - p_0.x);
            }
            
            
            // always counter clockwise
            // to make clockwise just reverse
            template<class P> // should have x and y inside
            std::vector<Index> ConvexHull(const std::vector<P>& ps) {
                std::vector<Index> inds;
                int n = (int)ps.size();
                std::vector<Index> order(n);
                for (int i = 0; i < n; i++) order[i] = i;
                // Sort points lexicographically
                sort(order.begin(), order.end(), 
                     [&ps](int i1, int i2){ return ps[i1].x < ps[i2].x || 
                         (ps[i1].x == ps[i2].x && ps[i1].y < ps[i2].y); });
                inds.resize(2*n);
                int k = 0;
                // Build lower hull
                for (int i = 0; i < n; i++) {
                    while (k >= 2 && CrossProduct(ps[inds[k-1]], ps[inds[k-2]], ps[order[i]]) <= 0) k--;
                    inds[k++] = order[i];
                }
                // Build upper hull
                for (int i = n-2, t = k+1; i >= 0; i--) {
                    while (k >= t && CrossProduct(ps[inds[k-1]], ps[inds[k-2]], ps[order[i]]) <= 0) k--;
                    inds[k++] = order[i];
                }
                inds.resize(k-1);
                return inds;
            }
            
            template<class P>
            double Perimeter(const std::vector<P>& ps, const std::vector<Index>& order, bool isClosed) {
                double s = 0.;
                for (Index i = 0; i < order.size()-1; i++) {
                    s += ps[order[i]].distance(ps[order[i+1]]);
                }
                if (isClosed) s += ps[order[0]].distance(ps[order[order.size()-1]]);
                return s;
            }
            
            // result is stored incide inds. inds - from which indices we choose
            template<class P>
            std::vector<Index> K_NearestPoints(const std::vector<P>& ps, const P& p, 
                                               std::vector<Index>& inds, int k) {
                k = std::min(k, (int)ps.size()-1);
                int n = (int)ps.size();
                inds.resize(n);
                for (int i = 0; i < n; i++) {
                    inds[i] = i;
                }
                partial_sort(inds.begin(), inds.begin()+k, inds.end(), 
                             [&ps, p](int i1, int i2){ return ps[i1].distance(p) < ps[i2].distance(p); });
                inds.resize(k);
            }
            
            template<class P>
            Index NearestPoint(const std::vector<P>& ps, const std::vector<Index>& indices, const P& p) {
                Index i_min = indices[0];
                for (Index i : indices) {
                    if (p.distance(ps[i]) < p.distance(ps[i_min])) {
                        i_min = i;
                    }
                }
                return i_min;
            }
            
            
            template <class ForwardIterator, class P>
            ForwardIterator NearestPoint(ForwardIterator first, ForwardIterator last, const P& p) {
                return std::min_element(first, last, [&p](const P& p_0, const P& p_1) {
                    return p_0.distance(p) < p_1.distance(p);
                });
            }
            
            template <class ForwardIterator, class P>
            ForwardIterator NearestPoint(const std::vector<P>& points, 
                                         ForwardIterator firstIndex, 
                                         ForwardIterator lastIndex, 
                                         const P& p) {
                return std::min_element(firstIndex, lastIndex, [&](size_t i_0, size_t i_1) {
                    return points[i_0].distance(p) < points[i_1].distance(p);
                });
            }
            
            template <class ForwardIterator, class P>
            ForwardIterator FarthestPoint(const std::vector<P>& points,
                                          ForwardIterator firstIndex,
                                          ForwardIterator lastIndex,
                                          const P& p) {
                return std::max_element(firstIndex, lastIndex, [&](Index i_0, Index i_1) {
                    return points[i_0].distance(p) < points[i_1].distance(p);
                });
            }
            
            
            
            
            
            
        } // namespace d2
        
    } // namespace geometry
    
} // namespace ant


namespace ant {
    
    namespace geometry {
        
        namespace d2 {
            
            namespace i {
                
                bool operator==(const Point& p_0, const Point& p_1) {
                    return p_0.x == p_1.x && p_0.y == p_1.y;
                }
                
                
                bool operator!=(const Point& p_0, const Point& p_1) {
                    return p_0.x != p_1.x || p_0.y != p_1.y;
                }
                
                
                
                int ShoelaceFormula(const std::vector<Point>& ps) {
                    double s = 0;
                    s += ps[ps.size()-1].x*ps[0].y - ps[0].x*ps[ps.size()-1].y;
                    for (auto i = 0; i < ps.size()-1; ++i) {
                        s += ps[i].x*ps[i+1].y - ps[i+1].x*ps[i].y;
                    }
                    return std::abs(s);
                }
                
                int ShoelaceFormula(const std::vector<Point>& ps, const std::vector<Index>& order) {
                    double s = 0;
                    s += ps[order[order.size()-1]].x*ps[order[0]].y - ps[order[0]].x*ps[order[order.size()-1]].y;
                    for (auto i = 0; i < order.size()-1; ++i) {
                        s += ps[order[i]].x*ps[order[i+1]].y - ps[order[i+1]].x*ps[order[i]].y;
                    }
                    return std::abs(s);
                }
                
                
                
            }
            
            
            
            namespace f {
                
                
                Point& operator+=(Point& p_0, const Point& p_1) {
                    p_0.x += p_1.x;
                    p_0.y += p_1.y;
                    return p_0;
                }
                
                Point& operator/=(Point& p_0, Float f) {
                    p_0.x /= f;
                    p_0.y /= f;
                    return p_0;
                } 
                
                Point operator+(Point p_0, const Point& p_1) {
                    return p_0 += p_1;
                }
                
                Point operator/(Point p_0, Float f) {
                    return p_0 /= f;
                }
                
                Indent& operator/=(Indent& i, Float f) {
                    i.dx /= f;
                    i.dy /= f;
                    return i;
                }
                
                Indent& operator*=(Indent& i, Float f) {
                    i.dx *= f;
                    i.dy *= f;
                    return i;
                }
                
                Indent operator/(Indent i, Float f) {
                    i /= f;
                    return i;
                }
                
                Indent operator*(Indent i, Float f) {
                    i *= f;
                    return i;
                }
                
                Indent operator+(Indent i_0, Indent i_1) {
                    i_0.dx += i_1.dx;
                    i_0.dy += i_1.dy;
                    return i_0;
                }
                
                Point& operator+=(Point& p, Indent i) {
                    p.x += i.dx;
                    p.y += i.dy;
                    return p;
                }
                
                Indent operator-(const Point& p_0, const Point& p_1) {
                    Indent i;
                    i.dx = p_0.x - p_1.x;
                    i.dy = p_0.y - p_1.y;
                    return i;
                }
                
                Indent operator*(Float d, Indent i) {
                    i *= d;
                    return i;
                }
                
                Point operator+(Indent i, Point p) {
                    p += i;
                    return p;
                }
                
                Point operator-(Indent i, Point p) {
                    return Point(i.dx - p.x, i.dy - p.y);
                }
                
                Point operator+(Point p, Indent i) {
                    return i + p;
                }
                
                Point operator-(Point p, Indent i) {
                    return Point(p.x - i.dx, p.y - i.dy);
                } 
                
                
                
            }
            
        }
        
    }
    
}




namespace ant {
    
    namespace opt {
        
        namespace tsp {
            
            using namespace ant::geometry::d2::f;
            using namespace ant::geometry::d2;
            
            
            typedef int City;
            
            
            struct TSP {
                typedef double Distance;
                typedef size_t Count;
                
                struct Edge : std::array<City, 2> {
                    Edge() {}
                    Edge(const std::array<City, 2>& arr) : array<City, 2>(arr) {}
                    Edge(City c_0, City c_1) : array<City, 2>({c_0, c_1}) {}
                    Edge(std::initializer_list<City> list) : array<City, 2>({*list.begin(), *(list.begin()+1)}) {}
                    
                    bool hasCity(City c) const {
                        return at(0) == c || at(1) == c; 
                    }
                    
                    City otherCity(City c) const {
                        return at(0) == c ? at(1) : at(0);
                    }
                };
                
                struct DisjointSet {
                    struct Record {
                        Record(ant::Index indexBeingChanged, bool didSizeIncrese) 
                        : indexBeingChanged(indexBeingChanged), didSizeIncrese(didSizeIncrese) {}
                        Index indexBeingChanged;
                        bool didSizeIncrese;
                    };
                    
                    std::vector<size_t> sz;
                    std::vector<Index> pt; 
                    std::stack<Record> records;
                    
                    DisjointSet(size_t count) : sz(count, 1), pt(count) {
                        iota(pt.begin(), pt.end(), 0);
                    }
                    
                    void unite(Index s_0, Index s_1) {
                        Index p_0 = find(s_0);
                        Index p_1 = find(s_1);
                        
                        if(sz[p_0] < sz[p_1]){
                            std::swap(p_0, p_1);
                        }
                        // if sz[0] > sz[1]  
                        pt[p_1] = p_0;
                        Record r(p_1, false);
                        if (sz[p_0] == sz[p_1]) {
                            sz[p_0]++;
                            r.didSizeIncrese = true;
                        } 
                        records.push(r);
                    }
                    
                    bool isUnited(Index s_0, Index s_1) const {
                        return find(s_0) == find(s_1);
                    }
                    
                    Index find(Index p) const {
                        while(p != pt[p]){
                            p = pt[p];
                        }
                        return p;
                    }
                    
                    Index operator[](Index i) const {
                        return find(i);
                    } 
                    
                    void undoUnite() {
                        assert(!records.empty());
                        auto &r = records.top();
                        Index index = pt[r.indexBeingChanged];
                        pt[r.indexBeingChanged] = r.indexBeingChanged;
                        if (r.didSizeIncrese) sz[index]--;
                        records.pop();
                    }
                };
                
                virtual std::vector<City> solve(const std::vector<Point>& points) = 0;
                
                static std::vector<Edge> tourToEdges(const std::vector<City>& cities) {
                    std::vector<Edge> edges;
                    for (Index i = 0; i < cities.size(); i++) {
                        edges.push_back(Edge(cities[i], cities[(i+1)%cities.size()]));   
                    }
                    return edges; 
                }
                
                static std::vector<City> edgesToTour(const std::vector<Edge>& edges) {
                    std::vector<City> cities(edges.size());
                    std::vector<Edge> es(edges);
                    cities[0] = es[0][0];
                    cities[1] = es[0][1];
                    swap(es[0], es.back());
                    es.pop_back();
                    for (size_t i = 2; i < cities.size(); i++) {
                        for (size_t e_i = 0; e_i < es.size(); e_i++) {
                            bool found = false;
                            const Edge &e = es[e_i];
                            if (e[0] == cities[i-1]) {
                                cities[i] = e[1]; 
                                found = true;
                            }
                            if (e[1] == cities[i-1]) {
                                cities[i] = e[0]; 
                                found = true;
                            }
                            if (found) {
                                swap(es[e_i], es.back());
                                es.pop_back();
                                break;
                            }
                        }
                    }
                    return cities;
                }
                
                // returns random city in [0,bound)
                static City random(City bound) {
                    std::default_random_engine rng(std::chrono::system_clock::now().time_since_epoch().count());
                    std::uniform_int_distribution<City> d(0, bound-1);
                    return d(rng);
                }
                
                static double random() {
                    std::default_random_engine rng(std::chrono::system_clock::now().time_since_epoch().count());
                    std::uniform_real_distribution<double> d(0, 1);
                    return d(rng);
                }
            };
            
            typedef TSP TSP_Solver;
            
            struct TSP_InsertionSolver : TSP_Solver {
                std::vector<City> startingTour(const std::vector<Point>& points) {
                    return ConvexHull(points);
                }
            };
            
            struct TSP_Improver {
                virtual std::vector<City> improve(const std::vector<Point>& points, 
                                                  const std::vector<City>& tour) = 0;
            };
            
            struct TSP_TwoOpt : TSP_Improver {
                typedef size_t Index;
                
                double swapEpsilon;
                
                TSP_TwoOpt() : swapEpsilon(1e-7) {}
                
                std::vector<City> improve(const std::vector<Point>& ps,
                                          const std::vector<City>& in_tour) override {
                    std::vector<City> tour = in_tour;
                    
#define next(i, t) (((i)+1)%(t).size())
#define distance(i_0, i_1) (ps[tour[i_0]].distance(ps[tour[i_1]])) 
                    //        function<Index(Index)> next = [&](Index i) { 
                    //            return (i+1)%tour.size(); 
                    //        };
                    // you need to store those distances
                    //        function<double(Index, Index)> distance = [&](Index i0, Index i1) {
                    //            return ps[tour[i0]].distance(ps[tour[i1]]);
                    //        };
                    //        
                    Index reverseCount = 0, printedReverseCount = 0;
                    bool start_again;
                    do {
                        start_again = false;
                        Index a1 = 0, a2, b1, b2;
                        while ((a2 = next(a1, tour)) != tour.size()-1) { 
                            b1 = next(a2, tour);
                            b2 = next(b1, tour);
                            while (b1 != 0) {
                                if (distance(a1, a2)+distance(b1, b2) > distance(a1, b1) + distance(a2, b2) + swapEpsilon) {
                                    reverse(tour.begin()+a2, tour.begin()+b1+1);
                                    reverseCount++;
                                    start_again = true;
                                }
                                b1 = b2;
                                b2 = next(b1, tour);
                            }
                            a1 = a2;
                        }
                        if (reverseCount > printedReverseCount+1000) {
                            std::cout << "reverce count: " << reverseCount << " obj: " << Perimeter(ps, tour, true)  << std::endl;
                            printedReverseCount = reverseCount;
                        }
                    } while (start_again);
#undef next
#undef distance
                    return tour;
                }
            }; 
            
            struct TSP_TwoOptAndHalf : TSP_Improver {
                typedef size_t Index;
                
                double swapEpsilon;
                
                TSP_TwoOptAndHalf() : swapEpsilon(1e-7) {}
                
                std::vector<City> improve(const std::vector<Point>& ps,
                                          const std::vector<City>& in_tour) override {
                    std::vector<City> tour = in_tour;
                    std::function<Index(Index)> next = [&](Index i) { 
                        return (i+1)%tour.size(); 
                    };
                    std::function<double(Index, Index)> distance = [&](Index i0, Index i1) {
                        return ps[tour[i0]].distance(ps[tour[i1]]);
                    };
                    
                    Index reverseCount = 0, printedReverseCount = 0;
                    bool start_again;
                    do {
                        start_again = false;
                        Index a0, a1 = 0, a2, b1, b2;
                        while ((a2 = next(a1)) != tour.size()-1) { 
                            b1 = next(a2);
                            b2 = next(b1);
                            while (b1 != 0) {
                                if (distance(a1, a2) + distance(b1, b2) > distance(a1, b1) + distance(a2, b2) + swapEpsilon) {
                                    reverse(tour.begin()+a2, tour.begin()+b1+1);
                                    reverseCount++;
                                    start_again = true;
                                }
                                else {
                                    a0 = (a1 + tour.size() - 1)%tour.size();
                                    if (distance(a0, a2) + distance(b1, a1) + distance(a1, b2) < 
                                        distance(a0, a1) + distance(a1, a2) + distance(b1, b2) - swapEpsilon) {
                                        //cout << "before: " << perimeter(ps, tour, true);
                                        City c = tour[a1];
                                        tour.erase(tour.begin() + a1);
                                        // cuz of erase shift
                                        tour.insert(tour.begin() + b2-1, c);
                                        //cout << " after: " << perimeter(ps, tour, true);
                                    }
                                } 
                                b1 = b2;
                                b2 = next(b1);
                            }
                            a1 = a2;
                        }
                        if (reverseCount > printedReverseCount+1000) {
                            std::cout << "reverce count: " << reverseCount << " obj: " << Perimeter(ps, tour, true)  << std::endl;
                            printedReverseCount = reverseCount;
                        }
                    } while (start_again);
                    return tour;
                }
                
            };
            
            struct TSP_ThreeOpt : TSP_Improver {
                
                double swapEpsilon;
                
                TSP_ThreeOpt() : swapEpsilon(1e-7) {}
                
                std::vector<City> improve(const std::vector<Point>& ps,
                                          const std::vector<City>& in_tour) override {
                    std::vector<City> tour = in_tour;
                    std::function<Index(Index)> next = [&](Index i) { 
                        return (i+1)%tour.size(); 
                    };
                    std::function<double(Index, Index)> distance = [&](Index i0, Index i1) {
                        return ps[tour[i0]].distance(ps[tour[i1]]);
                    };
                    bool start_again;
                    do {
                        start_again = false;
                        Index a1 = 0, a2, b1, b2, c1, c2;
                        while ((a2 = a1+1) != tour.size()-1-1) {
                            b1 = a2;
                            while ((b2 = b1+1) != tour.size()-1) {
                                c1 = b2;
                                // want to capture first-last edge
                                while ((c2 = (c1+1)%tour.size()) != 1) {
                                    bool did_transform = true;
                                    // 4 possible permutations
                                    double d_a = distance(a1, a2),
                                    d_b = distance(b1, b2),
                                    d_c = distance(c1, c2),
                                    total = d_a + d_b + d_c;
                                    if (total > distance(a1, c1) + distance(a2, b2) + distance(b1, c2) + swapEpsilon) {
                                        reverse(tour.begin()+a2, tour.begin()+b1+1);
                                        reverse(tour.begin()+a2, tour.begin()+c1+1);
                                        //cout << 0 << endl;
                                    } else
                                        if (total > distance(a1, b2) + distance(a2, c2) + distance(b1, c1) + swapEpsilon) {
                                            reverse(tour.begin()+b2, tour.begin()+c1+1);
                                            reverse(tour.begin()+a2, tour.begin()+c1+1);
                                            //cout << 1 << endl;
                                        } else
                                            if (total > distance(a1, b2) + distance(a2, c1) + distance(b1, c2) + swapEpsilon) {
                                                reverse(tour.begin()+a2, tour.begin()+b1+1);
                                                // now we have variant 2 just copy paste
                                                reverse(tour.begin()+b2, tour.begin()+c1+1);
                                                reverse(tour.begin()+a2, tour.begin()+c1+1);
                                                //cout << 2 << endl;
                                            } else
                                                if (total > distance(a1, b1) + distance(a2, c1) + distance(b2, c2) + swapEpsilon) {
                                                    reverse(tour.begin()+a2, tour.begin()+b1+1);
                                                    reverse(tour.begin()+b2, tour.begin()+c1+1);
                                                    //cout << 3 << endl;
                                                    
                                                } else did_transform = false;
                                    
                                    if (did_transform) {
                                        start_again = true;
                                        //    cout << "obj: " << perimeter(ps, tour, true)  << endl;
                                    }
                                    c1 = c2;
                                }
                                b1 = b2;
                            }
                            a1 = a2;
                        }
                    } while(start_again);
                    
                    return tour;
                }
            };
            
        } // end namespace tsp
        
    } // end namespace opt
    
} // end namespace ant


    

namespace ant {
    
    namespace opt {
        
        namespace tsp {
            
            struct TSP_RandomInsertion : TSP_InsertionSolver {
                
                struct Vector : std::vector<City> {
                    Vector(const std::vector<City>& cities) : std::vector<City>(cities) {}
                    
                    void insertBefore(Vector::iterator position, City city) {
                        insert(position, city);
                    }
                    
                    void insertBefore(Index i, City city) {
                        insertBefore(begin()+i, city);
                    }
                    
                    void insertAfter(Vector::iterator position, City city) {
                        insert(position+1, city);
                    }
                    
                    void insertAfter(Index i, City city) {
                        insertAfter(begin()+i, city);
                    }
                    
                    bool hasCity(City c) {
                        return find(begin(), end(), c) != end();
                    }
                };
                
                std::vector<City> solve(const std::vector<Point>& points) override {
                    Vector tour = startingTour(points);
                    std::vector<City> outCities(points.size());
                    iota(outCities.begin(), outCities.end(), 0);
                    // delete in tour cities
                    for (int i = (int)points.size()-1; i >= 0; i--) {
                        if (tour.hasCity(outCities[i])) {
                            std::swap(outCities[i], outCities.back());
                            outCities.pop_back();
                        }
                    } 
                    
                    for (size_t i = tour.size(); i < points.size(); i++) {
                        Index ind_out = random(outCities.size());
                        City city_out = outCities[ind_out];
                        
                        Index ind_in = NearestPoint(points, tour.begin(), tour.end(), points[city_out]) - tour.begin();
                        City city_in = tour[ind_in];
                        
                        outCities.erase(outCities.begin()+ind_out);
                        
                        Index ind_prev = (ind_in + tour.size() - 1)%tour.size();
                        Index ind_next = (ind_in + 1)%tour.size();
                        City city_prev = tour[ind_prev];
                        City city_next = tour[ind_next];
                        
                        
                        auto &ps = points;
                        // is (prev city - city out - city in) more optimum chain
                        if (-ps[city_prev].distance(ps[city_in]) 
                            +ps[city_prev].distance(ps[city_out]) <  
                            
                            -ps[city_next].distance(ps[city_in]) 
                            +ps[city_next].distance(ps[city_out])) {
                            
                            tour.insertBefore(ind_in, city_out);
                        } else {
                            tour.insertAfter(ind_in, city_out);
                        }
                    }
                    return tour;
                }
            };
        } // end namespace tsp
        
    } // end namespace opt
    
} // end namespace ant




using namespace std;
using namespace ant;
using namespace ant::geometry::d2::f;
using namespace ant::opt::tsp;
using ant::geometry::d2::i::Segment;
using namespace ant::geometry::d2;

class SmallPolygons {
public:
    // return zero index indices of vertices for each polygon
    // they space-separated
    vector<string> choosePolygons(vector<int> points, int N) {
        vector<Point> ps(points.size()/2);
        for (int i = 0; i < points.size()/2; ++i) {
            ps[i].set(points[2*i], points[2*i+1]);
        }
        vector<Index> s;
        while (true) {
            TSP_RandomInsertion rand;
            s = rand.solve(ps);
            TSP_TwoOpt two;
            two.swapEpsilon = 0.1;
            s = two.improve(ps, s);
            
            vector<TSP::Edge> edges = TSP::tourToEdges(s);
            vector<Segment> ss(edges.size());
            for (int i = 0; i < ss.size(); ++i) {
                ss[i].fst = i::Point(ps[edges[i][0]].x, ps[edges[i][0]].y);
                ss[i].snd = i::Point(ps[edges[i][1]].x, ps[edges[i][1]].y);
            }
            bool has_intersections = false;
            for (int i = 0; i < edges.size(); ++i) {
                for (int j = i+1; j < edges.size(); ++j) {
                    if (ss[i].Intersect(ss[j]) 
                        && ss[i].fst != ss[j].fst && ss[i].fst != ss[j].snd
                        && ss[i].snd != ss[j].fst && ss[i].snd != ss[j].snd) {
                        has_intersections = true;
                        goto out;
                        cout << ss[i].fst.x << " " << ss[i].fst.y << endl;
                        cout << ss[i].snd.x << " " << ss[i].snd.y << endl;
                        cout << ss[j].fst.x << " " << ss[j].fst.y << endl;
                        cout << ss[j].snd.x << " " << ss[j].snd.y << endl;
                        throw logic_error("segments interect each other");   
                    }
                }
            }
            out:
            if (!has_intersections) {
                break;
            }
        }
        string result = "";
        for (int k : s) {
            result += std::to_string(k) + " ";
        }
        result.erase(result.end()-1);
        if (result[result.size()-1] == ' ') invalid_argument("fck");
        return {result};
    }
};
