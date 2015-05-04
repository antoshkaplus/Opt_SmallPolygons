
#include <vector>
#include <tuple>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <cassert>
#include <iostream>
#include <unordered_set>
#include <queue>
#include <memory>
#include <array>
#include <chrono>
#include <algorithm>
#include <sys/time.h>

unsigned getTickCount()
{
#ifdef WINDOWS
    return GetTickCount();
#else
    struct timeval tv;
    gettimeofday(&tv, 0);
    return unsigned((tv.tv_sec * 1000) + (tv.tv_usec / 1000));
#endif
}


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
    
    
    enum struct Enabler {}; 
    constexpr Enabler enabler = {};
    
    // need to avoid usage with one or 0 parameters sometimes
    
    template<typename... Nothing> 
    struct All {
        static constexpr bool value = true;
    };   
    template<typename Condition, typename... OtherConditions>
    struct All<Condition, OtherConditions...> {
        static constexpr bool value = Condition::value && All<OtherConditions...>::value;
    };
    
    template<typename... Nothing> 
    struct Any {
        static constexpr bool value = true;
    };   
    template<typename Condition, typename... OtherConditions>
    struct Any<Condition, OtherConditions...> {
        static constexpr bool value = Condition::value || All<OtherConditions...>::value;
    };
    
    
    template<typename Condition>
    using EnableIf = typename std::enable_if<Condition::value, Enabler>::type;
    
    template<typename Condition, typename... OtherConditions>
    using EnableIfAll = EnableIf<All<Condition, OtherConditions...>>;
    
    template<typename Condition, typename... OtherConditions>
    using EnableIfAny = EnableIf<Any<Condition, OtherConditions...>>;
    
    
    template<typename... Nothing>
    struct IsAnySame {
        static constexpr bool value = false;
    };
    template<typename Type, typename Another, typename... Other>
    struct IsAnySame<Type, Another, Other...> {
        static constexpr bool value = std::is_same<Type, Another>::value || IsAnySame<Type, Other...>::value;
    };
    
    
    
    template<class Key, class Value>
    std::tuple<std::vector<Key>, std::vector<Value>> Zip(std::map<Key, Value>& m) {
        std::tuple<std::vector<Key>, std::vector<Value>> r;
        auto& keys = std::get<0>(r);
        auto& values = std::get<1>(r);  
        keys.reserve(m.size());
        values.reserve(m.size());
        for (auto& p : m) {
            keys.push_back(p.first);
            values.push_back(p.second);
        }
        return r;
    }
    
    // sometimes someone can use with Long, not just Int type
    template<class T>
    class Range {
    public:
        class Iterator : std::iterator<std::input_iterator_tag, T> {
        public:
            Iterator(const Range& range, T current) 
            : range_(range), current_(current) {
                // just copied that lol
                if (range_.step_*(current_-range_.last_) > 0) current_ = range_.last_;
            }
            const T operator*() const { return current_; }
            bool operator==(const Iterator& it) const {
                return current_ == *it;
            }
            bool operator!=(const Iterator& it) const {
                return current_ != *it;
            }
            Iterator& operator++() {
                current_ += range_._step;
                if (range_.step_*(current_-range_.last_) > 0) current_ = range_.last_;
                return *this;
            }
            Iterator operator++(int) { 
                Iterator it(*this); 
                operator++(); 
                return it;
            }
        private:
            const Range& range_;
            T current_;
        };
        
        friend class Iterator;
        
        Range(T last) : first_(0), last_(last), step_(1) {}
        Range(T first, T last, T step = 1)
        : first_(first), last_(last), step_(step) {}
        
        Iterator begin() const { return Iterator(*this, first_); }
        Iterator end()   const { return Iterator(*this, last_); }
        
    private:
        T first_, last_, step_;
    };
    
    
    
    class DisjointSet {
    public:
        DisjointSet() {}
        DisjointSet(Count element_count) {
            Init(element_count);
        }
        
        void Init(Count element_count) {
            element_count_ = element_count;
            set_count_ = element_count;
            data_.resize(element_count);
            size_.resize(element_count);
            Reset();
        }
        
        void Reset() {
            std::iota(data_.begin(), data_.end(), 0);
            fill(size_.begin(), size_.end(), 1);
        }
        
        void Unite(Index i_0, Index i_1) {
            --set_count_;
            Index
            r_0 = root(i_0),
            r_1 = root(i_1);
            // will join r_0 to r_1, so r_1 height should be bigger
            if (size_[r_0] > size_[r_1]) {
                std::swap(r_0, r_1);
            }
            data_[r_0] = r_1;
            size_[r_1] += size_[r_0];
            
        }
        
        bool is_separate(Index i_0, Index i_1) {
            return root(i_0) != root(i_1);
        }
        
        Index root(Index i) {
            while (i != data_[i]) {
                i = data_[i] = data_[data_[i]];
            }
            return i;
        }
        
        size_t size() {
            return element_count_;
        }
        
        Count set_count() {
            return set_count_;
        }
        
    private:
        Count element_count_;
        Count set_count_;
        std::vector<Index> data_;
        // how many elements in set with index, if index is root
        std::vector<size_t> size_;
    }; 
    
    
    
    
    //struct discrete_distribution {
    //
    //    template<class ForwardIterator>
    //    discrete_distribution(ForwardIterator first, ForwardIterator last) {
    //        weight_.assign(first, last);
    //        cumulative_.resize(weight_.size());
    //        std::partial_sum(weight_.begin(), weight_.end(), cumulative_.begin());
    //        uni_ = std::uniform_real_distribution<double>(0, cumulative_.back());
    //    }
    //    discrete_distribution(const std::initializer_list<double>& weights) 
    //    : discrete_distribution(weights.begin(), weights.end()) {}
    //
    //    void set_weight(Index i, double w) {
    //        assert(w >= 0);
    //        //            auto d = w-weight_[i];
    //        //            for (auto k = i; k < weight_.size(); ++k) {
    //        //                cumulative_[k] += d;
    //        //            }
    //        //            
    //        weight_[i] = w;
    //        std::fill(cumulative_.begin(), cumulative_.end(), 0);
    //        std::partial_sum(weight_.begin(), weight_.end(), cumulative_.begin());
    //        uni_ = std::uniform_real_distribution<double>(0, cumulative_.back());
    //    }
    //
    //    double get_weight(Index i) {
    //        return weight_[i];
    //    }
    //
    //    template<class RNG> 
    //    Index operator()(RNG& rng) {
    //        Index i = std::lower_bound(cumulative_.begin(), cumulative_.end(), uni_(rng))-cumulative_.begin();
    //        if (cumulative_.back() != 0.) while ( weight_[i] == 0.) --i;
    //        return i;
    //    }
    //
    //    std::uniform_real_distribution<double> uni_;
    //    std::vector<double> cumulative_;
    //    std::vector<double> weight_;
    //};
    
    
    
    // current stack supports iteration!
    template<class T>
    class Stack {
    public:
        // can't inherit from vector iterator
        // too open class
        class ConstIterator : std::iterator<std::input_iterator_tag, T> {
        public:
            ConstIterator(const Stack& stack, typename std::vector<T>::const_iterator current) 
            : stack_(stack), current_(current) {}
            const T& operator*() const { return *current_; }
            bool operator==(const ConstIterator& it) const {
                return current_ == it;
            }
            bool operator!=(const ConstIterator& it) const {
                return current_ != it.current_;
            }
            ConstIterator& operator++() {
                ++current_;
                return *this;
            }
            // post iterator
            ConstIterator operator++(int) { 
                ConstIterator it(*this); 
                operator++(); 
                return it;
            }
        private:
            const Stack& stack_;
            typename std::vector<T>::const_iterator current_;
        };
        
        friend class ConstIterator;
        
        ConstIterator begin() const { return ConstIterator(*this, data_.begin()); }
        ConstIterator end()   const { return ConstIterator(*this, data_.end()); }
        
        T& top() {
            return data_.back();
        }
        void push(const T& val) {
            data_.push_back(val);
        }
        void pop() {
            data_.pop_back();
        }
        bool empty() const {
            return data_.empty();
        }
        size_t size() const {
            return data_.size();
        }
        
    private:
        std::vector<T> data_;
    };
    
    
    
    // probably should be an inheritance
    template <class T>
    class CountMap : public std::map<T, Count> {
    public:
        void decrease(const T& t) { decrease(t, 1); }
        
        void decrease(const T& t, Count val) {
            auto it = this->find(t);
            if ((it->second-=val) == 0) {
                this->erase(it);
            }
        }
        
        void increase(const T& t) { increase(t, 1); }
        
        void increase(const T& t, Count val) {
            this->emplace(t, 0).first->second+=val;
        }
        
        std::set<T> keys() const {
            std::set<T> r;
            for (auto p : *this) {
                r.insert(p.first);
            }
            return r;
        }
        
        size_t get(const T& t) const {
            auto it = this->find(t);
            return it == this->end() ? 0 : it->second;
        }
    };
    
    template<class T>
    class CircularList {
    private:
        struct Node {
            T value;
            Node* prev;
            Node* next;
        };
        
        template<class V>
        struct BaseIterator : std::iterator<std::bidirectional_iterator_tag, V> {
            BaseIterator(Node* n) : node_(n) {} 
            
            V& operator*() const { return node_->value; }
            V* operator->() const { return node_->value; }
            
            bool operator==(const BaseIterator& it) const {
                return node_ == it.node_;
            }
            bool operator!=(const BaseIterator& it) const {
                return node_ != it.node_;
            }
            BaseIterator& operator++() {
                node_ = node_->next;
                return *this;
            }
            // post iterator
            BaseIterator operator++(int) { 
                BaseIterator it(node_); 
                node_ = node_->next; 
                return it;
            }
            
            BaseIterator operator--() {
                node_ = node_->prev;
                return *this;
            }
            BaseIterator operator--(int) {
                BaseIterator it(node_);
                node_ = node_->prev;
                return it;
            }
            
        private:
            Node* node_;
            friend struct circular_list;
        };
        
    public:
        using Iterator = BaseIterator<T>;
        using ConstIterator = BaseIterator<const T>;
        
        CircularList() : focus_(nullptr) {}
        
        void InitFocus(const T& value) {
            focus_ = new Node();
            focus_->value = value;
            focus_->prev = focus_->next = focus_;
        }
        
        template<typename It, EnableIf<IsAnySame<It, Iterator, ConstIterator>>>
        It InsertAfter(It it_pos, const T& value) {
            ++count_;
            Node *pos = it_pos.node_;
            if (pos == nullptr) {
                init_focus(value);
                return Iterator(focus_);
            }
            Node *prev = pos;
            Node *next = pos->next;
            Node *cur = new Node();
            cur->next = next;
            cur->prev = prev;
            cur->value = value;
            next->prev = cur;
            prev->next = cur;
            return It(cur);
        }
        
        template<typename It, EnableIf<IsAnySame<It, Iterator, ConstIterator>>>
        It InsertBefore(Iterator it_pos, const T& value) {
            ++count_;
            Node *pos = it_pos.node_;
            if (pos == nullptr) {
                init_focus(value);
                return Iterator(focus_);
            }
            Node *prev = pos->prev;
            Node *next = pos;
            Node *cur = new Node();
            cur->next = next;
            cur->prev = prev;
            cur->value = value;
            prev->next = cur;
            next->prev = cur;
            return It(cur);
        }
        
        template<typename It, EnableIf<IsAnySame<It, Iterator, ConstIterator>>>
        void Erase(It it_pos) {
            --count_;
            Node *pos = it_pos.node_;
            if (pos == focus_) focus_ = pos->next;
            pos->prev->next = pos->next;
            pos->next->prev = pos->prev;
            delete pos;
            if (count_ == 0) focus_ = nullptr;
        }
        
        Iterator focus() {
            return Iterator(focus_);
        }
        
        ConstIterator focus() const {
            return ConstIterator(focus_);
        }
        
        Count size() const {
            return count_;
        }
        
        bool is_empty() const {
            return focus_ == nullptr;
        }
        
    private:
        Node* focus_;
        Count count_ = 0;
    };
    
    
    
    std::map<std::string, std::string> command_line_options(const char* argv[], int argc);
    
    int atoi(char* str, Int n);
    int atoi(char* first, char *last);
    
    
    struct command_line_parser {
        command_line_parser(const char* argv[], int argc) {
            options_ = command_line_options(argv, argc);
        }
        
        bool exists(const std::string& option) const {
            return options_.find(option) != options_.end();
        }
        
        bool hasValue(const std::string& option) const {
            return options_.at(option) != "";
        }
        
        std::string getValue(std::string option) const {
            std::string value = options_.at(option);
            if (value == "") {
                throw std::logic_error("command line option has no value");
            }
            return value;
        }
    private:
        std::map<std::string, std::string> options_;
    };
    
    
    // let it be unsigned char, int or long
    template<class T>
    struct CombinationGenerator {
        // should be like iterator
        
        static const Int kByteBitCount = 8;
        static const Int kMaxElementCount = sizeof(T)*kByteBitCount;
        
        struct Tails : std::array<T, kMaxElementCount> {
            using std::array<T, kMaxElementCount>::at;
            
            Tails() {
                at(0) = 1;
                init(1);
            }
            
            void init(Int i) {
                if (i == kMaxElementCount) return;
                at(i) = at(i-1) << 1;
                at(i) += 1;
                init(i+1);
            } 
        };
        
        static constexpr Tails tails = Tails();
        
        CombinationGenerator(Int selection_count, Int element_count, T starting_combination = 0)
        : selection_count_(selection_count), element_count_(element_count), data_(starting_combination) {
            if (data_ == 0) {
                data_ = tails[selection_count-1];
            }
        }
        
        const T& next() {
            Int i = 0;
            Int k = 0; // how many elements behind
            while (i < element_count_) {
                
                if ((data_ >> i) && 1) {
                    if (!((data_ >> (i+1)) & 1)) {
                        //data_ |= (1 << i+1);
                        data_ >>= i+1;
                        data_ += 1;
                        data_ <<= i+1;
                        data_ |= tails[k-1];
                        break;
                    }
                    ++k;
                }
                // can shift father
            }
            return data_;
        } 
        
        bool hasNext() {
            return data_ != (tails[selection_count_-1] << (element_count_ - selection_count_)); 
        }
        
        T data_;
        Int element_count_;
        Int selection_count_;
    };
    
    template<class T>
    struct BinomialHeap {
        const T& min() {
            auto& x = data_[0];
            auto* min = data_;
            while (x < data_.size()) {
                if (data_[x] < min) {}
            } 
        }
        std::vector<T> data_;
    };
    
    
    
    template<class T>
    struct bst_set {
    public:
        
        bst_set() : size_(0), root_(0) {}
        
        virtual ~bst_set() {
            clear();
        }
        
        void clear() {
            node::clear(root_);
            root_ = nullptr;
            size_ = 0;
        }
        
        size_t size() const {
            return size_;
        }
        
        
    private:
        
        struct node {
            
            node(const T& t) : value_(t) {}
            
            // well if you call on null node it's your problem
            static node* next(node* n) {
                if (exists(right(n))) {
                    n = right(n);
                    while (exists(left(n))) {
                        n = left(n);
                    }
                    return n;
                }
                node* n_2 = parent(n);
                if (!exists(n_2)) return nullptr;
                
                if (left(n_2) == n) {
                    return n_2;
                } else {
                    // n right child
                    // will return nullptr if can't find anything
                    while (exists(left(n_2)) && left(n_2) != n) {
                        n = n_2;
                        n_2 = parent(n_2);
                    }
                    return n_2;
                }
            }
            
            static node* prev(node* n) {
                // like next but change left and right functions
                if (exists(left(n))) {
                    n = left(n);
                    while (exists(right(n))) {
                        n = right(n);
                    }
                    return n;
                }
                node* n_2 = parent(n);
                if (!exists(n_2)) return nullptr;
                
                if (right(n_2) == n) {
                    return n_2;
                } else {
                    // n right child
                    // will return nullptr if can't find anything
                    while (exists(right(n_2)) && right(n_2) != n) {
                        n = n_2;
                        n_2 = parent(n_2);
                    }
                    return n_2;
                }
                
            }
            
            // min element in subtree
            // can't be called will nullptr argument
            static node* min(node* n) {
                while (exists(left(n))) {
                    n = left(n);
                }
                return n;
            }
            
            static node* max(node* n) {
                while (exists(right(n))) {
                    n = right(n);
                }
                return n;
            }
            
            static void substitute_child(node* parent, node* child, node* substitution) {
                if (substitution != nullptr) substitution->parent_ = parent;
                if (left(parent) == child) {
                    parent->left_ = substitution;
                } else {
                    parent->right_ = substitution;
                }
            }
            
            
            
            static node* right(node* n) {
                return n->right_;
            }
            static node* parent(node* n) {
                return n->parent_;
            }
            static node* left(node* n) {
                return n->left_;
            }
            static bool exists(node* n) {
                return n != nullptr;
            }
            
            // clears whole subtree
            static void clear(node* n) {
                if (!exists(n)) return;
                // need to find out is right or left child then null everything
                if (exists(parent(n))) {
                    if (n->parent_->left == n) n->parent_->left = nullptr;
                    else n->parent_->right = nullptr;
                    n->parent_ = nullptr;
                }
                // now just clear separate tree with n
                while (n != nullptr) {
                    if (exists(left(n))) {
                        n = left(n);
                    } else if (exists(right(n))) {
                        n = right(n);
                    } else {
                        node* n_old = n;
                        // if nullptr will finish
                        n = parent(n);
                        delete n_old;
                    }
                }
            }
            
            node* right_    = nullptr;
            node* left_     = nullptr;
            node* parent_   = nullptr; 
            T value_;
        }; 
        
    public:    
        
        struct iterator : std::iterator<std::bidirectional_iterator_tag, node> {
            
            iterator(node* node) : current_(node) {}
            iterator() {}
            
            const T& operator*() const { 
                return current_->value;  
            } 
            
            bool operator==(const iterator& it) {
                return current_ == it.current_; 
            }
            bool operator!=(const iterator& it) {
                return current_ != it.current_;
            }
            
            // pred
            iterator& operator++() {
                current_ = node::next(current_);
                return *this;
            }
            // post
            iterator operator++(int) { 
                iterator it(current_);
                current_ = node::next(current_); 
                return it;
            }
            
        private:
            node* current_;
        };
        
        
        iterator begin() const {
            if (root_ == nullptr) return end();
            node* b = root_;
            while (node::exists(node::left(b))) {
                b = node::left(b);
            }
            return iterator(b);
        } 
        
        // probably should just return nullptr as node
        iterator end() const {
            return iterator(nullptr);
        }
        
        iterator find(const T& t) const {
            node* n = root_; 
            while (n != nullptr) {
                if (n->value_ < t) {
                    n = n->right_;
                    continue;
                }
                if (t < n->value_) {
                    n = n->left_;
                    continue;
                }
                return iterator(n);
            }
            return end();
        }
        
        bool exists(const T& t) const {
            return find(t) != end();
        }
        
        // like find but should keep track of element parent
        std::pair<iterator, bool> insert(const T& t) {
            node* n_new = new node(t); 
            if (root_ == nullptr) {
                root_ = n_new;
                size_ = 1;
                return {begin(), true};
            }
            node* n = root_;
            while (true) {
                if (n->value_ < t) {
                    if (n->right_ == nullptr) {
                        break;
                    }
                    n = n->right_;
                    continue;
                }
                if (t < n->value_) {
                    if (n->left_ == nullptr) {
                        break;
                    }
                    n = n->left_;
                    continue;
                }
                return {iterator(n), false};
            }
            // now should have not null parent
            n_new->parent_ = n;
            if (t < n->value_) {
                n->left_ = n_new;
            } else {
                // equality could not be
                n->right_ = n_new;
            }
            ++size_; 
            return {iterator(n_new), true};
        }
        
        // should return iterator on next element or end
        void erase(const T& t) {
            erase(find(t));
        }
        
        void erase(iterator it) {
            if (it == end()) return;
            --size_;
            node* n = it.current_;
            // what happens if n is root???
            if (n->right_ == nullptr && n->left_ == nullptr) {
                // ??? root_ goes to nullptr 
                
                node::substitute_child(n->parent_, n, nullptr);
                delete n;
                return;
            }
            if (n->right_ == nullptr && n->left_ != nullptr) {
                // ??? root_ goes to child
                
                node::substitute_child(n->parent_, n, n->left_);
                delete n;
                return;
            } 
            if (n->left_ == nullptr && n->right_ != nullptr) {
                // ??? root_ goes to child
                
                node::substitute_child(n->parent_, n, n->right_);
                delete n;
                return;
            }
            // ok, both children presents
            // ??? root_ goes to child
            if (std::uniform_int_distribution<>(0, 1)(rng_) == 0) {
                // first left to right
                node::substitute_child(n->parent_, n, n->right_);
                node* p = node::min(n->right_);
                n->left_->parent_ = p;
                p->left_ = n->left_;
            } else { 
                // second right to left
                node::substitute_child(n->parent_, n, n->left_);
                node* p = node::max(n->left_);
                n->right_->parent_ = p;
                p->right_ = n->right_;
            }
            delete n;
        }
        
        private :
        
        // can use in find and insert
        iterator find_closest(const T& t) {
            return end();
        }
        
        std::default_random_engine rng_{(unsigned)std::chrono::system_clock::now().time_since_epoch().count()};
        size_t size_;
        node* root_;
    };
    
    //
    //
    //P LogicBinarySearch(P a, P b, Func func, Distance dist, Condition cond, double eps) {
    //    
    //    bool v_a = func(a);
    //    bool v_b = func(b);
    //    while (dist(b, a) > eps) {
    //        if (cond())
    //    }
    //    
    //    
    //    
    //}
    
    template<class T>
    uint64_t Hash(T c_0, T c_1, T c_2, T c_3) {
        uint64_t r = 0;
        r += c_0;
        r <<= 16;
        r += c_1;
        r <<= 16;
        r += c_2;
        r <<= 16;
        r += c_3;
        return r;
    }
    
    
    
    
}

namespace ant { 
    
    std::map<std::string, std::string> command_line_options(const char* argv[], int argc) {
        std::map<std::string, std::string> res;
        for (Index i = 0; i < argc; ++i) {
            if (argv[i][0] == '-') {
                std::string key(argv[i]+1);
                res[key] = "";
                if (i+1 < argc && argv[i+1][0] != '-') {
                    res[key] = std::string(argv[i+1]);
                    ++i;
                } 
            }
        }
        return res;
    } 
    
    int atoi(char* str, Int n) {
        return atoi(str, str+n);
    }
    
    int atoi(char* first, char *last) {
        char ch = *last;
        *last = '\0';
        int r = std::atoi(first);
        *last = ch;
        return r;
    }
    
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
                    
                private:
                    // are listed in counterclockwise order
                    // later should put it outside
                    static bool CCW(const Point& A, const Point& B, const Point& C) {
                        return (C.y-A.y) * (B.x-A.x) > (B.y-A.y) * (C.x-A.x);
                    }
                    
                public:
                    bool Intersect(const Segment& s) const {
                        return CCW(fst, s.fst, s.snd) != CCW(snd, s.fst, s.snd) && 
                        CCW(fst, snd, s.fst) != CCW(fst, snd, s.snd);
                    }
                    
                    bool Lie(Point q) const
                    {
                        return (q.x <= std::max(fst.x, snd.x) && q.x >= std::min(fst.x, snd.x) &&
                                q.y <= std::max(fst.y, snd.y) && q.y >= std::min(fst.y, snd.y));
                    }
                    
                    bool IntersectOrLie(const Segment& s) const {
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
            
            
            //
            //template<class P> 
            //double HeronFormula(P p_0, P p_1, P p_2) {
            //    doi
            //    return 
            //}
            
            
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
            void K_NearestPoints(const std::vector<P>& ps, const P& p, 
                                 std::vector<Index>& inds, int k) {
                k = std::min(k, (int)inds.size());
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
                    return std::abs(s)/2;
                }
                
                int ShoelaceFormula(const std::vector<Point>& ps, const std::vector<Index>& order) {
                    double s = 0;
                    s += ps[order[order.size()-1]].x*ps[order[0]].y 
                    - ps[order[0]].x*ps[order[order.size()-1]].y;
                    for (auto i = 0; i < order.size()-1; ++i) {
                        s += ps[order[i]].x*ps[order[i+1]].y - ps[order[i+1]].x*ps[order[i]].y;
                    }
                    return std::abs(s)/2;
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
                
                virtual std::vector<City> solve(const std::vector<Point>& points) = 0;
                
                static std::vector<Edge> tourToEdges(const std::vector<City>& cities) {
                    std::vector<Edge> edges;
                    for (Index i = 0; i < cities.size(); i++) {
                        edges.push_back(Edge(cities[i], cities[(i+1)%cities.size()]));   
                    }
                    return edges; 
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
    
    namespace grid {
    
        
        template<class T>
        struct Grid {
            using value_type = T;
            
            typename std::vector<T>::const_iterator begin() const {
                return grid_.begin();
            }
            typename std::vector<T>::iterator begin() {
                return grid_.begin();
            }
            typename std::vector<T>::const_iterator end() const {
                return grid_.end();
            }
            typename std::vector<T>::iterator end() {
                return grid_.end();
            }
            
            
            Grid() : Grid(0, 0) {}
            Grid(Count row_count, Count col_count)
            :   row_count_(row_count), 
            col_count_(col_count), 
            grid_(row_count*col_count) {}
            
            Grid(Count row_count, Count col_count, const T& value_type) 
            :   row_count_(row_count),
            col_count_(col_count),
            grid_(row_count*col_count, value_type) {}
            
            Grid(std::initializer_list<std::vector<T>> list)
            :   Grid(list.size(), list.size() == 0 ? 0 : list.begin()->size()){
                auto b = list.begin();
                for (auto r = 0; r < row_count(); ++r, ++b) {
                    std::copy(b->begin(), b->end(), grid_.begin() + r*col_count());
                }
            }
            bool isInside(Int row, Int col) {
                return row >= 0 && row < row_count_ && 
                col >= 0 && col < col_count_;
            }
            
            void resize(Count row_count, Count col_count) {
                row_count_ = row_count;
                col_count_ = col_count;
                grid_.resize(row_count*col_count);
            }
            
            void fill(const T& t) {
                std::fill(grid_.begin(), grid_.end(), t);
            }
            
            Count row_count() const { return row_count_; }
            Count col_count() const { return col_count_; }
            Count cell_count() const { return row_count()*col_count(); } 
            
            
            T& operator()(Int row, Int col) {
                return grid_[row*col_count_ + col];
            }
            const T& operator()(Int row, Int col) const {
                return grid_[row*col_count_ + col];
            }
            
        private:
            Count row_count_, col_count_;
            std::vector<T> grid_;
            
            friend struct const_iterator;
            template<class U>
            friend bool operator==(const Grid<U>& g_0, const Grid<U>& g_1);
        };
        

    
    }
    
}







using namespace std;
using namespace ant;
using namespace ant::geometry::d2::f;
using namespace ant::opt::tsp;
using Segment = ant::geometry::d2::i::Segment;
using namespace ant::geometry::d2;


using AdjacentEdges = vector<pair<City, City>>;
using Polygons = vector<vector<City>>;
using Intersections = vector<pair<TSP::Edge, TSP::Edge>>;
using Tour = std::vector<City>;
using Edge = TSP::Edge;
using Seed = std::default_random_engine::result_type;



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
                        !Segment{ps[e[0]], ps[e[1]]}.Intersect(
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


vector<bool> ComputeVisitedCities(const vector<Tour>& tours, Count city_count) {
    vector<bool> visited(city_count, false);
    for (auto& t : tours) {
        for (auto c : t) {
            visited[c] = true;
        }
    }
    return visited;
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
                        
                        if (!Segment{ps[e[0]], ps[e[1]]}.Intersect(
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
                                   Count max_tour_count, Seed seed) {
    std::default_random_engine rng(seed);
    Count city_count = closest_cities.row_count();
    vector<City> cities(city_count);
    iota(cities.begin(), cities.end(), 0);
    vector<Tour> result;
    vector<bool> visited(city_count, false);
    for (auto t = 0; t < max_tour_count; ++t) {
        if (cities.empty()) break;
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
            continue;
        }
        // not a fact that k_i is inside cities
        // could be removed before
        visited[k_0] = visited[k_1] = visited[k_2] = true;
        result.push_back(vector<City>{k_0, k_1, k_2});
    }
    return result;
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
                      Tour& tour) {
    Count city_count = points.size();
    default_random_engine rng(std::chrono::system_clock::now().time_since_epoch().count());
    vector<bool> vs(points.size());
    unordered_set<int> inters_hashes;
    vector<int> ii;
    while (true) {
        InitVisitedCities({tour}, vs);
        // need to do it everytime
        Intersections inters = FindIntersectionsForTour(points, 
                                                        closest_cities, tour, adj_edges, vs);
        if (inters.empty()) break;
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
    return true;
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




struct SimplexInsertion_1 {
    
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
    
    const std::vector<Point> *points;
    const grid::Grid<City> *closest_cities; 
    const grid::Grid<double> *edge_distance;
    
    Count city_count;
    // get matrix for each point get 10-20 closest, sort them
    // row - each city, col - closest cities
    grid::Grid<char> edge_exists;
    std::vector<bool> visited;
    std::priority_queue<Item> insertion_queue;
    
    double Distance(City c_0, City c_1) {
        return (*edge_distance)(c_0, c_1);
    }
    
    bool Exists(Edge e) {
        return edge_exists(e[0], e[1]);
    }
    
    double Profit(Edge e, City v) {
        double b = Distance(e[0], v);
        double c = Distance(e[1], v);
        return b + c; 
        //        return HeronFormula(Distance(v, e[0]), Distance(v, e[1]), Distance(e[0], e[1]));
    }
    
    void Remove(Edge e) {
        edge_exists(e[0], e[1]) = edge_exists(e[1], e[0]) = false;
    }
    
    
    void Add(Edge e) {
        edge_exists(e[0], e[1]) = edge_exists(e[1], e[0]) = true; 
    }
    
    
    void UpdateQueue(City new_v, City e_0, City e_1) {
        auto& cc = *closest_cities; 
        for (int j = 0; j < cc.col_count(); ++j) {
            int k = cc(new_v, j);
            if (visited[k]) continue;
            // queue is max
            insertion_queue.emplace(Edge{e_0, new_v}, k,
                                    -Profit({e_0, new_v}, k));
            insertion_queue.emplace(Edge{e_1, new_v}, k, 
                                    -Profit({e_1, new_v}, k));
        }
    }
    
    vector<Tour> BuildTours() {
        for (int i = 0; i < city_count; ++i) {
            if (!visited[i]) throw logic_error("city not visited");
        }
        fill(visited.begin(), visited.end(), false);
        vector<Tour> tours;
        for (int k = 0; k < city_count; ++k) {
            if (visited[k]) continue;
            City first = k;
            City cur = k;
            City prev = -1;
            tours.resize(tours.size()+1);
            do {
                tours.back().push_back(cur);
                for (int j = 0; j < city_count; ++j) {
                    if (Exists(Edge{cur, j}) && j != prev) {
                        prev = cur;
                        cur = j;
                        break;
                    }
                }
            } while (cur != first);
            for (auto t : tours.back()) {
                visited[t] = true;
            }
        }
        if (find(visited.begin(), visited.end(), false) != visited.end()) {
            throw logic_error("tours are not correct");
        }
        return tours; 
    }
    
    void InitializeWithStartingTours(const vector<Tour>& starting_tours) {
        for (const Tour& t : starting_tours) {
            for (Index i = 0; i < t.size(); ++i) {
                Count sz = t.size();
                City c = t[i];
                City prev = t[(i-1+sz) % sz];
                City next = t[(i+1) % sz];
                Add({c, next});
                visited[c] = true;
                UpdateQueue(c, prev, next);
            }
        }
        
    }
    
    
    std::vector<Tour> Solve(const std::vector<Point>& points, 
                            const grid::Grid<City>& closest_cities,
                            const grid::Grid<double>& edge_distance, 
                            const vector<Tour>& starting_tours) {
        this->edge_distance = &edge_distance;
        this->points = &points;
        this->closest_cities = &closest_cities;
        city_count = points.size();
        // can be reused
        visited.resize(city_count);
        fill(visited.begin(), visited.end(), false);
        edge_exists.resize(city_count, city_count);
        edge_exists.fill(false);
        InitializeWithStartingTours(starting_tours);
        while (!insertion_queue.empty()) {
            Item t = insertion_queue.top();
            insertion_queue.pop();
            if (visited[t.city] || !Exists(t.edge)) continue;
            Remove(t.edge);
            Add(Edge{t.edge[0], t.city});
            Add(Edge{t.edge[1], t.city});
            visited[t.city] = true;
            UpdateQueue(t.city, t.edge[0], t.edge[1]);
        }
        return BuildTours();
    }
    
};




class SmallPolygons_2 {
public:
    
    grid::Grid<City> closest_cities;
    Count max_tour_count;
    Count city_count;
    vector<Point> points;
    vector<i::Point> points_int;
    // index of tour to which edge corresponds
    grid::Grid<double> edge_distance;
    grid::Grid<char> edge_valid;
    
    
    vector<Tour> GenerateStartingTours() {
        Seed seed = std::chrono::system_clock::now().time_since_epoch().count();
        auto starting_tours = ::GenerateStartingTours(closest_cities, max_tour_count, seed);
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
    
    
    // return zero index indices of vertices for each polygon
    // they space-separated
    vector<string> choosePolygons(vector<int> point_coordinates, int max_tour_count) {
        unsigned t = getTickCount();
        
        this->max_tour_count = max_tour_count;
        city_count = point_coordinates.size()/2;
        InitializePoints(point_coordinates);
        edge_distance = ComputeEdgeDistance(points);
        closest_cities = ComputeClosestCities(edge_distance, 30);
        edge_valid = ConstructValidEdges(closest_cities);
        // should assign only those that are closest
        SimplexInsertion_1 rand;
        vector<Tour> best_sol_tours;
        double best_area = numeric_limits<double>::max();
        Index iter = 0;
        while (true) {
            vector<Tour> sol_tours;
            while (true) {
                if (getTickCount() - t >= 8000) {
                    goto out;
                }
                
                auto starting_tours = GenerateStartingTours();
                // created initial solution probably infeasible
                sol_tours = rand.Solve(points, closest_cities, edge_distance, starting_tours);
                bool solo_success = true;
                auto adj_edges = ConstructAdjacentEdges(sol_tours, city_count);
                for (auto& t : sol_tours) {
                    auto vs = ComputeVisitedCities({t}, city_count);
                    // adj_edges aren't changing for the particular tour
                    auto inters = FindIntersectionsForTour(points_int, closest_cities, t, adj_edges, vs);
                    if (!inters.empty()) {
                        if (!MakeTourFeasible(points_int, closest_cities, edge_valid, adj_edges, t)) {
                            solo_success = false;
                            break;
                        } 
                    }
                }
                ++iter;
                if (solo_success) {
                    // need to update those dudes
                    adj_edges = ConstructAdjacentEdges(sol_tours, city_count);
                    auto city_groups = ConstructCityGroups(sol_tours, city_count);
                    auto inters = FindIntersectionsBetweenTours(points_int, closest_cities, 
                                                                city_groups, sol_tours, adj_edges);
                    if (inters.empty()) break;
                    
                    if (max_tour_count > 1) {
                        --max_tour_count;
                        cout << "iter #" << iter << " decreasing number of polygons to " 
                        << max_tour_count << endl; 
                    } 
                }
            }
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
};


