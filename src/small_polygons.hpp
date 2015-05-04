
#ifndef SMALL_POLYGONS
#define SMALL_POLYGONS


#include "small_polygons_solver.hpp"


class SmallPolygons {

    // we use slash for the right name of getters
    unsigned time_millis_ = 9800;
    shared_ptr<SimplexInsertion> insertion_{new SimplexInsertion()};
    SmallPolygonsSolver solver;
    
    
    
public:
    vector<string> choosePolygons(vector<int> point_coordinates, int max_tour_count) {
        Count city_count = point_coordinates.size()/2; 
        
        if (city_count < 100) {
            insertion_->set_param_exclude_edge(0.95);
            insertion_->set_param_area(1.45);
            //solver.set_should_improve(true);
        } else if (city_count < 500) {
            insertion_->set_param_exclude_edge(0.87);
            insertion_->set_param_area(1.23);
        } else {
            insertion_->set_param_exclude_edge(0.51);
            insertion_->set_param_area(0.7);
        }
        solver.set_simplex_insertion(*insertion_);
        solver.set_time_millis(time_millis_);
        return solver.choosePolygons(point_coordinates, max_tour_count);
    }
    
    void set_time_millis(unsigned time) {
        time_millis_ = time;
    }
    
    void set_insertion(shared_ptr<SimplexInsertion> insertion) {
        insertion_ = insertion; 
    }
    
    void set_should_improve(bool improve) {
        solver.set_should_improve(improve);
    }
    
    shared_ptr<SimplexInsertion> insertion() {
        return insertion_;
    }
    
    double area() const {
        return solver.area();
    }
};


#endif