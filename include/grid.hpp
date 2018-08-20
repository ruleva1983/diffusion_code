#ifndef GRID_H_
#define GRID_H_

#include <cassert>
#include <vector>

template <class DataType = std::vector<float>>
class Grid1D
{
public:
    
    Grid1D(): dim(0){
        values.resize(0);
        points.resize(0);
        dx = 0.0;
    }
    
    // Constuctor for equally spaced grid of N points
    Grid1D(int N, float x0, float xf, DataType& valuesInit=std::vector<float> (0)): dim(N){
        assert(valuesInit.size() == N);
        dx = (xf-x0)/static_cast<float>(N);
        points.resize(N);
        for (int i=0 ; i < N; ++i)
            points[i] = x0 + i*dx;
        values = valuesInit;
    }
    
    // Returns the value stored in the grid
    float value(int i){
        return values[i];
    }
    
    // Returns the position of the grid point
    float point(int i){
        return points[i];
    }
    
    // Sets the value in a grid point
    void set_value(int i, float value){
        values[i] = value;
    }
    
    // Applies boundary conditions
    void set_boundary(float xmin, float xmax){
        values[0] = xmin;
        values[values.size()-1] = xmax;
    }
    
    DataType& get_state(){
        return values;
    }
    
    float get_deltax(){
        return dx;
    }
    
private:
    int dim;
    float dx;
    DataType points;
    DataType values;
    
};


#endif
