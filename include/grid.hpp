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
    Grid1D(int N, float xL, float xR, DataType& valuesInit): dim(N) {
        assert(xR > xL);
        assert(valuesInit.size() == N);
        dx = (xR-xL)/static_cast<float>(N);
        points.resize(N);
        for (int i=0 ; i < N; ++i)
            points[i] = xL + i*dx;
        values = valuesInit;
    }

    float get_value(int i){
        return values[i];
    }

    float get_point(int i){
        return points[i];
    }

    void set_value(int i, float value){
        values[i] = value;
    }

    void set_boundary(float xmin, float xmax, int imin, int imax){
        assert (imin >= 0);
        assert (imax <= dim - 1);
        values[imin] = xmin;
        values[imax] = xmax;
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
