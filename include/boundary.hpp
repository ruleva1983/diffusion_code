#ifndef BOUNDARY_H_
#define BOUNDARY_H_

template <class BoundaryType=float>
class Boundary1D{
public:
    Boundary1D(float l=0.0, float h=7.0): low(l), high(h){}
    BoundaryType low;
    BoundaryType high;
};




#endif
