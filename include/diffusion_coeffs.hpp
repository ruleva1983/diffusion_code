#ifndef DIFFUSION_COEFFS_H_
#define DIFFUSION_COEFFS_H_

#include <cmath>
#include <functional>

class Coeff1D{
public:
    virtual float operator()(const float x=0.0, const float t=0.0) = 0;
};

class Constant: public Coeff1D{
public:
    Constant (float c=0): a(c)
    {
    }

    float operator()(const float x=0.0, const float t=0.0){
        return a;
    }
private:
    float a;
};

template <typename F=std::function<float(float, float)>>
class Functional: public Coeff1D{
public:
    Functional (F f): function(f)
    {
    }

    float operator()(const float x=0.0, const float t=0.0){
        return function(x, t);
    }
private:
    F function;

};


/*



float CoeffA(float L, float kp){
    return std::pow(10, 0.506*kp-9.325)*std::pow(L,10);
}

float CoeffB(float L, float kp){
    return 8*std::pow(10, 0.506*kp-9.325)*std::pow(L,9);
}
*/

#endif
