#ifndef DIFFUSION_COEFFS_H_
#define DIFFUSION_COEFFS_H_

#include <cmath>

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

template <typename F=float(*)(float, float)>
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






template <class T = Constant, class F = Functional>
class DiffusionCoefficients{
public:
    DiffusionCoefficients(T CA=Constant(1.0), T CB=Constant(1.0), T CC=Constant(1.0),
                          T CD=Constant(1.0)): coeffA(CA), coeffB(CB), coeffC(CC), coeffD(CD)
    {
    }

    DiffusionCoefficients(F CA, F CB, F CC, F CD): coeffA(CA), coeffB(CB), coeffC(CC), coeffD(CD)
    {
    }

    float evalA(float x, float t){
        return coeffA(x, t);
    }

    float evalB(float x, float t){
        return coeffB(x, t);
    }

    float evalC(float x, float t){
        return coeffC(x, t);
    }

    float evalD(float x, float t){
        return coeffD(x, t);
    }

private:
    T coeffA, coeffB, coeffC, coeffD;
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
