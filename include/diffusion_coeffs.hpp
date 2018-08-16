#ifndef DIFFUSION_COEFFS_H_
#define DIFFUSION_COEFFS_H_

#include <cmath>

float CoeffA(float L, float kp){
    return std::pow(10, 0.506*kp-9.325)*10*std::pow(L,10);
}

float CoeffB(float L, float kp){
    return std::pow(10, 0.506*kp-9.325)*(10*std::pow(L,9) - 2*std::pow(L,7));
}


#endif
