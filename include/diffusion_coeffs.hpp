#ifndef DIFFUSION_COEFFS_H_
#define DIFFUSION_COEFFS_H_

#include <cmath>

float CoeffA(float L, float kp){
    return std::pow(10, 0.506*kp-9.325)*std::pow(L,10);
}

float CoeffB(float L, float kp){
    return 8*std::pow(10, 0.506*kp-9.325)*std::pow(L,9);
}


#endif
