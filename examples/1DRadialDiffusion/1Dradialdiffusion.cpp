#include <iostream>
#include "integrate.hpp"
#include "solver.hpp"
#include "grid.hpp"
#include "boundary.hpp"
#include <vector>
#include <cmath>
#include <fstream>

# define M_PI  3.14159265358979323846



const int N = 50;
const float xL = 1, xR = 7;
const float dt = 0.02;
const float t0 = 0;
const float tF = 100.0;
const float upperboundary = 7.2;
const float lowerboundary = 0.0;
const float tau=1.0;

using FuncType = std::function<float(float, float)>;


float InitialState(float x){
    /*if (x < 2)
        return 1.0;
    else
        return 1.0;*/

    return std::exp(-(x-3)*(x-3)/(2*0.1));
}

float kp(float t){
    return 6.0;
}


int main(int argc, char** argv)
{

    FuncType a = [](float x, float t){return std::pow(x, 10)*std::pow(10,0.506*kp(t)-9.325);};
    FuncType b = [](float x, float t){return 8.0*std::pow(x, 9)*std::pow(10,0.506*kp(t)-9.325);};
    FuncType c = [tau](float x, float t){return -1/tau;};
    FuncType d = [](float x, float t){return 0.0;};
    FuncType diff = [](float x, float t){return 1.0;};
    Functional<FuncType> A(diff);
    Functional<FuncType> B(d);
    Functional<FuncType> C(c);
    Functional<FuncType> D(d);
    
    std::cout << A(1,1) << " " << A(2, 1);

    float dx = (xR - xL)/ static_cast<float> (N);
    std::vector<float> Xin (N);

    std::vector<std::vector<float>> Obs(0);

    for (int i=0 ; i < N ; ++i)
        Xin[i] = InitialState(xL + dx*i);
    Grid1D<> grid(N, xL, xR, Xin);

    Boundary1D<> boundary(lowerboundary, upperboundary);
    CrankNicholson scheme(grid.get_deltax(), A, B, C, D, 0.5);
    solve1D<CrankNicholson, Grid1D<>, Boundary1D<>,
            std::vector<std::vector<float>>>(grid, scheme, boundary, Obs, t0, tF, dt);

    std::ofstream file("1DRadialDiffusion.dat");
    for (int i=0 ; i < Obs.size() ; ++i){
        for (int j=0 ; j < Obs[i].size() - 1; j++)
            file << Obs[i][j] << ",";
        file << Obs[i][Obs[i].size() - 1]<< std::endl;
    }
    file.close();

    return 0;
}

