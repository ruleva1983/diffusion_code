#include <iostream>
#include "integrate.hpp"
#include "solver.hpp"
#include "grid.hpp"
#include "boundary.hpp"
#include <vector>
#include <cmath>
#include <fstream>

# define M_PI  3.14159265358979323846

using FuncType = std::function<float(float, float)>;


float Gaussian(float x, float mean=0.0, float sigma=1.0){
    return std::exp(-(x-mean)*(x-mean)/(2*sigma*sigma))/std::sqrt(2*M_PI*sigma*sigma);
}


int main()
{

    FuncType a = [](float x, float t){return 1.0;};
    FuncType b = [](float x, float t){return 0.0;};
    Functional<FuncType> A(a);
    Functional<FuncType> B(b);
    Functional<FuncType> C(b);
    Functional<FuncType> D(b);

    const int N = 200;
    const float xL = -20, xR = 20;
    float dx = (xR - xL)/ static_cast<float> (N);
    std::vector<float> Xin (N);

    std::vector<std::vector<float>> Obs(0);

    for (int i=0 ; i < N ; ++i)
        Xin[i] = Gaussian(xL + dx*i);
    Grid1D<> grid(N, xL, xR, Xin);

    Boundary1D<> boundary(0.0, 0.0);
    ExplicitScheme1D scheme(grid.get_deltax(), A, B, C, D);
    solve1D<ExplicitScheme1D, Grid1D<>, Boundary1D<>,
            std::vector<std::vector<float>>>(grid, scheme, boundary, Obs, 0.0, 10.0, 0.1);

    std::ofstream file("prova.dat");

    for (int i=0 ; i < Obs.size() ; ++i){
        for (int j=0 ; j <Obs[i].size(); j++)
            file << Obs[i][j] << ",";
        file << std::endl;
    }

    file.close();



    return 0;
}
