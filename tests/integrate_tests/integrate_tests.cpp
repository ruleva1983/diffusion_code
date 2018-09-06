#include <gtest/gtest.h>
#include <iostream>
#include <random>

#include <integrate.hpp>
#include <diffusion_coeffs.hpp>
#include <functional>
#include <Eigen/Dense>


using namespace Eigen;
using FuncType = std::function<float(float, float)>;


TEST(ExplicitScheme1D, SimpleDiffusion)
{
FuncType a = [](float x, float t){return 1.0;};
FuncType b = [](float x, float t){return 0.0;};
FuncType c = [](float x, float t){return 0.0;};
FuncType d = [](float x, float t){return 0.0;};
Functional<FuncType> A(a);
Functional<FuncType> B(b);
Functional<FuncType> C(c);
Functional<FuncType> D(d);

float deltax = 1.0, deltat = 1.0;
ExplicitScheme1D scheme(deltax, A, B, C, D);

VectorXf X(5);
X << 0.0,0.0,0.0,0.0,0.0;

scheme.evaluate_A(X, 0, deltat);
scheme.evaluate_b(X, 0, deltat);

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> dis(0, X.size()-2);
int i = dis(gen), j = dis(gen);
if (i == j)
    ASSERT_EQ(scheme.getMatA(i,j), 1.0);
else
    ASSERT_EQ(scheme.getMatA(i,j), 0.0);

ASSERT_EQ(scheme.getVecb(i), 0.0);

}
