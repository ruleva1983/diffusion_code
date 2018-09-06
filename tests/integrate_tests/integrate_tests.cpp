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

ExplicitScheme1D scheme(1.0, A, B, C, D);

VectorXf X = {0.0,1.0,2.0,1.0,0.0};
scheme.evaluate_A(X, 0.0, 1.0);

//ASSERT_EQ(grid.get_deltax(), 0.0);
//ASSERT_EQ(grid.get_state().size(), 0);

}
