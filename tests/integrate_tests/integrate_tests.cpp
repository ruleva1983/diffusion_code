#include <gtest/gtest.h>
#include <iostream>
#include <random>

#include <integrate.hpp>
#include <diffusion_coeffs.hpp>

#include <functional>
     
TEST(ExplicitScheme1D, SimpleDiffusion)
{


std::function<float, (float float)> a = [](float x, float t){return 1.0;};
    
Functional<> A(a);
//Functional B([](double x, double t) { return 0.0; });
//Functional C([](double x, double t) { return 0.0; });
//Functional D([](double x, double t) { return 0.0; });
//ExplicitScheme1D scheme();
//ASSERT_EQ(grid.get_deltax(), 0.0);
//ASSERT_EQ(grid.get_state().size(), 0);

}
