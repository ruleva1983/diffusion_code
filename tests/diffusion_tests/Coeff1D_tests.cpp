#include <gtest/gtest.h>
#include <iostream>
#include <random>

#include <diffusion_coeffs.hpp>

TEST(Coeff1D, Constant)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_real(0.0, 100.0);

    float a = dis_real(gen);
    Constant A(a);
    ASSERT_EQ(a,A());
}

float sum(float x, float t){
    return x+t;
}

float product(float x, float t){
    return x*t;
}



TEST(Coeff1D, Functional){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_real(0.0, 100.0);

    float a1 = dis_real(gen) , a2 = dis_real(gen);
    Functional<> A1(sum);
    ASSERT_EQ(sum(a1, a2), A1(a1, a2));

    Functional<> A2(product);
    ASSERT_EQ(product(a1, a2), A2(a1, a2));

}