#include <gtest/gtest.h>
#include <iostream>
#include <random>

#include <grid.hpp>

TEST(GridClass, Default_Constructor)
{
    Grid1D<std::vector<float>> grid;
    ASSERT_EQ(grid.get_deltax(), 0.0);
    ASSERT_EQ(grid.get_state().size(), 0);
}

TEST(GridClass, Constructor)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(5, 100);
    std::uniform_real_distribution<> dis_real(0.0, 100.0);

    float xL = dis_real(gen);
    float xR = xL + dis_real(gen);
    int N = dis(gen);
    std::vector<float> values(0);
    for (int i = 0 ; i < N; ++i)
        values.push_back(dis_real(gen));
    Grid1D<std::vector<float>> grid(N, xL, xR, values);
    ASSERT_EQ(grid.get_deltax(), (xR-xL)/static_cast<float>(N));
}

TEST(GridClass, GetSetValues){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(5, 100);
    std::uniform_real_distribution<> dis_real(0.0, 100.0);

    float xL = dis_real(gen);
    float xR = xL + dis_real(gen);
    int N = dis(gen);
    std::vector<float> values(0);
    for (int i = 0 ; i < N; ++i)
        values.push_back(dis_real(gen));
    Grid1D<std::vector<float>> grid(N, xL, xR, values);

    std::uniform_int_distribution<> dis2(0, N-1);
    int i = dis2(gen);
    float value = dis_real(gen);
    grid.set_value(i, value);
    ASSERT_EQ(value, grid.get_value(i));
}


TEST(GridClass, GetPoints){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(5, 100);
    std::uniform_real_distribution<> dis_real(0.0, 100.0);

    float xL = dis_real(gen);
    float xR = xL + dis_real(gen);
    int N = dis(gen);
    std::vector<float> values(0);
    for (int i = 0 ; i < N; ++i)
        values.push_back(dis_real(gen));
    Grid1D<std::vector<float>> grid(N, xL, xR, values);

    std::uniform_int_distribution<> dis2(0, N-1);
    int i = dis2(gen);
    float dx = (xR-xL)/static_cast<float>(N);
    ASSERT_EQ(grid.get_point(i), xL + i*dx);
}

TEST(GridClass, GetState){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(5, 100);
    std::uniform_real_distribution<> dis_real(0.0, 100.0);

    float xL = dis_real(gen);
    float xR = xL + dis_real(gen);
    int N = dis(gen);
    std::vector<float> values(0);
    for (int i = 0 ; i < N; ++i)
        values.push_back(dis_real(gen));
    Grid1D<std::vector<float>> grid(N, xL, xR, values);

    ASSERT_EQ(grid.get_state(), values);
}