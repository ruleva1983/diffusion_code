#include <gtest/gtest.h>
#include <iostream>
#include <random>

#include <grid.hpp>

TEST(ToolsTest, GridClass_Default_Constructor)
{
    Grid1D<std::vector<float>> grid;
    ASSERT_EQ(grid.get_deltax(), 0.0);
    ASSERT_EQ(grid.get_state().size(), 0);
}

TEST(ToolsTest, GridClass_Constructor)
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
    grid.set_value(0,1000.0);
    ASSERT_EQ(grid.get_value(0), values[0]);


}