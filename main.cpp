#include <iostream>

#include "integrate.hpp"
#include "solver.hpp"
#include "grid.hpp"
#include "boundary.hpp"
#include <vector>
#include <boost/date_time/gregorian/gregorian.hpp>
#include "io.hpp"

using namespace boost::gregorian;
using FuncType = std::function<float(float, float)>;




int main()
{

    FuncType a = [](float x, float t){return 1.0;};
    FuncType b = [](float x, float t){return 0.0;};
    FuncType c = [](float x, float t){return 0.0;};
    FuncType d = [](float x, float t){return 0.0;};
    Functional<FuncType> A(a);
    Functional<FuncType> B(b);
    Functional<FuncType> C(c);
    Functional<FuncType> D(d);
/*
    //date d(2002,1,10);
    const int N = 50;
    std::vector<float> Xin (N), Obs(0);
    for (int i=0 ; i < N ; ++i)
        Xin[i] = 1.0;
    Grid1D<> grid(N, 0.9, 6.1, Xin);
    Boundary1D<> boundary(0.0, 0.0);
    //ExplicitScheme1D scheme(grid.get_deltax(), A, B, C, D);
    //solve1D(grid, scheme, boundary,Obs , 0.0, 10.0, 0.1);
*/
	return 0;
}
