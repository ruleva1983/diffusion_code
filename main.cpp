#include <iostream>

#include "integrate.hpp"
#include "solver.hpp"
#include "grid.hpp"
#include "boundary.hpp"
#include <vector>
#include <boost/date_time/gregorian/gregorian.hpp>
#include "io.hpp"

using namespace boost::gregorian;

int main()
{
    date d(2002,1,10);
    const int N = 50;
    std::vector<float> Xin (N), Obs(0);
    for (int i=0 ; i < N ; ++i)
        Xin[i] = 1.0;
    Grid1D<> grid(N, 1.0, 6.0, Xin);
    Boundary1D<> boundary(0.0, 0.0);
    ImplicitScheme1D scheme(grid.get_deltax(),1.0);
    solve1D(grid, scheme, boundary,Obs , 0.0, 10.0, 0.1);
    
	return 0;
}
