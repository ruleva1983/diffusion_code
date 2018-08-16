#include <iostream>

#include "integrate.hpp"

using namespace Eigen;

int main()
{
    Eigen::VectorXf Xi (50), Xf;
    for (int i=0 ; i < Xi.size()  ; ++i)
        Xi(i) = 2;
    
    CrankNicholson scheme(50, 5, 1);
    Xf = scheme.make_step(Xi, 0.1);
   
    for (int j = 0; j< 500 ; ++j)
        Xf = scheme.make_step(Xf, 0.1);
    for (int i=0 ; i < Xi.size()  ; ++i)
        std::cout << Xf(i) << std::endl;
    
	return 0;
}
