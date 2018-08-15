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
    Xf = scheme.make_step(Xf, 0.1);
    for (int i=0 ; i < Xi.size()  ; ++i)
        std::cout << Xf(i) << std::endl;
    
    std::cout << "Euler implicit: " ;
    Eigen::VectorXf Xi2 (50);
    for (int i=0 ; i < Xi2.size()  ; ++i)
        Xi2(i) = 2;
    ImplicitScheme1D scheme2(50, 5, 1);
    Xf = scheme2.make_step(Xi2, 0.1);
    Xf = scheme2.make_step(Xf, 0.1);
    for (int i=0 ; i < Xi2.size()  ; ++i)
        std::cout << Xf(i) << std::endl;
    
	return 0;
}
