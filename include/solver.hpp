#ifndef SOLVER_H_
#define SOLVER_H_


#include <vector>
#include "integrate.hpp"
#include "grid.hpp"


template <typename Integrator, typename Grid, typename Boundary, typename Observer=std::vector<float>>
void solve1D(Grid& grid, Integrator& scheme, const Boundary& bound, Observer& obs, float t0, float tf, float dt){
    auto Xin = grid.get_state();
    while (tf > t0){
            Xin = scheme.make_step(Xin, dt);
            grid.set_boundary(bound.low, bound.high, 0, Xin.size() - 1);
            t0 += dt;
            std::cout << t0 << std::endl;
            //Add observer
        }
}


#endif
