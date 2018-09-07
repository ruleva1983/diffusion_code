#ifndef SOLVER_H_
#define SOLVER_H_


#include <vector>
#include "integrate.hpp"
#include "grid.hpp"


template <typename Integrator, typename Grid, typename Boundary, typename Observer>
void solve1D(Grid& grid, Integrator& scheme, const Boundary& bound, Observer& obs, float t,
             float tf, float dt){
    auto Xin = grid.get_state();
    grid.set_boundary(bound.low, bound.high, 0, Xin.size() - 1);
    obs.push_back(grid.get_state());
    while (tf > t){
            Xin = scheme.make_step(grid.get_state(), t, t + dt);
            grid.set_state(Xin);
            grid.set_boundary(bound.low, bound.high, 0, Xin.size() - 1);
            t += dt;
            obs.push_back(grid.get_state());
        }
}


#endif
