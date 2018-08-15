#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#include "Eigen/Dense"

using namespace Eigen;

class Scheme1D{
public:
    Scheme1D(int dim = 0): grid_dim(dim)
    {
    }
    
    Eigen::VectorXf make_step(const Eigen::VectorXf& X, float dt) {
        evaluate_A(X, dt);
        evaluate_b(X, dt);
        return A.colPivHouseholderQr().solve(b);
    }
    
    virtual void evaluate_A(const Eigen::VectorXf&, float) = 0;
    virtual void evaluate_b(const Eigen::VectorXf&, float) = 0;
    
    
protected:
    int grid_dim;       // The dimension of the spatial grid
    Eigen::MatrixXf A;  // Each scheme is characterized by a matrix A to be inverted
    Eigen::VectorXf b;
    
};


class ExplicitScheme1D: public Scheme1D{

public:
    ExplicitScheme1D(int dim, float dx, float d) : Scheme1D(dim), deltax(dx), D(d)
    {
    }
    
    void evaluate_A(const Eigen::VectorXf& X, float dt){
        A = Eigen::MatrixXf::Zero (grid_dim, grid_dim);
        for (int i = 0 ; i < grid_dim ; ++i )
            A(i,i) = 1;
    }
    
    //TODO Works with zero boundary conditions at the edges
    void evaluate_b(const Eigen::VectorXf& X, float dt){
        float coeff = dt*D/(dt*dt);
        b.resize(grid_dim);
        for (int i = 1 ; i < grid_dim - 1 ; ++i )
            b(i) = coeff*(X(i-1) + X(i+1))+(1-2*coeff)*X(i);
    }
    
private:
    float deltax;
    float D;
};


class ImplicitScheme1D: public Scheme1D{

public:
    ImplicitScheme1D(int dim, float dx, float d) : Scheme1D(dim), deltax(dx), D(d)
    {
    }
    
    void evaluate_A(const Eigen::VectorXf& X, float dt){
        float coeff = dt*D/(dt*dt);
        A = Eigen::MatrixXf::Zero (grid_dim, grid_dim);
        A(0,0) = 1 + 2*coeff;
        A(0,1) = -coeff;
        for (int i = 1 ; i < grid_dim-1 ; ++i )
        {
            A(i,i) = 1 + 2*coeff;
            A(i,i-1) = -coeff;
            A(i,i+1) = -coeff;
        }
        A(grid_dim-1, grid_dim-1) = 1+2*coeff;
        A(grid_dim-1, grid_dim-2) = -coeff;
    }
    
    //TODO Works with zero boundary conditions at the edges
    void evaluate_b(const Eigen::VectorXf& X, float dt){
        b.resize(grid_dim);
        for (int i = 0 ; i < grid_dim ; ++i )
            b(i) = X(i);
    }
    
    
private:
    float deltax;
    float D;
};


class CrankNicholson: public Scheme1D{

public:
    CrankNicholson(int dim, float dx, float d) : Scheme1D(dim), deltax(dx), D(d)
    {
    }
    
    void evaluate_A(const Eigen::VectorXf& X, float dt){
        float coeff = dt*D/(2*dt*dt);
        A = Eigen::MatrixXf::Zero (grid_dim, grid_dim);
        A(0,0) = 1 + 2*coeff;
        A(0,1) = -coeff;
        for (int i = 1 ; i < grid_dim-1 ; ++i )
        {
            A(i,i) = 1 + 2*coeff;
            A(i,i-1) = -coeff;
            A(i,i+1) = -coeff;
        }
        A(grid_dim-1, grid_dim-1) = 1+2*coeff;
        A(grid_dim-1, grid_dim-2) = -coeff;
    }
    
    //TODO Works with zero boundary conditions at the edges
    void evaluate_b(const Eigen::VectorXf& X, float dt){
        
        float coeff = dt*D/(2*dt*dt);
        b.resize(grid_dim);
        for (int i = 1 ; i < grid_dim - 1 ; ++i )
            b(i) = coeff*(X(i-1) + X(i+1))+(1-2*coeff)*X(i);

    }
    
    
private:
    float deltax;
    float D;
};








#endif
