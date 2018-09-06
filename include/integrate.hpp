#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#include "Eigen/Dense"
#include "diffusion_coeffs.hpp"
#include "grid.hpp"
#include <functional>

using namespace Eigen;

using Coeff = std::function<float(float, float)>;

class Scheme1D{
public:
    Scheme1D(){}
    
    Eigen::VectorXf make_step(Eigen::VectorXf& X, float t1, float t2) {
        evaluate_A(X, t1, t2-t1);
        evaluate_b(X, t1, t2-t1);
        return MatA.colPivHouseholderQr().solve(Vecb);
    }
    
    std::vector<float> make_step(std::vector<float>& Xvec, float t1, float t2) {
        Eigen::VectorXf X = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(Xvec.data(), Xvec.size());
        evaluate_A(X, t1, t2-t1);
        evaluate_b(X, t1, t2-t1);
        X = MatA.colPivHouseholderQr().solve(Vecb);
        std::vector<float> Xf(X.data(), X.data() + X.rows() * X.cols());
        return Xf;
    }
    
    virtual void evaluate_A(const Eigen::VectorXf&, float, float) = 0;
    virtual void evaluate_b(const Eigen::VectorXf&, float, float) = 0;
    
    float getMatA(int i, int j) const{
        return MatA(i,j);
    }

    float getVecb(int i) const{
        return Vecb(i);
    }

protected:
    Eigen::MatrixXf MatA;
    Eigen::VectorXf Vecb;
};


class ExplicitScheme1D: public Scheme1D{

    


public:
    ExplicitScheme1D(float deltax, Coeff a, Coeff b, Coeff c, Coeff d ) : Scheme1D(), dx(deltax), A(a), B(b), C(c), D(d)
    {
    }
    
    void evaluate_A(const Eigen::VectorXf& X, float t, float dt){
        MatA = Eigen::MatrixXf::Zero (X.size(), X.size());
        for (int i = 0 ; i < X.size() ; ++i )
            MatA(i,i) = 1;
    }
    
    void evaluate_b(const Eigen::VectorXf& X, float t, float dt){

        Vecb.resize(X.size());
        for (int i = 1 ; i < X.size() - 1 ; ++i ){
            float coeff_i = 1 + dt*C(X(i),t) - 2*dt/(dx*dx)*A(X(i),t);
            float coeff_ip1 = dt/(dx*dx)*A(X(i),t) + dt/(2*dx)*B(X(i),t);
            float coeff_im1 = dt/(dx*dx)*A(X(i),t) - dt/(2*dx)*B(X(i),t);
            float coeff = D(X(i),t);
            Vecb(i) = coeff + coeff_i*X(i) + coeff_im1*X(i-1) + coeff_ip1*X(i+1);
        }
    }
    
private:
    float dx;
    Coeff A, B, C, D;
};


class ImplicitScheme1D: public Scheme1D{

public:
    ImplicitScheme1D(float deltax, Coeff a, Coeff b, Coeff c, Coeff d ) : Scheme1D(), dx(deltax), A(a), B(b), C(c), D(d)
    {
    }
    
    void evaluate_A(const Eigen::VectorXf& X, float t, float dt){
        MatA = Eigen::MatrixXf::Zero (X.size(), X.size());
        MatA(0,0) = 1 + 2*dt/(dx*dx)*A(X(0),t+dt) - dt*C(X(0), t+dt);
        MatA(0,1) = -dt/(dx*dx)*A(X(0),t+dt) + dt*B(X(0),t+dt);
        for (int i = 1 ; i < X.size()-1 ; ++i )
        {
            MatA(i,i) = 1 + 2*dt/(dx*dx)*A(X(i),t+dt) - dt*C(X(i), t+dt);
            MatA(i,i-1) = -dt/(dx*dx)*A(X(i),t+dt) - dt*B(X(i),t+dt);
            MatA(i,i+1) = -dt/(dx*dx)*A(X(i),t+dt) + dt*B(X(i),t+dt);
        }
        MatA(X.size()-1, X.size()-1) = 1 + 2*dt/(dx*dx)*A(X(X.size()-1),t+dt) - dt*C(X(X.size()-1), t+dt);
        MatA(X.size()-1, X.size()-2) = -dt/(dx*dx)*A(X(0),t+dt) - dt*B(X(0),t+dt);
    }

    void evaluate_b(const Eigen::VectorXf& X, float t, float dt){
        Vecb.resize(X.size());
        for (int i = 0 ; i < X.size() ; ++i )
            Vecb(i) = X(i) + D(X(i), t+dt);
    }

private:
    float dx;
    Coeff A, B, C, D;
};

/*
class CrankNicholson: public Scheme1D{

public:
    CrankNicholson(float dx, float d) : Scheme1D(), deltax(dx), D(d)
    {
    }
    
    void evaluate_A(const Eigen::VectorXf& X, float dt){
        float coeff = dt*D/(2*dt*dt);
        A = Eigen::MatrixXf::Zero (X.size(), X.size());
        A(0,0) = 1 + 2*coeff;
        A(0,1) = -coeff;
        for (int i = 1 ; i < X.size()-1 ; ++i )
        {
            A(i,i) = 1 + 2*coeff;
            A(i,i-1) = -coeff;
            A(i,i+1) = -coeff;
        }
        A(X.size()-1, X.size()-1) = 1+2*coeff;
        A(X.size()-1, X.size()-2) = -coeff;
    }
    
    //TODO Works with zero boundary conditions at the edges
    void evaluate_b(const Eigen::VectorXf& X, float dt){
        
        float coeff = dt*D/(2*dt*dt);
        b.resize(X.size());
        for (int i = 1 ; i < X.size() - 1 ; ++i )
            b(i) = coeff*(X(i-1) + X(i+1))+(1-2*coeff)*X(i);

    }
    
private:
    float deltax;
    float D;
};
*/


#endif
