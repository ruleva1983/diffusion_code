#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#include "Eigen/Dense"
#include "diffusion_coeffs.hpp"
#include "grid.hpp"

using namespace Eigen;

class Scheme1D{
public:
    Scheme1D(){}
    
    Eigen::VectorXf make_step(const Eigen::VectorXf& X, float dt) {
        evaluate_A(X, dt);
        evaluate_b(X, dt);
        return A.colPivHouseholderQr().solve(b);
    }
    
    std::vector<float> make_step(std::vector<float>& Xvec, float dt) {
        Eigen::VectorXf X= Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(Xvec.data(), Xvec.size());
        evaluate_A(X, dt);
        evaluate_b(X, dt);
        X = A.colPivHouseholderQr().solve(b);
        std::vector<float> Xf(X.data(), X.data() + X.rows() * X.cols());
        return Xf;
    }
    
    virtual void evaluate_A(const Eigen::VectorXf&, float) = 0;
    virtual void evaluate_b(const Eigen::VectorXf&, float) = 0;
    
    
protected:
    Eigen::MatrixXf A;  // Each scheme is characterized by a matrix A to be inverted
    Eigen::VectorXf b;
};


class ExplicitScheme1D: public Scheme1D{

public:
    ExplicitScheme1D(float dx, float d) : Scheme1D(), deltax(dx), D(d)
    {
    }
    
    void evaluate_A(const Eigen::VectorXf& X, float dt){
        A = Eigen::MatrixXf::Zero (X.size(), X.size());
        for (int i = 0 ; i < X.size() ; ++i )
            A(i,i) = 1;
    }
    
    void evaluate_b(const Eigen::VectorXf& X, float dt){
        float coeff = dt*D/(dt*dt);
        b.resize(X.size());
        for (int i = 1 ; i < X.size() - 1 ; ++i )
            b(i) = coeff*(X(i-1) + X(i+1))+(1-2*coeff)*X(i);
    }
    
private:
    float deltax;
    float D;
};


class ImplicitScheme1D: public Scheme1D{

public:
    ImplicitScheme1D(float dx, float d) : Scheme1D(), deltax(dx), D(d)
    {
    }
    
    void evaluate_A(const Eigen::VectorXf& X, float dt){
        float coeff = dt*D/(dt*dt);
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
        b.resize(X.size());
        for (int i = 0 ; i < X.size() ; ++i )
            b(i) = X(i);
    }
    
    
private:
    float deltax;
    float D;
};


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

/*
class ImplicitRadialDiffusion: public Scheme1D{

public:
    ImplicitRadialDiffusion(float dx, float T) : Scheme1D(), deltax(dx), tau(T)
    {
    }
    
    
    std::vector<float> make_step(std::vector<float>& Xvec, std::vector<float>& Lvec, float kp, float dt) {
        Eigen::VectorXf X= Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(Xvec.data(), Xvec.size());
        Eigen::VectorXf L= Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(Lvec.data(), Lvec.size());
        evaluate_A(X, L, kp, dt);
        evaluate_b(X, dt);
        X = A.colPivHouseholderQr().solve(b);
        std::vector<float> Xf(X.data(), X.data() + X.rows() * X.cols());
        return Xf;
    }
    
    void evaluate_A(const Eigen::VectorXf& X, const Eigen::VectorXf& L, float kp, float dt){
        A = Eigen::MatrixXf::Zero (X.size(), X.size());
        A(0,0) = 1.0/dt + 1.0/tau + 2*CoeffA(L[0], kp)/std::pow(deltax, 2);
        A(0,1) = -CoeffA(L[0], kp)/std::pow(deltax, 2) - CoeffB(L[0], kp)/(2*deltax);
        for (int i = 1 ; i < X.size()-1 ; ++i ){
            A(i,i) = 1.0/dt + 1.0/tau + 2*CoeffA(L[i], kp)/std::pow(deltax, 2);
            A(i,i-1) = -CoeffA(L[i], kp)/std::pow(deltax, 2) + CoeffB(L[i], kp)/(2*deltax);
            A(i,i+1) = -CoeffA(L[i], kp)/std::pow(deltax, 2) - CoeffB(L[i], kp)/(2*deltax);
        }
        A(X.size()-1, X.size()-1) = 1.0/dt + 1.0/tau + 2*CoeffA(L[X.size()-1], kp)/std::pow(deltax, 2);
        A(X.size()-1, X.size()-2) = -CoeffA(L[X.size()-1], kp)/std::pow(deltax, 2) + CoeffB(L[X.size()-1], kp)/(2*deltax);
    }
    
    void evaluate_b(const Eigen::VectorXf& X, float dt){
        b = Eigen::VectorXf (X.size());
        for (int i = 1 ; i < X.size() - 1 ; ++i)
            b(i) = X(i)/dt;
    }
    
private:
    float deltax;
    float tau;        //Decay rate Model hyperparameter
};


/*
class CrankNicholsonRadial: public Scheme1D{

public:
    CrankNicholsonRadial(float dx, float T) : Scheme1D(), deltax(dx), tau(T)
    {
    }
    
    
    
    std::vector<float> make_step(std::vector<float>& Xvec, float dt) {
        Eigen::VectorXf X= Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(Xvec.data(), Xvec.size());
        evaluate_A(X, dt);
        evaluate_b(X, dt);
        X = A.colPivHouseholderQr().solve(b);
        std::vector<float> Xf(X.data(), X.data() + X.rows() * X.cols());
        return Xf;
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
    void evaluate_b(const Eigen::VectorXf& X, float dt, float L, float kp){
        b = Eigen::VectorXf (grid_dim);
        for (int i = 1 ; i < grid_dim - 1 ; ++i ){
            b(i) = X(i)*(1.0/dt - CoeffA(L, kp)/std::pow(deltax, 2) - 1.0/(2.0*tau));
            b(i) += X(i+1)*(CoeffA(L, kp)/(2*std::pow(deltax, 2))+CoeffB(L, kp)/(2*deltax));
            b(i) += X(i-1)*(CoeffA(L, kp)/(2*std::pow(deltax, 2))-CoeffB(L, kp)/(2*deltax));
        }
    }
    
private:
    float deltax;
    float tau;        //Decay rate Model hyperparameter
};
*/






#endif
