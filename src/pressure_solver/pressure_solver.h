#pragma once

#include <memory>

#include "discretization/1_discretization.h"

class PressureSolver
{
public:
    //! constructor
    PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);
 
    //! solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid
    virtual void solve()=0;
 
protected:
    //! set the boundary values to account for homogenous Neumann boundary conditions, this has to be called after every iteration
    void setBoundaryValues();

    //! compute the residual 
    double getResidual();

    //! object holding the needed field variables for rhs and p
    std::shared_ptr< Discretization > 	discretization_;
 
    double epsilon_;
 
    int maximumNumberOfIterations_;
};