#pragma once

#include "pressure_solver/pressure_solver_parallel.h"
#include "partitioning/partitioning.h"

class RedBlackSOR : public PressureSolverParallel
{
public:
    //! constructor
    RedBlackSOR(std::shared_ptr<Discretization> discretization,
                double epsilon, 
                int maximumNumberOfIterations, 
                double omega);

    //! solve the system of the Poisson equation for pressure
    void solve(); // TODO: override necessary?

private:
    double omega_;
};