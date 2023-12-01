#pragma once

#include <mpi.h>

#include "pressure_solver/pressure_solver.h"
#include "discretization/1_discretization.h"
#include "partitioning/partitioning.h"

class PressureSolverParallel : public PressureSolver
{
public:
    //! constructor
    PressureSolverParallel(std::shared_ptr<Discretization> discretization,
                           double epsilon,
                           int maximumNumberOfIterations,
                           std::shared_ptr<Partitioning> partitioning);
    
    virtual void solve() = 0;

protected:
    
    void sendAndBorrowValues(); 

    virtual double getResidual();

    std::shared_ptr<Partitioning> partitioning_;
};