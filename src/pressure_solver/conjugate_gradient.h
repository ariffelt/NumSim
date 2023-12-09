#pragma once

#include "pressure_solver/pressure_solver_parallel.h"

class ConjugateGradient : public PressureSolverParallel
{
public:
    //! constructor
    ConjugateGradient(std::shared_ptr<Discretization> discretization, //changed, war std::shared_ptr discretization,
                      double epsilon,
                      int maximumNumberOfIterations,
                      std::shared_ptr<Partitioning> partitioning); //changed, war std::shared_ptr partitioning);

    //! solve the system of the Poisson equation for pressure
    void solve() override;
protected:
    void sendAndBorrowValues_q();
};