#pragma once

#include <memory>

#include "pressure_solver/pressure_solver.h"

class SOR : public PressureSolver
{
public:
    //! constructor
    SOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega);

    //! solve the system of the Poisson equation for pressure
    void solve();

private:
    double omega_;
};