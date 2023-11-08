#pragma once

#include <memory>

#include "pressure_solver/pressure_solver.h"

class GaussSeidel : public PressureSolver
{
public:
    //! constructor TODO
    GaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

    //! solve the system of the Poisson equation for pressure
    void solve();
};