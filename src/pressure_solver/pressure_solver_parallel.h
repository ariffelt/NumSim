#pragma once

#include "pressure_solver/pressure_solver.h"

class PressureSolverParallel : public PressureSolver
{
public:
    //! constructor
    PressureSolverParallel(std::shared_ptr<Discretization> discretization,
                           double epsilon,
                           int maximumNumberOfIterations);

    //! TODO: what functions?
};