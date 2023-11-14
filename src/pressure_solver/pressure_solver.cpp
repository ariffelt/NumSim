#include "pressure_solver/pressure_solver.h"
#include <cassert>

// Constructor
PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization,
                               double epsilon,
                               int maximumNumberOfIterations) : discretization_(discretization), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations)
{
    assert(epsilon > 0);
    assert(maximumNumberOfIterations > 0);
}

double PressureSolver::getResidual()
{
    double res = 0;
    const double hx2 = discretization_->dx() * discretization_->dx();
    const double hy2 = discretization_->dy() * discretization_->dy();
    for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++)
    {
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++)
        {
            // Second order derivative of p in direction x
            double d2pdx2 = (discretization_->p(i - 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i + 1, j)) / (hx2);
            // Second order derivative of p in direction y
            double d2pdy2 = (discretization_->p(i, j - 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j + 1)) / (hy2);
            // compute residual
            res += (d2pdx2 + d2pdy2 - discretization_->rhs(i, j)) * (d2pdx2 + d2pdy2 - discretization_->rhs(i, j));
        }
    }

    return res / (((discretization_->pIEnd() + 1) - discretization_->pIBegin()) * ((discretization_->pJEnd() + 1) * discretization_->pJBegin()));
}

void PressureSolver::setBoundaryValues()
{
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
    {
        discretization_->p(i, discretization_->pJEnd() - 1) = discretization_->p(i, discretization_->pJEnd() - 2);
        discretization_->p(i, discretization_->pJBegin()) = discretization_->p(i, discretization_->pJBegin() + 1);
    }
    for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
    {
        discretization_->p(discretization_->pIEnd() - 1, j) = discretization_->p(discretization_->pIEnd() - 2, j);
        discretization_->p(discretization_->pIBegin(), j) = discretization_->p(discretization_->pIBegin() + 1, j);
    }
}