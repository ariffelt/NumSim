#include "pressure_solver/pressure_solver.h"

//Constructor
PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations):
    discretization_(discretization),
    epsilon_(epsilon),
    maximumNumberOfIterations_(maximumNumberOfIterations)
{
    assert(epsilon > 0);
    assert(maximumNumberOfIterations > 0);
}

double PressureSolver::getResidual(){
    double res = 0;
    for (i = discretisation_->pIbegin() + 1; i < discretisation_->pIEnd() - 1; i++)
    {
        for (j = discretisation_->pJbegin() + 1; j < discretisation_->pIEnd() - 1; j++)
        {
            // Second order derivative of p in direction x
            double d2pdx2 = (discretisation_->p(i - 1, j) - 2 * discretisation_->p(i, j) + discretisation_->p(i + 1, j)) / (hx2);
            // Second order derivative of p in direction y
            double d2pdy2 = (discretisation_->p(i, j - 1) - 2 * discretisation_->p(i, j) + discretisation_->p(i, j + 1)) / (hy2);
            // compute residual
            res02 += (d2pdx2 + d2pdy2 - discretization_->rhs(i, j)) * (d2pdx2 + d2pdy2 - discretization_->rhs(i, j));
        }
    }
    return res;
}

void PressureSolver::setBoundaryValues(){
    for (i = discretisation_->pIBegin(); i < discretisation_->pIEnd(); i++)
        {
            discretisation_->p(i, discretisation_->pJEnd() - 1) = discretisation_->p(i, discretisation_->pJEnd() - 2);
            discretisation_->p(i, discretisation_->pJBegin()) = discretisation_->p(i, discretisation_->pJBegin + 1);
        }
        for (j = discretisation_->pJBegin(); i < discretisation_->pJEnd(); j++)
        {
            discretisation_->p(discretisation_->pIEnd() - 1, j) = discretisation_->p(discretisation_->pIEnd() - 2, j);
            discretisation_->p(discretisation_->pIBegin(), j) = discretisation_->p(discretisation_->pIBegin + 1, j);
        }
}