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
    for (i = discretization_->pIbegin() + 1; i < discretization_->pIEnd() - 1; i++)
    {
        for (j = discretization_->pJbegin() + 1; j < discretization_->pJEnd() - 1; j++)
        {
            // Second order derivative of p in direction x
            double d2pdx2 = (discretization_->p(i - 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i + 1, j)) / (hx2);
            // Second order derivative of p in direction y
            double d2pdy2 = (discretization_->p(i, j - 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j + 1)) / (hy2);
            // compute residual
            res += (d2pdx2 + d2pdy2 - discretization_->rhs(i, j)) * (d2pdx2 + d2pdy2 - discretization_->rhs(i, j));
        }
    }
    return res/(discretization_ -> nCells[0]*nCells[1]);
}

void PressureSolver::setBoundaryValues(){
    for (i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
        {
            discretization_->p(i, discretization_->pJEnd() - 1) = discretization_->p(i, discretization_->pJEnd() - 2);_
            discretization_->p(i, discretization_->pJBegin()) = discretization_->p(i, discretization_->pJBegin + 1);
        }
        for (j = discretization_->pJBegin(); i < discretization_->pJEnd(); j++)
        {
            discretization_->p(discretization_->pIEnd() - 1, j) = discretization_->p(discretization_->pIEnd() - 2, j);
            discretization_->p(discretization_->pIBegin(), j) = discretization_->p(discretization_->pIBegin + 1, j);
        }
}