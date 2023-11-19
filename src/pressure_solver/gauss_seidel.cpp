#include "pressure_solver/gauss_seidel.h"

#include <iostream>
#include <cassert>

/**
 * Constructor
 * @param discretization discretization object
 * @param epsilon tolerance for the residual
 * @param maximumNumberOfIterations maximum number of iterations
 */
GaussSeidel::GaussSeidel(std::shared_ptr<Discretization> discretization,
                         double epsilon,
                         int maximumNumberOfIterations) : PressureSolver(discretization, epsilon, maximumNumberOfIterations)
{
}

/**
 * Solve the Poisson equation for pressure with the Gauss-Seidel method
 */
void GaussSeidel::solve()
{
    const double hx2 = discretization_->dx() * discretization_->dx();
    const double hy2 = discretization_->dy() * discretization_->dy();
    const double prefactor = hx2 * hy2 / (2.0 * (hx2 + hy2));

    // compute initial squared residual
    double res2 = getResidual();

    int iteration = 0;

    // tolerance for the residual, taken to the power of 2 to compare to squared residual
    const double epsilon2 = epsilon_ * epsilon_;

    // solve the system of the Poisson equation for pressure with the Gauss-Seidel method
    while (res2 >= epsilon2 && iteration < maximumNumberOfIterations_)
    { 
        // implement Gauss-Seidel for the inner points
        for (int i = discretization_->pIBegin() + 1; i <= discretization_->pIEnd() - 1; i++)
        {
            for (int j = discretization_->pJBegin() + 1; j <= discretization_->pJEnd() - 1; j++)
            {
                discretization_->p(i, j) = prefactor * ((discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / hx2 + (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / hy2 - discretization_->rhs(i, j));
            }
        }

        // update boundary values such that Neumann Boundary conditions are fulfilled
        setBoundaryValues();

        // increase iteration
        iteration++;

        // update residual
        res2 = getResidual();
    }
}