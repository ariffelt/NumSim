#include "pressure_solver/gauss_seidel.h"

#include <cassert>

GaussSeidel::GaussSeidel(std::shared_ptr<Discretization> discretization,
                           double epsilon,
                           int maximumNumberOfIterations) : PressureSolver(discretization, epsilon, maximumNumberOfIterations)
{
    
}

void GaussSeidel::solve()
{
    const double hx2 = discretization_->meshWidth_[0];
    const double hy2 = discretization_->meshWidth_[1];
    const double prefactor = hx2 * hy2 / (2 * (hx2 + hy2));

    // compute initial residual
    double res02 = getResidual();

    // solve the system of the Poisson equation for pressure with the Gauss-Seidel method
    int iteration = 0;
    const double epsilon2 = epsilon_ * epsilon_;
    double res2 = res02;

    while (res2 >= epsilon2 * res02 && iteration < maximumNumberOfIterations_)
    { // Implement GS for inner points
        for (i = discretisation_->pIBegin() + 1; i < discretisation_->pIEnd() - 1; i++)
        {
            for (j = discretisation_->pJBegin() + 1; j < discretisation_->pJEnd() - 1; j++)
            {
                    discretisation_ -> p(i,j) = prefactor * ((discretisation_ -> p(i-1,j) + discretisation_ -> p(p+1,j))/(hx2) + (discretisation_ -> p(i,j-1) + discretisation_ -> p(i,j+1))/(hy2)) - discretisation_ -> rhs(i,j));
            }
        }
        // Update boundary values such that Neumann Boundary conditions are fulfilled
        setBoundaryValues();

        // increase iteration
        iteration++;

        // update residual
        res2 = getResidual();
    }
