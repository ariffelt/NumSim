#include "pressure_solver/sor.h"

#include <cassert>
#include <cmath>

SOR::SOR(std::shared_ptr<Discretization> discretization,
         double epsilon,
         int maximumNumberOfIterations,
         double omega) : PressureSolver(discretization, epsilon, maximumNumberOfIterations),
                         omega_(omega)
{
}

void SOR::solve()
{
    // compute prefactor
    const double hx2 = discretization_->dx() * discretization_->dx();
    const double hy2 = discretization_->dy() * discretization_->dy();
    const double prefactor = hx2 * hy2 / (2 * (hx2 + hy2));

    // compute initial squared residual
    double res02 = getResidual();

    // set initial values
    int iteration = 0;
    const double epsilon2 = epsilon_ * epsilon_;
    double res2 = res02;

    // solve the system of the Poisson equation for pressure with the SOR method
    while (res2 >= epsilon2 * res02 && iteration < maximumNumberOfIterations_)
    {
        // Implement SOR for inner points
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++)
        {
            for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++)
            {
                double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / (hx2);
                double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / (hy2);
                discretization_->p(i, j) = (1 - omega_) * discretization_->p(i, j) + omega_ * prefactor * (px + py - discretization_->rhs(i, j));
            }
        }
        // Update boundary values such that Neumann Boundary conditions are fulfilled
        setBoundaryValues();

        // increase iteration
        iteration++;

        // update residual
        res2 = getResidual();
    }
}