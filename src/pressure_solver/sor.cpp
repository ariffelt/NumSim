#include "pressure_solver/sor.h"

#include <cassert>
#include <cmath>
#include <iostream>

/**
 * Constructor
 * @param discretization discretization object
 * @param epsilon tolerance for the residual
 * @param maximumNumberOfIterations maximum number of iterations
 * @param omega relaxation parameter
 */
SOR::SOR(std::shared_ptr<Discretization> discretization,
         double epsilon,
         int maximumNumberOfIterations,
         double omega) : PressureSolver(discretization, epsilon, maximumNumberOfIterations),
                         omega_(omega)
{
    assert(omega >= 1);
    assert(omega < 2);
    if (omega != 2 / (1 + sin(M_PI * discretization_->dx())) && omega != 2 / (1 + sin(M_PI * discretization_->dy())))
    {
        std::cout << "WARNING: omega is not the optimal value for SOR." << std::endl;
    }
}

/**
 * Solve the Poisson equation for pressure with the SOR method
 */
void SOR::solve()
{
    // compute prefactor
    const double hx2 = discretization_->dx() * discretization_->dx();
    const double hy2 = discretization_->dy() * discretization_->dy();
    const double prefactor = hx2 * hy2 / (2 * (hx2 + hy2));

    // compute initial squared residual
    double res2 = getResidual();

    // set initial values
    int iteration = 0;

    // tolerance for the residual, taken to the power of 2 to compare to squared residual
    const double epsilon2 = epsilon_ * epsilon_;

    // solve the system of the Poisson equation for pressure with the SOR method
    while (res2 >= epsilon2 && iteration < maximumNumberOfIterations_)
    {
        // implement SOR for the inner points

        // std::cout << "SOR iteration: " << iteration << " residual: " << res2 << std::endl;

        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++)
        {
            // std::cout << "i: " << i << std::endl;
            for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++)
            {
                // std::cout << "j: " << j << std::endl;
                if (discretization_->isInnerFluidCell(i, j))
                {
                    // std::cout << "i: " << i << " j: " << j << std::endl;
                    double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / (hx2);
                    double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / (hy2);

                    // std::cout << "px: " << px << " py: " << py << std::endl;

                    discretization_->p(i, j) = (1 - omega_) * discretization_->p(i, j) + omega_ * prefactor * (px + py - discretization_->rhs(i, j));

                    // std::cout << "p(i, j): " << discretization_->p(i, j) << std::endl;
                }
            }
        }

        // std::cout << "SOR iteration: " << iteration << " residual: " << res2 << std::endl;

        // update boundary values such that Neumann Boundary conditions are fulfilled
        setBoundaryValues();

        // std::cout << "SOR iteration: " << iteration << " residual: " << res2 << std::endl;

        // increase iteration
        iteration++;

        // update residual
        res2 = getResidual();
    }
}