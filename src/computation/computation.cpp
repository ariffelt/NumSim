#include "computation/computation.h"
#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/sor.h"
#include "discretization/2_donor_cell.h"
#include "discretization/2_central_differences.h"

#include <cassert>

#include <cmath>

// TODO: constructor
/* Computation::Computation()
{
    // nothing to do
} */

void Computation::initialize(int argc, char *argv[])
{
    //! load and print settings
    settings_ = Settings();
    settings_.loadFromFile(argv[1]);
    settings_.printSettings();

    //! compute the meshWidth from the physical size and the number of cells
    meshWidth_[0] = settings_.physicalSize[0] / settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1] / settings_.nCells[1];

    //! initialize discretization
    if (settings_.useDonorCell)
    {
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    }
    else
    {
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
    }

    //! initialize the solver
    if (settings_.pressureSolver == "GaussSeidel")
    {
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
    }
    else if (settings_.pressureSolver == "SOR")
    {
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    }
    else
    {
        std::cout << "Error: Unknown pressure solver " << settings_.pressureSolver << std::endl;
        exit(1);
    }

    //! initialize output writers
    outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
}
void Computation::testBC()
{
    applyBoundaryValues();
    outputWriterParaview_->writeFile(0);
    outputWriterText_->writeFile(0);
}
/**
 * Run the whole simulation until tend.
 */

void Computation::runSimulation()
{
    double t = 0.0;
    while (t < settings_.endTime)
    {   
        //std::cout << "starting t = " << t << std::endl;
        applyBoundaryValues();
        computeTimeStepWidth();
        // decrease time step width in last time step, s.t. the end time will be reached exactly
        if (t + dt_ > settings_.endTime)
        {
            dt_ = settings_.endTime - t;
        }
        t += dt_;
        computePreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        //std::cout << "dt = " << dt_ << std::endl;
        outputWriterParaview_->writeFile(t);
        outputWriterText_->writeFile(t);
    }
}
//! Compute the time step width dt from maximum velocities.
void Computation::computeTimeStepWidth()
{
    double dt_diffusion = (settings_.re / 2) * (pow(discretization_->dx(), 2) * pow(discretization_->dy(), 2)) / (pow(discretization_->dx(), 2) + pow(discretization_->dy(), 2));

    double maxU = 0;
    double maxV = 0;
    std::array<int, 2> uSize = discretization_->uSize();
    //! TODO: danger zone
    for (int i = 0; i < uSize[0]; i++)
    {
        for (int j = 0; j < uSize[1]; j++)
        {
            //! possible because the grid for u and v have the same dimensions in all directions
            maxU = std::max(maxU, std::fabs(discretization_->u(i, j)));
            maxV = std::max(maxV, std::fabs(discretization_->v(i, j)));
        }
    }
    double dt_convection = std::min(discretization_->dx() / maxU, discretization_->dy() / maxV);
    //! TODO: ensure somewhere that tau<1
    //! because of the scaling with the security factor tau < 1, a subtraction of a small value,
    // to ensure dt smaller and not smaller/equal than required, is not necessary
    assert(settings_.tau < 1);
    dt_ = settings_.tau * std::min(dt_diffusion, dt_convection);
}

//! Set velocity boundary values for u, v, F and G
void Computation::applyBoundaryValues()
{
    //! set the u boundary values for bottom and top first, as for corner cases, the left and right border should be used
    //changed
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        //! u boundary values bottom, assuming inhomogenous Dirichlet conditions
        discretization_->u(i, discretization_->uJBegin()) = 2.0 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin() + 1);
        //! u boundary values top, assuming inhomogenous Dirichlet conditions
        //! TODO: why do we have dirichlet BC at top? I thought we had flow at top? Or do we set the movement in the dirichletBcTop?
        discretization_->u(i, discretization_->uJEnd()) = 2.0 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() - 1);
    }
    //! set the u boundary values for left and right side now
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++)
    {
        //! u boundary values left, assuming inhomogenous Dirichlet conditions
        discretization_->u(discretization_->uIBegin(), j) = settings_.dirichletBcLeft[0];
        //! u boundary values right, assuming inhomogenous Dirichlet conditions
        //! changed (?)
        discretization_->u(discretization_->uIEnd() - 1, j) = settings_.dirichletBcRight[0];
    }

    //! set the v boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++)
    {
        //! v boundary values bottom, assuming inhomogenous Dirichlet conditions
        discretization_->v(i, discretization_->vJBegin()) = settings_.dirichletBcBottom[1];
        //! v boundary values top, assuming inhomogenous Dirichlet conditions
        //! TODO: find out whether we apply the BC on the very top, artificially added v, or one below
        //! changed (?)
        discretization_->v(i, discretization_->vJEnd() - 1) = settings_.dirichletBcTop[1];
    }

    //! set the v boundary values for left and right side now
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        //! v boundary values left, assuming inhomogenous Dirichlet conditions
        discretization_->v(discretization_->vIBegin(), j) = 2 * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin() + 1, j);
        //! v boundary values right, assuming inhomogenous Dirichlet conditions
        discretization_->v(discretization_->vIEnd(), j) = 2 * settings_.dirichletBcRight[1] - discretization_->v(discretization_->vIEnd() - 1, j);
    }

    //! set the F boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        //! F boundary values bottom
        discretization_->f(i, discretization_->uJBegin()) = discretization_->u(i, discretization_->uJBegin());
        //! F boundary values top
        //! TODO: Danger zone!
        discretization_->f(i, discretization_->uJEnd()) = discretization_->u(i, discretization_->uJEnd());
    }

    //! set the F boundary values for left and right side now
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++)
    {
        //! F boundary values left
        discretization_->f(discretization_->uIBegin(), j) = discretization_->u(discretization_->uIBegin(), j);
        //! F boundary values right
        //changed (?)
        discretization_->f(discretization_->uIEnd() - 1, j) = discretization_->u(discretization_->uIEnd() - 1, j);
    }

    //! set the G boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++)
    {
        //! G boundary values bottom
        discretization_->g(i, discretization_->vJBegin()) = discretization_->v(i, discretization_->vJBegin());
        //! G boundary values top
         discretization_->g(i, discretization_->vJEnd() - 1) = discretization_->v(i, discretization_->vJEnd() - 1);
    }

    //! set the G boundary values for left and right side now
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        //! G boundary values left
        discretization_->g(discretization_->vIBegin(), j) = discretization_->v(discretization_->vIBegin(), j);
        //! G boundary values right
        discretization_->g(discretization_->vIEnd(), j) = discretization_->v(discretization_->vIEnd(), j);
    }
}
    //! compute the preliminary velocities, F and G
    void Computation::computePreliminaryVelocities()
    {
        //! compute the preliminary velocities only for the inside of the area, as on the direct boundary, they are constantly 0

        // ! compute preliminary F
        //changed
        for (int i = discretization_->uIBegin() + 1; i < discretization_->uIEnd() - 1; i++)
        {
            //! TODO: evtl uJEnd()-1? Weil die oberen Werte zum wegwerfen sind? - changed
            for (int j = discretization_->uJBegin() + 1; j < discretization_->uJEnd(); j++)
            {
                double diffusionTerms = (1 / settings_.re) * (discretization_->computeD2uDx2(i, j) + discretization_->computeD2uDy2(i, j));
                double convectionTerms = discretization_->computeDu2Dx(i, j) + discretization_->computeDuvDy(i, j);
                //! TODO: check whether we update the g
                discretization_->f(i, j) = discretization_->u(i, j) + dt_ * (diffusionTerms - convectionTerms + settings_.g[0]);
            }
        }
        //! compute preliminary G
        //! TODO: evtl vIEnd()-1? Weil die rechten Werte zum wegwerfen sind? - changed
        //changed
        for (int i = discretization_->vIBegin() + 1; i < discretization_->vIEnd(); i++)
        {
            for (int j = discretization_->vJBegin() + 1; j < discretization_->vJEnd() - 1; j++)
            {
                double diffusionTerms = (1 / settings_.re) * (discretization_->computeD2vDx2(i, j) + discretization_->computeD2vDy2(i, j));
                double convectionTerms = discretization_->computeDuvDx(i, j) + discretization_->computeDv2Dy(i, j);

                discretization_->g(i, j) = discretization_->v(i, j) + dt_ * (diffusionTerms - convectionTerms + settings_.g[1]);
            }
        }
    }

    //! compute the right hand side of the Poisson equation for the pressure
    void Computation::computeRightHandSide()
    {
        //! TODO: < or <= or <-1? imaginary points
        //changed
        for (int i = 1; i < discretization_->rhsSize()[0] - 1; i++)
        {
            //changed
            for (int j = 1; j < discretization_->rhsSize()[1] - 1; j++)
            {
                double change_F = ((discretization_->f(i, j) - discretization_->f(i - 1, j)) / discretization_->dx());
                double change_G = ((discretization_->g(i, j) - discretization_->g(i, j - 1)) / discretization_->dy());
                discretization_->rhs(i, j) = (change_F + change_G) / dt_;
            }
        }
    }
    /**
     * Compute the pressure by solving the Poisson equation for the pressure.
     */
    void Computation::computePressure()
    {
        pressureSolver_->solve();
    }
    
    void Computation::computeVelocities()
    {
        // compute final u
        //xhanged
        for (int i = discretization_->uIBegin() + 1; i < discretization_->uIEnd() - 1; i++)
        {
            for (int j = discretization_->uJBegin() + 1; j < discretization_->uJEnd() ; j++)
            {
                discretization_->u(i, j) = discretization_->f(i, j) - dt_ * (discretization_->p(i + 1, j) - discretization_->p(i, j)) / discretization_->dx();
            }
        }
        // compute final v
        for (int i = discretization_->vIBegin() + 1; i < discretization_->vIEnd(); i++)
        {
            //changed
            for (int j = discretization_->vJBegin() + 1; j < discretization_->vJEnd() - 1; j++)
            {
                discretization_->v(i, j) = discretization_->g(i, j) - dt_ * (discretization_->p(i, j + 1) - discretization_->p(i, j)) / discretization_->dy();
            }
        }
    }
