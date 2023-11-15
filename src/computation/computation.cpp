#include "computation/computation.h"
#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/sor.h"

#include <cassert>

#include <cmath>

// TODO: constructor
/* Computation::Computation()
{
    // nothing to do
} */

void Computation::initialize(Settings settings)
{
    //! load and print settings
    settings_ = settings;
    settings_.printSettings();

    //! compute the meshWidth from the physical size and the number of cells
    meshWidth_[0] = settings_.physicalSize[0] / settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1] / settings_.nCells[1];

    //! initialize discretization
    if(settings_.useDonorCell) {
        discretization_ = std::make_shared<Discretization>(settings_.nCells, meshWidth_, settings_.alpha);
    } else {
        discretization_ = std::make_shared<Discretization>(settings_.nCells, meshWidth_);
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
/**
 * Run the whole simulation until tend. 
 */
//! TODO: to check @Louisa
void Computation::runSimulation()
{	double t = 0;
    while(t < settings_.endTime)
    {
        applyBoundaryValues();
        computeTimeStepWidth();
        computePreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        t += dt_;
        outputWriterParaview_ -> writeFile(t);
        outputWriterText_ -> writeFile(t);
    }
}
//! Compute the time step width dt from maximum velocities.
void Computation::computeTimeStepWidth()
{
    double dt_diffusion = ( settings_.re / 2 ) * ( pow(discretization_->dx(),2) * pow(discretization_->dy(),2) ) / (pow(discretization_->dx(),2) + pow(discretization_->dy(),2));

    double maxU = 0;
    double maxV = 0;
    std::array<int, 2> uSize = discretization_->uSize();
    //! TODO: danger zone
    for (int i = 0; i < uSize[0]; i++)
    {
        for (int j = 0; j < uSize[1]; j++)
        {
            //! possible because the grid for u and v have the same dimensions in all directions
            maxU = std::max(maxU, std::abs(discretization_->u(i,j)));
            maxV = std::max(maxV, std::abs(discretization_->v(i,j)));
        }
    }
    double dt_convection = std::min(discretization_->dx() / maxU, discretization_->dy() / maxV);

    //! TODO: ensure somewhere that tau<1
    //! because of the scaling with the security factor tau < 1, a subtraction of a small value,
    // to ensure dt smaller and not smaller/equal than required, is not necessary
    double dt = settings_.tau * std::min(dt_diffusion, dt_convection);
}

//! Set velocity boundary values
void Computation::applyBoundaryValues()
{
    //! set the u boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); i++)
    {
        //! u boundary values bottom, assuming inhomogenous Dirichlet conditions 
        discretization_->u(i, discretization_->uJBegin()) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin() + 1);
        //! u boundary values top, assuming inhomogenous Dirichlet conditions 
        //! TODO: why do we have dirichlet BC at top? I thought we had flow at top? Or do we set the movement in the dirichletBcTop?
        discretization_->u(i, discretization_->uJEnd()) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() - 1);
    }
    //! set the u boundary values for left and right side now
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++)
    {
        //! u boundary values left, assuming inhomogenous Dirichlet conditions 
        discretization_->u(discretization_->uIBegin(), j) = settings_.dirichletBcLeft[0];
        //! u boundary values right, assuming inhomogenous Dirichlet conditions
        discretization_->u(discretization_->uIEnd(), j) = settings_.dirichletBcRight[0];
    }

    //! set the v boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++)
    {
        //! v boundary values bottom, assuming inhomogenous Dirichlet conditions 
        discretization_->v(i, discretization_->vJBegin()) = settings_.dirichletBcBottom[1];
        //! v boundary values top, assuming inhomogenous Dirichlet conditions 
        //! TODO: find out whether we apply the BC on the very top, artificially added v, or one below
        discretization_->v(i, discretization_->vJEnd()) = settings_.dirichletBcTop[1];
    }

    //! set the v boundary values for left and right side now
    for (int j = discretization_->vJBegin(); j <= discretization_->vJEnd(); j++)
    {
        //! v boundary values left, assuming inhomogenous Dirichlet conditions 
        discretization_->v(discretization_->vIBegin(), j) = 2 * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin() + 1, j);
        //! v boundary values right, assuming inhomogenous Dirichlet conditions
        discretization_->v(discretization_->vIEnd(), j) = 2 * settings_.dirichletBcRight[1] - discretization_->v(discretization_->vIEnd() - 1, j);
    }
}

//! compute the preliminary velocities, F and G    
void Computation::computePreliminaryVelocities()
{
    //! TODO: danger zone, i<uIEnd() or i<=uIEnd()?
    //! compute the preliminary velocities only for the inside of the area, as on the direct boundary, they are constantly 0

    // ! compute preliminary F
    for (int i = discretization_->uIBegin() + 1; i < discretization_->uIEnd(); i++)
    {
        for (int j = discretization_->uJBegin() + 1; j < discretization_->uJEnd(); j++)
        {
            double diffusionTerms =  (1/settings_.re) * (discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j));
            double convectionTerms = discretization_->computeDu2Dx(i,j) + discretization_->computeDuvDy(i,j);
            //! TODO: check whether we update the g
            discretization_->f(i,j) = discretization_->u(i,j) + dt_ * (diffusionTerms - convectionTerms + settings_.g[0] );
        }
    }
    //! compute preliminary G
    for (int i = discretization_->vIBegin() + 1; i < discretization_->vIEnd(); i++)
    {
        for (int j = discretization_->vJBegin() + 1; j < discretization_->vJEnd(); j++)
        {
            double diffusionTerms =  (1/settings_.re) * (discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j));
            double convectionTerms = discretization_->computeDuvDx(i,j) + discretization_->computeDv2Dy(i,j);

            discretization_->g(i,j) = discretization_->v(i,j) + dt_ * (diffusionTerms - convectionTerms + settings_.g[1] );
}

void Computation::computeRightHandSide()
{
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
}
