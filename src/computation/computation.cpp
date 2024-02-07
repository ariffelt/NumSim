#include "computation/computation.h"
#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/sor.h"
#include "discretization/2_donor_cell.h"
#include "discretization/2_central_differences.h"
#include "freeflow/virtualparticle.h"

#include <cassert>
#include <cmath>

/**
 * Constructor for the Computation class.
 * Initializes the settings, discretization, pressure solver and output writers
 * @param argc number of command line arguments
 * @param argv command line arguments
 */
void Computation::initialize(int argc, char *argv[])
{
    // load and print settings
    settings_ = Settings();
    settings_.loadFromFile(argv[1]);
    settings_.printSettings();

    // compute the meshWidth from the physical size and the number of cells
    meshWidth_[0] = settings_.physicalSize[0] / settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1] / settings_.nCells[1];

    // initialize discretization
    if (settings_.useDonorCell)
    {
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha, settings_.gamma);
    }
    else
    {
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
    }

    // initialize the pressure solver
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

    // set markers for the fixed boundary conditions
    setBoundaryMarkers();

    // initialize output writers
    outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
}

/**
 * Set fixed boundary conditions markers to 2
 */
void Computation::setBoundaryMarkers()
{
    // set markers for the fixed boundary conditions
    for (int i = 0; i < discretization_->markerfieldSize()[0]; i++)
    {
        for (int j = 0; j < discretization_->markerfieldSize()[1]; j++)
        {
            if (i == 0 || i == discretization_->markerfieldSize()[0] - 1 || j == 0 || j == discretization_->markerfieldSize()[1] - 1)
            {
                discretization_->markerfield(i, j) = 2;
            }
        }
    }
}

/**
 * Test the boundary conditions by setting the boundary values and writing them out.
 */
void Computation::testBC()
{
    applyBoundaryValues();
    outputWriterParaview_->writeFile(0);
    outputWriterText_->writeFile(0);
}

/**
 * Run the simulation, starting at 0 until tend.
 * Prints 10 time steps to the console.
 */
void Computation::runSimulation()
{
    double t = 0.0;                                               // starting time
    double numberOfPrints = 10.0;                                 // number of prints to the console
    double timestepInterval = settings_.endTime / numberOfPrints; // time interval between two prints

    generateVirtualParticles();

    updateMarkerField();

    outputWriterParaview_->writeFile(t); // output initial state
    outputWriterText_->writeFile(t);

    // std::cout << "updated marker field" << std::endl;

    while (t < settings_.endTime)
    {
        setFountainVelocity(); // set the velocity of the fountain

        if (settings_.particelShape =="FOUNTAINWITHTEMPUP")
        {
            setFountainTemperatureUp(); // set the temperature of the fountain
        }

        applyBoundaryValues(); // set boundary values for u, v, F and G

        computeTimeStepWidth();
        if (t + dt_ > settings_.endTime) // decrease time step width in last time step, s.t. the end time will be reached exactly
        {
            dt_ = settings_.endTime - t;
        }
        t += dt_;
        std::cout << "t = " << t  << std::endl; // print time and time step width to the console

        computeTemperature(); // compute the temperature t

        // only compute preliminary velocities and rhs and solve pressure eq. on inner fluid cells
        
        computePreliminaryVelocities(); // compute preliminary velocities, F and G

        computeRightHandSide(); // compute rhs of the Poisson equation for the pressure
        
        computePressure(); // solve the Poisson equation for the pressure
        
        computeVelocities(); // compute the new velocities, u,v, from the preliminary velocities, F,G and the pressure, p
        
        outputWriterText_->writeFile(t);

        resetEmptyEdges();

        updateSurfacePs_ = false;
        updateSurfaceVelocities_ = true;

        freeflowBC(); // apply free flow boundary conditions

        //updateSurfacePs_ = true;
        //updateSurfaceVelocities_ = false;

        //freeflowBC(); // apply free flow boundary conditions

        outputWriterText_->writeFile(t);

        applyBoundaryValues(); // set boundary values for u, v, F and G

        std::cout << "vor updateParticleVelocities" << std::endl;

        computeParticleVelocities();

        std::cout << "nach updateParticleVelocities" << std::endl;

        updateMarkerField(); 

        std::cout << "nach updateMarkerfield" << std::endl;

        resetEmptyEdges();


        updateSurfacePs_ = false;
        updateSurfaceVelocities_ = true;

        freeflowBC(); // apply free flow boundary conditions

        updateSurfacePs_ = true;
        updateSurfaceVelocities_ = false;

        freeflowBC(); // apply free flow boundary conditions

        std::cout << "nach freeFlowBC" << std::endl;

        outputWriterParaview_->writeFile(t); // output simulation results
        outputWriterText_->writeFile(t);

        std::cout << "nach Outputwriter" << std::endl;

        if (t >= timestepInterval)
        {
            std::cout << "t = " << t << ", dt = " << dt_ << std::endl; // print time and time step width to the console
            timestepInterval += settings_.endTime / numberOfPrints;
        }
    }
}

void Computation::printParticles()
{
    for (int i = 0; i < particlesX_.size(); i++)
    {
        std::cout << "\t Particle " << i << " at position (" << particlesX_[i] << ", " << particlesY_[i] << ")" << std::endl;
    }
}

/**
 * Compute the time step width dt from maximum velocities.
 */
void Computation::computeTimeStepWidth()
{
    // compute maximal time step width from diffusion
    double dt_diffusion = (settings_.re / 2) * (pow(discretization_->dx(), 2) * pow(discretization_->dy(), 2)) / (pow(discretization_->dx(), 2) + pow(discretization_->dy(), 2));
    double dt_diffusion_temp = dt_diffusion * settings_.pr;

    double maxU = 0;
    double maxV = 0;
    std::array<int, 2> uSize = discretization_->uSize();

    for (int i = 0; i < uSize[0]; i++)
    {
        for (int j = 0; j < uSize[1]; j++)
        {
            // possible because the grid for u and v have the same dimensions in all directions
            maxU = std::max(maxU, std::fabs(discretization_->u(i, j)));
            maxV = std::max(maxV, std::fabs(discretization_->v(i, j)));
        }
    }

    double dt_convection = std::min(discretization_->dx() / maxU, discretization_->dy() / maxV);

    // subtraction of small value to ensure dt smaller and not smaller/equal than required not necessary since we scale with security factor tau < 1
    assert(settings_.tau < 1);
    //dt_ = std::min(settings_.tau * std::min(dt_diffusion, dt_convection), settings_.maximumDt);
    dt_ = std::min(settings_.tau * std::min(std::min(dt_diffusion, dt_diffusion_temp), dt_convection), settings_.maximumDt);
}

/**
 * Set velocity boundary values for u, v, F and G
 */
void Computation::applyBoundaryValues()
{
    // set u boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        //if (discretization_->markerfield(i,discretization_->uJBegin()) == 1)
        {
            // u boundary values bottom, assuming inhomogenous Dirichlet conditions
            discretization_->u(i, discretization_->uJBegin()) = 2.0 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin() + 1);
        }
        //if (discretization_->markerfield(i,discretization_->uJEnd()) == 1)
        {
            // u boundary values top, assuming inhomogenous Dirichlet conditions
            discretization_->u(i, discretization_->uJEnd()) = 2.0 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() - 1);
        }
    }

    // set u boundary values for left and right side
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++)
    {
        //if (discretization_->markerfield(discretization_->uIBegin(),j) == 1)
        {
            // u boundary values left, assuming inhomogenous Dirichlet conditions
            discretization_->u(discretization_->uIBegin(), j) = settings_.dirichletBcLeft[0];
        }
        //if (discretization_->markerfield(discretization_->uIEnd() - 1,j) == 1)
        {
            // u boundary values right, assuming inhomogenous Dirichlet conditions
            discretization_->u(discretization_->uIEnd() - 1, j) = settings_.dirichletBcRight[0];
        }
    }

    // set v boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++)
    {
        //if (discretization_->markerfield(i,discretization_->vJBegin()) == 1)
        {
            // v boundary values bottom, assuming inhomogenous Dirichlet conditions
            discretization_->v(i, discretization_->vJBegin()) = settings_.dirichletBcBottom[1];
        }
        //if (discretization_->markerfield(i,discretization_->vJEnd() - 1) == 1)
        {
            // v boundary values top, assuming inhomogenous Dirichlet conditions
            discretization_->v(i, discretization_->vJEnd() - 1) = settings_.dirichletBcTop[1];
        }
    }

    // set v boundary values for left and right side
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        //if (discretization_->markerfield(discretization_->vIBegin(),j) == 1)
        {
            // v boundary values left, assuming inhomogenous Dirichlet conditions
            discretization_->v(discretization_->vIBegin(), j) = 2 * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin() + 1, j);
        }
        //if (discretization_->markerfield(discretization_->vIEnd(),j) == 1)
        {
            // v boundary values right, assuming inhomogenous Dirichlet conditions
            discretization_->v(discretization_->vIEnd(), j) = 2 * settings_.dirichletBcRight[1] - discretization_->v(discretization_->vIEnd() - 1, j);
        }
    }

    // set F boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        //if (discretization_->markerfield(i,discretization_->uJBegin()) == 1)
        {
            // F boundary values bottom
            discretization_->f(i, discretization_->uJBegin()) = discretization_->u(i, discretization_->uJBegin());
        }
        //if (discretization_->markerfield(i,discretization_->uJEnd()) == 1)
        {
            // F boundary values top
            discretization_->f(i, discretization_->uJEnd()) = discretization_->u(i, discretization_->uJEnd());
        }
    }

    // set F boundary values for left and right side
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++)
    {
        //if (discretization_->markerfield(discretization_->uIBegin(),j) == 1)
        {
            // F boundary values left
            discretization_->f(discretization_->uIBegin(), j) = discretization_->u(discretization_->uIBegin(), j);
        }
        //if (discretization_->markerfield(discretization_->uIEnd() - 1,j) == 1)
        {
            // F boundary values right
            discretization_->f(discretization_->uIEnd() - 1, j) = discretization_->u(discretization_->uIEnd() - 1, j);
        }
    }

    // set G boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++)
    {
        //if (discretization_->markerfield(i,discretization_->vJBegin()) == 1)
        {
            // G boundary values bottom
            discretization_->g(i, discretization_->vJBegin()) = discretization_->v(i, discretization_->vJBegin());
        }
        //if (discretization_->markerfield(i,discretization_->vJEnd() - 1) == 1)
        {
            // G boundary values top
            discretization_->g(i, discretization_->vJEnd() - 1) = discretization_->v(i, discretization_->vJEnd() - 1);
        }
    }

    // set G boundary values for left and right side
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        //if (discretization_->markerfield(discretization_->vIBegin(),j) == 1)
        {
            // G boundary values left
            discretization_->g(discretization_->vIBegin(), j) = discretization_->v(discretization_->vIBegin(), j);
        }
        //if (discretization_->markerfield(discretization_->vIEnd(),j) == 1)
        {
            // G boundary values right
            discretization_->g(discretization_->vIEnd(), j) = discretization_->v(discretization_->vIEnd(), j);
        }
    }

    // set T boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->tIBegin(); i <= discretization_->tIEnd(); i++)
    {
        // T boundary values bottom, assuming inhomogenous Dirichlet conditions
        discretization_->t(i, discretization_->tJBegin()) = 2.0 * settings_.dirichletBcBottomT - discretization_->t(i, discretization_->tJBegin() + 1);
        // T boundary values top, assuming inhomogenous Dirichlet conditions
        discretization_->t(i, discretization_->tJEnd()) = 2.0 * settings_.dirichletBcTopT - discretization_->t(i, discretization_->tJEnd() - 1);
    }

    // set T boundary values for left and right side
    for (int j = discretization_->tJBegin(); j <= discretization_->tJEnd(); j++)
    {
        // T boundary values left, assuming inhomogenous Dirichlet conditions
        discretization_->t(discretization_->tIBegin(), j) = 2.0 * settings_.dirichletBcLeftT - discretization_->t(discretization_->tIBegin() + 1, j);
        // T boundary values right, assuming inhomogenous Dirichlet conditions
        discretization_->t(discretization_->tIEnd(), j) = 2.0 * settings_.dirichletBcRightT - discretization_->t(discretization_->tIEnd() - 1, j);
    }
}

/**
 * Compute the temperature, T, from the velocities, u,v
 */
void Computation::computeTemperature()
{
    for (int i = discretization_->tIBegin() + 1; i < discretization_->tIEnd(); i++)
    {
        for (int j = discretization_->tJBegin() + 1; j < discretization_->tJEnd(); j++)
        {
            if (discretization_->isInnerFluidCell(i,j))
            {
                double diffusionTerms = (1 / (settings_.re * settings_.pr)) * (discretization_->computeD2tDx2(i, j) + discretization_->computeD2tDy2(i, j));
                double convectionTerms = discretization_->computeDutDx(i, j) + discretization_->computeDvtDy(i, j);

                discretization_->t(i, j) = discretization_->t(i, j) + dt_ * (diffusionTerms - convectionTerms + discretization_->q(i, j));
            }
        }
    }
}

/**
 * Compute the preliminary velocities, F and G
 * F and G are computed only for the inside of the area, as on the direct boundary, they are constantly 0
 */
void Computation::computePreliminaryVelocities()
{
    // compute preliminary F
    for (int i = discretization_->uIBegin() + 1; i < discretization_->uIEnd() - 1; i++)
    {
        for (int j = discretization_->uJBegin() + 1; j < discretization_->uJEnd(); j++)
        {
            if (discretization_->isFluidCell(i,j) && discretization_->isFluidCell(i+1,j))
            {
                double diffusionTerms = (1 / settings_.re) * (discretization_->computeD2uDx2(i, j) + discretization_->computeD2uDy2(i, j));
                double convectionTerms = discretization_->computeDu2Dx(i, j) + discretization_->computeDuvDy(i, j);

                discretization_->f(i, j) = discretization_->u(i, j) + dt_ * (diffusionTerms - convectionTerms + settings_.g[0]) - (dt_ * settings_.beta *settings_.g[0] * (discretization_->t(i, j) + discretization_->t(i + 1, j)) / 2);
            }
        }
    }

    // compute preliminary G
    for (int i = discretization_->vIBegin() + 1; i < discretization_->vIEnd(); i++)
    {
        for (int j = discretization_->vJBegin() + 1; j < discretization_->vJEnd() - 1; j++)
        {
            if (discretization_->isFluidCell(i,j) && discretization_->isFluidCell(i,j+1))
            {
                double diffusionTerms = (1 / settings_.re) * (discretization_->computeD2vDx2(i, j) + discretization_->computeD2vDy2(i, j));
                double convectionTerms = discretization_->computeDuvDx(i, j) + discretization_->computeDv2Dy(i, j);

                discretization_->g(i, j) = discretization_->v(i, j) + dt_ * (diffusionTerms - convectionTerms + settings_.g[1]) - (dt_ * settings_.beta *settings_.g[1] * (discretization_->t(i, j) + discretization_->t(i, j + 1)) / 2);
                // std::cout << "computed G for (" << i << ", " << j <<"): " <<discretization_->g(i, j) << std::endl;
            }
        }
    }
}

/**
 * Compute the right hand side of the Poisson equation for the pressure.
 */
void Computation::computeRightHandSide()
{
    for (int i = 1; i < discretization_->rhsSize()[0] - 1; i++)
    {
        for (int j = 1; j < discretization_->rhsSize()[1] - 1; j++)
        {
            if (discretization_->isInnerFluidCell(i,j))
            {
                double change_F = ((discretization_->f(i, j) - discretization_->f(i - 1, j)) / discretization_->dx());
                double change_G = ((discretization_->g(i, j) - discretization_->g(i, j - 1)) / discretization_->dy());

                discretization_->rhs(i, j) = (change_F + change_G) / dt_;
            }
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

/**
 * Compute the new velocities, u,v, from the preliminary velocities, F,G and the pressure, p.
 */
void Computation::computeVelocities()
{
    // compute final u
    for (int i = discretization_->uIBegin() + 1; i < discretization_->uIEnd() - 1; i++)
    {
        for (int j = discretization_->uJBegin() + 1; j < discretization_->uJEnd(); j++)
        {
            if (discretization_->isFluidCell(i,j) && discretization_->isFluidCell(i+1,j))
            {
                discretization_->u(i, j) = discretization_->f(i, j) - dt_ * (discretization_->p(i + 1, j) - discretization_->p(i, j)) / discretization_->dx();
            }
        }
    }

    // compute final v
    for (int i = discretization_->vIBegin() + 1; i < discretization_->vIEnd(); i++)
    {
        for (int j = discretization_->vJBegin() + 1; j < discretization_->vJEnd() - 1; j++)
        {
            if (discretization_->isFluidCell(i,j) && discretization_->isFluidCell(i,j+1))
            {
                discretization_->v(i, j) = discretization_->g(i, j) - dt_ * (discretization_->p(i, j + 1) - discretization_->p(i, j)) / discretization_->dy();
            }
        }
    }
}

/**
 * Introduce virtual particles in initial state.
 * Equally distribute them in the whole domain on a finer mesh.
 * Only placed in volume filled with fluid, not obstacles or air.
 */
void Computation::generateVirtualParticles()
{
    // todo: remove hardcoding
    if (settings_.particelShape == "DAM")
    {
        generateDam(16);
    }
    else if (settings_.particelShape == "FULL")
    {
        generateFull(20);
    }
    else if (settings_.particelShape == "BOX")
    {
        generateBox(20);
    }
    else if (settings_.particelShape == "DROP")
    {
        generateDrop(1);
    }
    else if (settings_.particelShape == "BIGDROP")
    {
        generateBigDrop(4);
    }
    else if (settings_.particelShape == "BAR")
    {
        generateBar(10);
    }
    else if(settings_.particelShape =="DROPWATER")
    {
        generateDropInWater(10);
    }
    else if(settings_.particelShape =="FOUNTAIN")
    {
        generateFountain(10, 100);
    }
    else if(settings_.particelShape =="FOUNTAINWITHTEMP")
    {
        generateFountainWithTemp(10, 100);
    }
    else if(settings_.particelShape =="FOUNTAINWITHTEMPUP")
    {
        generateFountainWithTempUp(10, 100);
    }
    else
    {
    particlesX_ = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0,1.0, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1,1.1,1.1};	
    particlesY_ = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0,1.0, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1,1.1,1.1};	
    }
}

void Computation::generateDam(int noParticles)
{
    // Distribute the noParticles equally in a box in the left lower corner of the domain
    double dx = discretization_->dx();
    double dy = discretization_->dy();

    particlesX_ = {};
    particlesY_ = {};

    for (int i=int(settings_.nCells[1]/4); i<int(3*settings_.nCells[1]/4); i++)
    {
        for (int j=0; j<int(settings_.nCells[1]) - 3; j++)
        {
            for (int k=0; k<noParticles; k++)
            {
                particlesX_.push_back(i*dx + k*dx/noParticles);
                particlesY_.push_back(j*dy + k*dy/noParticles);
            }
        }
    }
}

void Computation::generateFountain(int noParticles, int noParticlesFountain)
{
    // Distribute the noParticles equally in a box in the left lower corner of the domain
    double dx = discretization_->dx();
    double dy = discretization_->dy();

    particlesX_ = {};
    particlesY_ = {};

    // for (int i=int(settings_.nCells[1]/4); i<int(3*settings_.nCells[1]/4); i++)
    for (int i=int(0); i<int(settings_.nCells[0]); i++)

    {
        //for (int j=0; j<int(settings_.nCells[1]); j++)
        for (int j=0; j<int(settings_.nCells[1])/5; j++)
        {
            if ((i <= int(settings_.nCells[0]/2 + 5)) && (i >= int(settings_.nCells[0]/2 - 5)) && (j <= int(settings_.nCells[1]/5)) && (j >= int(settings_.nCells[1]/5 - 6))){

                    for (int k=0; k<noParticlesFountain; k++)
                    {
                        particlesX_.push_back(i*dx + k*dx/noParticlesFountain);
                        particlesY_.push_back(j*dy + k*dy/noParticlesFountain);
                    }
                }
            else {
                for (int k=0; k<noParticles; k++)
                {
                    particlesX_.push_back(i*dx + k*dx/noParticles);
                    particlesY_.push_back(j*dy + k*dy/noParticles);
                }
            }
        }
    }
}

void Computation::generateFountainWithTemp(int noParticles, int noParticlesFountain)
{
    setFountainTemperature();
    // Distribute the noParticles equally in a box in the left lower corner of the domain
    double dx = discretization_->dx();
    double dy = discretization_->dy();

    particlesX_ = {};
    particlesY_ = {};

    // for (int i=int(settings_.nCells[1]/4); i<int(3*settings_.nCells[1]/4); i++)
    for (int i=int(0); i<int(settings_.nCells[0]); i++)

    {
        //for (int j=0; j<int(settings_.nCells[1]); j++)
        for (int j=0; j<int(settings_.nCells[1])/5; j++)
        {
            if ((i <= int(settings_.nCells[0]/2 + 5)) && (i >= int(settings_.nCells[0]/2 - 5)) && (j <= int(settings_.nCells[1]/5)) && (j >= int(settings_.nCells[1]/5 - 6))){

                    for (int k=0; k<noParticlesFountain; k++)
                    {
                        particlesX_.push_back(i*dx + k*dx/noParticlesFountain);
                        particlesY_.push_back(j*dy + k*dy/noParticlesFountain);
                    }
                }
            else {
                for (int k=0; k<noParticles; k++)
                {
                    particlesX_.push_back(i*dx + k*dx/noParticles);
                    particlesY_.push_back(j*dy + k*dy/noParticles);
                }
            }
        }
    }
}

void Computation::generateFountainWithTempUp(int noParticles, int noParticlesFountain)
{
    // Distribute the noParticles equally in a box in the left lower corner of the domain
    double dx = discretization_->dx();
    double dy = discretization_->dy();

    particlesX_ = {};
    particlesY_ = {};

    // for (int i=int(settings_.nCells[1]/4); i<int(3*settings_.nCells[1]/4); i++)
    for (int i=int(0); i<int(settings_.nCells[0]); i++)

    {
        //for (int j=0; j<int(settings_.nCells[1]); j++)
        for (int j=0; j<int(settings_.nCells[1])/5; j++)
        {
            if ((i <= int(settings_.nCells[0]/2 + 5)) && (i >= int(settings_.nCells[0]/2 - 5)) && (j <= int(settings_.nCells[1]/5)) && (j >= int(settings_.nCells[1]/5 - 6))){

                    for (int k=0; k<noParticlesFountain; k++)
                    {
                        particlesX_.push_back(i*dx + k*dx/noParticlesFountain);
                        particlesY_.push_back(j*dy + k*dy/noParticlesFountain);
                    }
                }
            else {
                for (int k=0; k<noParticles; k++)
                {
                    particlesX_.push_back(i*dx + k*dx/noParticles);
                    particlesY_.push_back(j*dy + k*dy/noParticles);
                }
            }
        }
    }
}

void Computation::setFountainVelocity()
{
    discretization_->v(int(settings_.nCells[0]/2 - 1), int(settings_.nCells[1]/5 - 2)) = 2;
    discretization_->v(int(settings_.nCells[0]/2 + 1), int(settings_.nCells[1]/5 - 2)) = 2;

    discretization_->v(int(settings_.nCells[0]/2), int(settings_.nCells[1]/5 - 3)) = 2;

    discretization_->v(int(settings_.nCells[0]/2 - 1), int(settings_.nCells[1]/5 - 4)) = 2;
    discretization_->v(int(settings_.nCells[0]/2 + 1), int(settings_.nCells[1]/5 - 4)) = 2;

    discretization_->v(int(settings_.nCells[0]/2), int(settings_.nCells[1]/5 - 5)) = 2;

    discretization_->v(int(settings_.nCells[0]/2 - 1), int(settings_.nCells[1]/5 - 6)) = 2;
    discretization_->v(int(settings_.nCells[0]/2 + 1), int(settings_.nCells[1]/5 - 6)) = 2;
}

void Computation::setFountainTemperature()
{
    discretization_->q(int(settings_.nCells[0]/2 - 1), int(settings_.nCells[1]/5 - 2)) = 2;
    discretization_->q(int(settings_.nCells[0]/2 + 1), int(settings_.nCells[1]/5 - 2)) = 2;

    discretization_->q(int(settings_.nCells[0]/2), int(settings_.nCells[1]/5 - 3)) = 2;

    discretization_->q(int(settings_.nCells[0]/2 - 1), int(settings_.nCells[1]/5 - 4)) = 2;
    discretization_->q(int(settings_.nCells[0]/2 + 1), int(settings_.nCells[1]/5 - 4)) = 2;

    discretization_->q(int(settings_.nCells[0]/2), int(settings_.nCells[1]/5 - 5)) = 2;

    discretization_->q(int(settings_.nCells[0]/2 - 1), int(settings_.nCells[1]/5 - 6)) = 2;
    discretization_->q(int(settings_.nCells[0]/2 + 1), int(settings_.nCells[1]/5 - 6)) = 2;
}

void Computation::setFountainTemperatureUp()
{
    int maxHeight = 0;
    for (int j = 1; j < int(settings_.nCells[1]) - 1; j++)
    {
       if (discretization_->isInnerFluidCell(int(settings_.nCells[0]/2), j))
       {
           maxHeight = j;
       }
    }
    // discretization_->q(int(settings_.nCells[0]/2 - 1), int(settings_.nCells[1]/5 - 2)) = 2;
    // discretization_->q(int(settings_.nCells[0]/2 + 1), int(settings_.nCells[1]/5 - 2)) = 2;

    discretization_->q(int(settings_.nCells[0]/2), maxHeight) = 2;
    discretization_->q(int(settings_.nCells[0]/2), maxHeight+1) = 2;
    discretization_->q(int(settings_.nCells[0]/2), maxHeight-1) = 2;
    discretization_->q(int(settings_.nCells[0]/2 + 1), maxHeight) = 2;
    discretization_->q(int(settings_.nCells[0]/2 - 1), maxHeight) = 2;
    

    // discretization_->q(int(settings_.nCells[0]/2 - 1), int(settings_.nCells[1]/5 - 4)) = 2;
    // discretization_->q(int(settings_.nCells[0]/2 + 1), int(settings_.nCells[1]/5 - 4)) = 2;

    // discretization_->q(int(settings_.nCells[0]/2), int(settings_.nCells[1]/5 - 5)) = 2;

    // discretization_->q(int(settings_.nCells[0]/2 - 1), int(settings_.nCells[1]/5 - 6)) = 2;
    // discretization_->q(int(settings_.nCells[0]/2 + 1), int(settings_.nCells[1]/5 - 6)) = 2;
}

void Computation::generateBox(int noParticles)
{
    // Generate a box in the middle of the domain
    double dx = discretization_->dx();
    double dy = discretization_->dy();

    particlesX_ = {};
    particlesY_ = {};

    for (int i=int(3*settings_.nCells[0]/8); i<int(5*settings_.nCells[0]/8); i++)
    {
        for (int j=int(3*settings_.nCells[1]/8); j<int(5*settings_.nCells[1]/8); j++)
        {
            for (int k=0; k<noParticles; k++)
            {
                particlesX_.push_back(i*dx + k*dx/noParticles);
                particlesY_.push_back(j*dy + k*dy/noParticles);
            }
        }
    }
}

void Computation::generateBigDrop(int noParticles)
{
    particlesX_ = {};
    particlesY_ = {};

    for (int i = 0; i<2; i++)
    {
        for (int j = 0; j<2;j++)
        {
            particlesX_.push_back(settings_.nCells[0]/ 2 * discretization_->dx()+discretization_->dx()/2 - i*discretization_->dx());
            particlesY_.push_back(settings_.nCells[1]/ 2 * discretization_->dy()+discretization_->dy()/2 - j* discretization_->dy());  
        }
    }
}

void Computation::generateDrop(int noParticles)
{
    // Generate a single particle in the middle of the domain
    particlesX_ = {};
    particlesY_ = {};

    particlesX_.push_back(settings_.nCells[0]/ 2 * discretization_->dx()+discretization_->dx()/2);
    particlesY_.push_back(settings_.nCells[1]/ 2 * discretization_->dy()+discretization_->dy()/2);
}

void Computation::generateFull(int noParticles)
{
    // Distribute the noParticles equally in a box in the left lower corner of the domain
    double dx = discretization_->dx();
    double dy = discretization_->dy();

    particlesX_ = {};
    particlesY_ = {};

    for (int i=0; i<int(settings_.nCells[0]); i++)
    {
        for (int j=0; j<int(settings_.nCells[1]); j++)
        {
            particlesX_.push_back(i*dx);
            particlesY_.push_back(j*dy);
        }
    }
}

void Computation::generateBar(int noParticles)
{

    double dx = discretization_->dx();
    double dy = discretization_->dy();
    particlesX_ = {};
    particlesY_ = {};

    for (int i=0; i<int(settings_.nCells[0]); i++)
    {
        particlesX_.push_back(i*dx);
        particlesX_.push_back(i*dx);
        particlesY_.push_back((settings_.nCells[1]-1)*dy);
        particlesY_.push_back((settings_.nCells[1]-2)*dy);
    }
}

void Computation::generateDropInWater(int noParticles)
{
    double dx = discretization_->dx();
    double dy = discretization_->dy();

    for (int i = 0; i< settings_.nCells[0]; i++)
    {
        for (int j = 0; j< int(settings_.nCells[1]/16); j++)
        {
            for (int k=0; k<noParticles; k++)
            {
                particlesX_.push_back(i*dx + k*dx/noParticles);
                particlesY_.push_back(j*dy + k*dy/noParticles);
            }
        }
    }
    for (int i=int(3*settings_.nCells[0]/8); i<int(5*settings_.nCells[0]/8); i++)
    {
        for (int j=int(3*settings_.nCells[1]/8); j<int(5*settings_.nCells[1]/8); j++)
        {
            for (int k=0; k<noParticles; k++)
            {
                particlesX_.push_back(i*dx + k*dx/noParticles);
                particlesY_.push_back(j*dy + k*dy/noParticles);
            }
        }
    }
}
/**
 * Compute the new particle velocities.
 * And move particles according to these.
 */
void Computation::computeParticleVelocities()
{
    // double dx = discretization_->dx();
    // double dy = discretization_->dy();

    // interpolate velocities to the particle positions (do not coincide with velocity grid points)
    for (int k = 0; k < particlesX_.size(); k++)
    {
    //     // compute particle velocity in x direction

    //     // index of upper right corner
    //     int iUpperRight = int(particlesX_[k] / dx + 1);
    //     int jUpperRight = int((particlesY_[k] + 1 / 2) / dy + 1);

    //     // position of 4 neighbouring grid points with values u
    //     double x1 = (iUpperRight - 1) * dx;
    //     double x2 = iUpperRight * dx;
    //     double y1 = (jUpperRight - 1) * dy - dy / 2;
    //     double y2 = jUpperRight * dy - dy / 2;



    //     // bilinear interpolation
    //     double u = 1 / (dx * dy) * ((x2 - particlesX_[k]) * (y2 - particlesY_[k]) * discretization_->u(iUpperRight - 1, jUpperRight - 1) 
    //                                 + (particlesX_[k] - x1) * (y2 - particlesY_[k]) * discretization_->u(iUpperRight, jUpperRight - 1) 
    //                                 + (x2 - particlesX_[k]) * (particlesY_[k] - y1) * discretization_->u(iUpperRight - 1, jUpperRight) 
    //                                 + (particlesX_[k] - x1) * (particlesY_[k] - y1) * discretization_->u(iUpperRight, jUpperRight));

    //     // compute particle velocity in y direction
    //     // index of upper right corner
    //     iUpperRight = int((particlesX_[k] + 1 / 2) / dx + 1);
    //     jUpperRight = int(particlesY_[k] / dy + 1);

    //     // position of 4 neighbouring grid points with values v
    //     x1 = (iUpperRight - 1) * dx - dx / 2;
    //     x2 = iUpperRight * dx - dx / 2;
    //     y1 = (jUpperRight - 1) * dy;
    //     y2 = jUpperRight * dy;

    //     // bilinear interpolation
    //     double v = 1 / (dx * dy) * ((x2 - particlesX_[k]) * (y2 - particlesY_[k]) * discretization_->v(iUpperRight - 1, jUpperRight - 1) 
    //                                 + (particlesX_[k] - x1) * (y2 - particlesY_[k]) * discretization_->v(iUpperRight, jUpperRight - 1) 
    //                                 + (x2 - particlesX_[k]) * (particlesY_[k] - y1) * discretization_->v(iUpperRight - 1, jUpperRight) 
    //                                 + (particlesX_[k] - x1) * (particlesY_[k] - y1) * discretization_->v(iUpperRight, jUpperRight));

    //     particlesX_[k] += dt_ * u;
    //     particlesY_[k] += dt_ * v;

    // move particle
    particlesX_[k] += dt_ * discretization_->u().interpolateAt(particlesX_[k], particlesY_[k]);
    particlesY_[k] += dt_ * discretization_->v().interpolateAt(particlesX_[k], particlesY_[k]);
    }

    // std::cout << "Particle velocities computed" << std::endl;
}

/**
 * Update the cell types.
 */
void Computation::updateMarkerField()
{
    // update marker field
    for (int i = 0; i < discretization_->markerfieldSize()[0]; i++)
    {
        for (int j = 0; j < discretization_->markerfieldSize()[1]; j++)
        {   
            // check if it is a fixed boundary cell
            if (discretization_->markerfield(i, j) == 1)
            {   
                discretization_->markerfield(i, j) = 0; // assume cell is empty
            }
        }
    }

    for (int k = 0; k < particlesX_.size(); k++)
    {
        int i = int(particlesX_[k] / discretization_->dx() + 1);
        int j = int(particlesY_[k] / discretization_->dy() + 1);
        discretization_->markerfield(i, j) = 1; // assume cell is fluid    
    }
}

/**
 * Apply free flow boundary conditions.
 */
void Computation::freeflowBC()
{
    // for loops over all p cells
    // check if cell is fluid cell
    // check if it has 
    for (int i = 1; i <= discretization_->pSize()[0] - 1; i++)
    {
        for (int j = 1; j <= discretization_->pSize()[1] - 1; j++)
        {
            if (!discretization_->isInnerFluidCell(i,j) && discretization_->markerfield(i, j) == 1)
            {
                // check if the right, top, left or bottom cell is a fluid cell
                if (discretization_->markerfield(i + 1, j) >= 1)
                {
                    if (discretization_->markerfield(i, j + 1) >= 1)
                    {
                        if (discretization_->markerfield(i - 1, j) >= 1)
                        {
                            // cell is surrounded by 3 fluid cells (bottom wall)
                            bottomWallBC(i,j);
                            continue;
                        }
                        else
                        {
                            if (discretization_->markerfield(i, j - 1) >= 1)
                            {
                                // cell is surrounded by 3 fluid cells (left wall)
                                leftWallBC(i,j);
                                continue;
                            }
                            else
                            {
                                // cell is surrounded by 2 fluid cells (left, bottom corner)
                                bottomLeftCornerBC(i,j);
                                continue;
                            }
                        }
                    }
                    else
                    {
                        if (discretization_->markerfield(i - 1, j) >= 1)
                        {
                            if (discretization_->markerfield(i, j - 1) >= 1)
                            {
                                // cell is surrounded by 3 fluid cells (top wall)
                                topWallBC(i,j);
                                continue;
                            }
                            else
                            {
                                // cell is surrounded by 2 fluid cells (horizontal pipe)
                                horizontalPipeBC(i,j);
                                continue;
                            }
                        }
                        else
                        {
                            if (discretization_->markerfield(i, j - 1) >= 1)
                            {
                                // cell is surrounded by 2 fluid cells (top, left corner)
                                topLeftCornerBC(i,j);
                                continue;
                            }
                            else
                            {
                                // cell is surrounded by 1 fluid cell (tip from right)
                                tipFromRightBC(i,j);
                                continue;
                            }
                        }
                    }
                }
                else 
                {
                    if (discretization_->markerfield(i, j + 1) >= 1)
                    {
                        if (discretization_->markerfield(i - 1, j) >= 1)
                        {
                            if (discretization_->markerfield(i, j - 1) >= 1)
                            {
                                // cell is surrounded by 3 fluid cells (right wall)
                                rightWallBC(i,j);
                                continue;
                            }
                            else
                            {
                                // cell is surrounded by 2 fluid cells (bottom, right corner)
                                bottomRightCornerBC(i,j);
                                continue;
                            }
                        }
                        else
                        {
                            if (discretization_->markerfield(i, j - 1) >= 1)
                            {
                                // cell is surrounded by 2 fluid cells (vertical pipe)
                                verticalPipeBC(i,j);
                                continue;
                            }
                            else
                            {
                                // cell is surrounded by 1 fluid cell (tip from top)
                                tipFromTopBC(i,j);
                                continue;
                            }
                        }
                    }
                    else
                    {
                        if (discretization_->markerfield(i - 1, j) >= 1)
                        {
                            if (discretization_->markerfield(i, j - 1) >= 1)
                            {
                                // cell is surrounded by 2 fluid cells (top, right corner)
                                topRightCornerBC(i,j);
                                continue;
                            }
                            else
                            {
                                // cell is surrounded by 1 fluid cell (tip from left)
                                tipFromLeftBC(i,j);
                                continue;
                            }
                        }
                        else
                        {
                            if (discretization_->markerfield(i, j - 1) >= 1)
                            {
                                // cell is surrounded by 1 fluid cell (tip from bottom)
                                tipFromBottomBC(i,j);
                                continue;
                            }
                            else
                            {
                                // cell is surrounded by 0 fluid cells (drop)
                                dropBC(i,j);
                                continue;
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * compute the boundary conditions for the bottom wall
 */
void Computation::bottomWallBC(int i, int j)
{   
    // std::cout << "bottom wall ( "<<i<<", "<<j<<"); ";

    // mass balance
    if (updateSurfaceVelocities_)
    {
        discretization_->v(i,j-1) = discretization_->v(i,j) + discretization_->dy() / discretization_->dx() * (discretization_->u(i,j) - discretization_->u(i-1,j));

        // tangential stress
        if (discretization_->markerfield(i-1,j-1) == 0)
        {
            discretization_->u(i-1,j-1) = discretization_->u(i-1,j) + discretization_->dy() / discretization_->dx() * (discretization_->v(i,j-1) - discretization_->v(i-1,j-1));
        }
    }

    if (updateSurfacePs_)
    {
        // normal stress
        discretization_->p(i,j) = 2.0 / settings_.re * (discretization_->v(i,j) - discretization_->v(i,j-1)) / discretization_->dy();

        //temperature
        discretization_->t(i,j) = 0.0;
    }
}

/**
 * compute the boundary conditions for the left wall
 */
void Computation::leftWallBC(int i, int j)
{
    // std::cout << "left wall ( "<<i<<", "<<j<<"); ";
    if (updateSurfaceVelocities_)
    {
        // mass balance
        discretization_->u(i-1,j) = discretization_->u(i,j) + discretization_->dx() / discretization_->dy() * (discretization_->v(i,j) - discretization_->v(i,j-1));

        // tangential stress
        if (discretization_->markerfield(i-1,j-1) == 0)
        {
            discretization_->v(i-1,j-1) = discretization_->v(i,j-1) + discretization_->dx() / discretization_->dy() * (discretization_->u(i-1,j) - discretization_->u(i-1,j-1));
        }
    }

    // normal stress
    if (updateSurfacePs_)
    {
        discretization_->p(i,j) = 2.0 / settings_.re * (discretization_->u(i,j) - discretization_->u(i-1,j)) / discretization_->dx();

        //temperature
        discretization_->t(i,j) = 0.0;
    }
}   

/**
 * compute the boundary conditions for the bottom left corner
 */
void Computation::bottomLeftCornerBC(int i, int j)
{
    // std::cout << "bottom left corner ( "<<i<<", "<<j<<"); ";

    if (updateSurfaceVelocities_)
    {
        // mass balance + tangential stress
        discretization_->u(i-1,j) = discretization_->u(i,j);
        discretization_->v(i,j-1) = discretization_->v(i,j);

        if (discretization_->markerfield(i-1,j-1) == 0)
        {
            discretization_->u(i-1,j-1) = discretization_->u(i-1,j);
            discretization_->v(i-1,j-1) = discretization_->v(i,j-1);
        }
    }
                                
    // normal stress
    if (updateSurfacePs_)
    {
        discretization_->p(i,j) = 1.0 / (2.0 * settings_.re) * ((discretization_->u(i,j+1) + discretization_->u(i-1,j+1) - discretization_->u(i,j) - discretization_->u(i-1,j)) / (discretization_->dy()) +
                                                             (discretization_->v(i+1,j) + discretization_->v(i+1,j-1) - discretization_->v(i,j) - discretization_->v(i,j-1)) / (discretization_->dx()));
        
        //temperature   
        discretization_->t(i,j) = 0.0;
    }
    // std::cout << "\n USING VALUES: v(i,j)-,v(i+1,j)+,v(i+1,j-1)+,v(i,j-1)-:"<< discretization_->v(i,j) << ", " << discretization_->v(i+1,j) << ", " << discretization_->v(i+1,j-1) << ", " << discretization_->v(i,j-1) << std::endl;
}

/**
 * compute the boundary conditions for the top wall
 */
void Computation::topWallBC(int i, int j)
{

    // std::cout << "top wall ( "<<i<<", "<<j<<"); ";
    if (updateSurfaceVelocities_)
    {
        // mass balance
        discretization_->v(i,j) = discretization_->v(i,j-1) - discretization_->dy() / discretization_->dx() * (discretization_->u(i,j) - discretization_->u(i-1,j));

        // tangential stress
        if (discretization_->markerfield(i-1,j+1) == 0)
        {
            discretization_->u(i-1,j+1) = discretization_->u(i-1,j) - discretization_->dy() / discretization_->dx() * (discretization_->v(i,j) - discretization_->v(i-1,j));
        }
    }

    // normal stress
    if (updateSurfacePs_)
    {
        discretization_->p(i,j) = 2.0 / settings_.re * (discretization_->v(i,j) - discretization_->v(i,j-1)) / discretization_->dy();

        //temperature
        discretization_->t(i,j) = 0.0;
    }
}

/**
 * compute the boundary conditions for the horizontal pipe
 */
void Computation::horizontalPipeBC(int i, int j)
{

    // std::cout << "horizontal pipe ( "<<i<<", "<<j<<"); ";
    if (updateSurfaceVelocities_)
    {
        // external forces
        discretization_->v(i,j) = discretization_->v(i,j) + dt_ * settings_.g[1];
        discretization_->v(i,j-1) = discretization_->v(i,j-1) + dt_ * settings_.g[1];

        // tangential stress
        if (discretization_->markerfield(i-1,j+1) == 0)
        {
            discretization_->u(i-1,j+1) = discretization_->u(i-1,j) - discretization_->dy() / discretization_->dx() * (discretization_->v(i,j) - discretization_->v(i-1,j));
        }
        if (discretization_->markerfield(i-1,j-1) == 0)
        {
            discretization_->u(i-1,j-1) = discretization_->u(i-1,j) + discretization_->dy() / discretization_->dx() * (discretization_->v(i,j-1) - discretization_->v(i-1,j-1));
        }
    }
                                
    // normal stress
    if (updateSurfacePs_)
    {
        discretization_->p(i,j) = 2.0 / settings_.re * (discretization_->v(i,j) - discretization_->v(i-1,j)) / discretization_->dy();

        //temperature
        discretization_->t(i,j) = 0.0;
    }
}

/**
 * compute the boundary conditions for the top left corner
 */
void Computation::topLeftCornerBC(int i, int j)
{

    // // std::cout << "top left corner ( "<<i<<", "<<j<<"); ";
    if (updateSurfaceVelocities_)
    {
        // mass balance + tangential stress
        discretization_->v(i,j) = discretization_->v(i,j-1);
        discretization_->u(i-1,j) = discretization_->u(i,j);

        if (discretization_->markerfield(i-1,j+1) == 0)
        {
            discretization_->v(i-1,j) = discretization_->v(i,j);
            discretization_->u(i-1,j+1) = discretization_->u(i-1,j);
        }

        // tangential stress
        if (discretization_->markerfield(i-1,j-1) == 0)
        {
            discretization_->v(i-1,j-1) = discretization_->v(i,j-1) + discretization_->dx() / discretization_->dy() * (discretization_->u(i-1,j) - discretization_->u(i-1,j-1));
        }
    }

    // normal stress
    if (updateSurfacePs_)
    {
        discretization_->p(i,j) = - 1.0 / (2.0 * settings_.re) * ((discretization_->u(i,j) + discretization_->u(i-1,j) - discretization_->u(i,j-1) - discretization_->u(i-1,j-1)) / (discretization_->dy()) +
                                                                (discretization_->v(i+1,j) + discretization_->v(i+1,j-1) - discretization_->v(i,j) - discretization_->v(i,j-1)) / (discretization_->dx()));
        
        //temperature
        discretization_->t(i,j) = 0.0;
    }
}

/**
 * compute the boundary conditions for the tip from right
 */
void Computation::tipFromRightBC(int i, int j)
{

    // std::cout << "tip from right ( "<<i<<", "<<j<<"); ";
    if (updateSurfaceVelocities_)
    {
        // mass balance
        discretization_->u(i-1,j) = discretization_->u(i,j) + discretization_->dx() / discretization_->dy() * (discretization_->v(i,j) - discretization_->v(i,j-1));

        // external forces
        discretization_->v(i,j) = discretization_->v(i,j) + dt_ * settings_.g[1];
        discretization_->v(i,j-1) = discretization_->v(i,j-1) + dt_ * settings_.g[1];
                                    
        // mass balance + tangential stress
        if (discretization_->markerfield(i-1,j+1) == 0)
        {
            discretization_->u(i-1,j+1) = discretization_->u(i-1,j);
            discretization_->v(i-1,j) = discretization_->v(i,j);
        }
        if (discretization_->markerfield(i-1,j-1) == 0)
        {
            discretization_->v(i-1,j-1) = discretization_->v(i,j-1);
            discretization_->u(i-1,j-1) = discretization_->u(i-1,j);
        }
    }

    // normal stress
    if (updateSurfacePs_)
    {
        discretization_->p(i,j) = 0.0;

        //temperature
        discretization_->t(i,j) = 0.0;
    }
}

/**
 * compute the boundary conditions for the right wall
 */
void Computation::rightWallBC(int i, int j)
{

    // std::cout << "right wall ( "<<i<<", "<<j<<"); ";
    if (updateSurfaceVelocities_)
    {
        // mass balance
        discretization_->u(i,j) = discretization_->u(i-1,j) - discretization_->dx() / discretization_->dy() * (discretization_->v(i,j) - discretization_->v(i,j-1));

        // tangential stress
        if (discretization_->markerfield(i+1,j-1) == 0)
        {
            discretization_->v(i+1,j-1) = discretization_->v(i,j-1) - discretization_->dx() / discretization_->dy() * (discretization_->u(i,j) - discretization_->u(i,j-1));
        }
    }
    if (updateSurfacePs_)
    {
        // normal stress
        discretization_->p(i,j) = 2.0 / settings_.re * (discretization_->u(i,j) - discretization_->u(i-1,j)) / discretization_->dx();

        //temperature
        discretization_->t(i,j) = 0.0;
    }
}

/**
 * compute the boundary conditions for the bottom right corner
 */
void Computation::bottomRightCornerBC(int i, int j)
{

    // std::cout << "bottom right corner ( "<<i<<", "<<j<<"); ";
    if (updateSurfaceVelocities_)
    {
        // mass balance + tangential stress
        discretization_->u(i,j) = discretization_->u(i-1,j);
        discretization_->v(i,j-1) = discretization_->v(i,j);

        if (discretization_->markerfield(i+1,j-1) == 0)
        {
            discretization_->u(i,j-1) = discretization_->u(i,j);
            discretization_->v(i+1,j-1) = discretization_->v(i,j-1);
        }

        // tangential stress
        if (discretization_->markerfield(i-1,j-1) == 0)
        {
            discretization_->u(i-1,j-1) = discretization_->u(i-1,j) + discretization_->dy() / discretization_->dx() * (discretization_->v(i,j-1) - discretization_->v(i-1,j-1));
        }
    }
                                
    // normal stress
    if (updateSurfacePs_)
    {
        discretization_->p(i,j) = - 1.0 / (2.0 * settings_.re) * ((discretization_->u(i,j+1) + discretization_->u(i-1,j+1) - discretization_->u(i,j) - discretization_->u(i-1,j)) / (discretization_->dy()) + (discretization_->v(i,j) + discretization_->v(i,j-1) - discretization_->v(i-1,j) - discretization_->v(i-1,j-1)) / (discretization_->dx()));

        //temperature
        discretization_->t(i,j) = 0.0;
    }
}


/**
 * compute the boundary conditions for the vertical pipe
 */
void Computation::verticalPipeBC(int i, int j)
{

    // std::cout << "vertical pipe ( "<<i<<", "<<j<<"); ";
    if (updateSurfaceVelocities_)
    {
        // external forces
        discretization_->u(i,j) = discretization_->u(i,j) + dt_ * settings_.g[0];
        discretization_->u(i-1,j) = discretization_->u(i-1,j) + dt_ * settings_.g[0];

        // tangential stress
        if (discretization_->markerfield(i-1,j-1) == 0)
        {
            discretization_->v(i-1, j-1) = discretization_->v(i,j-1) + discretization_->dx() / discretization_->dy() * (discretization_->u(i-1,j) - discretization_->u(i-1,j-1));
        }
        if (discretization_->markerfield(i+1,j-1) == 0)
        {
            discretization_->v(i+1,j-1) = discretization_->v(i,j-1) - discretization_->dx() / discretization_->dy() * (discretization_->u(i,j) - discretization_->u(i,j-1));
        }
    }

    // normal stress
    //TODO: check, should maybe not be zero
    if (updateSurfacePs_)
    {
        discretization_->p(i,j) = 2.0 / settings_.re * (discretization_->u(i,j) - discretization_->u(i-1,j)) / discretization_->dx();

        //temperature
        discretization_->t(i,j) = 0.0;
    }
}

/**
 * compute the boundary conditions for the tip from top
 */
void Computation::tipFromTopBC(int i, int j)
{

    // std::cout << "tip from top ( "<<i<<", "<<j<<"); ";
    if (updateSurfaceVelocities_)
    {
        // mass balance
        discretization_->v(i,j-1) = discretization_->v(i,j) + discretization_->dy() / discretization_->dx() * (discretization_->u(i,j) - discretization_->u(i-1,j));

        // external forces
        discretization_->u(i,j) = discretization_->u(i,j) + dt_ * settings_.g[0];
        discretization_->u(i-1,j) = discretization_->u(i-1,j) + dt_ * settings_.g[0];

        // mass balance + tangential stress
        if (discretization_->markerfield(i-1,j-1) == 0)
        {
            discretization_->u(i-1,j-1) = discretization_->u(i-1,j);
            discretization_->v(i-1,j-1) = discretization_->v(i,j-1);
        }
        if (discretization_->markerfield(i+1,j-1) == 0)
        {
            discretization_->u(i,j-1) = discretization_->u(i,j);
            discretization_->v(i+1,j-1) = discretization_->v(i,j-1);
        }
    }

    if (updateSurfacePs_)
    {
        // normal stress
        discretization_->p(i,j) = 0.0;

        //temperature
        discretization_->t(i,j) = 0.0;
    }
}

/**
 * compute the boundary conditions for the top right corner
 */
void Computation::topRightCornerBC(int i, int j)
{

    // std::cout << "top right corner ( "<<i<<", "<<j<<"); ";
    if (updateSurfaceVelocities_)
    {
        // mass balance + tangential stress
        discretization_->u(i,j) = discretization_->u(i-1,j);
        discretization_->v(i,j) = discretization_->v(i,j-1);

        if (discretization_->markerfield(i+1,j+1) == 0)
        {
            discretization_->u(i,j+1) = discretization_->u(i,j);
            discretization_->v(i+1,j) = discretization_->v(i,j);
        }

        // tangential stress
        if (discretization_->markerfield(i-1,j+1) == 0)
        {
            discretization_->u(i-1,j+1) = discretization_->u(i-1,j) - discretization_->dy() / discretization_->dx() * (discretization_->v(i,j) - discretization_->v(i-1,j));
        }
        if (discretization_->markerfield(i+1,j-1) == 0)
        {
            discretization_->v(i+1,j-1) = discretization_->v(i,j-1) - discretization_->dx() / discretization_->dy() * (discretization_->u(i,j) - discretization_->u(i,j-1));
        }
    }

    
    // normal stress
    if (updateSurfacePs_)
    {
        discretization_->p(i,j) = 1.0 / (2.0 * settings_.re) * ((discretization_->u(i,j) + discretization_->u(i-1,j) - discretization_->u(i,j-1) - discretization_->u(i-1,j-1)) / (discretization_->dy()) + (discretization_->v(i,j) + discretization_->v(i,j-1) - discretization_->v(i-1,j) - discretization_->v(i-1,j-1)) / (discretization_->dx()));

        //temperature
        discretization_->t(i,j) = 0.0;
    }
}

/**
 * compute the boundary conditions for the tip from left
*/
void Computation::tipFromLeftBC(int i, int j)
{

    // std::cout << "tip from left ( "<<i<<", "<<j<<"); ";
    if (updateSurfaceVelocities_)
    {
        // mass balance
        discretization_->u(i,j) = discretization_->u(i-1,j) - discretization_->dx() / discretization_->dy() * (discretization_->v(i,j) - discretization_->v(i,j-1));

        // external forces
        discretization_->v(i,j) = discretization_->v(i,j) + dt_ * settings_.g[1];
        discretization_->v(i,j-1) = discretization_->v(i,j-1) + dt_ * settings_.g[1];

        // mass balance + tangential stress
        if (discretization_->markerfield(i+1,j+1) == 0)
        {
            discretization_->u(i,j+1) = discretization_->u(i,j);
            discretization_->v(i+1,j) = discretization_->v(i,j);
        }
        if (discretization_->markerfield(i+1,j-1) == 0)
        {
            discretization_->v(i+1,j-1) = discretization_->v(i,j-1);
            discretization_->u(i,j-1) = discretization_->u(i,j);
        }

        // tangential stress
        if (discretization_->markerfield(i-1,j+1) == 0)
        {
            discretization_->u(i-1,j+1) = discretization_->u(i-1,j) - discretization_->dy() / discretization_->dx() * (discretization_->v(i,j) - discretization_->v(i-1,j));
        }
        if (discretization_->markerfield(i-1,j-1) == 0)
        {
            discretization_->u(i-1,j-1) = discretization_->u(i-1,j) + discretization_->dy() / discretization_->dx() * (discretization_->v(i,j-1) - discretization_->v(i-1,j-1));
        }
    
    }
    // pressure is set to 0
    if (updateSurfacePs_)
    {
        discretization_->p(i,j) = 0.0;

        //temperature
        discretization_->t(i,j) = 0.0;
    }
}

/*
 * compute the boundary conditions for the tip from bottom
*/
void Computation::tipFromBottomBC(int i, int j)
{

    //std::cout << "tip from bottom ( "<<i<<", "<<j<<"); ";
    if (updateSurfaceVelocities_)
    {
        // mass balance
        discretization_->v(i,j) = discretization_->v(i,j-1) - discretization_->dy() / discretization_->dx() * (discretization_->u(i,j) - discretization_->u(i-1,j));
        
        // external forces
        discretization_->u(i,j) = discretization_->u(i,j) + dt_ * settings_.g[0];
        discretization_->u(i-1,j) = discretization_->u(i-1,j) + dt_ * settings_.g[0];

        // mass balance + tangential stress
        if (discretization_->markerfield(i+1,j+1) == 0)
        {
            discretization_->v(i+1,j) = discretization_->v(i,j);
            discretization_->u(i,j+1) = discretization_->u(i,j);
        }
        if (discretization_->markerfield(i-1,j+1) == 0)
        {
            discretization_->u(i-1,j+1) = discretization_->u(i-1,j);
            discretization_->v(i-1,j) = discretization_->v(i,j);
        }

        // tangential stress
        if (discretization_->markerfield(i-1,j-1) == 0)
        {
            discretization_->v(i-1,j-1) = discretization_->v(i,j-1) + discretization_->dx() / discretization_->dy() * (discretization_->u(i-1,j) - discretization_->u(i-1,j-1));
        }
        if (discretization_->markerfield(i+1,j-1) == 0)
        {
            discretization_->v(i+1,j-1) = discretization_->v(i,j-1) - discretization_->dx() / discretization_->dy() * (discretization_->u(i,j) - discretization_->u(i,j-1));
        }
    }

    // pressure is set to 0
    if (updateSurfacePs_)
    {
        discretization_->p(i,j) = 0.0;

        //temperature
        discretization_->t(i,j) = 0.0;
    }
}

/*
 * compute the boundary conditions for the drop
*/
void Computation::dropBC(int i, int j)
{

    //std::cout << "drop ( "<<i<<", "<<j<<"); ";
    if (updateSurfaceVelocities_)
    {
        // external forces
        discretization_->u(i,j) = discretization_->u(i,j) + dt_ * settings_.g[0];
        discretization_->u(i-1,j) = discretization_->u(i-1,j) + dt_ * settings_.g[0];
        discretization_->v(i,j) = discretization_->v(i,j) + dt_ * settings_.g[1];
        discretization_->v(i,j-1) = discretization_->v(i,j-1) + dt_ * settings_.g[1];
                                    


        // mass balance + tangential stress
        if (discretization_->markerfield(i+1,j+1) == 0)
        {
            discretization_->u(i,j+1) = discretization_->u(i,j);
            discretization_->v(i+1,j) = discretization_->v(i,j);
        }
        if (discretization_->markerfield(i+1,j-1) == 0)
        {
            discretization_->u(i,j-1) = discretization_->u(i,j);
            discretization_->v(i+1,j-1) = discretization_->v(i,j-1);
        }
        if (discretization_->markerfield(i-1,j-1) == 0)
        {
            discretization_->v(i-1,j-1) = discretization_->v(i,j-1);
            discretization_->u(i-1,j-1) = discretization_->u(i-1,j);
        }
        if (discretization_->markerfield(i-1,j+1) == 0)
        {
            discretization_->u(i-1,j+1) = discretization_->u(i-1,j);
            discretization_->v(i-1,j) = discretization_->v(i,j);
        }
    }
    
    // pressure is set to 0
    if (updateSurfacePs_)
    {
        discretization_->p(i,j) = 0.0;

        //temperature
        discretization_->t(i,j) = 0.0;
    }
}

void Computation::resetEmptyEdges()
{
    for (int i = 1; i < discretization_->rhsSize()[0] - 1; i++)
    {
        for (int j = 1; j < discretization_->rhsSize()[1] - 1; j++)
        {
            if (discretization_->markerfield(i, j) == 0)
            {
                if (discretization_->markerfield(i + 1, j) == 0)
                {
                    discretization_->u(i, j) = 0;
                }
                if (discretization_->markerfield(i, j + 1) == 0)
                {
                    discretization_->v(i, j) = 0;
                }
            }
        }
    }
}