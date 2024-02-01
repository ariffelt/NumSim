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
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
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

    // initialize output writers
    outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
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

    while (t < settings_.endTime)
    {
        applyBoundaryValues(); // set boundary values for u, v, F and G

        computeTimeStepWidth();
        if (t + dt_ > settings_.endTime) // decrease time step width in last time step, s.t. the end time will be reached exactly
        {
            dt_ = settings_.endTime - t;
        }
        t += dt_;

        // only compute preliminary velocities and rhs and solve pressure eq. on inner fluid cells
        
        computePreliminaryVelocities(); // compute preliminary velocities, F and G

        computeRightHandSide(); // compute rhs of the Poisson equation for the pressure

        computePressure(); // solve the Poisson equation for the pressure

        computeVelocities(); // compute the new velocities, u,v, from the preliminary velocities, F,G and the pressure, p

        freeflowBC(); // apply free flow boundary conditions

        computeParticleVelocities();

        updateMarkerField();

        freeflowBC();

        outputWriterParaview_->writeFile(t); // output simulation results
        outputWriterText_->writeFile(t);

        if (t >= timestepInterval)
        {
            std::cout << "t = " << t << ", dt = " << dt_ << std::endl; // print time and time step width to the console
            timestepInterval += settings_.endTime / numberOfPrints;
        }
    }
}

/**
 * Compute the time step width dt from maximum velocities.
 */
void Computation::computeTimeStepWidth()
{
    // compute maximal time step width from diffusion
    double dt_diffusion = (settings_.re / 2) * (pow(discretization_->dx(), 2) * pow(discretization_->dy(), 2)) / (pow(discretization_->dx(), 2) + pow(discretization_->dy(), 2));

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
    dt_ = std::min(settings_.tau * std::min(dt_diffusion, dt_convection), settings_.maximumDt);
}

/**
 * Set velocity boundary values for u, v, F and G
 */
void Computation::applyBoundaryValues()
{
    // set u boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        if (discretization_->markerfield(i,discretization_->uJBegin()) == 1)
        {
            // u boundary values bottom, assuming inhomogenous Dirichlet conditions
            discretization_->u(i, discretization_->uJBegin()) = 2.0 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin() + 1);
        }
        if (discretization_->markerfield(i,discretization_->uJEnd()) == 1)
        {
            // u boundary values top, assuming inhomogenous Dirichlet conditions
            discretization_->u(i, discretization_->uJEnd()) = 2.0 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() - 1);
        }
    }

    // set u boundary values for left and right side
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++)
    {
        if (discretization_->markerfield(discretization_->uIBegin(),j) == 1)
        {
            // u boundary values left, assuming inhomogenous Dirichlet conditions
            discretization_->u(discretization_->uIBegin(), j) = settings_.dirichletBcLeft[0];
        }
        if (discretization_->markerfield(discretization_->uIEnd() - 1,j) == 1)
        {
            // u boundary values right, assuming inhomogenous Dirichlet conditions
            discretization_->u(discretization_->uIEnd() - 1, j) = settings_.dirichletBcRight[0];
        }
    }

    // set v boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++)
    {
        if (discretization_->markerfield(i,discretization_->vJBegin()) == 1)
        {
            // v boundary values bottom, assuming inhomogenous Dirichlet conditions
            discretization_->v(i, discretization_->vJBegin()) = settings_.dirichletBcBottom[1];
        }
        if (discretization_->markerfield(i,discretization_->vJEnd() - 1) == 1)
        {
            // v boundary values top, assuming inhomogenous Dirichlet conditions
            discretization_->v(i, discretization_->vJEnd() - 1) = settings_.dirichletBcTop[1];
        }
    }

    // set v boundary values for left and right side
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        if (discretization_->markerfield(discretization_->vIBegin(),j) == 1)
        {
            // v boundary values left, assuming inhomogenous Dirichlet conditions
            discretization_->v(discretization_->vIBegin(), j) = 2 * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin() + 1, j);
        }
        if (discretization_->markerfield(discretization_->vIEnd(),j) == 1)
        {
            // v boundary values right, assuming inhomogenous Dirichlet conditions
            discretization_->v(discretization_->vIEnd(), j) = 2 * settings_.dirichletBcRight[1] - discretization_->v(discretization_->vIEnd() - 1, j);
        }
    }

    // set F boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        if (discretization_->markerfield(i,discretization_->uJBegin()) == 1)
        {
            // F boundary values bottom
            discretization_->f(i, discretization_->uJBegin()) = discretization_->u(i, discretization_->uJBegin());
        }
        if (discretization_->markerfield(i,discretization_->uJEnd()) == 1)
        {
            // F boundary values top
            discretization_->f(i, discretization_->uJEnd()) = discretization_->u(i, discretization_->uJEnd());
        }
    }

    // set F boundary values for left and right side
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++)
    {
        if (discretization_->markerfield(discretization_->uIBegin(),j) == 1)
        {
            // F boundary values left
            discretization_->f(discretization_->uIBegin(), j) = discretization_->u(discretization_->uIBegin(), j);
        }
        if (discretization_->markerfield(discretization_->uIEnd() - 1,j) == 1)
        {
            // F boundary values right
            discretization_->f(discretization_->uIEnd() - 1, j) = discretization_->u(discretization_->uIEnd() - 1, j);
        }
    }

    // set G boundary values for bottom and top first, as for corner cases, the left and right border should be used
    for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++)
    {
        if (discretization_->markerfield(i,discretization_->vJBegin()) == 1)
        {
            // G boundary values bottom
            discretization_->g(i, discretization_->vJBegin()) = discretization_->v(i, discretization_->vJBegin());
        }
        if (discretization_->markerfield(i,discretization_->vJEnd() - 1) == 1)
        {
            // G boundary values top
            discretization_->g(i, discretization_->vJEnd() - 1) = discretization_->v(i, discretization_->vJEnd() - 1);
        }
    }

    // set G boundary values for left and right side
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        if (discretization_->markerfield(discretization_->vIBegin(),j) == 1)
        {
            // G boundary values left
            discretization_->g(discretization_->vIBegin(), j) = discretization_->v(discretization_->vIBegin(), j);
        }
        if (discretization_->markerfield(discretization_->vIEnd(),j) == 1)
        {
            // G boundary values right
            discretization_->g(discretization_->vIEnd(), j) = discretization_->v(discretization_->vIEnd(), j);
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
            if (isInnerFluidCell(i,j))
            {
                double diffusionTerms = (1 / settings_.re) * (discretization_->computeD2uDx2(i, j) + discretization_->computeD2uDy2(i, j));
                double convectionTerms = discretization_->computeDu2Dx(i, j) + discretization_->computeDuvDy(i, j);

                discretization_->f(i, j) = discretization_->u(i, j) + dt_ * (diffusionTerms - convectionTerms + settings_.g[0]);
            }
        }
    }

    // compute preliminary G
    for (int i = discretization_->vIBegin() + 1; i < discretization_->vIEnd(); i++)
    {
        for (int j = discretization_->vJBegin() + 1; j < discretization_->vJEnd() - 1; j++)
        {
            if (isInnerFluidCell(i,j))
            {
                double diffusionTerms = (1 / settings_.re) * (discretization_->computeD2vDx2(i, j) + discretization_->computeD2vDy2(i, j));
                double convectionTerms = discretization_->computeDuvDx(i, j) + discretization_->computeDv2Dy(i, j);

                discretization_->g(i, j) = discretization_->v(i, j) + dt_ * (diffusionTerms - convectionTerms + settings_.g[1]);
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
            if (isInnerFluidCell(i,j))
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
            if (isInnerFluidCell(i, j))
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
            if (isInnerFluidCell(i, j))
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
    int numParticles = 10;
    // std::vector<Particle> particles_(numParticles);
    for (int i = 0; i < numParticles; ++i)
    {
        particles_[i] = Particle(settings_.inlet[0], settings_.inlet[1]);
    }
}

/**
 * Compute the new particle velocities.
 * And move particles according to these.
 */
void Computation::computeParticleVelocities()
{
    double dx = discretization_->dx();
    double dy = discretization_->dy();

    // interpolate velocities to the particle positions (do not coincide with velocity grid points)
    for (int k = 0; k < particles_.size(); k++)
    {
        // compute particle velocity in x direction

        // index of upper right corner
        int iUpperRight = int(particles_[k].x / dx + 1);
        int jUpperRight = int((particles_[k].y + 1 / 2) / dy + 1);

        // position of 4 neighbouring grid points with values u
        double x1 = (iUpperRight - 1) * dx;
        double x2 = iUpperRight * dx;
        double y1 = (jUpperRight - 1) * dy - dy / 2;
        double y2 = jUpperRight * dy - dy / 2;

        // bilinear interpolation
        double u = 1 / (dx * dy) * ((x2 - particles_[k].x) * (y2 - particles_[k].y) * discretization_->u(iUpperRight - 1, jUpperRight - 1) 
                                    + (particles_[k].x - x1) * (y2 - particles_[k].y) * discretization_->u(iUpperRight, jUpperRight - 1) 
                                    + (x2 - particles_[k].x) * (particles_[k].y - y1) * discretization_->u(iUpperRight - 1, jUpperRight) 
                                    + (particles_[k].x - x1) * (particles_[k].y - y1) * discretization_->u(iUpperRight, jUpperRight));

        // compute particle velocity in y direction
        // index of upper right corner
        iUpperRight = int((particles_[k].x + 1 / 2) / dx + 1);
        jUpperRight = int(particles_[k].y / dy + 1);

        // position of 4 neighbouring grid points with values v
        x1 = (iUpperRight - 1) * dx - dx / 2;
        x2 = iUpperRight * dx - dx / 2;
        y1 = (jUpperRight - 1) * dy;
        y2 = jUpperRight * dy;

        // bilinear interpolation
        double v = 1 / (dx * dy) * ((x2 - particles_[k].x) * (y2 - particles_[k].y) * discretization_->v(iUpperRight - 1, jUpperRight - 1) 
                                    + (particles_[k].x - x1) * (y2 - particles_[k].y) * discretization_->v(iUpperRight, jUpperRight - 1) 
                                    + (x2 - particles_[k].x) * (particles_[k].y - y1) * discretization_->v(iUpperRight - 1, jUpperRight) 
                                    + (particles_[k].x - x1) * (particles_[k].y - y1) * discretization_->v(iUpperRight, jUpperRight));


        // move particle
        particles_[k].x += dt_ * u;
        particles_[k].y += dt_ * v;
    }
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
            discretization_->markerfield(i, j) = 0; // assume cell is empty
            for (int k = 0; k < particles_.size(); k++)
            {
                if (particles_[k].x == i && particles_[k].y == j)
                {
                    discretization_->markerfield(i, j) = 1; // fluid cell
                    break;
                }
            }
        }
    }
}

/**
 * Check if cell has any neighbouring empty cell, if it has none then it is an inner fluid cell.
 */
bool Computation::isInnerFluidCell(int i, int j)
{
    // check if cell has any neighbouring empty cell, if it has none then it is an inner fluid cell
    // watch out for boundary cells
    // return true if the current cell is a fluid cell and does not have any empty cell or any boundaries as neighbour
    bool isInnerFluidCell = false;
    if (discretization_->markerfield(i, j) == 1)
    {
        if (discretization_->markerfield(i - 1, j) == 1 && discretization_->markerfield(i + 1, j) == 1 && discretization_->markerfield(i, j - 1) == 1 && discretization_->markerfield(i, j + 1) == 1)
        {
            isInnerFluidCell = true;
        }
    }
    return isInnerFluidCell;
}

/**
 * Apply free flow boundary conditions.
 */
void Computation::freeflowBC()
{
    // for loops over all p cells
    // check if cell is fluid cell
    // check if it has 
    for (int i = 1; i < discretization_->pSize()[0] - 1; i++)
    {
        for (int j = 1; j < discretization_->pSize()[1] - 1; j++)
        {
            if (!isInnerFluidCell(i,j))
            {
                // check if the right, top, left or bottom cell is a fluid cell
                if (discretization_->markerfield(i + 1, j) == 1)
                {
                    if (discretization_->markerfield(i, j + 1) == 1)
                    {
                        if (discretization_->markerfield(i - 1, j) == 1)
                        {
                            // cell is surrounded by 3 fluid cells (bottom wall)
                            break;
                        }
                        else
                        {
                            if (discretization_->markerfield(i, j - 1) == 1)
                            {
                                // cell is surrounded by 3 fluid cells (right wall)
                                break;
                            }
                            else
                            {
                                // cell is surrounded by 2 fluid cells (right, bottom corner)
                                break;
                            }
                        }
                    }
                    else
                    {
                        if (discretization_->markerfield(i - 1, j) == 1)
                        {
                            if (discretization_->markerfield(i, j - 1) == 1)
                            {
                                // cell is surrounded by 3 fluid cells (top wall)
                                break;
                            }
                            else
                            {
                                // cell is surrounded by 2 fluid cells (horizontal pipe)
                                break;
                            }
                        }
                        else
                        {
                            if (discretization_->markerfield(i, j - 1) == 1)
                            {
                                // cell is surrounded by 2 fluid cells (top, right corner)
                                break;
                            }
                            else
                            {
                                // cell is surrounded by 1 fluid cell (tip from right)
                                break;
                            }
                        }
                    }
                }
                else 
                {
                    if (discretization_->markerfield(i, j + 1) == 1)
                    {
                        if (discretization_->markerfield(i - 1, j) == 1)
                        {
                            if (discretization_->markerfield(i, j - 1) == 1)
                            {
                                // cell is surrounded by 3 fluid cells (right wall)
                                break;
                            }
                            else
                            {
                                // cell is surrounded by 2 fluid cells (bottom, left corner)
                                break;
                            }
                        }
                        else
                        {
                            if (discretization_->markerfield(i, j - 1) == 1)
                            {
                                // cell is surrounded by 2 fluid cells (vertikal pipe)
                                break;
                            }
                            else
                            {
                                // cell is surrounded by 1 fluid cell (tip from top)
                                break;
                            }
                        }
                    }
                    else
                    {
                        if (discretization_->markerfield(i - 1, j) == 1)
                        {
                            if (discretization_->markerfield(i, j - 1) == 1)
                            {
                                // cell is surrounded by 2 fluid cells (top, right corner)
                                break;
                            }
                            else
                            {
                                // cell is surrounded by 1 fluid cell (tip from left)
                                break;
                            }
                        }
                        else
                        {
                            if (discretization_->markerfield(i, j - 1) == 1)
                            {
                                // cell is surrounded by 1 fluid cell (tip from bottom)
                                break;
                            }
                            else
                            {
                                // cell is surrounded by 0 fluid cells (drop)
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}
