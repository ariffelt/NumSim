#include "computation/computation.h"

// TODO: constructor
Computation::Computation()
{
    // nothing to do
}

void Computation::initialize(Settings settings)
{
    if (settings_.pressureSolver == "GaussSeidel")
    {
        pressureSolver_ = std::make_unique<Gauss_Seidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
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
}
/**
 * Run the whole simulation until tend. 
 */
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
void Computation::computeTimeStepWidth()
{
}
void Computation::applyBoundaryValues()
{
}
void Computation::computePreliminaryVelocities()
{
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
