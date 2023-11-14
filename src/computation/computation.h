#pragma once

#include <memory>

#include "settings/settings.h"
#include "discretization/1_discretization.h"
#include "pressure_solver/pressure_solver.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"

class Computation
{
public:
    //! constructor TODO
    Computation();

    void initialize(Settings settings);
 
    //! run the whole simulation until tend
    void runSimulation();

private:    
    //! compute the time step width dt from maximum velocities
    void computeTimeStepWidth();

    //! set boundary values of u and v to correct values   
    void applyBoundaryValues();

    //! compute the preliminary velocities, F and G    
    void computePreliminaryVelocities();

    //! compute the right hand side of the Poisson equation for the pressure    
    void computeRightHandSide();

    //! solve the Poisson equation for the pressure    
    void computePressure();

    //! compute the new velocities, u,v, from the preliminary velocities, F,G and the pressure, p    
    void computeVelocities();
    
    Settings settings_;
    
    std::shared_ptr<Discretization> discretization_;
    
    std::unique_ptr<PressureSolver> pressureSolver_;
    
    std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
    
    std::unique_ptr<OutputWriterText> outputWriterText_;
    
    std::array<double, 2> meshWidth_;
    
    double 	dt_;
};