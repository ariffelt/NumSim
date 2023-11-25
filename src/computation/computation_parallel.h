#pragma once

#include <mpi.h>
#include "partitioning/partitioning.h"
#include "computation.h"

class ComputationParallel : public Computation
{
public:
    virtual void initialize(std::string filename);

    virtual void runSimulation();

protected:
    //! TODO: which of these functions do we actually need to override?


    //! compute the time step width dt from maximum velocities
    // void computeTimeStepWidth();

    // //! set boundary values of u and v to correct values
    // void applyBoundaryValues();

    // //! compute the preliminary velocities, F and G
    // void computePreliminaryVelocities();

    // //! compute the right hand side of the Poisson equation for the pressure
    // void computeRightHandSide();

    // //! solve the Poisson equation for the pressure
    // void computePressure();

    // //! compute the new velocities, u,v, from the preliminary velocities, F,G and the pressure, p
    // void computeVelocities();
};
