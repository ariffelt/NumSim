#pragma once

#include <mpi.h>

#include "computation.h"

#include "partitioning/partitioning.h"

#include "discretization/2_central_differences.h"
#include "discretization/2_donor_cell.h"

#include "pressure_solver/red_black_sor.h"

#include "output_writer/output_writer_text_parallel.h"
#include "output_writer/output_writer_paraview_parallel.h"


class ComputationParallel : public Computation
{
public:
    virtual void initialize(int argc, char *argv[]);

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
    std::shared_ptr<Partitioning> partitioning_;
};
