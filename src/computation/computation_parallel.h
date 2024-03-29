#pragma once

#include <mpi.h>

#include "computation/computation.h"

#include "partitioning/partitioning.h"

#include "discretization/1_discretization.h"
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
    //! compute the time step width dt from maximum velocities
    void computeTimeStepWidthParallel();

    //! communicate the preliminary velocities btw subdomains
    void communicatePreliminaryVelocities();

    //! set boundary values for u, v, F and G and exchange values at borders btw subdomains
    void applyBoundaryValues();

    void applyBoundaryValuesBottom();
    void applyBoundaryValuesTop();
    void applyBoundaryValuesLeft();
    void applyBoundaryValuesRight();

    void exchangeVelocitiesBottom();
    void exchangeVelocitiesTop();
    void exchangeVelocitiesLeft();
    void exchangeVelocitiesRight();

    std::unique_ptr<OutputWriterParaviewParallel> outputWriterParaviewParallel_;

    std::unique_ptr<OutputWriterTextParallel> outputWriterTextParallel_;

    std::shared_ptr<Partitioning> partitioning_;
};
