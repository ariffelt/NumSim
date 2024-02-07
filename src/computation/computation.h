#pragma once

#include <memory>

#include "settings/settings.h"
#include "discretization/1_discretization.h"
#include "pressure_solver/pressure_solver.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"
// #include "freeflow/virtualparticle.h"

class Computation
{
public:

    //! constructor, initialize settings, discretization, pressure solver and output writers
    void initialize(int argc, char *argv[]);

    //! run the whole simulation until tend
    void runSimulation();

    //! test the implementation of boundary conditions
    void testBC();
    //! set the boundary markers for the cells
    void setBoundaryMarkers();

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

    void generateVirtualParticles();

    void computeParticleVelocities();

    void generateDam(int noParticles);

    void generateFountain(int noParticles, int noParticlesFountain);

    void setFountainVelocity();

    void setFountainTemperature();

    void generateFull(int noParticles);

    void generateBox(int noParticles);

    void generateDrop(int noParticles);

    void generateBigDrop(int noParticles);

    void generateBar(int noParticles);

    void generateDropInWater(int noParticles);  

    // void updateCellTypes();

    void updateMarkerField();

    void freeflowBC();

    void bottomWallBC(int i, int j);

    void leftWallBC(int i, int j);

    void bottomLeftCornerBC(int i, int j);

    void topWallBC(int i, int j);

    void horizontalPipeBC(int i, int j);

    void topLeftCornerBC(int i, int j);

    void tipFromRightBC(int i, int j);  

    void rightWallBC(int i, int j);

    void bottomRightCornerBC(int i, int j); 

    void verticalPipeBC(int i, int j);

    void tipFromTopBC(int i, int j);

    void topRightCornerBC(int i, int j);

    void tipFromLeftBC(int i, int j);

    void tipFromBottomBC(int i, int j);

    void dropBC(int i, int j);

    void printParticles();

    void resetEmptyEdges();

    //! compute new temperature t from the old temperature t_old and the velocities u and v
    void computeTemperature();

    Settings settings_;

    std::shared_ptr<Discretization> discretization_;

    std::unique_ptr<PressureSolver> pressureSolver_;

    std::unique_ptr<OutputWriterParaview> outputWriterParaview_;

    std::unique_ptr<OutputWriterText> outputWriterText_;

    std::array<double, 2> meshWidth_;

    double dt_;

    std::vector<double> particlesX_;
    std::vector<double> particlesY_;

    bool updateSurfacePs_ = false;
    bool updateSurfaceVelocities_ = false;
};