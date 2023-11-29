#include "computation_parallel.h"
#include "pressure_solver/red_black_sor.h"

/**
 * Constructor for the Computation Parallel class.
 * Initializes the settings, discretization, pressure solver, partitioning and output writers
 * @param argc number of command line arguments
 * @param argv command line arguments
 */
void ComputationParallel::initialize(int argc, char *argv[])
{
    // load and print settings
    settings_ = Settings();
    settings_.loadFromFile(argv[1]);
    settings_.printSettings();

    partitioning_ = std::make_shared<Partitioning>(settings_.nCells);

    meshWidth_[0] = settings_.physicalSize[0] / settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1] / settings_.nCells[1];

    // initialize discretization
    if (settings_.useDonorCell)
    {
        discretization_ = std::make_shared<DonorCell>(partitioning_, meshWidth_, settings_.alpha);
    }
    else
    {
        discretization_ = std::make_shared<CentralDifferences>(partitioning_, meshWidth_);
    }


    if (settings_.pressureSolver == "SOR")
    {
        pressureSolver_ = std::make_unique<RedBlackSOR>(discretization_, settings_.epsilon,
                                                        settings_.maximumNumberOfIterations, settings_.omega, partitioning_);
    }
    else
    {
        std::cout << "Error: Unknown pressure solver for parallel computing: " << settings_.pressureSolver << std::endl;
        exit(1);
    }
    outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discretization_, partitioning_);
    outputWriterParaview_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, partitioning_);
}

/**
 * Run the simulation, starting at 0 until tend.
 * Prints 10 time steps to the console.
 */
void ComputationParallel::runSimulation()
{

}