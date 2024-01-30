#include "settings/settings.h"
#include "computation/computation_parallel.h"

/**
 * Entry point of simulation
*/
int main(int argc, char *argv[])
{
    // if the number of given command line arguments is only 1 (= the program name),
    // print out usage information and exit
    if (argc == 1)
    {
        std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

        return EXIT_FAILURE;
    }

    // read in the first argument
    std::string filename = argv[1];

    // print message
    // std::cout << "Filename: \"" << filename << "\"" << std::endl;

    //MPI_Init(NULL, NULL);

    // create computation object and run simulation
    Computation computation = Computation();
    computation.initialize(argc, argv);
    computation.runSimulation();

    //MPI_Finalize();

    return EXIT_SUCCESS;
}