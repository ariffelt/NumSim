#include "settings/settings.h"
#include "computation/computation.h"

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
    std::cout << "Filename: \"" << filename << "\"" << std::endl;

    Settings settings;
    // load settings from file
    settings.loadFromFile(filename);

    // display all settings on console
    settings.printSettings();

    // create computation object and run simulation
    //TODO: in initialize settings wiederverwenden
    Computation computation = Computation();
    computation.initialize(settings);
    computation.runSimulation();

    return EXIT_SUCCESS;
}
