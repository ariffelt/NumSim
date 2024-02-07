#include "settings.h"

#include <fstream>
#include <iomanip>

/**
 * Loads the settings from a parameter file
 * @param filename name of the parameter file
*/
void Settings::loadFromFile(std::string filename)
{
    // open file
    std::ifstream file(filename.c_str(), std::ios::in);

    // check if file is open
    if (!file.is_open())
    {
        std::cout << "Could not open parameter file \"" << filename << "\"." << std::endl;
        return;
    }

    // loop over lines of file
    for (int lineNo = 0;; lineNo++)
    {
        // read line
        std::string line;
        getline(file, line);

        // at the end of the file break for loop
        if (file.eof())
            break;

        // remove whitespace at beginning of line (if there is any)
        while (line[0] == ' ' || line[0] == '\t')
        {
            line.erase(0, 1);
        }

        // if first character is a '#', skip line (line[0] == '#')
        if (line[0] == '#')
            continue;

        // if line does not contain a '=' sign, skip line
        if (line.find_first_of("=") == std::string::npos)
            continue;

        // parse parameter name
        std::string parameterName = line.substr(0, line.find_first_of("="));

        // remove trailing spaces from parameterName
        if (parameterName.find_first_of(" \t") != std::string::npos)
        {
            parameterName.erase(parameterName.find_first_of(" \t"));
        }

        // parse value
        std::string value = line.substr(line.find_first_of("=") + 1);

        // remove whitespace at beginning of value (if there is any)
        while (value[0] == ' ' || value[0] == '\t')
        {
            value.erase(0, 1);
        }

        // remove comments at end of value
        if (value.find_first_of("#") != std::string::npos)
        {
            value.erase(value.find_first_of("#"));
        }

        // remove whitespace at end of value
        while (value[value.size() - 1] == ' ' || value[value.size() - 1] == '\t')
        {
            value.erase(value.size() - 1, 1);
        }

        // parse actual value and set corresponding parameter
        if (parameterName == "nCellsX")
        {
            nCells[0] = atoi(value.c_str());
        }
        else if (parameterName == "nCellsY")
        {
            nCells[1] = atoi(value.c_str());
        }
        else if (parameterName == "physicalSizeX")
        {
            physicalSize[0] = atof(value.c_str());
        }
        else if (parameterName == "physicalSizeY")
        {
            physicalSize[1] = atof(value.c_str());
        }
        else if (parameterName == "re")
        {
            re = atof(value.c_str());
        }
        else if (parameterName == "endTime")
        {
            endTime = atof(value.c_str());
        }
        else if (parameterName == "tau")
        {
            tau = atof(value.c_str());
        }
        else if (parameterName == "maximumDt")
        {
            maximumDt = atof(value.c_str());
        }
        else if (parameterName == "gX")
        {
            g[0] = atof(value.c_str());
        }
        else if (parameterName == "gY")
        {
            g[1] = atof(value.c_str());
        }
        else if (parameterName == "useDonorCell")
        {
            useDonorCell = (value == "true" ? true : false);
        }
        else if (parameterName == "alpha")
        {
            alpha = atof(value.c_str());
        }
        else if (parameterName == "dirichletBcBottomU" || parameterName == "dirichletBottomX")
        {
            dirichletBcBottom[0] = atof(value.c_str());
        }
        else if (parameterName == "dirichletBcBottomV" || parameterName == "dirichletBottomY")
        {
            dirichletBcBottom[1] = atof(value.c_str());
        }
        else if (parameterName == "dirichletBcTopU" || parameterName == "dirichletTopX")
        {
            dirichletBcTop[0] = atof(value.c_str());
        }
        else if (parameterName == "dirichletBcTopV" || parameterName == "dirichletTopY")
        {
            dirichletBcTop[1] = atof(value.c_str());
        }
        else if (parameterName == "dirichletBcLeftU" || parameterName == "dirichletLeftX")
        {
            dirichletBcLeft[0] = atof(value.c_str());
        }
        else if (parameterName == "dirichletBcLeftV" || parameterName == "dirichletLeftY")
        {
            dirichletBcLeft[1] = atof(value.c_str());
        }
        else if (parameterName == "dirichletBcRightU" || parameterName == "dirichletRightX")
        {
            dirichletBcRight[0] = atof(value.c_str());
        }
        else if (parameterName == "dirichletBcRightV" || parameterName == "dirichletRightY")
        {
            dirichletBcRight[1] = atof(value.c_str());
        }
        else if (parameterName == "dirichletBcRightT" || parameterName == "dirichletRightT")
        {
            dirichletBcRightT = atof(value.c_str());
        }
        else if (parameterName == "dirichletBcLeftT" || parameterName == "dirichletLeftT")
        {
            dirichletBcLeftT = atof(value.c_str());
        }
        else if (parameterName == "dirichletBcTopT" || parameterName == "dirichletTopT")
        {
            dirichletBcTopT = atof(value.c_str());
        }
        else if (parameterName == "dirichletBcBottomT" || parameterName == "dirichletBottomT")
        {
            dirichletBcBottomT = atof(value.c_str());
        }
        else if (parameterName == "pr")
        {
            pr = atof(value.c_str());
        }
        else if (parameterName == "gamma")
        {
            gamma = atof(value.c_str());
        }
        else if (parameterName == "beta")
        {
            beta = atof(value.c_str());
        }
        else if (parameterName == "pressureSolver")
        {
            pressureSolver = value;
        }
        else if (parameterName == "omega")
        {
            omega = atof(value.c_str());
        }
        else if (parameterName == "epsilon")
        {
            epsilon = atof(value.c_str());
        }
        else if (parameterName == "maximumNumberOfIterations")
        {
            maximumNumberOfIterations = static_cast<int>(atof(value.c_str()));
        }
        else if (parameterName == "inletX")
        {
            inlet[0] = atof(value.c_str());
        }
        else if (parameterName == "inletY")
        {
            inlet[1] = atof(value.c_str());
        }
        else if (parameterName == "particleShape")
        {
            particelShape = value;
        }

        else
        {
            std::cout << "unknown parameter \"" << parameterName << "\"." << std::endl;
        }
    }
}

void Settings::printSettings()
{
    std::cout << "Settings:" << std::endl
              << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl
              << "  re: " << re << ", endTime: " << endTime << ", tau: " << tau << ", maximumDt: " << maximumDt << ", g: (" << g[0] << "," << g[1] << ")" << std::endl
              << "  useDonorCell: " << useDonorCell << ", alpha: " << alpha << std::endl
              << "  dirichletBc: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1] << ")"
              << ", top: (" << dirichletBcTop[0] << "," << dirichletBcTop[1] << ")"
              << ", left: (" << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"
              << ", right: (" << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << std::endl
              << "  pressureSolver: " << pressureSolver << ", omega: " << omega << ", epsilon: " << epsilon << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
}