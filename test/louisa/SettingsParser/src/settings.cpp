#include "settings.h"
#include <fstream>
#include <iomanip>

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

    // print line
    std::cout << "line " << lineNo << ": " << line << std::endl;

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
    //missing: also removing the letters behind the comment

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
        useDonorCell = (value.c_str() == "true");
    }
    else if (parameterName == "alpha")
    {
        alpha = atof(value.c_str());
    }
    else if (parameterName == "dirichletBcBottomU")
    {
        dirichletBcBottom[0] = atof(value.c_str());
    }
    else if (parameterName == "dirichletBcBottomV")
    {
        dirichletBcBottom[1] = atof(value.c_str());
    }
    else if (parameterName == "dirichletBcTopU")
    {
        dirichletBcTop[0] = atof(value.c_str());
    }
    else if (parameterName == "dirichletBcTopV")
    {
        dirichletBcTop[1] = atof(value.c_str());
    }
    else if (parameterName == "dirichletBcLeftU")
    {
        dirichletBcLeft[0] = atof(value.c_str());
    }
    else if (parameterName == "dirichletBcLeftV")
    {
        dirichletBcLeft[1] = atof(value.c_str());
    }
    else if (parameterName == "dirichletBcRightU")
    {
        dirichletBcRight[0] = atof(value.c_str());
    }
    else if (parameterName == "dirichletBcRightV")
    {
        dirichletBcRight[1] = atof(value.c_str());
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
        maximumNumberOfIterations = atoi(value.c_str());
    }
    else
    {
        std::cout << "unknown parameter \"" << parameterName << "\"." << std::endl;
    }
  }


}

void Settings::printSettings()
{
  std::cout << "Settings: " << std::endl
    << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl
    << "  endTime: " << endTime << " s, re: " << re << ", g: (" << g[0] << "," << g[1] << "), tau: " << tau << ", maximum dt: " << maximumDt << std::endl
    << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1]  << ")"
    << ", top: ("  << dirichletBcTop[0] << "," << dirichletBcTop[1]  << ")"
    << ", left: ("  << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"
    << ", right: ("  << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << std::endl
    << "  useDonorCell: " << std::boolalpha << useDonorCell << ", alpha: " << alpha << std::endl
    << "  pressureSolver: " << pressureSolver << ", omega: " << omega << ", epsilon: " << epsilon << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
}