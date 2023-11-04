#include <fstream> // for file operations
#include <iostream> // for cout

void loadFromFile(std::string filename)
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

        // remove whitespace at end of value
        while (value[value.size() - 1] == ' ' || value[value.size() - 1] == '\t')
        {
            value.erase(value.size() - 1, 1);
        }

        // parse actual value and set corresponding parameter
        if (parameterName == "endTime")
        {
            // set parameter endTime
            std::cout << "parameter endTime = " << atof(value.c_str()) << std::endl;
        }
        else if (parameterName == "re")
        {
            // set parameter re
            std::cout << "parameter re = " << atoi(value.c_str()) << std::endl;
        }
        else
        {
            // unknown parameter
            std::cout << "unknown parameter \"" << parameterName << "\"." << std::endl;
        }
    }

    #ifndef NDEBUG
        // only run this code in debug target
        std::cout << "lots of inefficient but informative output . . .";
    #endif
}