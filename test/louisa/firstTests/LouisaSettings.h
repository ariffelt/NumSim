#pragma once // This is an include guard to prevent multiple inclusion

#include <iostream>

//! parse a txt file with settings, each line contains "<parameterName> = <value>"
void loadFromFile(std::string filename);