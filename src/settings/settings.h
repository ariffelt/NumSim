#pragma once

#include <iostream>
#include <array>

/** 
 * All settings that parametrize a simulation run
 */
struct Settings
{
    std::array<int, 2> nCells;          //< number of cells in x and y direction
    std::array<double, 2> physicalSize; //< physical size of the domain
    double re = 1000;                   //< Reynolds number
    double pr = 7.0;                    //< Prandtl number
    double endTime = 10.0;              //< end time of the simulation
    double tau = 0.5;                   //< safety factor for time step width
    double maximumDt = 0.1;             //< maximum time step width

    std::array<double, 2> g{0., 0.}; //< external forces

    bool useDonorCell = false; //< if the donor cell scheme should be used
    double alpha = 0.5;        //< factor for donor-cell scheme
    double gamma = 0.5;        //< factor for donor-cell scheme for temperature

    double beta = 0.04;        //< dimensionless volume expansion coefficient

    std::array<double, 2> dirichletBcBottom; //< prescribed values of u,v at bottom of domain
    std::array<double, 2> dirichletBcTop;    //< prescribed values of u,v at top of domain
    std::array<double, 2> dirichletBcLeft;   //< prescribed values of u,v at left of domain
    std::array<double, 2> dirichletBcRight;  //< prescribed values of u,v at right of domain
    std::array<double, 2> inlet;             //< cell index of waterfountain inlet

    double dirichletBcBottomT; //< prescribed value of T at bottom of domain
    double dirichletBcTopT;    //< prescribed value of T at top of domain
    double dirichletBcLeftT;   //< prescribed value of T at left of domain
    double dirichletBcRightT;  //< prescribed value of T at right of domain

    std::string pressureSolver = "SOR";  //< which pressure solver to use, "GaussSeidel" or "SOR"
    double omega = 1.0;                  //< overrelaxation factor
    double epsilon = 1e-5;               //< tolerance for the residual in the pressure solver
    int maximumNumberOfIterations = 1e5; //< maximum number of iterations in the solver

    std::string particelShape = "DEFAULT";

    //! parse a text file with settings, each line contains "<parameterName> = <value>"
    void loadFromFile(std::string filename);

    //! output all settings to console
    void printSettings();
};