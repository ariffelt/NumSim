// class VirtualParticle {
// public:
//     double x;
//     double y;

//     void setX(double x) {
//         this->x = x;
//     }

//     void setY(double y) {
//         this->y = y;
//     }
// };
#pragma once

#include "storage/array2d.h"

#include <iostream>
#include <array>

// Define a class for a particle
class Particles
{
public:
    //! Constructor to initialize a particle
    Particles(int noParticles);

    //! get the position of particle i
    std::array<double,2> particleI(int i);

protected:
    int noParticles_; // Number of particles
    Array2D particlePositions_; // Position of the particles 
};