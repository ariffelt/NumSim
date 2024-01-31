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

#include <iostream>
#include <vector>

// Define a struct for a particle
struct Particle {
    double x, y; // Position

    // Default constructor
    Particle();

    // Constructor to initialize a particle
    Particle(double initX, double initY);
};