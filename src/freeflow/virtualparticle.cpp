#include "virtualparticle.h"

// VirtualParticle::VirtualParticle() : x(x), y(y)
// {
// }

#include <iostream>
#include <vector>

// Default constructor
Particle::Particle() : x(0.0), y(0.0) {}

// Constructor to initialize a particle
Particle::Particle(double initX, double initY)
    : x(initX), y(initY) {}

// int main() {
//     // Create an array of particles
//     const int numParticles = 100; // Adjust the number of particles as needed
//     // Particle particles[numParticles];

//     // Alternatively, you can use a vector for dynamic sizing
//     std::vector<Particle> particles(numParticles);

//     // Initialize particles with random positions for demonstration
//     for (int i = 0; i < numParticles; ++i) {
//         particles[i] = Particle(
//             static_cast<double>(rand()) / RAND_MAX * 10.0, // Random x position between 0 and 10
//             static_cast<double>(rand()) / RAND_MAX * 10.0 // Random y position between 0 and 10
//         );
//     }

//     // Print the positions of the particles
//     for (int i = 0; i < numParticles; ++i) {
//         std::cout << "Particle " << i << ": x = " << particles[i].x << ", y = " << particles[i].y << std::endl;
//     }

//     return 0;
// }
