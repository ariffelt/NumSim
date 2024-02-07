#include "virtualparticle.h"

// VirtualParticle::VirtualParticle() : x(x), y(y)
// {
// }

#include <iostream>

// Constructor to initialize a particle
Particles::Particles(int noParticles) : noParticles_(noParticles), 
                                        particlePositions_({noParticles,2})
{
    
}

std::array<double,2> Particles::particleI(int i)
{   
    double posX = particlePositions_(i,0);
    double posY = particlePositions_(i,1);
    return {posX, posY};
}