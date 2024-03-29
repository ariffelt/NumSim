#pragma once

#include <array>

#include "discretization/1_discretization.h"

class CentralDifferences : public Discretization
{
public:
    //! constructor
    CentralDifferences(const std::shared_ptr<Partitioning> partitioning, std::array<double, 2> meshWidth);

    //! TODO: do I need override { behind const? Because I am using the pure virtual function computeDuDx declared in discretization?
    //! compute the 1st derivative ∂ u^2 / ∂x
    virtual double computeDu2Dx(int i, int j) const;

    //! compute the 1st derivative ∂ v^2 / ∂x
    virtual double computeDv2Dy(int i, int j) const;

    //! compute the 1st derivative ∂ (uv) / ∂x
    virtual double computeDuvDx(int i, int j) const;

    //! compute the 1st derivative ∂ (uv) / ∂y
    virtual double computeDuvDy(int i, int j) const;

    //! compute the 1st derivative ∂ (ut) / ∂x
    virtual double computeDutDx(int i, int j) const;

    //! compute the 1st derivative ∂ (vt) / ∂y
    virtual double computeDvtDy(int i, int j) const;
};