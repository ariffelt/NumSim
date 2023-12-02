#pragma once

#include <array>

#include "discretization/1_discretization.h"

class DonorCell : public Discretization
{
public:
    //! constructor
    DonorCell(Partitioning partitioning, std::array<double, 2> meshWidth, double alpha);

    //! compute the 1st derivative ∂ u^2 / ∂x
    virtual double computeDu2Dx(int i, int j) const;

    //! compute the 1st derivative ∂ v^2 / ∂x
    virtual double computeDv2Dy(int i, int j) const;

    //! compute the 1st derivative ∂ (uv) / ∂x
    virtual double computeDuvDx(int i, int j) const;

    //! compute the 1st derivative ∂ (uv) / ∂y
    virtual double computeDuvDy(int i, int j) const;

private:
    double alpha_;
};