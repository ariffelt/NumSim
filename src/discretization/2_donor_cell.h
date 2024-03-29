#pragma once

#include <array>

#include "discretization/1_discretization.h"

class DonorCell : public Discretization
{
public:
    //! constructor
    DonorCell(const std::shared_ptr<Partitioning> partitioning, std::array<double, 2> meshWidth, double alpha, double gamma);

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

private:
    double alpha_;
    double gamma_;
};