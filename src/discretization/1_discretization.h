#pragma once

#include <array>

#include "discretization/0_staggered_grid.h"

class Discretization : public StaggeredGrid
{
public:
    //! constructor
    Discretization(std::array<int, 2> nCells, std::array<double, 2> meshWidth);

    //! compute the 1st derivative ∂ u^2 / ∂x
    virtual double computeDu2Dx(int i, int j) const = 0;

    //! compute the 1st derivative ∂ v^2 / ∂x
    virtual double computeDv2Dy(int i, int j) const = 0;

    //! compute the 1st derivative ∂ (uv) / ∂x
    virtual double computeDuvDx(int i, int j) const = 0;

    //! compute the 1st derivative ∂ (uv) / ∂y
    virtual double computeDuvDy(int i, int j) const = 0;

    //! compute the 2nd derivative ∂^2 u / ∂x^2
    virtual double computeD2uDx2(int i, int j) const;

    //! compute the 2nd derivative ∂^2 u / ∂y^2
    virtual double computeD2uDy2(int i, int j) const;

    //! compute the 2nd derivative ∂^2 v / ∂x^2
    virtual double computeD2vDx2(int i, int j) const;

    //! compute the 2nd derivative ∂^2 v / ∂y^2
    virtual double computeD2vDy2(int i, int j) const;

    //! compute the 1st derivative ∂p / ∂x
    virtual double computeDpDx(int i, int j) const;

    //! compute the 1st derivative ∂p / ∂y
    virtual double computeDpDy(int i, int j) const;

    //! compute the 2nd derivative ∂^2 t / ∂x^2
    virtual double computeD2tDx2(int i, int j) const;

    //! compute the 2nd derivative ∂^2 t / ∂y^2
    virtual double computeD2tDy2(int i, int j) const;

    //! compute the 1st derivative ∂ (ut) / ∂x
    virtual double computeDutDx(int i, int j) const = 0;

    //! compute the 1st derivative ∂ (vt) / ∂y
    virtual double computeDvtDy(int i, int j) const = 0;
};