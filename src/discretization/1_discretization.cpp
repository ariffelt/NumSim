#include "discretization/1_discretization.h"

#include <cassert>

/**
 * Constructor
 * @param nCells number of cells in x and y direction
 * @param meshWidth mesh width in x and y direction
 */
Discretization::Discretization(std::array<int, 2> nCells, std::array<double, 2> meshWidth) : StaggeredGrid(nCells, meshWidth)
{
}

/**
 * Compute the 2nd derivative ∂^2 u / ∂x^2
 */
double Discretization::computeD2uDx2(int i, int j) const
{
    return (u(i + 1, j) - 2.0 * u(i, j) + u(i - 1, j)) / (dx() * dx());
}

/**
 * Compute the 2nd derivative ∂^2 u / ∂y^2
 */
double Discretization::computeD2uDy2(int i, int j) const
{
    return (u(i, j + 1) - 2.0 * u(i, j) + u(i, j - 1)) / (dy() * dy());
}

/**
 * Compute the 2nd derivative ∂^2 v / ∂x^2
 */
double Discretization::computeD2vDx2(int i, int j) const
{
    return (v(i + 1, j) - 2.0 * v(i, j) + v(i - 1, j)) / (dx() * dx());
}

/**
 * Compute the 2nd derivative ∂^2 v / ∂y^2
 */
double Discretization::computeD2vDy2(int i, int j) const
{
    return (v(i, j + 1) - 2.0 * v(i, j) + v(i, j - 1)) / (dy() * dy());
}

/**
 * Compute the 1st derivative ∂p / ∂x
 */
double Discretization::computeDpDx(int i, int j) const
{
    return (p(i + 1, j) - p(i, j)) / dx();
}

/**
 * Compute the 1st derivative ∂p / ∂y
 */
double Discretization::computeDpDy(int i, int j) const
{
    return (p(i, j + 1) - p(i, j)) / dy();
}

/**
 * Compute the 2nd derivative ∂^2 t / ∂x^2
 */
double Discretization::computeD2tDx2(int i, int j) const
{
    return (t(i + 1, j) - 2.0 * t(i, j) + t(i - 1, j)) / (dx() * dx());
}

/**
 * Compute the 2nd derivative ∂^2 t / ∂y^2
 */
double Discretization::computeD2tDy2(int i, int j) const
{
    return (t(i, j + 1) - 2.0 * t(i, j) + t(i, j - 1)) / (dy() * dy());
}