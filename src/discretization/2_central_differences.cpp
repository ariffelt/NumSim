#include "discretization/2_central_differences.h"

#include <cassert>
#include <cmath>

/**
 * Constructor
 * @param partitioning partitioning of the grid
 * @param nCells number of cells in x and y direction
 * @param meshWidth mesh width in x and y direction
 */
CentralDifferences::CentralDifferences(std::shared_ptr<Partitioning> partitioning, std::array<int, 2> nCells, std::array<double, 2> meshWidth) : Discretization(partitioning, nCells, meshWidth)
{
}

/**
 * Compute the 1st derivative ∂ u^2 / ∂x
 */
double CentralDifferences::computeDu2Dx(int i, int j) const
{
    const double u_half_step_left = (u(i - 1, j) + u(i, j)) / 2.0;
    const double u_half_step_right = (u(i, j) + u(i + 1, j)) / 2.0;

    return (pow(u_half_step_right, 2) - pow(u_half_step_left, 2)) / dx();
}

/**
 * Compute the 1st derivative ∂ v^2 / ∂y
 */
double CentralDifferences::computeDv2Dy(int i, int j) const
{
    const double v_half_step_bottom = (v(i, j - 1) + v(i, j)) / 2.0;
    const double v_half_step_top = (v(i, j) + v(i, j + 1)) / 2.0;

    return (pow(v_half_step_top, 2) - pow(v_half_step_bottom, 2)) / dy();
}

/**
 * Compute the 1st derivative ∂ (uv) / ∂x
 */
double CentralDifferences::computeDuvDx(int i, int j) const
{
    const double u_half_step_top = (u(i, j + 1) + u(i, j)) / 2.0;
    const double u_half_step_up_step_left = (u(i - 1, j + 1) + u(i - 1, j)) / 2.0;
    const double v_half_step_right = (v(i + 1, j) + v(i, j)) / 2.0;
    const double v_half_step_left = (v(i - 1, j) + v(i, j)) / 2.0;

    return ((v_half_step_right * u_half_step_top) - (v_half_step_left * u_half_step_up_step_left)) / dx();
}

/**
 * Compute the 1st derivative ∂ (uv) / ∂y
 */
double CentralDifferences::computeDuvDy(int i, int j) const
{
    const double u_half_step_top = (u(i, j + 1) + u(i, j)) / 2.0;
    const double u_half_step_down = (u(i, j) + u(i, j - 1)) / 2.0;
    const double v_half_step_right = (v(i + 1, j) + v(i, j)) / 2.0;
    const double v_half_step_right_step_down = (v(i + 1, j - 1) + v(i, j - 1)) / 2.0;

    return ((v_half_step_right * u_half_step_top) - (v_half_step_right_step_down * u_half_step_down)) / dy();
}