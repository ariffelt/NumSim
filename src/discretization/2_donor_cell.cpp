#include "discretization/2_donor_cell.h"

#include <cassert>
#include <cmath>

/**
 * Constructor
 * @param partitioning partitioning of the grid
 * @param nCells number of cells in x and y direction
 * @param meshWidth mesh width in x and y direction
 * @param alpha donor cell weight parameter
 * @param gamma donor cell weight parameter for temperature
*/
DonorCell::DonorCell(const std::shared_ptr<Partitioning> partitioning, std::array<double, 2> meshWidth, double alpha, double gamma) : 
Discretization(partitioning, meshWidth), alpha_(alpha), gamma_(gamma)
{
}

/**
 *  compute the 1st derivative ∂ u^2 / ∂x
*/
double DonorCell::computeDu2Dx(int i, int j) const
{
    const double u_half_step_right = (u(i, j) + u(i + 1, j)) / 2.0;
    const double u_half_step_left = (u(i - 1, j) + u(i, j)) / 2.0;

    // compute original central difference part of the equation
    const double central_difference = (pow(u_half_step_right, 2) - pow(u_half_step_left, 2)) / dx();

    // if u_speed_diff_to_right > 0 : flow gets slower from current cell to right cell
    const double u_speed_diff_to_right = (u(i, j) - u(i + 1, j)) / 2.0;
    // if u_speed_diff_to_left > 0: flow gets slower from left cell to current cell
    const double u_speed_diff_to_left = (u(i - 1, j) - u(i, j)) / 2.0;

    // compute donor cell part of the equation
    const double donor_cell = (fabs(u_half_step_right) * u_speed_diff_to_right - fabs(u_half_step_left) * u_speed_diff_to_left)/ dx();

    return central_difference + alpha_ * donor_cell;
}

/**
 * compute the 1st derivative ∂ v^2 / ∂x
*/
double DonorCell::computeDv2Dy(int i, int j) const
{
    const double v_half_step_top = (v(i, j) + v(i, j + 1)) / 2.0;
    const double v_half_step_bottom = (v(i, j - 1) + v(i, j)) / 2.0;

    // compute original central difference part of the equation
    const double central_difference = (pow(v_half_step_top, 2) - pow(v_half_step_bottom, 2)) / dy();

    // if v_speed_diff_to_top > 0 : flow gets slower from current cell to top cell
    const double v_speed_diff_to_top = (v(i, j) - v(i, j + 1)) / 2.0;
    // if v_speed_diff_to_bottom > 0: flow gets slower from bottom cell to current cell
    const double v_speed_diff_to_bottom = (v(i, j - 1) - v(i, j)) / 2.0;

    // compute donor cell part of the equation
    const double donor_cell = (fabs(v_half_step_top) * v_speed_diff_to_top - fabs(v_half_step_bottom) * v_speed_diff_to_bottom)/ dy();

    return central_difference + alpha_ * donor_cell;
}

/**
 * compute the 1st derivative ∂ (uv) / ∂x
*/
double DonorCell::computeDuvDx(int i, int j) const
{
    const double u_half_step_top = (u(i, j + 1) + u(i, j)) / 2.0;
    const double v_half_step_right = (v(i + 1, j) + v(i, j)) / 2.0;
    const double u_half_step_up_step_left = (u(i - 1, j + 1) + u(i - 1, j)) / 2.0;
    const double v_half_step_left = (v(i - 1, j) + v(i, j)) / 2.0;

    // compute original central difference part of the equation
    const double central_difference = (v_half_step_right * u_half_step_top - v_half_step_left * u_half_step_up_step_left) / dx();

    // if v_speed_diff_to_right > 0: flow gets slower from current cell to right cell
    const double v_speed_diff_to_right = (v(i, j) - v(i + 1, j)) / 2.0;
    // if v_speed_diff_to_left > 0: flow gets slower from left cell to current cell
    const double v_speed_diff_to_left = (v(i - 1, j) - v(i, j)) / 2.0;

    // compute donor cell part of the equation
    const double donor_cell = (fabs(u_half_step_top) * v_speed_diff_to_right - fabs(u_half_step_up_step_left) * v_speed_diff_to_left)/ dx();

    return central_difference + alpha_ * donor_cell;
}

/**
 * compute the 1st derivative ∂ (uv) / ∂y
*/
double DonorCell::computeDuvDy(int i, int j) const
{
    const double v_half_step_right = (v(i + 1, j) + v(i, j)) / 2.0;
    const double u_half_step_top = (u(i, j + 1) + u(i, j)) / 2.0;
    const double v_half_step_right_step_down = (v(i, j - 1) + v(i + 1, j - 1)) / 2.0;
    const double u_half_step_down = (u(i, j) + u(i, j - 1)) / 2.0;

    // compute original central difference part of the equation
    const double central_difference = (v_half_step_right * u_half_step_top - v_half_step_right_step_down * u_half_step_down) / dy();

    // if u_speed_diff_to_top > 0: flow gets slower from current cell to top cell
    const double u_speed_diff_to_top = (u(i, j) - u(i, j + 1)) / 2.0;
    // if u_speed_diff_to_down > 0: flow gets slower from cell below to current cell
    const double u_speed_diff_to_down = (u(i, j - 1) - u(i, j)) / 2.0;

    // compute donor cell part of the equation
    const double donor_cell = (fabs(v_half_step_right) * u_speed_diff_to_top - fabs(v_half_step_right_step_down) * u_speed_diff_to_down)/ dy();

    return central_difference + alpha_ * donor_cell;
}

/**
 * compute the 1st derivative ∂ (ut) / ∂x
*/
double DonorCell::computeDutDx(int i, int j) const
{
    const double t_half_step_right = (t(i + 1, j) + t(i, j)) / 2.0;
    const double t_half_step_left = (t(i - 1, j) + t(i, j)) / 2.0;

    // compute original central difference part of the equation
    const double central_difference = (u(i,j) * t_half_step_right - u(i-1,j) * t_half_step_left) / dx();

    // if t_diff_to_right > 0: flow gets slower from current cell to right cell
    const double t_diff_to_right = (t(i, j) - t(i + 1, j)) / 2.0;
    // if t_diff_to_left > 0: flow gets slower from left cell to current cell
    const double t_diff_to_left = (t(i - 1, j) - t(i, j)) / 2.0;

    // compute donor cell part of the equation
    const double donor_cell = (fabs(u(i,j)) * t_diff_to_right - fabs(u(i-1,j)) * t_diff_to_left)/ dx();

    return central_difference + gamma_ * donor_cell;
}

/**
 * compute the 1st derivative ∂ (vt) / ∂y
*/
double DonorCell::computeDvtDy(int i, int j) const
{
    const double t_half_step_top = (t(i, j + 1) + t(i, j)) / 2.0;
    const double t_half_step_bottom = (t(i, j - 1) + t(i, j)) / 2.0;

    // compute original central difference part of the equation
    const double central_difference = (v(i,j) * t_half_step_top - v(i,j-1) * t_half_step_bottom) / dy();

    // if t_diff_to_top > 0: flow gets slower from current cell to top cell
    const double t_diff_to_top = (t(i, j) - t(i, j + 1)) / 2.0;
    // if t_diff_to_bottom > 0: flow gets slower from bottom cell to current cell
    const double t_diff_to_bottom = (t(i, j - 1) - t(i, j)) / 2.0;

    // compute donor cell part of the equation
    const double donor_cell = (fabs(v(i,j)) * t_diff_to_top - fabs(v(i,j-1)) * t_diff_to_bottom)/ dy();

    return central_difference + gamma_ * donor_cell;
}