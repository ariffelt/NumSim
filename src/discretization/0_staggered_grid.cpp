#include "discretization/0_staggered_grid.h"

#include <cassert>
#include <iostream>

/**
 * Constructor of staggered grid.
 * Provides several parameters for the staggered grid.
 * @param nCells number of cells in each coordinate direction
 * @param meshWidth mesh width in each coordinate direction
*/
StaggeredGrid::StaggeredGrid(std::array<int, 2> nCells, std::array<double, 2> meshWidth) : nCells_(nCells),
                                                                                           meshWidth_(meshWidth),
                                                                                           u_(uSize(), {0.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                           v_(vSize(), {-meshWidth[0] / 2.0, 0.0}, meshWidth),
                                                                                           p_(pSize(), {-meshWidth[0] / 2.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                           t_(tSize(), {-meshWidth[0] / 2.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                           markerfield_(markerfieldSize(), {-meshWidth[0] / 2.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                           rhs_(rhsSize(), {-meshWidth[0] / 2.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                           f_(uSize(), {0.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                           g_(vSize(), {-meshWidth[0] / 2.0, 0.0}, meshWidth),
                                                                                           q_(tSize(), {-meshWidth[0] / 2.0, -meshWidth[1] / 2.0}, meshWidth)
{
    assert(nCells[0] > 0);
    assert(nCells[1] > 0);
    assert(meshWidth[0] > 0);
    assert(meshWidth[1] > 0);
}

/**
 * get the mesh width, i.e. the length of a single cell in x and y direction
*/
const std::array<double, 2> StaggeredGrid::meshWidth() const
{
    return meshWidth_;
}

/**
 * get the number of cells in each coordinate direction
*/
const std::array<int, 2> StaggeredGrid::nCells() const
{
    return nCells_;
}

/**
 * get a reference to field variable u
*/
const FieldVariable &StaggeredGrid::u() const
{
    return u_;
}

/**
 * get a reference to field variable v
*/
const FieldVariable &StaggeredGrid::v() const
{
    return v_;
}

/**
 * get a reference to field variable p
*/
const FieldVariable &StaggeredGrid::p() const
{
    return p_;
}

/**
 * get a reference to field variable t
 */
const FieldVariable &StaggeredGrid::t() const
{
    return t_;
}

/**
 * access value of u in element (i,j)
*/
double StaggeredGrid::u(int i, int j) const
{
    return u_(i, j);
}

/**
 * access value of u in element (x,y)
*/
double &StaggeredGrid::u(int i, int j)
{
    return u_(i, j);
}

/**
 * access value of v in element (i,j)
*/
double StaggeredGrid::v(int i, int j) const
{
    return v_(i, j);
}

/**
 * access value of v in element (x,y)
*/
double &StaggeredGrid::v(int i, int j)
{
    return v_(i, j);
}

/**
 * access value of p in element (i,j)
*/
double StaggeredGrid::p(int i, int j) const
{
    return p_(i, j);
}

/**
 * access value of p in element (x,y)
*/
double &StaggeredGrid::p(int i, int j)
{
    return p_(i, j);
}

/**
 * access value of t in element (i,j)
 */
double StaggeredGrid::t(int i, int j) const
{
    return t_(i, j);
}

/**
 * access value of t in element (x,y)
 */
double &StaggeredGrid::t(int i, int j)
{
    return t_(i, j);
}

/**
 * access value of rhs in element (i,j)
*/
double &StaggeredGrid::rhs(int i, int j)
{
    return rhs_(i, j);
}

/**
 * access value of F in element (i,j)
*/
double &StaggeredGrid::f(int i, int j)
{
    return f_(i, j);
}

/**
 * access value of G in element (i,j)
*/
double &StaggeredGrid::g(int i, int j)
{
    return g_(i, j);
}

/**
 * access value of markerfield in element (i,j)
*/
double StaggeredGrid::markerfield(int i, int j) const
{
    return markerfield_(i, j);
}

/**
 * access value of markerfield in element (x,y)
*/
double &StaggeredGrid::markerfield(int i, int j)
{
    return markerfield_(i, j);
}

/**
 * access value of Q in element (i,j)
*/
double &StaggeredGrid::q(int i, int j)
{
    return q_(i, j);
}

/**
 * get the mesh width in x direction
*/
double StaggeredGrid::dx() const
{
    return meshWidth_[0];
}

/**
 * get the mesh width in y direction
*/
double StaggeredGrid::dy() const
{
    return meshWidth_[1];
}

/**
 * get first valid index for u in x direction, x-index of the left boundary cells
*/
int StaggeredGrid::uIBegin() const
{
    return 0;
}

/**
 * get last valid index for u in x direction
*/
int StaggeredGrid::uIEnd() const
{
    return nCells_[0] + 1;
}

/**
 * get first valid index for u in y direction
*/
int StaggeredGrid::uJBegin() const
{
    return 0;
}

/**
 * get last valid index for u in y direction
*/
int StaggeredGrid::uJEnd() const
{
    return nCells_[1] + 1;
}

/**
 * get the size of u
*/
std::array<int, 2> StaggeredGrid::uSize() const
{
    return {nCells_[0] + 2, nCells_[1] + 2};
}

/**
 * get first valid index for v in x direction
*/
int StaggeredGrid::vIBegin() const
{
    return 0;
}

/**
 * get last valid index for v in x direction
*/
int StaggeredGrid::vIEnd() const
{
    return nCells_[0] + 1;
}

/**
 * get first valid index for v in y direction
*/
int StaggeredGrid::vJBegin() const
{
    return 0;
}

/**
 * get last valid index for v in y direction
*/
int StaggeredGrid::vJEnd() const
{
    return nCells_[1] + 1;
}

/**
 * get the size of v
*/
std::array<int, 2> StaggeredGrid::vSize() const
{
    return {nCells_[0] + 2, nCells_[1] + 2};
}

/**
 * get first valid index for p in x direction
*/
int StaggeredGrid::pIBegin() const
{
    return 0;
}

/**
 * get last valid index for p in x direction
*/
int StaggeredGrid::pIEnd() const
{
    return nCells_[0] + 1;
}

/**
 * get first valid index for p in y direction
*/
int StaggeredGrid::pJBegin() const
{
    return 0;
}

/**
 * get last valid index for p in y direction
*/
int StaggeredGrid::pJEnd() const
{
    return nCells_[1] + 1;
}

/**
 * get the size of p
*/
std::array<int, 2> StaggeredGrid::pSize() const
{
    return {nCells_[0] + 2, nCells_[1] + 2};
}

/**
 * get first valid index for markerfield in x direction
*/
int StaggeredGrid::markerfieldIBegin() const
{
    return 0;
}

/**
 * get last valid index for markerfield in x direction
*/
int StaggeredGrid::markerfieldIEnd() const
{
    return nCells_[0] + 1;
}

/**
 * get first valid index for markerfield in y direction
*/
int StaggeredGrid::markerfieldJBegin() const
{
    return 0;
}

/**
 * get last valid index for markerfield in y direction
*/
int StaggeredGrid::markerfieldJEnd() const
{
    return nCells_[1] + 1;
}

/**
 * get the size of markerfield
*/
std::array<int, 2> StaggeredGrid::markerfieldSize() const
{
    return {nCells_[0] + 2, nCells_[1] + 2};
}

/**
 * get first valid index for t in x direction
 */
int StaggeredGrid::tIBegin() const
{
    return 0;
}

/**
 * get last valid index for t in x direction
 */
int StaggeredGrid::tIEnd() const
{
    return nCells_[0] + 1;
}

/**
 * get first valid index for t in y direction
 */
int StaggeredGrid::tJBegin() const
{
    return 0;
}

/**
 * get last valid index for t in y direction
 */
int StaggeredGrid::tJEnd() const
{
    return nCells_[1] + 1;
}

/**
 * get the size of t
 */
std::array<int, 2> StaggeredGrid::tSize() const
{
    return {nCells_[0] + 2, nCells_[1] + 2};
}

/**
 * get the size of rhs, same as size of u,v,p
*/
std::array<int, 2> StaggeredGrid::rhsSize() const
{
    return {nCells_[0] + 2, nCells_[1] + 2};
}

bool StaggeredGrid::isFluidCell(int i, int j)
{
    return markerfield(i, j) == 1;
}

/**
 * Check if cell has any neighbouring empty cell, if it has none then it is an inner fluid cell.
 */
bool StaggeredGrid::isInnerFluidCell(int i, int j)
{
    // check if cell has any neighbouring empty cell, if it has none then it is an inner fluid cell
    // watch out for boundary cells
    // return true if the current cell is a fluid cell and does not have any empty cell or any boundaries as neighbour
    bool isInnerFluidCell = false;
    if (markerfield(i, j) == 1)
    {
        if (markerfield(i - 1, j) >= 1 && markerfield(i + 1, j) >= 1 && markerfield(i, j - 1) >= 1 && markerfield(i, j + 1) >= 1)
        {
            isInnerFluidCell = true;
            // std::cout << "Inner fluid cell at (" << i << "," << j << ")"<<std::endl;
        }
    }
    return isInnerFluidCell;
}