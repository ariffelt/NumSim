#include "discretization/0_staggered_grid.h"

#include <cassert>
#include <memory>
/**
 * Constructor of staggered grid.
 * Provides several parameters for the staggered grid.
 * @param partitioning partitioning of the grid
 * @param nCells number of cells in each coordinate direction
 * @param meshWidth mesh width in each coordinate direction
 * TODO: implement ghost layers
 */
StaggeredGrid::StaggeredGrid(const std::shared_ptr<Partitioning> partitioning, std::array<double, 2> meshWidth) : partitioning_(partitioning),
                                                                                                                  nCells_(partitioning->nCellsLocal()),
                                                                                                                  meshWidth_(meshWidth),
                                                                                                                  u_(uSize(), {0.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                                                  v_(vSize(), {-meshWidth[0] / 2.0, 0.0}, meshWidth),
                                                                                                                  p_(pSize(), {-meshWidth[0] / 2.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                                                  rhs_(rhsSize(), {-meshWidth[0] / 2.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                                                  f_(uSize(), {0.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                                                  g_(vSize(), {-meshWidth[0] / 2.0, 0.0}, meshWidth)
{
    assert(nCells_[0] > 0);
    assert(nCells_[1] > 0);
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
 * get the size of rhs, same as size of u,v,p
 */
std::array<int, 2> StaggeredGrid::rhsSize() const
{
    return {nCells_[0] + 2, nCells_[1] + 2};
}

/**
 * get offset used in parallel computePreliminaryVelocities and computeVelocities
 */
int StaggeredGrid::getOffsetRight() const
{
    if (partitioning_->ownPartitionContainsRightBoundary())
    {
        return 0;
    }
    return 1;
}

/**
 * get offset used in parallel computePreliminaryVelocities and computeVelocities
 */
int StaggeredGrid::getOffsetTop() const
{
    if (partitioning_->ownPartitionContainsTopBoundary())
    {
        return 0;
    }
    return 1;
}

/**
 * get offset used in parallel red_black_sor
 */
int StaggeredGrid::sor_offset() const
{
    int sum_offsets = partitioning_->nodeOffset()[0] + partitioning_->nodeOffset()[1];
    if (sum_offsets % 2 == 0)
    {
        return 0;
    }
    return 1;
}