#include "discretization/0_staggered_grid.h"
//! TODO: Check wheter cassert counts as third party library, which is not allowed according to rule 1
#include <cassert>

//! Constructor
StaggeredGrid::StaggeredGrid(std::array<int, 2> nCells, std::array<double, 2> meshWidth) : nCells_(nCells),
                                                                                           meshWidth_(meshWidth),
                                                                                           //! initialize field variables
                                                                                           //! CHECK: Is origin (the 2nd argument) the physical origin, or index origin?
                                                                                           //! CHECK: Do we need meshwidth/2 as rhs's last argument?
                                                                                           //! CHECK: do we need as origin {meshwidth[0],meshwidh[1]/2} - VAR1
                                                                                           // or VAR2: {0,-meshwidth[1]/2}? So are we considering all u's or just the u's in the computation region?
                                                                                           //! CHECK: same for v,p,rhs,f,g

                                                                                           //! VAR1
                                                                                           // u_(uSize(), {meshWidth[0], meshWidth[1]/2.0}, meshWidth),
                                                                                           // v_(vSize(), {meshWidth[0]/2.0, meshWidth[1]}, meshWidth),
                                                                                           // p_(pSize(), {meshWidth[0]/2.0, meshWidth[1]/2.0}, meshWidth),
                                                                                           // rhs_(rhsSize(), {meshWidth[0]/2.0, meshWidth[1]/2.0}, meshWidth),
                                                                                           // f_(uSize(), {meshWidth[0], meshWidth[1]/2.0}, meshWidth),
                                                                                           // g_(vSize(), {meshWidth[0]/2.0, meshWidth[1]}, meshWidth)

                                                                                           //! VAR2
                                                                                           u_(uSize(), {0.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                           v_(vSize(), {-meshWidth[0] / 2.0, 0.0}, meshWidth),
                                                                                           p_(pSize(), {-meshWidth[0] / 2.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                           rhs_(rhsSize(), {-meshWidth[0] / 2.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                           f_(uSize(), {0.0, -meshWidth[1] / 2.0}, meshWidth),
                                                                                           g_(vSize(), {-meshWidth[0] / 2.0, 0.0}, meshWidth)
{
    assert(nCells[0] > 0);
    assert(nCells[1] > 0);
    assert(meshWidth[0] > 0);
    assert(meshWidth[1] > 0);
}

//! get the mesh width, i.e. the length of a single cell in x and y direction
const std::array<double, 2> StaggeredGrid::meshWidth() const
{
    return meshWidth_;
}

//! get number of cells in each coordinate direction
const std::array<int, 2> StaggeredGrid::nCells() const
{
    return nCells_;
}

//! get a reference to field variable u
const FieldVariable &StaggeredGrid::u() const
{
    return u_;
}
//! get a reference to field variable v
const FieldVariable &StaggeredGrid::v() const
{
    return v_;
}
//! get a reference to field variable p
const FieldVariable &StaggeredGrid::p() const
{
    return p_;
}

//! access value of u in element (i,j)
double StaggeredGrid::u(int i, int j) const
{
    return u_(i, j);
}
//! access value of u in element (x,y)
double &StaggeredGrid::u(int i, int j)
{
    return u_(i, j);
}
//! access value of v in element (i,j)
double StaggeredGrid::v(int i, int j) const
{
    return v_(i, j);
}
//! access value of v in element (x,y)
double &StaggeredGrid::v(int i, int j)
{
    return v_(i, j);
}
//! access value of p in element (i,j)
double StaggeredGrid::p(int i, int j) const
{
    return p_(i, j);
}
//! access value of p in element (x,y)
double &StaggeredGrid::p(int i, int j)
{
    return p_(i, j);
}
//! access value of rhs in element (i,j)
double &StaggeredGrid::rhs(int i, int j)
{
    return rhs_(i, j);
}
//! access value of F in element (i,j)
double &StaggeredGrid::f(int i, int j)
{
    return f_(i, j);
}
//! access value of G in element (i,j)
double &StaggeredGrid::g(int i, int j)
{
    return g_(i, j);
}

//! get the mesh width in x direction
double StaggeredGrid::dx() const
{
    return meshWidth_[0];
}
//! get the mesh width in y direction
double StaggeredGrid::dy() const
{
    return meshWidth_[1];
}

//! TODO: check if -1 is necessary - if we need to allow for additional column of ghost cells
//! first valid index for u in x direction, is 0 because this is the x-index of the left boundary cells
int StaggeredGrid::uIBegin() const
{
    return 0;
}
//! TODO: also danger zone - +1?+2?
//! last valid index for u in x direction
int StaggeredGrid::uIEnd() const
{
    return nCells_[0] + 1;
}
//! TODO: also danger zone
//! first valid index for u in y direction
int StaggeredGrid::uJBegin() const
{
    return 0;
}
//! TODO: also danger zone
//! last valid index for u in y direction
int StaggeredGrid::uJEnd() const
{
    return nCells_[1] + 1;
}
//! TODO: also danger zone
//! size of u
std::array<int, 2> StaggeredGrid::uSize() const
{
    return {nCells_[0] + 2, nCells_[1] + 2};
}

//! TODO: also danger zone
//! first valid index for v in x direction
int StaggeredGrid::vIBegin() const
{
    return 0;
}
//! TODO: also danger zone
//! last valid index for v in x direction
int StaggeredGrid::vIEnd() const
{
    return nCells_[0] + 1;
}
//!  TODO: also danger zone
//! first valid index for v in y direction
int StaggeredGrid::vJBegin() const
{
    return 0;
}
//! TODO: also danger zone
//! last valid index for v in y direction
int StaggeredGrid::vJEnd() const
{
    return nCells_[1] + 1;
}
//! TODO: also danger zone
//! size of v
std::array<int, 2> StaggeredGrid::vSize() const
{
    return {nCells_[0] + 2, nCells_[1] + 2};
}

//! TODO: also danger zone
//!  first valid index for p in x direction
int StaggeredGrid::pIBegin() const
{
    return 0;
}
//! TODO: also danger zone
//! last valid index for p in x direction
int StaggeredGrid::pIEnd() const
{
    return nCells_[0] + 1;
}
//! TODO: also danger zone
//! first valid index for p in y direction
int StaggeredGrid::pJBegin() const
{
    return 0;
}
//! TODO: also danger zone
//! last valid index for p in y direction
int StaggeredGrid::pJEnd() const
{
    return nCells_[1] + 1;
}
//! TODO: also danger zone
//! size of p
std::array<int, 2> StaggeredGrid::pSize() const
{
    return {nCells_[0] + 2, nCells_[1] + 2};
}

//! TODO: check rhsSize - since we use the same mesh for u,v,p and rhs,
// the size of rhs is the same as the size of u,v,p (?)
//! size of rhs
std::array<int, 2> StaggeredGrid::rhsSize() const
{
    return {nCells_[0] + 2, nCells_[1] + 2};
}