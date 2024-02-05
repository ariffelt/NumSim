#pragma once

#include <array>

#include "storage/fieldvariable.h"

class StaggeredGrid
{
public:
    //! constructor
    StaggeredGrid(std::array<int, 2> nCells, std::array<double, 2> meshWidth);

    //! get the mesh width, i.e. the length of a single cell in x and y direction
    const std::array<double, 2> meshWidth() const;

    //! get number of cells in each coordinate direction
    const std::array<int, 2> nCells() const;

    //! get a reference to field variable u
    const FieldVariable &u() const;

    //! get a reference to field variable v
    const FieldVariable &v() const;

    //! get a reference to field variable p
    const FieldVariable &p() const;

    //! get a reference to field variable markerfield
    const FieldVariable &markerfield() const;

    //! access value of u in element (i,j)
    double u(int i, int j) const;

    //! access value of u in element (x,y)
    double &u(int i, int j);

    //! access value of v in element (i,j)
    double v(int i, int j) const;

    //! access value of v in element (x,y)
    double &v(int i, int j);

    //! access value of p in element (i,j)
    double p(int i, int j) const;

    //! access value of p in element (x,y)
    double &p(int i, int j);

    //! access value of rhs in element (i,j)
    double &rhs(int i, int j);

    //! access value of F in element (i,j)
    double &f(int i, int j);

    //! access value of G in element (i,j)
    double &g(int i, int j);

    //! access value of markerfield in element (i,j)
    double markerfield(int i, int j) const;

    //! access value of markerfield in element (x,y)
    double &markerfield(int i, int j);

    //! get the mesh width in x direction
    double dx() const;

    //! get the mesh width in y direction
    double dy() const;

    //! first valid index for u in x direction
    int uIBegin() const;

    //! last valid index for u in x direction
    int uIEnd() const;

    //! first valid index for u in y direction
    int uJBegin() const;

    //! last valid index for u in y direction
    int uJEnd() const;

    //! size of u
    std::array<int, 2> uSize() const;

    //! first valid index for v in x direction
    int vIBegin() const;

    //! last valid index for v in x direction
    int vIEnd() const;

    //! first valid index for v in y direction
    int vJBegin() const;

    //! last valid index for v in y direction
    int vJEnd() const;

    //! size of v
    std::array<int, 2> vSize() const;

    //! first valid index for p in x direction
    int pIBegin() const;

    //! last valid index for p in x direction
    int pIEnd() const;

    //! first valid index for p in y direction
    int pJBegin() const;

    //! last valid index for p in y direction
    int pJEnd() const;

    //! size of p
    std::array<int, 2> pSize() const;

    //! first valid index for markerfield in x direction
    int markerfieldIBegin() const;

    //! last valid index for markerfield in x direction
    int markerfieldIEnd() const;

    //! first valid index for markerfield in y direction
    int markerfieldJBegin() const;

    //! last valid index for markerfield in y direction
    int markerfieldJEnd() const;

    //! size of markerfield
    std::array<int, 2> markerfieldSize() const;

    //! size of rhs
    std::array<int, 2> rhsSize() const;

    bool isFluidCell(int i, int j);

    bool isInnerFluidCell(int i, int j);



protected:
    const std::array<int, 2> nCells_;
    const std::array<double, 2> meshWidth_;

    FieldVariable u_;
    FieldVariable v_;
    FieldVariable p_;
    FieldVariable rhs_;
    FieldVariable f_;
    FieldVariable g_;
    FieldVariable markerfield_;
};