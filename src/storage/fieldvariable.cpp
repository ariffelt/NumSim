#include "storage/fieldvariable.h"

#include <cassert>

FieldVariable::FieldVariable(std::array<int, 2> size,
                             std::array<double, 2> origin,
                             std::array<double, 2> meshWidth) : Array2D(size),
                                                                origin_(origin),
                                                                meshWidth_(meshWidth)
{
    assert(size[0] > 0);
    assert(size[1] > 0);
    assert(meshWidth[0] > 0);
    assert(meshWidth[1] > 0);
}

double FieldVariable::interpolateAt(double x, double y) const
{
    // get the indices of the cell that contains the point (x,y), + 1 because of cell indices definition
    int i = (int)((x - origin_[0]) / meshWidth_[0]);
    int j = (int)((y - origin_[1]) / meshWidth_[1]);

    // check for indices on top or right boundary
    if (i == size_[0] - 1)
        i--;
    if (j == size_[1] - 1)
        j--;

    //! TODO: add eps?

    // get the relative position of the point (x,y) within the cell, + 1 because of cell indices definition
    double relativeX = (x - origin_[0]) / meshWidth_[0] - i;
    double relativeY = (y - origin_[1]) / meshWidth_[1] - j;

    // check if indices are in range
    assert(i >= 0 && i <= size_[0] - 1);
    assert(j >= 0 && j <= size_[1] - 1);

    // get the values at the four corners of the cell
    //! TODO: check if dereferencing is necessary
    double bottomLeft = (*this)(i, j);
    double bottomRight = (*this)(i + 1, j);
    double topLeft = (*this)(i, j + 1);
    double topRight = (*this)(i + 1, j + 1);

    // interpolate linearly in x direction
    double bottom = (1.0 - relativeX) * bottomLeft + relativeX * bottomRight;
    double top = (1.0 - relativeX) * topLeft + relativeX * topRight;

    // interpolate linearly in y direction
    double interpolation = (1.0 - relativeY) * bottom + relativeY * top;

    return interpolation;
}
