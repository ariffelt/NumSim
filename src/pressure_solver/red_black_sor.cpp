#include "red_black_sor.h"

#include <cassert>
#include <cmath>

RedBlackSOR::RedBlackSOR(std::shared_ptr<Discretization> discretization,
                         double epsilon, int maximumNumberOfIterations, double omega,
                         std::shared_ptr<Partitioning> partitioning)
    : PressureSolverParallel(discretization, epsilon, maximumNumberOfIterations, partitioning),
      omega_(omega)
{
}
/**
 * Solve the Poisson equation for pressure with the SOR method
 */
void RedBlackSOR::solve()
{
  // compute prefactor
  const double hx2 = discretization_->dx() * discretization_->dx();
  const double hy2 = discretization_->dy() * discretization_->dy();
  const double prefactor = hx2 * hy2 / (2 * (hx2 + hy2));

  double res2 = getResidual();

  // set initial values
  int iteration = 0;

  // tolerance for the residual, taken to the power of 2 to compare to squared residual
  const double epsilon2 = epsilon_ * epsilon_;
  while (res2 >= epsilon2 && iteration < maximumNumberOfIterations_)
  {
    //* wnn wir ungerade anzahl zellen haben, dann ...?
    // first half, black step
    // exclude the boundary/ghost layers
    for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++)
    {
      // use herefore that j is 1 on the first step, such that the first cell we look at is pIBegin()+1
      for (int i = discretization_->pIBegin() + (j % 2); i < discretization_->pIEnd(); i += 2)
      {
        double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / (hx2);
        double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / (hy2);

        discretization_->p(i, j) = (1 - omega_) * discretization_->p(i, j) + omega_ * prefactor * (px + py - discretization_->rhs(i, j));
      }
    }
    sendAndBorrowValues();

    // second half, red step
    // exclude the boundary/ghost layers
    for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++)
    {
      // one further on the right than before
      for (int i = discretization_->pIBegin() + 1 + (j % 2); i < discretization_->pIEnd(); i += 2)
      {
        double px = (discretization_->p(i - 1, j) + discretization_->p(i + 1, j)) / (hx2);
        double py = (discretization_->p(i, j - 1) + discretization_->p(i, j + 1)) / (hy2);

        discretization_->p(i, j) = (1 - omega_) * discretization_->p(i, j) + omega_ * prefactor * (px + py - discretization_->rhs(i, j));
      }
    }
    sendAndBorrowValues();

    // increase iteration
    iteration++;

    // update residual
    res2 = getResidual();
  }
}