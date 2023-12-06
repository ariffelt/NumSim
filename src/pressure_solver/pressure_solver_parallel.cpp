#include "pressure_solver/pressure_solver_parallel.h"

#include <cassert>

PressureSolverParallel::PressureSolverParallel(std::shared_ptr<Discretization> discretization,
                                               double epsilon,
                                               int maximumNumberOfIterations,
                                               std::shared_ptr<Partitioning> partitioning)
    : PressureSolver(discretization, epsilon, maximumNumberOfIterations),
      partitioning_(partitioning)
{
}

void PressureSolverParallel::sendAndBorrowValues()
{
    MPI_Request req_p_right;
    MPI_Request req_p_left;
    MPI_Request req_p_top;
    MPI_Request req_p_bottom;
    

    int length_p_leftright = discretization_->pJEnd() - discretization_->pJBegin();
    int length_p_bottomtop = discretization_->pIEnd() - discretization_->pIBegin();
    std::vector<double> p_bottom(length_p_bottomtop, 0);
    std::vector<double> p_left(length_p_leftright, 0);
    std::vector<double> p_top(length_p_bottomtop, 0);
    std::vector<double> p_right(length_p_leftright, 0);

    // one partition always borrows the own values to the lower and left neighbour and borrows
    // the values from the upper and right neighbour

    // BOTTOM:
    // first: lower neigbour, send and receive:
    if (partitioning_->ownPartitionContainsBottomBoundary())
    {
        for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++)
        {
            discretization_->p(i, discretization_->pJBegin()) = discretization_->p(i, discretization_->pJBegin() + 1);
        }
    }
    else
    {
        // send bottom cells
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++)
        { // dont overwrite the first and last value
            p_bottom[i] = discretization_->p(i, discretization_->pJBegin() + 1);
        }
        partitioning_->MPI_isend(partitioning_->bottomNeighbourRankNo(), p_bottom, req_p_bottom);


        // receive call bottom cells
        partitioning_->MPI_irecv(partitioning_->bottomNeighbourRankNo(), p_bottom, length_p_bottomtop, req_p_bottom);
        // set bottom cells
        MPI_Wait(&req_p_bottom, MPI_STATUS_IGNORE); // wait for receive to finish
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++)
        {
            discretization_->p(i, discretization_->pJBegin()) = p_bottom[i];
        }


    }

    // LEFT:
    // left neighbour: send and receive
    if (partitioning_->ownPartitionContainsLeftBoundary())
    {
        for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++)
        {
            discretization_->p(discretization_->pIBegin(), j) = discretization_->p(discretization_->pIBegin() + 1, j);
        }
    }
    else
    {
        // send left cells
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++)
        { // dont overwrite the first and last value
            p_left[j] = discretization_->p(discretization_->pIBegin() + 1, j);
        }
        partitioning_->MPI_isend(partitioning_->leftNeighbourRankNo(), p_left, req_p_left);

        // receive call left cells
        partitioning_->MPI_irecv(partitioning_->leftNeighbourRankNo(), p_left, length_p_leftright, req_p_left);
        // set left cells
        MPI_Wait(&req_p_left, MPI_STATUS_IGNORE); // wait for receive to finish
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++)
        {
            discretization_->p(discretization_->pIBegin(), j) = p_left[j];
        }


    }

    // TOP
    //  now top neighbor: careful! this time we borrow values from the top!
    //  I think it does not even make a difference.
    //  Even though we would want to first receive, then compute, then send, the current values that we have in the ghost layer are the ones received
    //  one black/red iteration before. The values I would update are not the ones I would use/compute anyways in the next step (until ) (?)
    //  Folge: Ich schicke und receive einfach immmer blind alle randwerte.
    //  - das was ich hier mache ist das, was vor und nach Schritt 2 passieren soll: dh eigentlihc mache ich hier 3-1-2 (weil 2=berechnen immer aufgerufen wird wenn diese Methode hier fertig ist)
    //-> bei bottom und left: reveive, dann send
    //-> bei top und right: send, dann receive
    if (partitioning_->ownPartitionContainsTopBoundary())
    {
        for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++)
        {
            discretization_->p(i, discretization_->pJEnd()) = discretization_->p(i, discretization_->pJEnd() - 1);
        }
    }
    else
    {
        // send top cells
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++)
        { // dont overwrite the first and last value
            p_top[i] = discretization_->p(i, discretization_->pJEnd() - 1);
        }
        partitioning_->MPI_isend(partitioning_->topNeighbourRankNo(), p_top, req_p_top);

        // receive top cells
        partitioning_->MPI_irecv(partitioning_->topNeighbourRankNo(), p_top, length_p_bottomtop, req_p_top);
        // set top cells
        MPI_Wait(&req_p_top, MPI_STATUS_IGNORE); // wait for receive to finish
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++)
        {
            discretization_->p(i, discretization_->pJEnd()) = p_top[i];
        }
    }

    // RIGHT
    if (partitioning_->ownPartitionContainsRightBoundary())
    {
        for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++)
        {
            discretization_->p(discretization_->pIEnd(), j) = discretization_->p(discretization_->pIEnd() - 1, j);
        }
    }
    else
    {
        // send right cells
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++)
        { // dont overwrite the first and last value
            p_right[j] = discretization_->p(discretization_->pIEnd() - 1, j);
        }
        partitioning_->MPI_isend(partitioning_->rightNeighbourRankNo(), p_right, req_p_right);

        // receive right cells
        partitioning_->MPI_irecv(partitioning_->rightNeighbourRankNo(), p_right, length_p_leftright, req_p_right);
        // set right cells
        MPI_Wait(&req_p_right, MPI_STATUS_IGNORE); // wait for receive to finish
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++)
        {
            discretization_->p(discretization_->pIEnd(), j) = p_right[j];
        }
    }
}

// Only change compared to PressureSolver::getResidual() is the denominator in the return statement (the totalAmountCells)
double PressureSolverParallel::getResidual()
{
    double res = 0.0;
    const double hx2 = discretization_->dx() * discretization_->dx();
    const double hy2 = discretization_->dy() * discretization_->dy();

    // loop over all cells in this partition, compute residual
    for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++)
    {
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++)
        {
            // second order derivative of p in direction x
            double d2pdx2 = (discretization_->p(i - 1, j) - 2.0 * discretization_->p(i, j) + discretization_->p(i + 1, j)) / (hx2);
            // second order derivative of p in direction y
            double d2pdy2 = (discretization_->p(i, j - 1) - 2.0 * discretization_->p(i, j) + discretization_->p(i, j + 1)) / (hy2);
            // compute squared residual
            res += (d2pdx2 + d2pdy2 - discretization_->rhs(i, j)) * (d2pdx2 + d2pdy2 - discretization_->rhs(i, j));
        }
    }

    const int totalAmountCells = partitioning_->nCellsGlobal()[0] * partitioning_->nCellsGlobal()[1];
    return res / totalAmountCells;
}