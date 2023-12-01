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
        // receive call bottom cells
        int length_p_bottom = discretization_->pIEnd() - discretization_->pIBegin(); // make the list of p_bottom values 2 cells, such that value at position 1 will be at the new position 1
        std::vector<double> p_bottom_rcv(length_p_bottom, 0);
        partitioning_->MPI_irecv(partitioning_->bottomNeighbourRankNo(), p_bottom_rcv, length_p_bottom, req_p_bottom);
        // set bottom cells
        MPI_Wait(&req_p_bottom, MPI_STATUS_IGNORE); // wait for receive to finish
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++)
        {
            discretization_->p(i, discretization_->pJBegin()) = p_bottom_rcv[i];
        }

        // send bottom cells
        std::vector<double> p_bottom_send(length_p_bottom, 0);
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++)
        { // dont overwrite the first and last value
            p_bottom_send[i] = discretization_->p(i, discretization_->pJBegin() + 1);
        }
        partitioning_->MPI_isend(partitioning_->bottomNeighbourRankNo(), p_bottom_send, req_p_bottom);
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
        // receive call left cells
        int length_p_left = discretization_->pJEnd() - discretization_->pJBegin(); // make the list of p_bottom values 2 cells, such that value at position 1 will be at the new position 1
        std::vector<double> p_left_rcv(length_p_left, 0);
        partitioning_->MPI_irecv(partitioning_->leftNeighbourRankNo(), p_left_rcv, length_p_left, req_p_left);
        // set left cells
        MPI_Wait(&req_p_left, MPI_STATUS_IGNORE); // wait for receive to finish
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++)
        {
            discretization_->p(discretization_->pIBegin(), j) = p_left_rcv[j];
        }

        // send left cells
        std::vector<double> p_left_send(length_p_left, 0);
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++)
        { // dont overwrite the first and last value
            p_left_send[j] = discretization_->p(discretization_->pIBegin() + 1, j);
        }
        partitioning_->MPI_isend(partitioning_->leftNeighbourRankNo(), p_left_send, req_p_left);
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
        int length_p_top = discretization_->pIEnd() - discretization_->pIBegin(); // make the list of p_bottom values 2 cells, such that value at position 1 will be at the new position 1
        std::vector<double> p_top_send(length_p_top, 0);
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++)
        { // dont overwrite the first and last value
            p_top_send[i] = discretization_->p(i, discretization_->pJEnd() - 1);
        }
        partitioning_->MPI_isend(partitioning_->topNeighbourRankNo(), p_top_send, req_p_top);

        // receive top cells
        std::vector<double> p_top_rcv(length_p_top, 0);
        partitioning_->MPI_irecv(partitioning_->topNeighbourRankNo(), p_top_rcv, length_p_top, req_p_top);
        // set top cells
        MPI_Wait(&req_p_top, MPI_STATUS_IGNORE); // wait for receive to finish
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++)
        {
            discretization_->p(i, discretization_->pJEnd()) = p_top_rcv[i];
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
        int length_p_right = discretization_->pJEnd() - discretization_->pJBegin(); // make the list of p_bottom values 2 cells, such that value at position 1 will be at the new position 1
        std::vector<double> p_right_send(length_p_right, 0);
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++)
        { // dont overwrite the first and last value
            p_right_send[j] = discretization_->p(discretization_->pIEnd() - 1, j);
        }
        partitioning_->MPI_isend(partitioning_->rightNeighbourRankNo(), p_right_send, req_p_right);

        // receive right cells
        std::vector<double> p_right_rcv(length_p_right, 0);
        partitioning_->MPI_irecv(partitioning_->rightNeighbourRankNo(), p_right_rcv, length_p_right, req_p_right);
        // set right cells
        MPI_Wait(&req_p_right, MPI_STATUS_IGNORE); // wait for receive to finish
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++)
        {
            discretization_->p(discretization_->pIEnd(), j) = p_right_rcv[j];
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