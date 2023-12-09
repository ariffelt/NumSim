#include "pressure_solver/conjugate_gradient.h"

// ConjguateGradient::ConjguateGradient(std::shared_ptr discretization,
//                                      double epsilon,
//                                      int maximumNumberOfIterations,
//                                      std::shared_ptr partitioning)
//     : PressureSolverParallel(discretization, epsilon, maximumNumberOfIterations, partitioning)
//changed to:
ConjugateGradient::ConjugateGradient(std::shared_ptr<Discretization> discretization,
                                     double epsilon,
                                     int maximumNumberOfIterations,
                                     std::shared_ptr<Partitioning> partitioning)
    : PressureSolverParallel(discretization, epsilon, maximumNumberOfIterations, partitioning)
{
    direction_ = std::make_unique<Array2D>(discretization_->pSize());
    // Constructor implementation goes here
}

/**
 * Solve the Poisson equation for pressure with the preconditioned conjugate gradient method
 */
void ConjugateGradient::solve()
{
    // initialize variables
    const double hx2 = discretization_->dx() * discretization_->dx();
    const double hy2 = discretization_->dy() * discretization_->dy();
    const double prefactor = hx2 * hy2 / (2 * (hx2 + hy2));
    const double epsilon2 = epsilon_ * epsilon_;
    Array2D Aq_ = Array2D(discretization_->pSize());
    //Array2D direction_ = Array2D(discretization_->pSize());
    Array2D residual_ = Array2D(discretization_->pSize());
    int iteration = 0;
    double alpha_local = 0.0;
    double alpha_global = 0.0;

    sendAndBorrowValues();

    // Compute the initial residual
    double res2 = getResidual();

    for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd(); i++)
    {
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd(); j++)
        {
            double D2pDx2 = (discretization_->p(i + 1, j) - 2 * discretization_->p(i, j) + discretization_->p(i - 1, j)) / hx2;
            double D2pDy2 = (discretization_->p(i, j + 1) - 2 * discretization_->p(i, j) + discretization_->p(i, j - 1)) / hy2;

            residual_(i, j) = discretization_->rhs(i, j) - (D2pDx2 + D2pDy2);

            (*direction_)(i, j) = residual_(i, j);

            alpha_local += residual_(i, j) * residual_(i, j);
        }
    }

    MPI_Reduce(&alpha_local, &alpha_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    do {
        sendAndBorrowValues_q();

        double lambda_local = 0.0;

        for (int i = 1; i < discretization_->pIEnd(); i++)
        {
            for (int j = 1; j < discretization_->pJEnd(); j++)
            {
                double D2qDx2 = ((*direction_)(i + 1, j) - 2 * (*direction_)(i, j) + (*direction_)(i - 1, j)) / hx2;
                double D2qDy2 = ((*direction_)(i, j + 1) - 2 * (*direction_)(i, j) + (*direction_)(i, j - 1)) / hy2;

                Aq_(i, j) = D2qDx2 + D2qDy2;

                lambda_local += (*direction_)(i, j) * Aq_(i, j);
            }
        }

        double lambda_global = 0.0;
        MPI_Reduce(&lambda_local, &lambda_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        lambda_global = alpha_global / lambda_global;

        // Increase iteration
        iteration++;

        alpha_local = 0.0;
        double alpha_prev = alpha_global;

        for (int i = 1; i < discretization_->pIEnd(); i++)
        {
            for (int j = 1; j < discretization_->pJEnd(); j++)
            {
                discretization_->p(i-1, j-1) += lambda_global * (*direction_)(i, j); //maybe change, seltsam verschiben

                residual_(i, j) -= lambda_global * Aq_(i, j);

                alpha_local += residual_(i, j) * residual_(i, j);
            }
        }

        MPI_Reduce(&alpha_local, &alpha_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        double beta = alpha_global / alpha_prev;

        for (int i = 1; i < discretization_->pIEnd(); i++)
        {
            for (int j = 1; j < discretization_->pJEnd(); j++)
            {
                (*direction_)(i, j) = residual_(i, j) + beta * (*direction_)(i, j);
            }
        }

        res2 = alpha_global / (partitioning_->nCellsGlobal()[0] * partitioning_->nCellsGlobal()[1]); //wenn nicht mal nur damit probieren
        // alpha_global is already the sum over all partitions, so res2 doen't need to be reduced

    } while (res2 >= epsilon2 && iteration < maximumNumberOfIterations_);

    sendAndBorrowValues();
}

void ConjugateGradient::sendAndBorrowValues_q()
{
    // send and borrow values for q

    // Initialize vectors for sending and receiving
    std::vector<double> p_bottom(discretization_->pIEnd() - discretization_->pIBegin()); //die haben pIEnd-pIBegin-2
    std::vector<double> p_top(discretization_->pIEnd() - discretization_->pIBegin());
    std::vector<double> p_left(discretization_->pJEnd() - discretization_->pJBegin());
    std::vector<double> p_right(discretization_->pJEnd() - discretization_->pJBegin());

    // Initialize lengths
    const int length_p_topbottom = p_top.size();
    const int length_p_leftright = p_left.size();

    // Initialize MPI request variables
    MPI_Request req_p_bottom, req_p_top, req_p_left, req_p_right;

    // BOTTOM
    if (partitioning_->ownPartitionContainsBottomBoundary())
    {
        for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++)
        {
            (*direction_)(i,discretization_->pJBegin()) = (*direction_)(i,discretization_->pJBegin() + 1);
        }
    }
    else
    {
        // send bottom cells
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++)
        { // dont overwrite the first and last value
            p_bottom[i] = (*direction_)(i + 1, 1);
        }
        partitioning_->MPI_isend(partitioning_->bottomNeighbourRankNo(), p_bottom, req_p_bottom);

        // receive bottom cells
        partitioning_->MPI_irecv(partitioning_->bottomNeighbourRankNo(), p_bottom, length_p_topbottom, req_p_bottom);
        // set bottom cells
        MPI_Wait(&req_p_bottom, MPI_STATUS_IGNORE); // wait for receive to finish
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++)
        {
            (*direction_)(i + 1,0) = p_bottom[i]; //maybe change
        }
    }

    // TOP
    if (partitioning_->ownPartitionContainsTopBoundary())
    {
        for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++)
        {
            (*direction_)(i, discretization_->pJEnd()) = (*direction_)(i,discretization_->pJEnd() - 1);
        }
    }
    else
    {
        // send top cells
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++)
        { 
            p_top[i] = (*direction_)(i + 1, discretization_->pJEnd() - 1);
        }
        partitioning_->MPI_isend(partitioning_->topNeighbourRankNo(), p_top, req_p_top);

        // receive top cells
        partitioning_->MPI_irecv(partitioning_->topNeighbourRankNo(), p_top, length_p_topbottom, req_p_top);
        // set top cells
        MPI_Wait(&req_p_top, MPI_STATUS_IGNORE); // wait for receive to finish
        for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++)
        {
            (*direction_)(i + 1, discretization_->pJEnd()) = p_top[i];
        }
    }

    // LEFT
    if (partitioning_->ownPartitionContainsLeftBoundary())
    {
        for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++)
        {
            (*direction_)(discretization_->pIBegin(), j) = (*direction_)(discretization_->pIBegin() + 1, j);
        }
    }
    else
    {
        // send left cells
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++)
        { 
            p_left[j] = (*direction_)(1, j + 1);
        }
        partitioning_->MPI_isend(partitioning_->leftNeighbourRankNo(), p_left, req_p_left);

        // receive left cells
        partitioning_->MPI_irecv(partitioning_->leftNeighbourRankNo(), p_left, length_p_leftright, req_p_left);
        // set left cells
        MPI_Wait(&req_p_left, MPI_STATUS_IGNORE); // wait for receive to finish
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++)
        {
            (*direction_)(discretization_->pIBegin(), j + 1) = p_left[j];
        }
    }

    // RIGHT
    if (partitioning_->ownPartitionContainsRightBoundary())
    {
        for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++)
        {
            (*direction_)(discretization_->pIEnd(), j) = (*direction_)(discretization_->pIEnd() - 1, j);
        }
    }
    else
    {
        // send right cells
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++)
        { // dont overwrite the first and last value
            p_right[j] = (*direction_)(discretization_->pIEnd() - 1, j + 1);
        }
        partitioning_->MPI_isend(partitioning_->rightNeighbourRankNo(), p_right, req_p_right);

        // receive right cells
        partitioning_->MPI_irecv(partitioning_->rightNeighbourRankNo(), p_right, length_p_leftright, req_p_right);
        // set right cells
        MPI_Wait(&req_p_right, MPI_STATUS_IGNORE); // wait for receive to finish
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++)
        {
            (*direction_)(discretization_->pIEnd(), j + 1) = p_right[j];
        }
    }
}