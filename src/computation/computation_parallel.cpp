#include "computation/computation_parallel.h"

/**
 * Constructor for the Computation Parallel class.
 * Initializes the settings, discretization, pressure solver, partitioning and output writers
 * @param argc number of command line arguments
 * @param argv command line arguments
 */
void ComputationParallel::initialize(int argc, char *argv[])
{
    // load settings
    settings_ = Settings();
    settings_.loadFromFile(argv[1]);

    // initialize partitioning
    partitioning_ = std::make_shared<Partitioning>();
    partitioning_->initialize(settings_.nCells);

    // print settings
    if (partitioning_->ownRankNo() == 0)
    {
        settings_.printSettings();
    }

    // initialize mesh width
    meshWidth_[0] = settings_.physicalSize[0] / settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1] / settings_.nCells[1];

    // initialize discretization
    if (settings_.useDonorCell)
    {
        discretization_ = std::make_shared<DonorCell>(partitioning_, meshWidth_, settings_.alpha, settings_.gamma);
    }
    else
    {
        discretization_ = std::make_shared<CentralDifferences>(partitioning_, meshWidth_);
    }

    // initialize pressure solver
    if (settings_.pressureSolver == "SOR")
    {
        pressureSolver_ = std::make_unique<RedBlackSOR>(discretization_, settings_.epsilon,
                                                        settings_.maximumNumberOfIterations, settings_.omega, partitioning_);
    }
    else
    {
        std::cout << "Error: Unknown pressure solver for parallel computing: " << settings_.pressureSolver << std::endl;
        exit(1);
    }

    // initialize output writers
    outputWriterTextParallel_ = std::make_unique<OutputWriterTextParallel>(discretization_, partitioning_);
    outputWriterParaviewParallel_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, partitioning_);
}

/**
 * Run the simulation, starting at 0 until tend.
 * Write out information every whole second.
 */
void ComputationParallel::runSimulation()
{
    // initialize variables
    double t = 0.0;
    double last_printed_time = -1.0;

    while (t < settings_.endTime)
    {
        // set boundary values for u, v, F and G and exchange values at borders btw subdomains
        // exchange velocities at boundaries btw subdomains
        applyBoundaryValues();

        // compute time step size dt (contributions from all processes)
        computeTimeStepWidthParallel();

        // decrease time step width in last time step, s.t. the end time will be reached exactly
        if (t + dt_ > settings_.endTime)
        {
            dt_ = settings_.endTime - t;
        }
        t += dt_;

        // compute preliminary velocities (each process independent from the others)
        computePreliminaryVelocities();

        // communicate preliminary velocities btw subdomains
        communicatePreliminaryVelocities();

        // compute rhs of pressure equation (each process independent from the others)
        computeRightHandSide();

        // solve pressure equation
        // exchange pressure values at boundaries btw subdomain after each red/black iteration
        // collect the global residual from all processes after each iteration
        computePressure();

        // compute final velocities (each process independent from the others)
        computeVelocities();

// write outputs every second, only write text output in debug mode
// debug mode
#ifndef NDEBUG
        if (t - last_printed_time >= 1.0 || t == settings_.endTime)
        {
            outputWriterParaviewParallel_->writeFile(t);
            outputWriterTextParallel_->writeFile(t);
            if (partitioning_->ownRankNo() == 0)
            {
                std::cout << "t = " << t << std::endl;
            }
            last_printed_time = t;
        }
// release mode
#else
        if (t - last_printed_time >= 1.0 || t == settings_.endTime)
        {
            outputWriterParaviewParallel_->writeFile(t);
            last_printed_time = t;
        }
#endif
    }
}

/**
 * Compute the time step width dt from maximum velocities
 * Each process computes its own dt
 * We then compare the local dt to find the global minimum, which is chosen as dt
 */
void ComputationParallel::computeTimeStepWidthParallel()
{
    // compute time step width dt from maximum velocities for each subdomain
    // use old method from Computation
    computeTimeStepWidth();

    // initialize global time step width and dt_ as local
    double dtGlobal;
    double dtLocal = dt_;

    // reduce dtLocal to dtGlobal by taking the minimum over all subdomains
    MPI_Allreduce(&dtLocal, &dtGlobal, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    // set time step width to global minimum
    dt_ = dtGlobal;
}

/**
 * Set velocity boundary values for u, v, F and G
 * Exchange values at borders btw subdomains
 * Start at bottom and top and then left and right
 */
void ComputationParallel::applyBoundaryValues()
{
    // set bottom boundary values
    if (partitioning_->ownPartitionContainsBottomBoundary())
    {
        applyBoundaryValuesBottom();
    }
    else // exchange velocities at bottom subdomain-boundary
    {
        exchangeVelocitiesBottom();
    }

    // set top boundary values
    if (partitioning_->ownPartitionContainsTopBoundary())
    {
        applyBoundaryValuesTop();
    }
    else // exchange velocities at top subdomain-boundary
    {
        exchangeVelocitiesTop();
    }

    // set left boundary values
    if (partitioning_->ownPartitionContainsLeftBoundary())
    {
        applyBoundaryValuesLeft();
    }
    else // exchange velocities at left subdomain-boundary
    {
        exchangeVelocitiesLeft();
    }

    // set right boundary values
    if (partitioning_->ownPartitionContainsRightBoundary())
    {
        applyBoundaryValuesRight();
    }
    else // exchange velocities at right subdomain-boundary
    {
        exchangeVelocitiesRight();
    }
}

/**
 * Apply boundary values for u, v, F and G at bottom domain boundary
 */
void ComputationParallel::applyBoundaryValuesBottom()
{
    // set bottom boundary values for u and F
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        discretization_->u(i, discretization_->uJBegin()) = 2.0 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin() + 1);
        discretization_->f(i, discretization_->uJBegin()) = discretization_->u(i, discretization_->uJBegin());
    }
    // set bottom boundary values for v and G
    for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++)
    {
        discretization_->v(i, discretization_->vJBegin()) = settings_.dirichletBcBottom[1];
        discretization_->g(i, discretization_->vJBegin()) = discretization_->v(i, discretization_->vJBegin());
    }
}

/**
 * Apply boundary values for u, v, F and G at top domain boundary
 */
void ComputationParallel::applyBoundaryValuesTop()
{
    // set top boundary values for u and F
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        discretization_->u(i, discretization_->uJEnd()) = 2.0 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() - 1);
        discretization_->f(i, discretization_->uJEnd()) = discretization_->u(i, discretization_->uJEnd());
    }
    // set top boundary values for v and G
    for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++)
    {
        discretization_->v(i, discretization_->vJEnd() - 1) = settings_.dirichletBcTop[1];
        discretization_->g(i, discretization_->vJEnd() - 1) = discretization_->v(i, discretization_->vJEnd() - 1);
    }
}

/**
 * Apply boundary values for u, v, F and G at left domain boundary
 */
void ComputationParallel::applyBoundaryValuesLeft()
{
    // set left boundary values for u and F
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++)
    {
        discretization_->u(discretization_->uIBegin(), j) = settings_.dirichletBcLeft[0];
        discretization_->f(discretization_->uIBegin(), j) = discretization_->u(discretization_->uIBegin(), j);
    }
    // set left boundary values for v and G
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        discretization_->v(discretization_->vIBegin(), j) = 2 * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin() + 1, j);
        discretization_->g(discretization_->vIBegin(), j) = discretization_->v(discretization_->vIBegin(), j);
    }
}

/**
 * Apply boundary values for u, v, F and G at right domain boundary
 */
void ComputationParallel::applyBoundaryValuesRight()
{
    // set right boundary values for u and F
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++)
    {
        discretization_->u(discretization_->uIEnd() - 1, j) = settings_.dirichletBcRight[0];
        discretization_->f(discretization_->uIEnd() - 1, j) = discretization_->u(discretization_->uIEnd() - 1, j);
    }
    // set right boundary values for v and G
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        discretization_->v(discretization_->vIEnd(), j) = 2 * settings_.dirichletBcRight[1] - discretization_->v(discretization_->vIEnd() - 1, j);
        discretization_->g(discretization_->vIEnd(), j) = discretization_->v(discretization_->vIEnd(), j);
    }
}

/**
 * Exchange velocities at bottom boundary btw subdomains
 * Send first inner row of u and v to bottom neighbouring subdomain
 * Receive ghost layer row of u and v from bottom neighbouring subdomain
 */
void ComputationParallel::exchangeVelocitiesBottom()
{
    // get rank no of bottom neighbouring subdomain
    int bottomNeigbhourRank = partitioning_->bottomNeighbourRankNo();

    // initialize send and receive requests
    MPI_Request send_bottom_u;
    MPI_Request send_bottom_v;
    MPI_Request recv_bottom_u;
    MPI_Request recv_bottom_v;

    // send first inner row of u to bottom neighbouring subdomain
    // first input: pointer to first element to send
    // second input: number of elements to send -> this works bc we are sending a row and thats how u is stored
    MPI_Isend(&discretization_->u(discretization_->uIBegin(), discretization_->uJBegin() + 1), discretization_->uIEnd() - discretization_->uIBegin(),
              MPI_DOUBLE, bottomNeigbhourRank, 0, MPI_COMM_WORLD, &send_bottom_u);

    // send first inner row of v to bottom neighbouring subdomain
    MPI_Isend(&discretization_->v(discretization_->vIBegin() + 1, discretization_->vJBegin() + 1), discretization_->vIEnd() - discretization_->vIBegin() - 1,
              MPI_DOUBLE, bottomNeigbhourRank, 0, MPI_COMM_WORLD, &send_bottom_v);

    // receive ghost layer row of u from bottom neighbouring subdomain
    MPI_Irecv(&discretization_->u(discretization_->uIBegin() + 1, discretization_->uJBegin()), discretization_->uIEnd() - discretization_->uIBegin() - 1,
              MPI_DOUBLE, bottomNeigbhourRank, 0, MPI_COMM_WORLD, &recv_bottom_u);

    // receive ghost layer row of v from bottom neighbouring subdomain
    MPI_Irecv(&discretization_->v(discretization_->vIBegin() + 1, discretization_->vJBegin()), discretization_->vIEnd() - discretization_->vIBegin() - 1,
              MPI_DOUBLE, bottomNeigbhourRank, 0, MPI_COMM_WORLD, &recv_bottom_v);

    // wait for receive to finish
    MPI_Wait(&recv_bottom_u, MPI_STATUS_IGNORE);
    MPI_Wait(&recv_bottom_v, MPI_STATUS_IGNORE);
}

/**
 * Exchange velocities at top boundary btw subdomains
 * Receive first inner row of u and v from top neighbouring subdomain into last row
 * Send last inner row of u and v to ghost layer row of top neighbouring subdomain
 */
void ComputationParallel::exchangeVelocitiesTop()
{
    // get rank no of top neighbouring subdomain
    int topNeigbhourRank = partitioning_->topNeighbourRankNo();

    // initialize send and receive requests
    MPI_Request recv_top_u;
    MPI_Request recv_top_v;
    MPI_Request send_top_u;
    MPI_Request send_top_v;

    // receive first inner row of u from top neighbouring subdomain into last row of the current subdomain
    MPI_Irecv(&discretization_->u(discretization_->uIBegin(), discretization_->uJEnd()), discretization_->uIEnd() - discretization_->uIBegin(),
              MPI_DOUBLE, topNeigbhourRank, 0, MPI_COMM_WORLD, &recv_top_u);

    // receive first inner row of v from top neighbouring subdomain into last row of the current subdomain
    MPI_Irecv(&discretization_->v(discretization_->vIBegin() + 1, discretization_->vJEnd()), discretization_->vIEnd() - discretization_->vIBegin() - 1,
              MPI_DOUBLE, topNeigbhourRank, 0, MPI_COMM_WORLD, &recv_top_v);

    // send last inner row of u to ghost layer row of top neighbouring subdomain
    MPI_Isend(&discretization_->u(discretization_->uIBegin() + 1, discretization_->uJEnd() - 1), discretization_->uIEnd() - discretization_->uIBegin() - 1,
              MPI_DOUBLE, topNeigbhourRank, 0, MPI_COMM_WORLD, &send_top_u);

    // send last inner row of v to ghost layer row of top neighbouring subdomain
    MPI_Isend(&discretization_->v(discretization_->vIBegin() + 1, discretization_->vJEnd() - 1), discretization_->vIEnd() - discretization_->vIBegin() - 1,
              MPI_DOUBLE, topNeigbhourRank, 0, MPI_COMM_WORLD, &send_top_v);

    // wait for receive to finish
    MPI_Wait(&recv_top_u, MPI_STATUS_IGNORE);
    MPI_Wait(&recv_top_v, MPI_STATUS_IGNORE);
}

/**
 * Exchange velocities at left boundary btw subdomains
 * Send first inner column of u and v to left neighbouring subdomain
 * Receive ghost layer column of u and v from left neighbouring subdomain
 */
void ComputationParallel::exchangeVelocitiesLeft()
{
    // get rank no of left neighbouring subdomain
    int leftNeigbhourRank = partitioning_->leftNeighbourRankNo();

    // initialize send and receive requests
    MPI_Request send_left_u;
    MPI_Request send_left_v;
    MPI_Request recv_left_u;
    MPI_Request recv_left_v;

    // create vectors for column data of u and v
    int num_rows_u = discretization_->uJEnd() - discretization_->uJBegin() + 1;
    int num_rows_v = discretization_->vJEnd() - discretization_->vJBegin() + 1;
    double *column_u = new double[num_rows_u];
    double *column_v = new double[num_rows_v];

    // save column data for sending
    for (int j = 0; j < num_rows_u; j++)
    {
        column_u[j] = discretization_->u(discretization_->uIBegin() + 1, discretization_->uJBegin() + j);
    }
    for (int j = 0; j < num_rows_v; j++)
    {
        column_v[j] = discretization_->v(discretization_->vIBegin() + 1, discretization_->vJBegin() + j);
    }

    // send first inner column of u to left neighbouring subdomain
    MPI_Isend(column_u, num_rows_u, MPI_DOUBLE, leftNeigbhourRank, 0, MPI_COMM_WORLD, &send_left_u);

    // send first inner column of v to left neighbouring subdomain
    MPI_Isend(column_v, num_rows_v, MPI_DOUBLE, leftNeigbhourRank, 0, MPI_COMM_WORLD, &send_left_v);

    // overwrite column vectors for ghost layer column data of u and v
    // receive ghost layer column of u from left neighbouring subdomain
    MPI_Irecv(column_u, num_rows_u, MPI_DOUBLE, leftNeigbhourRank, 0, MPI_COMM_WORLD, &recv_left_u);

    // receive ghost layer column of v from left neighbouring subdomain
    MPI_Irecv(column_v, num_rows_v, MPI_DOUBLE, leftNeigbhourRank, 0, MPI_COMM_WORLD, &recv_left_v);

    // wait for receive to finish before overwriting ghost layer column data
    MPI_Wait(&recv_left_u, MPI_STATUS_IGNORE);
    MPI_Wait(&recv_left_v, MPI_STATUS_IGNORE);

    // overwrite ghost layer column data of u and v
    for (int j = 0; j < num_rows_u; j++)
    {
        discretization_->u(discretization_->uIBegin(), discretization_->uJBegin() + j) = column_u[j];
    }
    for (int j = 0; j < num_rows_v; j++)
    {
        discretization_->v(discretization_->vIBegin(), discretization_->vJBegin() + j) = column_v[j];
    }
}

/**
 * Exchange velocities at right boundary btw subdomains
 * Receive first inner column of u and v from right neighbouring subdomain into right columns of current subdomain
 * Send last inner column of u and v to ghost layer column of right neighbouring subdomain
 */
void ComputationParallel::exchangeVelocitiesRight()
{
    // get rank no of right neighbouring subdomain
    int rightNeigbhourRank = partitioning_->rightNeighbourRankNo();

    // initialize send and receive requests
    MPI_Request recv_right_u;
    MPI_Request recv_right_v;
    MPI_Request send_right_u;
    MPI_Request send_right_v;

    // create vectors for column data of u and v
    int num_rows_u = discretization_->uJEnd() - discretization_->uJBegin() + 1;
    int num_rows_v = discretization_->vJEnd() - discretization_->vJBegin() + 1;
    double *column_u = new double[num_rows_u];
    double *column_v = new double[num_rows_v];

    // receive first inner column of u from right neighbouring subdomain into right columns of current subdomain
    MPI_Irecv(column_u, num_rows_u, MPI_DOUBLE, rightNeigbhourRank, 0, MPI_COMM_WORLD, &recv_right_u);

    // receive first inner column of v from right neighbouring subdomain into right columns of current subdomain
    MPI_Irecv(column_v, num_rows_v, MPI_DOUBLE, rightNeigbhourRank, 0, MPI_COMM_WORLD, &recv_right_v);

    // wait for receive to finish before overwriting right columns
    MPI_Wait(&recv_right_u, MPI_STATUS_IGNORE);
    MPI_Wait(&recv_right_v, MPI_STATUS_IGNORE);

    for (int j = 0; j < num_rows_u; j++)
    {
        // overwrite right columns of current subdomain
        discretization_->u(discretization_->uIEnd(), discretization_->uJBegin() + j) = column_u[j];
        // save column data for sending
        column_u[j] = discretization_->u(discretization_->uIEnd() - 1, discretization_->uJBegin() + j);
    }
    for (int j = 0; j < num_rows_v; j++)
    {
        // overwrite right columns of current subdomain
        discretization_->v(discretization_->vIEnd(), discretization_->vJBegin() + j) = column_v[j];
        // save column data for sending
        column_v[j] = discretization_->v(discretization_->vIEnd() - 1, discretization_->vJBegin() + j);
    }

    // send last inner column of u to ghost layer column of right neighbouring subdomain
    MPI_Isend(column_u, num_rows_u, MPI_DOUBLE, rightNeigbhourRank, 0, MPI_COMM_WORLD, &send_right_u);

    // send last inner column of v to ghost layer column of right neighbouring subdomain
    MPI_Isend(column_v, num_rows_v, MPI_DOUBLE, rightNeigbhourRank, 0, MPI_COMM_WORLD, &send_right_v);
}

/**
 * Communicate preliminary velocities F and G btw subdomains, s.t. we have values at the boundaries and can compute the rhs over whole inner domain
 */
void ComputationParallel::communicatePreliminaryVelocities()
{
    // initialize helper vectors
    int num_rows_F = discretization_->uJEnd() - discretization_->uJBegin();
    int num_rows_G = discretization_->vJEnd() - discretization_->vJBegin();
    int num_cols_F = discretization_->uIEnd() - discretization_->uIBegin();
    int num_cols_G = discretization_->vIEnd() - discretization_->vIBegin();

    // receive G values at bottom boundary
    if (!partitioning_->ownPartitionContainsBottomBoundary())
    {
        // initialize receive request
        MPI_Request request_G_bottom;

        // receive row from bottom neighbouring subdomain into first inner row
        MPI_Irecv(&discretization_->g(discretization_->vIBegin() + 1, discretization_->vJBegin()), num_cols_G, MPI_DOUBLE,
                  partitioning_->bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &request_G_bottom);

        // wait for receive to finish
        partitioning_->MPI_wait(request_G_bottom);
    }

    // send G values at top boundary
    if (!partitioning_->ownPartitionContainsTopBoundary())
    {
        // initialize send request
        MPI_Request request_G_top;

        // send last inner row of G to top neighbouring subdomain
        MPI_Isend(&discretization_->g(discretization_->vIBegin() + 1, discretization_->vJEnd() - 1), num_cols_G, MPI_DOUBLE,
                  partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD, &request_G_top);
    }

    // receive F values at left boundary
    if (!partitioning_->ownPartitionContainsLeftBoundary())
    {
        // initialize receive request
        MPI_Request request_F_left;

        // create vector for left column of F
        std::vector<double> F_left(num_rows_F, 0);

        // receive ghost layer column from left neighbouring subdomain
        partitioning_->MPI_irecv(partitioning_->leftNeighbourRankNo(), F_left, num_rows_F, request_F_left);

        // wait for receive to finish
        partitioning_->MPI_wait(request_F_left);

        // set first column of F to received values
        for (int j = discretization_->uJBegin() + 1; j < discretization_->uJEnd(); j++)
        {
            discretization_->f(discretization_->uIBegin(), j) = F_left[j - 1];
        }
    }

    // send and receive F values at right boundary
    if (!partitioning_->ownPartitionContainsRightBoundary())
    {
        // initialize send request
        MPI_Request request_F_right;

        // send last inner column of u to right neighbouring subdomain
        std::vector<double> F_right(num_rows_F, 0);

        // save last inner column of F to send
        for (int j = discretization_->uJBegin() + 1; j < discretization_->uJEnd(); j++)
        {
            F_right[j - 1] = discretization_->f(discretization_->uIEnd() - 1, j);
        }

        // send last inner column of F to right neighbouring subdomain
        partitioning_->MPI_isend(partitioning_->rightNeighbourRankNo(), F_right, request_F_right);
    }
}