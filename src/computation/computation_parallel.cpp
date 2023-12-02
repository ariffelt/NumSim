#include "computation/computation_parallel.h"

/**
 * Constructor for the Computation Parallel class.
 * Initializes the settings, discretization, pressure solver, partitioning and output writers
 * @param argc number of command line arguments
 * @param argv command line arguments
 */
void ComputationParallel::initialize(int argc, char *argv[])
{
    // load and print settings
    settings_ = Settings();
    settings_.loadFromFile(argv[1]);
    settings_.printSettings();

    partitioning_ = std::make_shared<Partitioning>(settings_.nCells);
    // TODO: necessary?
    partitioning_->initialize(settings_.nCells);

    // resize mpi_requests
    sendRequests_.resize(16, MPI_REQUEST_NULL);
    receiveRequests_.resize(16, MPI_REQUEST_NULL);

    meshWidth_[0] = settings_.physicalSize[0] / settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1] / settings_.nCells[1];

    // initialize discretization
    if (settings_.useDonorCell)
    {
        discretization_ = std::make_shared<DonorCell>(partitioning_, meshWidth_, settings_.alpha);
    }
    else
    {
        discretization_ = std::make_shared<CentralDifferences>(partitioning_, meshWidth_);
    }


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
    outputWriterTextParallel_ = std::make_unique<OutputWriterTextParallel>(discretization_, partitioning_);
    outputWriterParaviewParallel_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, partitioning_);
}

/**
 * Run the simulation, starting at 0 until tend.
 * Prints 10 time steps to the console.
 */
void ComputationParallel::runSimulation()
{
    // initialize variables
    double t = 0.0;

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

        // compute rhs of pressure equation (each process independent from the others)
        computeRightHandSide();

        // solve pressure equation
        // exchange pressure values at boundaries btw subdomain after each red/black iteration
        // collect the global residual from all processes after each iteration
        computePressure();

        // compute final velocities (each process independent from the others)
        computeVelocities();

        // write output, only write text output in debug mode
        // TODO: write output only every n-th time step
        #ifndef NDEBUG
        outputWriterTextParallel_->writeFile(t);
        outputWriterParaviewParallel_->writeFile(t);
        #else
        outputWriterParaviewParallel_->writeFile(t);
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
 * TODO: when do I want to wait for MPI_requests?
 */
void ComputationParallel::applyBoundaryValues()
{
    // set bottom boundary values
    if (partitioning_->ownPartitionContainsBottomBoundary())
    {
        applyBoundaryValuesBottom();
    }
    else // exchange velocities at bottom boundary
    {
        exchangeVelocitiesBottom();
    }

    // set top boundary values
    if (partitioning_->ownPartitionContainsTopBoundary())
    {
        applyBoundaryValuesTop();
    }
    else // exchange velocities at top boundary
    {
        exchangeVelocitiesTop();
    }

    // set left boundary values
    if (partitioning_->ownPartitionContainsLeftBoundary())
    {
        applyBoundaryValuesLeft();
    }
    else // exchange velocities at left boundary
    {
        exchangeVelocitiesLeft();
    }

    // set right boundary values
    if (partitioning_->ownPartitionContainsRightBoundary())
    {
        applyBoundaryValuesRight();
    }
    else // exchange velocities at right boundary
    {
        exchangeVelocitiesRight();
    }
}

void ComputationParallel::applyBoundaryValuesBottom()
{
    // set bottom boundary values for u and F
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        discretization_->u(i, discretization_->uJBegin()) = 2.0 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin() + 1);
        discretization_->f(i, discretization_->uJBegin()) = discretization_->u(i, discretization_->uJBegin());
    }
    // set bottom boundary values for v and G
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
    {
        discretization_->v(i, discretization_->vJBegin()) = settings_.dirichletBcBottom[1];
        discretization_->g(i, discretization_->vJBegin()) = discretization_->v(i, discretization_->vJBegin());
    }
}

void ComputationParallel::applyBoundaryValuesTop()
{
    // set top boundary values for u and F
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        discretization_->u(i, discretization_->uJEnd()) = 2.0 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() - 1);
        discretization_->f(i, discretization_->uJEnd()) = discretization_->u(i, discretization_->uJEnd());
    }
    // set top boundary values for v and G
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
    {
        discretization_->v(i, discretization_->vJEnd() - 1) = settings_.dirichletBcTop[1];
        discretization_->g(i, discretization_->vJEnd() - 1) = discretization_->v(i, discretization_->vJEnd() - 1);
    }
}

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
 * TODO: do we need to use MPI_Isend or MPI_send or MPI_isend (from partitioning) here?
 */
void ComputationParallel::exchangeVelocitiesBottom()
{
    // get rank no of bottom neighbouring subdomain
    int bottomNeigbhourRank = partitioning_->bottomNeighbourRankNo();

    // send first inner row of u to bottom neighbouring subdomain
    // first input: pointer to first element to send
    // second input: number of elements to send -> this works bc we are sending a row and thats how u is stored
    MPI_Isend(&discretization_->u(discretization_->uIBegin(), discretization_->uJBegin() + 1), discretization_->uIEnd() - discretization_->uIBegin(), 
        MPI_DOUBLE, bottomNeigbhourRank, 0, MPI_COMM_WORLD, &sendRequests_[sendRequestCounter_++]);

    // send first inner row of v to bottom neighbouring subdomain
    MPI_Isend(&discretization_->v(discretization_->vIBegin(), discretization_->vJBegin() + 1), discretization_->vIEnd() - discretization_->vIBegin(), 
        MPI_DOUBLE, bottomNeigbhourRank, 0, MPI_COMM_WORLD, &sendRequests_[sendRequestCounter_++]);

    // receive ghost layer row of u from bottom neighbouring subdomain
    MPI_Irecv(&discretization_->u(discretization_->uIBegin(), discretization_->uJBegin()), discretization_->uIEnd() - discretization_->uIBegin(), 
        MPI_DOUBLE, bottomNeigbhourRank, 0, MPI_COMM_WORLD, &receiveRequests_[receiveRequestCounter_++]);

    // receive ghost layer row of v from bottom neighbouring subdomain
    MPI_Irecv(&discretization_->v(discretization_->vIBegin(), discretization_->vJBegin()), discretization_->vIEnd() - discretization_->vIBegin(), 
        MPI_DOUBLE, bottomNeigbhourRank, 0, MPI_COMM_WORLD, &receiveRequests_[receiveRequestCounter_++]);
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

    // receive first inner row of u from top neighbouring subdomain into last row of the current subdomain
    MPI_Irecv(&discretization_->u(discretization_->uIBegin(), discretization_->uJEnd()), discretization_->uIEnd() - discretization_->uIBegin(), 
        MPI_DOUBLE, topNeigbhourRank, 0, MPI_COMM_WORLD, &receiveRequests_[receiveRequestCounter_++]);

    // receive first inner row of v from top neighbouring subdomain into last row of the current subdomain
    MPI_Irecv(&discretization_->v(discretization_->vIBegin(), discretization_->vJEnd()), discretization_->vIEnd() - discretization_->vIBegin(), 
        MPI_DOUBLE, topNeigbhourRank, 0, MPI_COMM_WORLD, &receiveRequests_[receiveRequestCounter_++]);

    // send last inner row of u to ghost layer row of top neighbouring subdomain
    MPI_Isend(&discretization_->u(discretization_->uIBegin(), discretization_->uJEnd() - 1), discretization_->uIEnd() - discretization_->uIBegin(), 
        MPI_DOUBLE, topNeigbhourRank, 0, MPI_COMM_WORLD, &sendRequests_[sendRequestCounter_++]);

    // send last inner row of v to ghost layer row of top neighbouring subdomain
    MPI_Isend(&discretization_->v(discretization_->vIBegin(), discretization_->vJEnd() - 1), discretization_->vIEnd() - discretization_->vIBegin(), 
        MPI_DOUBLE, topNeigbhourRank, 0, MPI_COMM_WORLD, &sendRequests_[sendRequestCounter_++]);
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

    // create vectors for column data of u and v
    int num_rows_u = discretization_->uJEnd() - discretization_->uJBegin();
    int num_rows_v = discretization_->vJEnd() - discretization_->vJBegin();
    double* column_u = new double[num_rows_u];
    double* column_v = new double[num_rows_v];

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
    MPI_Isend(&column_u, num_rows_u, MPI_DOUBLE, leftNeigbhourRank, 0, MPI_COMM_WORLD, &sendRequests_[sendRequestCounter_++]);

    // send first inner column of v to left neighbouring subdomain
    MPI_Isend(&column_v, num_rows_v, MPI_DOUBLE, leftNeigbhourRank, 0, MPI_COMM_WORLD, &sendRequests_[sendRequestCounter_++]); 


    // overwrite column vectors for ghost layer column data of u and v
    // receive ghost layer column of u from left neighbouring subdomain
    MPI_Irecv(&column_u, num_rows_u, MPI_DOUBLE, leftNeigbhourRank, 0, MPI_COMM_WORLD, &receiveRequests_[receiveRequestCounter_++]);

    // receive ghost layer column of v from left neighbouring subdomain
    MPI_Irecv(&column_v, num_rows_v, MPI_DOUBLE, leftNeigbhourRank, 0, MPI_COMM_WORLD, &receiveRequests_[receiveRequestCounter_++]);

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

    // create vectors for column data of u and v
    int num_rows_u = discretization_->uJEnd() - discretization_->uJBegin();
    int num_rows_v = discretization_->vJEnd() - discretization_->vJBegin();
    double* column_u = new double[num_rows_u];
    double* column_v = new double[num_rows_v];

    // receive first inner column of u from right neighbouring subdomain into right columns of current subdomain
    MPI_Irecv(&column_u, num_rows_u, MPI_DOUBLE, rightNeigbhourRank, 0, MPI_COMM_WORLD, &receiveRequests_[receiveRequestCounter_++]);

    // receive first inner column of v from right neighbouring subdomain into right columns of current subdomain
    MPI_Irecv(&column_v, num_rows_v, MPI_DOUBLE, rightNeigbhourRank, 0, MPI_COMM_WORLD, &receiveRequests_[receiveRequestCounter_++]);

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
    MPI_Isend(&column_u, num_rows_u, MPI_DOUBLE, rightNeigbhourRank, 0, MPI_COMM_WORLD, &sendRequests_[sendRequestCounter_++]);

    // send last inner column of v to ghost layer column of right neighbouring subdomain
    MPI_Isend(&column_v, num_rows_v, MPI_DOUBLE, rightNeigbhourRank, 0, MPI_COMM_WORLD, &sendRequests_[sendRequestCounter_++]);
}