#include "computation/computation_parallel.h"

#include <cmath>


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

    partitioning_ = std::make_shared<Partitioning>();
    partitioning_->initialize(settings_.nCells);

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
    outputWriterTextParallel_ = std::make_unique<OutputWriterTextParallel>(discretization_, *partitioning_);
    outputWriterParaviewParallel_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, *partitioning_);
}

/**
 * Run the simulation, starting at 0 until tend.
 * Prints 10 time steps to the console.
 */
void ComputationParallel::runSimulation()
{
    // initialize variables
    double t = 0.0;
    double last_printed_time = 0.0;

    while (t < settings_.endTime)
    {
        // set boundary values for u, v, F and G and exchange values at borders btw subdomains
        // exchange velocities at boundaries btw subdomains
        applyBoundaryValues();
        //std::cout << "Process " << partitioning_->ownRankNo() << ": t = " << t << std::endl;
        // compute time step size dt (contributions from all processes)
        computeTimeStepWidthParallel();
        //computeTimeStepWidthAlt();
        //std::cout << "finished computeTimeStepWidth, Process" << partitioning_->ownRankNo() << std::endl;
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

        if (t - last_printed_time >=1.0){
            outputWriterParaviewParallel_->writeFile(t);
            last_printed_time = t;
        }
        #ifdef DEBUG
        if (partitioning_->ownRankNo() == 0) {
            // Create progress bar
            std::cout << "\r" << std::flush;
            int barWidth = 70;
            double progress = t / settings_.endTime;
            std::cout << "[";
            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i)
            {
                if (i < pos)
                    std::cout << "=";
                else if (i == pos)
                    std::cout << ">";
                else
                    std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0) << " %\r" << std::flush;
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
    // std::cout << "Before Allreduce in computeTimeStepWidth, Process " << partitioning_->ownRankNo() << ": dtLocal = " << dtLocal << std::endl;
    // reduce dtLocal to dtGlobal by taking the minimum over all subdomains
    MPI_Allreduce(&dtLocal, &dtGlobal, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    // std::cout << "After Allreduce in computeTimeStepWidth, Process " << partitioning_->ownRankNo() << ": dtLocal = " << dtLocal << std::endl;

    // set time step width to global minimum
    dt_ = dtGlobal;
}

void ComputationParallel::computeTimeStepWidthAlt() {
    const double dx =  discretization_->dx();
    const double dy =  discretization_->dy();

    // Compute maximal time step width regarding the diffusion
    double dt_diff = settings_.re / 2 / (1 / (dx * dx) + 1 / (dy * dy) );

    double maxU = 0.0;
    double maxV = 0.0;
    // Compute maximal time step width regarding the convection u
    std::array<int, 2> uSize = discretization_->uSize();
    for (int i = 0; i < uSize[0]; i++)
    {
        for (int j = 0; j < uSize[1]; j++)
        {
            // possible because the grid for u and v have the same dimensions in all directions
            maxU = std::max(maxU, std::fabs(discretization_->u(i, j)));
            maxV = std::max(maxV, std::fabs(discretization_->v(i, j)));
        }
    }

    
    double u_absMax_local = maxU;

    double u_absMax = 0.0;
    
    MPI_Allreduce(&u_absMax_local,
                  &u_absMax,
                  1,
                  MPI_DOUBLE,
                  MPI_MAX,
                  MPI_COMM_WORLD
    );
    
    double dt_conv_u = std::numeric_limits<double>::max();
    if (u_absMax > 0.0)
        dt_conv_u = dx / u_absMax;

    
    // Compute maximal time step width regarding the convection v
    double v_absMax_local = maxV;

    double v_absMax = 0.0;

    MPI_Allreduce(&v_absMax_local,
                  &v_absMax,
                  1,
                  MPI_DOUBLE,
                  MPI_MAX,
                  MPI_COMM_WORLD
    );
    //std::cout << "u_absMax= "<< u_absMax << std::endl;
    //std::cout << "hier1, Process" << partitioning_->ownRankNo() << std::endl;
    double dt_conv_v = std::numeric_limits<double>::max();
    if (v_absMax > 0.0)
        dt_conv_v = dy / v_absMax;
    
    // Set the appropriate time step width by using a security factor tau
    //std::cout << "hier2, Process" << partitioning_->ownRankNo() << std::endl;
    dt_ = std::min(settings_.tau * std::min(dt_diff, std::min(dt_conv_u,dt_conv_v)), settings_.maximumDt);
    //std::cout << "hier3, Process" << partitioning_->ownRankNo() << std::endl;
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

    // TODO: do we need to wait for all requests to finish?

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
    MPI_Isend(&discretization_->v(discretization_->vIBegin()+1, discretization_->vJBegin() + 1), discretization_->vIEnd() - discretization_->vIBegin()-1, 
        MPI_DOUBLE, bottomNeigbhourRank, 0, MPI_COMM_WORLD, &send_bottom_v);

    // receive ghost layer row of u from bottom neighbouring subdomain
    MPI_Irecv(&discretization_->u(discretization_->uIBegin()+1, discretization_->uJBegin()), discretization_->uIEnd() - discretization_->uIBegin()-1, 
        MPI_DOUBLE, bottomNeigbhourRank, 0, MPI_COMM_WORLD, &recv_bottom_u);

    // receive ghost layer row of v from bottom neighbouring subdomain
    MPI_Irecv(&discretization_->v(discretization_->vIBegin()+1, discretization_->vJBegin()), discretization_->vIEnd() - discretization_->vIBegin()-1, 
        MPI_DOUBLE, bottomNeigbhourRank, 0, MPI_COMM_WORLD, &recv_bottom_v);

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
    MPI_Irecv(&discretization_->v(discretization_->vIBegin()+1, discretization_->vJEnd()), discretization_->vIEnd() - discretization_->vIBegin()-1, 
        MPI_DOUBLE, topNeigbhourRank, 0, MPI_COMM_WORLD, &recv_top_v);

    // send last inner row of u to ghost layer row of top neighbouring subdomain
    MPI_Isend(&discretization_->u(discretization_->uIBegin()+1, discretization_->uJEnd() - 1), discretization_->uIEnd() - discretization_->uIBegin()-1, 
        MPI_DOUBLE, topNeigbhourRank, 0, MPI_COMM_WORLD, &send_top_u);

    // send last inner row of v to ghost layer row of top neighbouring subdomain
    MPI_Isend(&discretization_->v(discretization_->vIBegin()+1, discretization_->vJEnd() - 1), discretization_->vIEnd() - discretization_->vIBegin()-1, 
        MPI_DOUBLE, topNeigbhourRank, 0, MPI_COMM_WORLD, &send_top_v);

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
    int num_rows_v = discretization_->vJEnd() - discretization_->vJBegin() + 1;
    int num_rows_u = discretization_->uJEnd() - discretization_->uJBegin() + 1;
    double* column_u = new double[num_rows_u];
    double* column_v = new double[num_rows_v];
    // std::fill_n(column_v, num_rows_v, 0);
    // std::fill_n(column_u, num_rows_u, 0);

    // std::cout << "Rank " << partitioning_->ownRankNo() << ": num_rows_u = " << num_rows_u << std::endl;
    // std::cout << "Rank " << partitioning_->ownRankNo() << ": num_rows_v = " << num_rows_v << std::endl;
    // std::cout << "1. column_v = [";
    // for (int i = 0; i < num_rows_v; i++) {
    //     std::cout << column_v[i];
    //     if (i < num_rows_v - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;
    // std::cout << "1. column_u = [";
    // for (int i = 0; i < num_rows_u; i++) {
    //     std::cout << column_u[i];
    //     if (i < num_rows_u - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;

    // save column data for sending
    for (int j = 0; j < num_rows_u; j++)
    {
        column_u[j] = discretization_->u(discretization_->uIBegin() + 1, discretization_->uJBegin() + j);
    }
    for (int j = 0; j < num_rows_v; j++)
    {
        column_v[j] = discretization_->v(discretization_->vIBegin() + 1, discretization_->vJBegin() + j);
    }

    // std::cout << "2. column_v = [";
    // for (int i = 0; i < num_rows_v; i++) {
    //     std::cout << column_v[i];
    //     if (i < num_rows_v - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;
    // std::cout << "2. column_u = [";
    // for (int i = 0; i < num_rows_u; i++) {
    //     std::cout << column_u[i];
    //     if (i < num_rows_u - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;

    // send first inner column of u to left neighbouring subdomain
    MPI_Isend(column_u, num_rows_u, MPI_DOUBLE, leftNeigbhourRank, 0, MPI_COMM_WORLD, &send_left_u);

    // send first inner column of v to left neighbouring subdomain
    MPI_Isend(column_v, num_rows_v, MPI_DOUBLE, leftNeigbhourRank, 0, MPI_COMM_WORLD, &send_left_v); 


    // overwrite column vectors for ghost layer column data of u and v
    // receive ghost layer column of u from left neighbouring subdomain
    MPI_Irecv(column_u, num_rows_u, MPI_DOUBLE, leftNeigbhourRank, 0, MPI_COMM_WORLD, &recv_left_u);

    // receive ghost layer column of v from left neighbouring subdomain
    MPI_Irecv(column_v, num_rows_v, MPI_DOUBLE, leftNeigbhourRank, 0, MPI_COMM_WORLD, &recv_left_v);

    // std::cout << "3. column_v = [";
    // for (int i = 0; i < num_rows_v; i++) {
    //     std::cout << column_v[i];
    //     if (i < num_rows_v - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;
    // std::cout << "3. column_u = [";
    // for (int i = 0; i < num_rows_u; i++) {
    //     std::cout << column_u[i];
    //     if (i < num_rows_u - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;

    MPI_Wait(&recv_left_u, MPI_STATUS_IGNORE);
    MPI_Wait(&recv_left_v, MPI_STATUS_IGNORE);

    // std::cout << "4. column_v = [";
    // for (int i = 0; i < num_rows_v; i++) {
    //     std::cout << column_v[i];
    //     if (i < num_rows_v - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;
    // std::cout << "4. column_u = [";
    // for (int i = 0; i < num_rows_u; i++) {
    //     std::cout << column_u[i];
    //     if (i < num_rows_u - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;

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
    double* column_u = new double[num_rows_u];
    double* column_v = new double[num_rows_v];

    //MPI_Wait(&recv_right_u, MPI_STATUS_IGNORE);
    //MPI_Wait(&recv_right_v, MPI_STATUS_IGNORE);

    // receive first inner column of u from right neighbouring subdomain into right columns of current subdomain
    MPI_Irecv(column_u, num_rows_u, MPI_DOUBLE, rightNeigbhourRank, 0, MPI_COMM_WORLD, &recv_right_u);

    // receive first inner column of v from right neighbouring subdomain into right columns of current subdomain
    MPI_Irecv(column_v, num_rows_v, MPI_DOUBLE, rightNeigbhourRank, 0, MPI_COMM_WORLD, &recv_right_v);

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