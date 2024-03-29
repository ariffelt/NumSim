#include "partitioning/partitioning.h"

#include <cassert>
#include <vector>
#include <cmath>
#include <iostream>

// compute partitioning, set internal variables
void Partitioning::initialize(std::array<int, 2> nCellsGlobal)
{
    nCellsGlobal_ = nCellsGlobal;
    //MPI_Comm_size(MPI_COMM_WORLD, &nRanks_);
    //MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo_);

    nRanks_ = 1;
    ownRankNo_ = 0;
    nDomains_ = {1, 1};

    // compute optimal partitioning
    computePartitioning(nRanks_);

    // compute local number of cells
    nCellsLocal_ = std::array<int, 2>{nCellsGlobal_[0] / nDomains_[0], nCellsGlobal_[1] / nDomains_[1]};
    std::array<int, 2> remainder = std::array<int, 2>{nCellsGlobal_[0] % nDomains_[0], nCellsGlobal_[1] % nDomains_[1]};

    computeNodeOffset();

    // add remainder to the local number of calls in the last rank in each direction
    if (nodeColumnIndex() == nDomains_[0])
    {
        nCellsLocal_[0] += remainder[0];
    }
    if (nodeRowIndex() == nDomains_[1])
    {
        nCellsLocal_[1] += remainder[1];
    }
}

// get the row index of the own subdomain
int Partitioning::nodeRowIndex() const
{
    return (int)(ownRankNo_ / nDomains_[0]) + 1;
}

// get the column index of the own subdomain
int Partitioning::nodeColumnIndex() const
{
    return ownRankNo_ % nDomains_[0] + 1;
}

// compute node offset
// (i_local,j_local) + nodeOffset = (i_global,j_global)
// used in OutputWriterParaviewParallel
void Partitioning::computeNodeOffset()
{
    // get number of columns in previous ranks within one row
    int nColumn = nCellsLocal_[0] * (nodeColumnIndex() - 1);
    // get number of rows in previous ranks within one column
    int nRow = nCellsLocal_[1] * (nodeRowIndex() - 1);

    nodeOffset_ = {nColumn, nRow};
}

// decompose computational domain into partitions
// works also for rectangular domains
void Partitioning::computePartitioning(int nRanks)
{
    assert(nRanks % 2 == 0 || nRanks == 1);
    assert(nRanks > 0);

    // compute optimal partitioning
    int xOpt = 1;
    int yOpt = nRanks;
    for (int x = 1; x <= nRanks; x++)
    {
        if (nRanks % x == 0)
        {
            int y = (int)(nRanks / x);

            // check if the x, y ratio is closer to 1 than the current optimal ratio
            if (nCellsGlobal_[1] / y - nCellsGlobal_[0] / x < fabs((nCellsGlobal_[1] / yOpt - nCellsGlobal_[0] / xOpt)))
            {
                xOpt = x;
                yOpt = y;
            }
        }
    }
    if (nCellsGlobal_[0] >= nCellsGlobal_[1])
    {
        std::swap(xOpt, yOpt);
    }
    nDomains_ = {xOpt, yOpt};
}

// get the local number of cells in the own subdomain
std::array<int, 2> Partitioning::nCellsLocal() const
{
    return nCellsLocal_;
}

// get the global number of cells in the whole computational domain
// used in OutputWriterParaviewParallel
std::array<int, 2> Partitioning::nCellsGlobal() const
{
    return nCellsGlobal_;
}

// get the own MPI rank no
// used in OutputWriterParaviewParallel and OutputWriterTextParallel
int Partitioning::ownRankNo() const
{
    return ownRankNo_;
}

// number of MPI ranks
int Partitioning::nRanks() const
{
    return nRanks_;
}

// MPI-send command
void Partitioning::MPI_send(int destinationRank, std::vector<double> &data)
{
    MPI_Send(data.data(), data.size(), MPI_DOUBLE, destinationRank, 0, MPI_COMM_WORLD);
}

// MPI-isend command
void Partitioning::MPI_isend(int destinationRank, std::vector<double> &data, MPI_Request &request)
{
    MPI_Isend(data.data(), data.size(), MPI_DOUBLE, destinationRank, 0, MPI_COMM_WORLD, &request);
}

// MPI-recv command
void Partitioning::MPI_recv(int sourceRank, std::vector<double> &data, int count)
{
    MPI_Recv(data.data(), count, MPI_DOUBLE, sourceRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

// MPI-irecv command
void Partitioning::MPI_irecv(int sourceRank, std::vector<double> &data, int count, MPI_Request &request)
{
    MPI_Irecv(data.data(), count, MPI_DOUBLE, sourceRank, 0, MPI_COMM_WORLD, &request);
}

// MPI-wait command
void Partitioning::MPI_wait(MPI_Request &request)
{
    MPI_Wait(&request, MPI_STATUS_IGNORE);
}

// MPI-waitall command
void Partitioning::MPI_waitall(std::vector<MPI_Request> &requests)
{
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}

// MPI-allreduce command
double Partitioning::MPI_allreduce(double &localVal, MPI_Op op)
{
    double globalVal = 0.0;
    MPI_Allreduce(&localVal, &globalVal, 1, MPI_DOUBLE, op, MPI_COMM_WORLD);
    return globalVal;
}

// checks if the own partition has part of the bottom boundary of the whole domain
bool Partitioning::ownPartitionContainsBottomBoundary() const
{
    return nodeRowIndex() == 1;
}

// if the own partition has part of the top boundary of the whole domain
// used in OutputWriterParaviewParallel
bool Partitioning::ownPartitionContainsTopBoundary() const
{
    return nodeRowIndex() == nDomains_[1];
}

// if the own partition has part of the left boundary of the whole domain
bool Partitioning::ownPartitionContainsLeftBoundary() const
{
    return nodeColumnIndex() == 1;
}

// if the own partition has part of the right boundary of the whole domain
// used in OutputWriterParaviewParallel
bool Partitioning::ownPartitionContainsRightBoundary() const
{
    return nodeColumnIndex() == nDomains_[0];
}

// get the rank no of the left neighbouring rank
int Partitioning::leftNeighbourRankNo()
{
    if (ownPartitionContainsLeftBoundary())
    {
        return ownRankNo_;
    }
    else
    {
        return ownRankNo_ - 1;
    }
}

// get the rank no of the right neighbouring rank
int Partitioning::rightNeighbourRankNo()
{
    if (ownPartitionContainsRightBoundary())
    {
        return ownRankNo_;
    }
    else
    {
        return ownRankNo_ + 1;
    }
}

// get the rank no of the top neighbouring rank
int Partitioning::topNeighbourRankNo()
{
    if (ownPartitionContainsTopBoundary())
    {
        return ownRankNo_;
    }
    else
    {
        return ownRankNo_ + nDomains_[0];
    }
}

// get the rank no of the bottom neighbouring rank
int Partitioning::bottomNeighbourRankNo()
{
    if (ownPartitionContainsBottomBoundary())
    {
        return ownRankNo_;
    }
    else
    {
        return ownRankNo_ - nDomains_[0];
    }
}

// get the offset values for counting local nodes in x and y direction.
// (i_local,j_local) + nodeOffset = (i_global,j_global)
// used in OutputWriterParaviewParallel
std::array<int, 2> Partitioning::nodeOffset() const
{
    return nodeOffset_;
}