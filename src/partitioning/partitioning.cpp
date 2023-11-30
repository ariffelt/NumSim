#include "partitioning/partitioning.h"

#include <cassert>
#include <vector>

//! compute partitioning, set internal variables
void Partitioning::initialize(std::array<int, 2> nCellsGlobal)
{
    nCellsGlobal_ = nCellsGlobal;
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks_);
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo_);

    // compute optimal partitioning
    //TODO: Dritten Input überprüfen
    MPI_Dims_create(nRanks_, 2, nCellsGlobal_.data());
}

//! decompose computational domain into partitions
void computePartitioning()
{
    //!TODO: override nPartitions_
}

//! get the local number of cells in the own subdomain
std::array<int, 2> Partitioning::nCellsLocal() const
{
}

//! get the global number of cells in the whole computational domain
//! used in OutputWriterParaviewParallel
std::array<int, 2> Partitioning::nCellsGlobal() const
{
}

//! get the own MPI rank no
//! used in OutputWriterParaviewParallel and OutputWriterTextParallel
int Partitioning::ownRankNo() const
{
}

//! number of MPI ranks
int Partitioning::nRanks() const
{
}

void Partitioning::MPI_isend(int destinationRank, std::vector<double> &data, MPI_Request &request)
{
}

void Partitioning::MPI_irecv(int sourceRank, std::vector<double> &data, int count, MPI_Request &request)
{
}

//! if the own partition has part of the bottom boundary of the whole domain
bool Partitioning::ownPartitionContainsBottomBoundary() const
{
}

//! if the own partition has part of the top boundary of the whole domain
//! used in OutputWriterParaviewParallel
bool Partitioning::ownPartitionContainsTopBoundary() const
{
}

//! if the own partition has part of the left boundary of the whole domain
bool Partitioning::ownPartitionContainsLeftBoundary() const
{
}

//! if the own partition has part of the right boundary of the whole domain
//! used in OutputWriterParaviewParallel
bool Partitioning::ownPartitionContainsRightBoundary() const
{
}

//! get the rank no of the left neighbouring rank
int Partitioning::leftNeighbourRankNo() const
{
}

//! get the rank no of the right neighbouring rank
int Partitioning::rightNeighbourRankNo() const
{
}

//! get the rank no of the top neighbouring rank
int Partitioning::topNeighbourRankNo() const
{
}

//! get the rank no of the bottom neighbouring rank
int Partitioning::bottomNeighbourRankNo() const
{
}

//! get the offset values for counting local nodes in x and y direction.
//! (i_local,j_local) + nodeOffset = (i_global,j_global)
//! used in OutputWriterParaviewParallel
std::array<int, 2> Partitioning::nodeOffset() const
{
}


