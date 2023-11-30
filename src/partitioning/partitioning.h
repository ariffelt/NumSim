#pragma once

#include <array>
#include <mpi.h>

class Partitioning
{
public:

  //! compute partitioning, set internal variables
  void initialize(std::array<int,2> nCellsGlobal);

  //! get the local number of cells in the own subdomain
  std::array<int,2> nCellsLocal() const;

  //! get the global number of cells in the whole computational domain
  //! used in OutputWriterParaviewParallel
  std::array<int,2> nCellsGlobal() const;

  //! get the own MPI rank no
  //! used in OutputWriterParaviewParallel and OutputWriterTextParallel
  int ownRankNo() const;

  //! number of MPI ranks
  int nRanks() const;

  void MPI_isend(int destinationRank, std::vector<double> &data, MPI_Request &request);

  void MPI_irecv(int sourceRank, std::vector<double> &data, int count, MPI_Request &request);

  //! if the own partition has part of the bottom boundary of the whole domain
  bool ownPartitionContainsBottomBoundary() const;

  //! if the own partition has part of the top boundary of the whole domain
  //! used in OutputWriterParaviewParallel
  bool ownPartitionContainsTopBoundary() const;

  //! if the own partition has part of the left boundary of the whole domain
  bool ownPartitionContainsLeftBoundary() const;

  //! if the own partition has part of the right boundary of the whole domain
  //! used in OutputWriterParaviewParallel
  bool ownPartitionContainsRightBoundary() const;

  //! get the rank no of the left neighbouring rank
  int leftNeighbourRankNo() const;

  //! get the rank no of the right neighbouring rank
  int rightNeighbourRankNo() const;

  //! get the rank no of the top neighbouring rank
  int topNeighbourRankNo() const;

  //! get the rank no of the bottom neighbouring rank
  int bottomNeighbourRankNo() const;

  //! get the offset values for counting local nodes in x and y direction. 
  //! (i_local,j_local) + nodeOffset = (i_global,j_global)
  //! used in OutputWriterParaviewParallel
  std::array<int,2> nodeOffset() const;

protected:

  //! number of cells in x and y direction
  std::array<int,2> nCellsGlobal_;

  //! number of cells in x and y direction in the own partition
  std::array<int,2> nCellsLocal_;

  //! number of MPI ranks in x and y direction
  std::array<int,2> nRanks_;

  //! number of partitions in x and y direction
  std::array<int,2> nPartitions_;

  //! rank no of the own MPI rank
  int ownRankNo_;

  //! rank no of the left neighbouring rank
  int leftNeighbourRankNo_;

  //! rank no of the right neighbouring rank
  int rightNeighbourRankNo_;

  //! rank no of the top neighbouring rank
  int topNeighbourRankNo_;

  //! rank no of the bottom neighbouring rank
  int bottomNeighbourRankNo_;

  //! if the own partition has part of the bottom boundary of the whole domain
  bool ownPartitionContainsBottomBoundary_;

  //! if the own partition has part of the top boundary of the whole domain
  bool ownPartitionContainsTopBoundary_;

  //! if the own partition has part of the left boundary of the whole domain
  bool ownPartitionContainsLeftBoundary_;

  //! if the own partition has part of the right boundary of the whole domain
  bool ownPartitionContainsRightBoundary_;

  //! offset values for counting local nodes in x and y direction. 
  //! (i_local,j_local) + nodeOffset = (i_global,j_global)
  std::array<int,2> nodeOffset_;

private:

  //! decompose computational domain into partitions
  void computePartitioning();
};
