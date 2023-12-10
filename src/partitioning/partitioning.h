#pragma once

#include <memory>
#include <array>
#include <vector>
#include <mpi.h>

/**
 * The class Partitioning computes the partitioning of the computational domain into subdomains.
 */

class Partitioning
{
public:
  //! compute partitioning, set internal variables
  void initialize(std::array<int, 2> nCellsGlobal);

  //! compute node offset
  //! (i_local,j_local) + nodeOffset = (i_global,j_global)
  void computeNodeOffset();

  //! get the local number of cells in the own subdomain
  std::array<int, 2> nCellsLocal() const;

  //! get the global number of cells in the whole computational domain
  //! used in OutputWriterParaviewParallel
  std::array<int, 2> nCellsGlobal() const;

  //! get the own MPI rank no
  //! used in OutputWriterParaviewParallel and OutputWriterTextParallel
  int ownRankNo() const;

  //! number of MPI ranks
  int nRanks() const;

  //! MPI-send command
  void MPI_send(int destinationRank, std::vector<double> &data);

  //! MPI-isend command
  void MPI_isend(int destinationRank, std::vector<double> &data, MPI_Request &request);

  //! MPI-recv command
  void MPI_recv(int sourceRank, std::vector<double> &data, int count);

  //! MPI-irecv command
  void MPI_irecv(int sourceRank, std::vector<double> &data, int count, MPI_Request &request);

  //! MPI-wait command
  void MPI_wait(MPI_Request &request);

  //! MPI-waitall command
  void MPI_waitall(std::vector<MPI_Request> &requests);

  //! MPI-allreduce command
  double MPI_allreduce(double &localVal, MPI_Op op);

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
  int leftNeighbourRankNo();

  //! get the rank no of the right neighbouring rank
  int rightNeighbourRankNo();

  //! get the rank no of the top neighbouring rank
  int topNeighbourRankNo();

  //! get the rank no of the bottom neighbouring rank
  int bottomNeighbourRankNo();

  //! get the offset values for counting local nodes in x and y direction.
  //! (i_local,j_local) + nodeOffset = (i_global,j_global)
  //! used in OutputWriterParaviewParallel
  std::array<int, 2> nodeOffset() const;

protected:
  //! get the row index of the node
  int nodeRowIndex() const;

  //! get the column index of the node
  int nodeColumnIndex() const;

  //! number of cells in x and y direction
  std::array<int, 2> nCellsGlobal_;

  //! number of cells in x and y direction in the own partition
  std::array<int, 2> nCellsLocal_;

  //! number of MPI ranks
  int nRanks_;

  //! number of domains in x and y direction
  //! 2D array
  std::array<int, 2> nDomains_;

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
  std::array<int, 2> nodeOffset_;

private:
  //! decompose computational domain into partitions
  void computePartitioning(int nRanks);
};
