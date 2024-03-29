#pragma once

#include "output_writer/output_writer_text.h"
#include "partitioning/partitioning.h"

/** Write *.txt files that are useful for debugging.
 *  All values are written to the file as they are stored in the field variables,
 *  no interpolation takes place.
 */
class OutputWriterTextParallel : 
  public OutputWriterText
{
public:
  //Constructor
  OutputWriterTextParallel(std::shared_ptr<Discretization> discretization, std::shared_ptr<Partitioning> partitioning);

  //! write current velocities to file, filename is output_<count>.<rankNo>.txt
  void writeFile(double currentTime);

  //! write only current values of pressure to file, filename is pressure_<count>.<rankNo>.txt
  void writePressureFile();

  //! include partitioning
  std::shared_ptr<Partitioning> partitioning_;
};
