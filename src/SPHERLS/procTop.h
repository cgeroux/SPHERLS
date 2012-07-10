/**
  @file
  
  Header file for the ProcTop class
  
*/

#ifndef PROCTOP_H
#define PROCTOP_H

class ProcTop{
  public:
    int nNumProcs;/**<
      Number of processors in global communicator MPI::COMM_WORLD. The value of this variable is
      independent of processor \ref ProcTop::nRank.
      */
    int *nProcDims;/**<
      Dimensions of the processor topology. It is an array of size 3 to hold the
      size of the processor grid in each dimension. The value of this variable is set in the 
      configuration file "config.xml" which is parsed by the function \ref init. The values of this 
      variable are independent of processor \ref ProcTop::nRank.
      */
    int *nPeriodic;/**<
      Periodic boundary conditions. It is an array of size 3 to tell if a 
      dimension is periodic (wraps) or not. It contains an interger value of 0 or 1. 0, the 
      boundary condition is not periodic, 1 the boundary condition is periodic. The value of this 
      variable is set in the configuration file "config.xml" which is parsed by the function 
      \ref init. The values of this variable are independent of processor \ref ProcTop::nRank.
      */
    int **nCoords;/**<
      Coordinates of the processors. It is of size \ref ProcTop::nNumProcs by 3.
      The values of this variable are independent of processor \ref ProcTop::nRank.
      */
    int nRank;/**<
      Is a unique integer which identifies the processor. The values of \c ProcTop::nRank range
      from 0 to \ref ProcTop::nNumProcs-1 depending on the processor.
      */
    int nNumNeighbors;/**<
      The number of neighbors surrounding the current processor.
      The maximum number of neighbors possible is 27, 3x3x3 don't forget the current processor
      itself can be its own neighbor because of periodic boundary conditions. The value of this
      variable is dependent on processor \ref ProcTop::nRank.
      */
    int *nNeighborRanks;/**<
      \ref ProcTop::nRank s of the neighboring processors.
      An array of size \ref nNumNeighbors to hold ranks of neighbouring processors.
      */
    int nNumRadialNeighbors;/**<
      The number of neighbors in the radial direction.
      Can range from 1 to 2 depending on weather there is a processor beneath or above the current 
      preccessor.
      */
    int *nRadialNeighborRanks;/**<
      \ref ProcTop::nRank s of the neighboring radial processors.
      It is an array of size \ref ProcTop::nNumRadialNeighbors.
      */
    int *nRadialNeighborNeighborIDs;/**<
      Holds the ID of a radialial neighbor, to be used to
      obtain their \ref ProcTop::nRank from \ref ProcTop::nNeighborRanks
      */
    ProcTop();/**<
      Constructor for class \ref ProcTop.
      */
};/**@class ProcTop
  This class manages information which pertains to the processor topology.
  */
#endif
