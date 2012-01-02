#ifndef DATAMANIPULATION_H
#define DATAMANIPULATION_H

/** @file
  
  Header file for \ref dataManipulation.cpp
*/

#include <mpi.h>
#include "global.h"


void init(ProcTop &procTop,Grid &grid,Output &output,Time &time,Parameters &parameters
  ,PeakKETracking &peakKETracking,MessPass &messPass,Performance& performance,Implicit &implicit
  ,int argc,char* argv[]);/**<
  Initializes the program. It does this by reading a number of configuration options from the config
  file "SPHERLS.xml". It also reads in the starting model, as specified in the "SPHERLS.xml" file,
  using the function \ref modelRead. During the reading of the initial model the \ref modelRead
  function also calls \ref setupLocalGrid to determine the sizes of the local grids and allocate
  memory for them. 
  
  Other things of note that are done in this function are:
    - the calulation timer is started, \ref Performance::dStartTimer
    - It also reads in the equation of state table if using a tabulated equation of state 
      (\ref Parameters::bEOSGammaLaw = false) by calling \ref eos::readBin
    - Initilize the peak kenetic energy tracking
    - Initilizes the watchZones, i.e. figure out which processors have which watch zones, opens the
      files and prints headers.
  
  @param[out] procTop all parts of this stucture are set, and do not change thoughout the
    rest of the calculation.
  @param[out] grid through the function \ref modelRead the function \ref setupLocalGrid is called to
    allocate memory for the grid, and set sizes of it.
  @param[out] output
  @param[out] time
  @param[out] parameters
  @param[out] peakKETracking
  @param[out] messPass
  @param[out] performance
  @param[out] implicit
  @param[in] argc
  @param[in] argv
  */
void setupLocalGrid(ProcTop &procTop, Grid &grid);/**<
  Determins size of local grids (\ref Grid::nLocalGridDims) based on processor topology, and 
  allocates memory for the local grids (\ref Grid::dLocalGridNew, \ref Grid::dLocalGridOld).
  It sets various other quantities aswell such as,
    - the coordinates of all processors (\ref ProcTop::nCoords)
    - the offset for interface centered quantities (\ref Grid::nCenIntOffset, which depends on zoning
      and boundary conditions
    - the position the local grid is in relative to the global grid 
      (\ref Grid::nGlobalGridPositionLocalGrid).
  
  @param[in,out] procTop contains information about the processor topology
  @param[in,out] grid contains information about gird
  */
void fin(bool bWriteCurrentStateToFile,Time &time, Output &output,ProcTop &procTop
  , Grid& grid, PeakKETracking& peakKETracking, Parameters &parameters, Functions &functions
  , Performance& performance,Implicit& implicit);/**<
  Finishes program execution by writing out last grid state, closing output files, and writting out
  run time.
  
  @param[in] bWriteCurrentStateToFile is a bool value which indicates wheather or not to write out
    current model state.
  @param[in] time
  @param[in] output
  @param[in] procTop
  @param[in] grid
  @param[in] peakKETracking
  @param[in] parameters
  @param[in] functions
  @param[in] performance
  @param[in] implicit
  */
void modelWrite_GL(std::string sFileName,ProcTop &procTop, Grid &grid, Time &time
  , Parameters &parameters, PeakKETracking &peakKETracking);/**<
  Writes out a model in distrubuted model format, meaning that each processor writes it's own local
  grid to a file in binary format. They can be combined, and or converted to ascii format using 
  SPHERLSanal. This is for a gamma-law gas model.
  
  @param[in] sFileName base name of the output files
  @param[in] procTop
  @param[in] grid
  @param[in] time
  @param[in] parameters
  */
void modelWrite_TEOS(std::string sFileName,ProcTop &procTop, Grid &grid, Time &time
  , Parameters &parameters, PeakKETracking &peakKETracking);/**<
  Writes out a model in distrubuted model format, meaning that each processor writes it's own local
  grid to a file in binary format. They can be combined, and or converted to ascii format using 
  SPHERLSanal. This is for a tabulated equation of state model.
  
  @param[in] sFileName base name of the output files
  @param[in] procTop
  @param[in] grid
  @param[in] time
  @param[in] parameters
  */
void modelRead(std::string sFileName,ProcTop &procTop, Grid &grid, Time &time
  , Parameters &parameters, PeakKETracking& peakKETracking);/**<
  Reads in a collected binary file into the local grid and calls \ref setupLocalGrid to allocate
  memory and set various parameters of the model. Works for both gamma-law gas, and tabulated
  equation of state models.
  
  @param[in] sFileName name of the file containing the model to be read in
  @param[out] procTop 
  @param[out] grid
  @param[out] time
  @param[out] parameters
  @param[out] peakKETracking
  */
void initUpdateLocalBoundaries(ProcTop &procTop, Grid &grid, MessPass &messPass
  ,Implicit &implicit);/**<
  Sets up MPI derived data types used for updating the local grid boundaries
  between processors. It sets where the local grids should start/stop updating the local grids
  (\ref Grid::nStartUpdateExplicit, \ref Grid::nEndUpdateExplicit, \ref Grid::nStartUpdateImplicit,
  \ref Grid::nEndUpdateImplicit, \ref Grid::nStartGhostUpdateExplicit,
  \ref Grid::nEndGhostUpdateExplicit, \ref Grid::nStartGhostUpdateImplicit,
  \ref Grid::nEndGhostUpdateImplicit). It sets the radial processor neighbors 
  (\ref ProcTop::nNumRadialNeighbors ). 
  
  It also allocates memeory for:
    - \ref MessPass::requestSend
    - \ref MessPass::requestRecv
    - \ref MessPass:statusSend
    - \ref MessPass:statusRecv
  
  @param[in,out] procTop
  @param[in,out] grid
  @param[in,out] messPass
  @param[in,out] implicit
  */
void updateLocalBoundaries(ProcTop &procTop, MessPass &messPass, Grid &grid);/**<
  Updates the boundaries of the local grids from the data in the local grids of other processors. It
  does this for all variables and updates to the old grid. It also has processor 
  \ref ProcTop::nRank=0 call \ref average3DTo1DBoundariesOld which averages the 3D information into
  the 1D boundaries.
  
  @param[in] procTop
  @param[in] messPass
  @param[in,out] grid
  */
void updateLocalBoundariesNewGrid(int nVar, ProcTop &procTop, MessPass &messPass,Grid &grid);/**<
  Updates the boundaries of the local grids from the data in the local grids of other processors. It
  does this for a specific variable specified by \c nVar and updates to the new grid. It also has
  processor \ref ProcTop::nRank=0 call \ref average3DTo1DBoundariesNew which averages the 3D 
  information into the 1D boundaries for that specific variable.
  
  @param[in] procTop
  @param[in] messPass
  @param[in,out] grid
  */
void updateOldGrid(ProcTop &procTop, Grid &grid);/**<
  Updates the old grid with the new grid, not including boundaries.
  
  @param[in] procTop
  @param[in,out] grid
  */
void updateNewGridWithOld(Grid &grid, ProcTop &procTop);/**<
  Copies the contents of the old grid to the new grid including ghost cells.

  @param[in,out] grid
  @param[in] procTop
  */
void average3DTo1DBoundariesOld(Grid &grid);/**<
  This function averages the 3D boundary recieved by the 1D processor (\ref ProcTop::nRank ==0)
  into 1D. This average is volume weighted. This function only needs to be called by the 1D 
  processor, and if called by other processors may have unexpected results. This function calculates
  the average from the old grid, and places the average into the old grid. It does so for all 
  variables external and internal. This function is used every time the grid boundaries are updated
  with \ref updateLocalBoundaries.
  
  @param[in,out] grid supplies the information for calculating the averages and recieves the 
  averages.
  */
void average3DTo1DBoundariesNew(Grid &grid, int nVar);/**<
  This function averages the 3D boundary recieved by the 1D processor (\ref ProcTop::nRank ==0)
  into 1D. This average is volume weighted. This function only needs to be called by the 1D 
  processor, and if called by other processors may have unexpected results. This function calculates
  the average from the new grid, and places the average into new old grid. It does so for only the
  specified variable. This function is used every time the grid boundaries are updated with 
  \ref updateLocalBoundariesNewGrid.
  
  @param[in,out] grid supplies the information for calculating the averages and recieves the 
  averages.
  @param[in] nVar index of the variable to be averaged with in the grid.
  */
void updateLocalBoundaryVelocitiesNewGrid_R(ProcTop &procTop,MessPass &messPass,Grid &grid);/**<
  Updates velocity boundaries of the new grid in a 1D calculations after the velocities have been
  newly calculated.
  */
void updateLocalBoundaryVelocitiesNewGrid_RT(ProcTop &procTop,MessPass &messPass,Grid &grid);/**<
  Updates velocity boundaries of the new grid in a 2D calculations after the velocities have been
  newly calculated.
  */
void updateLocalBoundaryVelocitiesNewGrid_RTP(ProcTop &procTop,MessPass &messPass,Grid &grid);/**<
  Updates velocity boundaries of the new grid in a 3D calculations after the velocities have been
  newly calculated.
  */
void initImplicitCalculation(Implicit &implicit, Grid &grid, ProcTop &procTop, int nNumArgs
  , char* cArgs[]);/**<
  This function initilizes data structures and defines indixes of non-zero elements in the 
  coeffecient matrix. It also sets up pathways for collection of the temperature corrections back to 
  the processors which need them for their local grids.
  
  @param[in,out] implicit
  @param[in] grid size information of the grid is used
  @param[in] procTop
  @param[in] nNumArgs number of command line arguments, PETSc wants them
  @param[in] cArgs a list of command line arguments, PETSc wants them
  */
#endif
