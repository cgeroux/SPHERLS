#ifndef GLOBAL_H
#define GLOBAL_H

/** @file
  
  Header file for \ref global.cpp.
  
  This file contains definitions which are required throughout the program. The classes defined 
  herein are used through out the program.
*/

#include <vector>
#include <mpi.h>
#include "watchzone.h"
#include "eos.h"
#include "petscksp.h"
#include <csignal>
#include <limits>
#include "profileData.h"
#include "procTop.h"
#include "time.h"

//Debugging flags
#define SIGNEGDEN 0/**<
  Raise signal on calculation of negative density if set to 1. Useful when debugging, it will stop 
  the debugger at the location of the calculation of the negative density. If not 1, it will speed 
  up calculation slightly and generate more useful output upon detection of negative densities. If 1
  and not being run in the debugger, it likely won't generate any usefull output upon negative
  density, and wil simply abort the program.
  */
#define SIGNEGENG 0/**<
  Raise signal on calculation of negative energy if set to 1, else don't rais a signal. Otherwise it
  will be handled through the normal exception method. This is useful when debugging, it will stop
  the debugger at the location of the calculation of the negative energy. If not 1, it will speed up
  calculation slightly and generate more useful output upon detection of negative energy. If 1 and
  not being run in the debugger, it likely won't generate any usefull output upon negative energies,
  and wil simply abort the program.
  */
#define SIGNEGTEMP 0/**<
  Raise signal on calculation of negative temperature if set to 1, else don't rais a signal. 
  Otherwise it will be handled through the normal exception method. This is useful when debugging, 
  it will stop the debugger at the location of the calculation of the negative energy. If not 1, it 
  will speed up calculation slightly and generate more useful output upon detection of negative 
  energy. If 1 and not being run in the debugger, it likely won't generate any usefull output upon 
  negative energies, and wil simply abort the program.
  */
#define TRACKMAXSOLVERERROR 0/**<
  Report the error of the linear equation solver if set to 1, else don't. Not tracking the error
  reduces the calculations per iteration and will speed up running, however if there is question of
  weather the solver is working accurately this is very handy to turn on.
  */
#define SEDOV 0/**<
  If 1 we are preforming the sedov test, which sets special boundary conditions, if 0 we use normal 
  boundary conditions. It also handles artificial viscosity, and timestep slightly differently.
  */
#define VISCOUS_ENERGY_EQ 1/**<
  If 1 will include viscosity in the energy equation. If 0 it won't. This normally should be set to 
  1
  */
#define DUMP_VERSION 1/**<
  Sets the version of the dump file. Should be incremented if changes are made to the information that
  is printed out in a dump.
  */
#define DEBUG_EQUATIONS 0/**<
  If 1 will write out in the form of a profile file, all the horizontal maximum values of all terms
  in all equations.
  */
#define DEDEM_CLAMP 1/**<
  If 1 a clamp on the DEDM gradient will be used to limit how large DE/DM becomes in the advection
  term in the energy equation.
  */

//classes
class MessPass{
  public:
    MPI::Datatype *typeSendNewGrid;/**<
      Send data types for entire grid. It is of size \ref ProcTop::nNumNeighbors.
      */
    MPI::Datatype *typeRecvOldGrid;/**<
      Recv data types for entire grid. It is of sizee \ref ProcTop::nNumNeighbors.
      */
    MPI::Datatype **typeSendNewVar;/**<
      Send data types for variables. It is of size \ref ProcTop::nNumNeighbors by
      \ref Grid::nNumVars.
      */
    MPI::Datatype **typeRecvNewVar;/**<
      Recieve data types for variables. It is of size \ref ProcTop::nNumNeighbors by
      \ref Grid::nNumVars.
      */
    MPI::Request *requestSend;/**<
      Message handles.
      */
    MPI::Request *requestRecv;/**<
      Message handles.
      */
    MPI::Status *statusSend;/**<
      Message status.
      */
    MPI::Status *statusRecv;/**<
      Message status.
      */
    MessPass();/**<
      Constructor for class \ref MessPass.
      */
};/**@class MessPass
  This class manages information which pertains to message passing between processors.
  */
class Grid{
  public:
    //Grid indices of variables
    int nM;/**<
      Index of \f$M_r\f$ independent variable in grid \ref Grid::dLocalGridOld and 
      \ref Grid::dLocalGridNew. This is an external grid variable included in the count 
      \ref Grid::nNumVars. This is an independent grid variable.
      */
    int nTheta;/**<
      Index of \f$\theta\f$ independent variable in grid \ref Grid::dLocalGridOld and
      \ref Grid::dLocalGridNew. This is an external grid variable included in the count
      \ref Grid::nNumVars. This is an independent grid variable.
      */
    int nPhi;/**<
      Index of \f$\phi\f$ independent variable in grid \ref Grid::dLocalGridOld and
      \ref Grid::dLocalGridNew. This is an external grid variable included in the count
      \ref Grid::nNumVars. This is an independent grid variable.
      */
    int nDM;/**<
      Index of \f$\delta M\f$ in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew. This
      is an external grid variable included in the count \ref Grid::nNumVars
      */
    int nR;/**<
      Index of \f$r\f$ in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew. This is an
      external grid variable included in the count \ref Grid::nNumVars
      */
    int nD;/**<
      Index of \f$\rho\f$ in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew. This is
      an external grid variable included in the count \ref Grid::nNumVars
      */
    int nU;/**<
      Index of \f$u\f$ in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew. This is an
      external grid variable included in the count \ref Grid::nNumVars
      */
    int nU0;/**<
      Index of \f$u_0\f$ in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew. This is an
      external grid variable included in the count \ref Grid::nNumVars
      */
    int nV;/**<
      Index of \f$v\f$ in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew. This is an
      external grid variable included in the count \ref Grid::nNumVars
      */
    int nW;/**<
      Index of \f$w\f$ in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew. This is an
      external grid variable included in the count \ref Grid::nNumVars
      */
    int nT;/**<
      Index of \f$T\f$ in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew. This is an
      external grid variable included in the count \ref Grid::nNumVars.  This variable is defined 
      at cell centers.
      */
    int nE;/**<
      Index of \f$E\f$ in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew. This is an
      internal grid variable included in the count \ref Grid::nNumIntVars, unless the calculation is
      adiabatic in which case it is an external grid variable. This variable is defined at cell
      centers.
      */
    int nP;/**<
      Index of Pressure in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew. This is an
      internal grid variable and is included in the count of \ref Grid::nNumIntVars. This variable 
      is defined at cell centers.
      */
    int nKappa;/**<
      Index of Opacity in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew. This is an
      internal grid variable and is included in the count of \ref Grid::nNumIntVars. This variable 
      is defined at cell centers.
      */
    int nGamma;/**<
      Index of adiabatic index in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew. This
      is an internal grid variable and is included in the count of \ref Grid::nNumIntVars. This 
      variable is defined at cell centers.
      */
    int nDenAve;/**<
      Index of \f$\langle\rho\rangle\f$ in grids, \ref Grid::dLocalGridOld and 
      \ref Grid::dLocalGridNew. This is an internal grid variable and is included in the count of 
      \ref Grid::nNumIntVars. This variable is defined at cell centers only in the radial direction.
      */
    int nQ0;/**<
      Index of the radial artificial viscosity in grids, \ref Grid::dLocalGridOld and
      \ref Grid::dLocalGridNew. This is an internal grid variable and is included in the count of
      \ref Grid::nNumIntVars. This variable is defined at cell centers.
      */
    int nQ1;/**<
      Index of the theta artificial viscosity in grids, \ref Grid::dLocalGridOld and
      \ref Grid::dLocalGridNew. This is an internal grid variable and is included in the count of
      \ref Grid::nNumIntVars. This variable is defined at cell centers.
      */
    int nQ2;/**<
      Index of the phi artificial viscosity in grids, \ref Grid::dLocalGridOld and 
      \ref Grid::dLocalGridNew. This is an internal grid variable and is included in the count of
      \ref Grid::nNumIntVars. This variable is defined at cell centers.
      */
    int nDTheta;/**<
      Index of \f$\Delta \theta\f$ in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew.
      This is an internal grid variable and is included in the count of \ref Grid::nNumIntVars. This
      variable is defined at cell centers.
      */
    int nDPhi;/**<
      Index of \f$\Delta \phi\f$ in grids, \ref Grid::dLocalGridOld and \ref Grid::dLocalGridNew. 
      This is an internal grid variable and is included in the count of \ref Grid::nNumIntVars. This
      variable is defined at cell centers.
      */
    int nSinThetaIJK;/**<
      Index of \f$\sin\theta\f$ defined at zone center in grids, \ref Grid::dLocalGridOld and
      \ref Grid::dLocalGridNew. This is an internal grid variable and is included in the count of
      \ref Grid::nNumIntVars.
      */
    int nSinThetaIJp1halfK;/**<
      Index of \f$\sin\theta\f$ at \f$\theta\f$ interfaces in grids. This is an internal grid 
      variable and is included in the count of \ref Grid::nNumIntVars.
      */
    int nCotThetaIJp1halfK;/**<
      Index of \f$\cot\theta\f$ at \f$\theta\f$ interfaces in grids. This is an internal grid 
      variable and is included in the count of \ref Grid::nNumIntVars.
      */
    int nCotThetaIJK;/**<
      Index of \f$\cot\theta\f$ at cell centeres of grids. This is an internal grid variable and is
      included in the count of \ref Grid::nNumIntVars.
      */
    int nDCosThetaIJK;/**<
      Index of \f$\Delta\cos\theta\f$ defined at zone center in grids. This is an internal grid 
      variable and is included in the count of \ref Grid::nNumIntVars.
      */
    int nEddyVisc;/**<
      Index of the eddy viscosity in the grid, it is defined at zone centers in the grids. This is 
      an internal grid  variable and is included in the count of \ref Grid::nNumIntVars.
    */
    int nDonorCellFrac;/**<
      Index of the amount of donor cell to use at that particular radial zone. It is defined at zone
      centers, and is an internal grid variable and is included in the count of
      \ref Grid::nNumIntVars.
      */
    int nNumDims; /**<
      Number of dimensions of the grid. It is used to chose the appropriate conservation equations.
      The value of this variable is independent of processor \ref ProcTop::nRank.
      */
    int nNumVars; /**<
      Number of grid variables. 
      This is set when reading in the model input file in the function \ref modelRead. It is the 
      number of variables that are printed and read from a file. The total number of variables also
      includes \ref Grid::nNumIntVars. The value of this variable is independent of processor 
      \ref ProcTop::nRank.
      */
    int nNumIntVars; /**<
      Number of internal variables.
      Internal variables are variables which are not reported in model dumps, and are not required
      to fully specify a starting model. They are used to save important information required during
      computation, an example is \f$\sin\theta\f$. The value of this variable is independent of 
      processor \ref ProcTop::nRank. This variable is set depending on the model read in 
      (adiabatic/non-adiabatic/number of dimensions) in the function \ref modelRead located in the 
      file \ref dataManipulation.cpp.
      */
    int nNum1DZones; /**<
      Number of zones in 1D region of grid.
      The number of zones in 3D region is (\ref Grid::nGlobalGridDims[0]- Grid::nNum1DZones ). 
      This is set when reading in the model input file in the function \ref modelRead. The value of 
      this variable is independent of processor \ref ProcTop::nRank.
      */
    int nNumGhostCells; /**<
      Number of cells which are not included in local grid updating. 
      This number is used in all dimensions to add to each local grid. When variables are not 
      defined in a given direction ghost cells are not included in that direction. This is set when
      reading in the model input file in the function \ref modelRead. The value of this variable is 
      independent of processor \ref ProcTop::nRank.
      */
    int *nGlobalGridDims; /**<
      Size of the entire global grid. It is an array of size 3 to hold size of each dimension of 
      global grid. This size does not include \ref Grid::nNumGhostCells or the extra size required 
      for interface centered quantities. The values of this variable are independent of processor 
      \ref ProcTop::nRank. In the case of 1D or 2D calculations the \f$\theta\f$ and \f$\phi\f$ 
      dimensions are set to 1 or just the \f$\phi\f$ dimensions is set to 1 depending on the 
      number of dimensions. The \f$r\f$, \f$\theta\f$ and \f$\phi\f$ dimensions are in the 0, 1 
      and 2 indices of the array respectively.
      */
    int **nVariables; /**<
      Provides information on grid variables.
      A 2D array of size \ref Grid::nNumVars+\ref Grid::nNumIntVars by 3+1.
      <tt> nVariables[n][l]</tt> has values:
      
        - -1: indicating that variable <tt>n</tt> is not defined 
        - 0: indicating that variable <tt>n</tt> is zone centered quantity
        - 1: indicating that variable <tt>n</tt> is an interface centered quantity
      
      in directions <tt>l=0,1,2</tt> which corresponding to \f$\hat{r}\f$, \f$\hat{\theta}\f$, and
      \f$\hat{\phi}\f$ respectively. <tt> nVariables[n][l]</tt> with \c l=3 is used to 
      indicate if a variable is dependent on time (1) or not(0). The values of this variable are
      independent of processor \ref ProcTop::nRank.
      */
    int ***nLocalGridDims; /**<
      Local grid dimensions. It is An array of size \ref ProcTop::nNumProcs by \ref 
      Grid::nNumVars+\ref Grid::nNumIntVars by 3. \c nLocalGridDims[p][n][l] gives the dimension 
      of the local grid on processor \c p for variable \c n in direction \c l. This variable does 
      not include \ref Grid::nNumGhostCells. The values of this variable are independent of 
      processor \ref ProcTop::nRank.
      */
    double ****dLocalGridNew; /**<
      Updated local grid values.
      An array of size \ref Grid::nNumVars+\ref Grid::nNumIntVars by \ref Grid::nLocalGridDims[0]
      +2*\ref Grid::nNumGhostCells by \ref Grid::nLocalGridDims[1]+2*\ref Grid::nNumGhostCells by
      \ref Grid::nLocalGridDims[2]+2*\ref Grid::nNumGhostCells provided that the variable is defined
      in all 3 directions. Variables that are not defined in all 3 directions will have the 
      additional two ghost cells left out in that direction and will also have a dimension of size 1
      in that direction. This array contains the current grid state as it is being updated through 
      calculations. This is a processor dependent variable and contains only the local grid for the
      current processor plus ghost cells.
      */
    double ****dLocalGridOld; /**<
      Grid values from previous time step.
      An array the same size as \ref Grid::dLocalGridNew but instead of containing the current grid 
      state, it contains the last complete grid state. This is a processor dependent variable and
      contains only the local grid for the current processor plus ghost cells.
      */
    int **nStartUpdateExplicit; /**<
      Positions to begin updating grid with explicit calculations. It is an array of size 
      \ref nNumVars+\ref nNumIntVars by 3. The start positions are defined in 
      \ref initUpdateLocalBoundaries(). These start values are dependent on processor
      \ref ProcTop::nRank.
      */
    int **nEndUpdateExplicit; /**<
      Positions to stop updating grid with explicit calculations. It is an array of size 
      \ref nNumVars+\ref nNumIntVars by 3. The end positions are defined in 
      \ref initUpdateLocalBoundaries(). These start values are dependent on processor
      \ref ProcTop::nRank.
      */
    int **nStartUpdateImplicit; /**<
      Positions to begin updating grid with implicit calculations. It is an array of size 
      \ref nNumVars+\ref nNumIntVars by 3. The start positions are defined in 
      \ref initUpdateLocalBoundaries(). These start values are dependent on processor
      \ref ProcTop::nRank.
      */
    int **nEndUpdateImplicit; /**<
      Positions to stop updating grid with implicit calculations. It is an array of size 
      \ref nNumVars+\ref nNumIntVars by 3. The end positions are defined in 
      \ref initUpdateLocalBoundaries(). These start values are dependent on processor
      \ref ProcTop::nRank.
      */
    int ***nStartGhostUpdateExplicit; /**<
      Positions to begin updating ghost cells with explicit calculations. It is an array of size 
      \ref Grid::nNumVars+\ref Grid::nNumIntVars by 2*3 by 3. The second dimension indicates a 
      particular ghost region. There are 2*3 since each direction can have two ghost regions. The 
      ghost region 0, is the outter ghost region in direction 0, 1 is the inner ghost region in 
      direction 0, etc.
      */
    int ***nEndGhostUpdateExplicit; /**<
      Positions to end updating ghost cells with explicit calculatiosn. Is an array of size \ref 
      Grid::nNumVars+\ref Grid::nNumIntVars by 2*3 by 3. The second dimension corresponds to which 
      ghost region, since each dimension can have two ghost regions. The ghost region 0, is the 
      outter ghost region in direction 0, 1 is the inner ghost region in direction 0, etc.
      */
    int ***nStartGhostUpdateImplicit; /**<
      Positions to begin updating ghost cells with implicit calculations. It is an array of size 
      \ref Grid::nNumVars+\ref Grid::nNumIntVars by 2*3 by 3. The second dimension indicates a 
      particular ghost region. There are 2*3 since each direction can have two ghost regions. The 
      ghost region 0, is the outter ghost region in direction 0, 1 is the inner ghost region in 
      direction 0, etc.
      */
    int ***nEndGhostUpdateImplicit; /**<
      Positions to end updating ghost cells with implicit calculations. Is an array of size \ref 
      Grid::nNumVars+\ref Grid::nNumIntVars by 2*3 by 3. The second dimension corresponds to which 
      ghost region, since each dimension can have two ghost regions. The ghost region 0, is the 
      outter ghost region in direction 0, 1 is the inner ghost region in direction 0, etc.
      */
    int *nCenIntOffset;/**<
      Indicates the offset between interface and center quantities.
      If <tt>nCenIntOffset[l]=0</tt> then the outter interface quantities have the same index as 
      zone centered quantities in direction <tt>l</tt>. If <tt>nCenIntOffset[l]=1</tt> then the 
      outter interface quantities are given by the index for the zone centered quantities +1, 
      in direction <tt>l</tt>. The values are dependent on \ref ProcTop::nRank and 
      \ref ProcTop::nPeriodic.
      */
    int nGlobalGridPositionLocalGrid[3];/**<
      The location at which the local grid starts in the global grid. This starts at 0, for the 
      inner most cell, including ghost zones.
      */
    int nNumZones1DBoundaryZeroHorizontalVelocity;/**
      sets how many zones out from the 1D-multi-D boundary that theta/phi velocities are not updated
      and thus kept at zero.
      */
    Grid(); /**<
      Constructor for the class \ref Grid.
      */
};/**@class Grid
  This class manages information which pertains to grid data.
  
  External variables used with Gamma Law (GL) gas equaiton of state and their array indexes:
  <table>
    <tr><th>1D (\ref nNumVars=7)</th><th>2D (\ref nNumVars=9)</th><th>3D (\ref nNumVars=11)</th></tr>
    <tr><td>
      <table>
      <tr><td>Variable</td><td>Index</td></tr>
      <tr><td>\ref nM</td><td>0</td></tr>
      <tr><td>\ref nDM</td><td>1</td></tr>
      <tr><td>\ref nR</td><td>2</td></tr>
      <tr><td>\ref nD</td><td>3</td></tr>
      <tr><td>\ref nU</td><td>4</td></tr>
      <tr><td>\ref nU0</td><td>5</td></tr>
      <tr><td>\ref nE</td><td>6</td></tr>
      </table>
    </td><td>
      <table>
      <tr><td>Variable</td><td>Index</td></tr>
      <tr><td>\ref nM</td><td>0</td></tr>
      <tr><td>\ref nTheta</td><td>1</td></tr>
      <tr><td>\ref nDM</td><td>2</td></tr>
      <tr><td>\ref nR</td><td>3</td></tr>
      <tr><td>\ref nD</td><td>4</td></tr>
      <tr><td>\ref nU</td><td>5</td></tr>
      <tr><td>\ref nU0</td><td>6</td></tr>
      <tr><td>\ref nV</td><td>7</td></tr>
      <tr><td>\ref nE</td><td>8</td></tr>
      </table>
    </td><td>
      <table>
      <tr><td>Variable</td><td>Index</td></tr>
      <tr><td>\ref nM </td><td>0</td></tr>
      <tr><td>\ref nTheta</td><td>1</td></tr>
      <tr><td>\ref nPhi</td><td>2</td></tr>
      <tr><td>\ref nDM</td><td>3</td></tr>
      <tr><td>\ref nR</td><td>4</td></tr>
      <tr><td>\ref nD</td><td>5</td></tr>
      <tr><td>\ref nU</td><td>6</td></tr>
      <tr><td>\ref nU0</td><td>7</td></tr>
      <tr><td>\ref nV</td><td>8</td></tr>
      <tr><td>\ref nW</td><td>9</td></tr>
      <tr><td>\ref nE</td><td>10</td></tr>
      </table>
    </td>
    </tr>
  </table>
  
  External variables used with Tabulated Equation Of State (TEOS) and their array indexes:
  <table>
    <tr><th>1D (\ref nNumVars=7)</th><th>2D (\ref nNumVars=9)</th><th>3D (\ref nNumVars=11)</th></tr>
    <tr><td>
      <table>
      <tr><td>Variable</td><td>Index</td></tr>
      <tr><td>\ref nM </td><td>0</td></tr>
      <tr><td>\ref nDM </td><td>1</td></tr>
      <tr><td>\ref nR </td><td>2</td></tr>
      <tr><td>\ref nD </td><td>3</td></tr>
      <tr><td>\ref nU </td><td>4</td></tr>
      <tr><td>\ref nU0 </td><td>5</td></tr>
      <tr><td>\ref nT </td><td>6</td></tr>
      </table>
    </td><td>
      <table>
      <tr><td>Variable</td><td>Index</td></tr>
      <tr><td>\ref nM </td><td>0</td></tr>
      <tr><td>\ref nTheta </td><td>1</td></tr>
      <tr><td>\ref nDM </td><td>2</td></tr>
      <tr><td>\ref nR </td><td>3</td></tr>
      <tr><td>\ref nD </td><td>4</td></tr>
      <tr><td>\ref nU </td><td>5</td></tr>
      <tr><td>\ref nU0 </td><td>6</td></tr>
      <tr><td>\ref nV </td><td>7</td></tr>
      <tr><td>\ref nT </td><td>8</td></tr>
      </table>
    </td><td>
      <table>
      <tr><td>Variable</td><td>Index</td></tr>
      <tr><td>\ref nM </td><td>0</td></tr>
      <tr><td>\ref nTheta </td><td>1</td></tr>
      <tr><td>\ref nPhi </td><td>2</td></tr>
      <tr><td>\ref nDM </td><td>3</td></tr>
      <tr><td>\ref nR </td><td>4</td></tr>
      <tr><td>\ref nD </td><td>5</td></tr>
      <tr><td>\ref nU </td><td>6</td></tr>
      <tr><td>\ref nU0 </td><td>7</td></tr>
      <tr><td>\ref nV </td><td>8</td></tr>
      <tr><td>\ref nW </td><td>9</td></tr>
      <tr><td>\ref nT </td><td>10</td></tr>
      </table>
    </td>
    </tr>
  </table>
  
  Internal variables  with GL gas equation of state:
  <table>
    <tr><th>1D (\ref nNumIntVars=2)</th><th>2D (\ref nNumIntVars=8)</th></tr>
    <tr><td>
      <table>
      <tr><td>Variable</td><td>Index</td></tr>
      <tr><td>\ref nP </td><td>\ref nNumVars+0</td></tr>
      <tr><td>\ref nQ0 </td><td>\ref nNumVars+1</td></tr>
      </table>
    </td><td>
      <table>
      <tr><td>Variable</td><td>Index</td></tr>
      <tr><td>\ref nP </td><td>\ref nNumVars+0</td></tr>
      <tr><td>\ref nQ0 </td><td>\ref nNumVars+1</td></tr>
      <tr><td>\ref nDenAve </td><td>\ref nNumVars+2</td></tr>
      <tr><td>\ref nDCosThetaIJK </td><td>\ref nNumVars+3</td></tr>
      <tr><td>\ref nQ1</td><td>\ref nNumVars+4</td></tr>
      <tr><td>\ref nDTheta</td><td>\ref nNumVars+5</td></tr>
      <tr><td>\ref nSinThetaIJK</td><td>\ref nNumVars+6</td></tr>
      <tr><td>\ref nSinThetaIJp1halfK</td><td>\ref nNumVars+7</td></tr>
      </table></td>
    </tr>
    <tr><th>3D (\ref nNumIntVars=12)</th></tr>
    <tr>
    <td>
      <table>
      <tr><td>Variable</td><td>Index</td></tr>
      <tr><td>\ref nP </td><td>\ref nNumVars+0</td></tr>
      <tr><td>\ref nQ0 </td><td>\ref nNumVars+1</td></tr>
      <tr><td>\ref nDenAve </td><td>\ref nNumVars+2</td></tr>
      <tr><td>\ref nDPhi </td><td>\ref nNumVars+3</td></tr>
      <tr><td>\ref nDCosThetaIJK </td><td>\ref nNumVars+4</td></tr>
      <tr><td>\ref nQ1</td><td>\ref nNumVars+5</td></tr>
      <tr><td>\ref nDTheta</td><td>\ref nNumVars+6</td></tr>
      <tr><td>\ref nSinThetaIJK</td><td>\ref nNumVars+7</td></tr>
      <tr><td>\ref nSinThetaIJp1halfK</td><td>\ref nNumVars+8</td></tr>
      <tr><td>\ref nCotThetaIJK</td><td>\ref nNumVars+9</td></tr>
      <tr><td>\ref nCotThetaIJp1halfK</td><td>\ref nNumVars+10</td></tr>
      <tr><td>\ref nQ2</td><td>\ref nNumVars+11</td></tr>
      </table>
    </td>
    </tr>
  </table>
  
  Internal variables  with TEOS:
  <table celpadding="0">
    <tr><th>1D (\ref nNumIntVars=5)</th><th>2D (\ref nNumIntVars=11)</th></tr>
    <tr><td>
      <table>
      <tr><td>Variable</td><td>Index</td></tr>
      <tr><td>\ref nP </td><td>\ref nNumVars+0</td></tr>
      <tr><td>\ref nQ0 </td><td>\ref nNumVars+1</td></tr>
      <tr><td>\ref nE </td><td>\ref nNumVars+2</td></tr>
      <tr><td>\ref nKappa </td><td>\ref nNumVars+3</td></tr>
      <tr><td>\ref nGamma </td><td>\ref nNumVars+4</td></tr>
      </table>
    </td><td>
      <table>
      <tr><td>Variable</td><td>Index</td></tr>
      <tr><td>\ref nP </td><td>\ref nNumVars+0</td></tr>
      <tr><td>\ref nQ0 </td><td>\ref nNumVars+1</td></tr>
      <tr><td>\ref nDenAve </td><td>\ref nNumVars+2</td></tr>
      <tr><td>\ref nDCosThetaIJK </td><td>\ref nNumVars+3</td></tr>
      <tr><td>\ref nE </td><td>\ref nNumVars+4</td></tr>
      <tr><td>\ref nKappa </td><td>\ref nNumVars+5</td></tr>
      <tr><td>\ref nGamma </td><td>\ref nNumVars+6</td></tr>
      <tr><td>\ref nQ1</td><td>\ref nNumVars+7</td></tr>
      <tr><td>\ref nDTheta</td><td>\ref nNumVars+8</td></tr>
      <tr><td>\ref nSinThetaIJK</td><td>\ref nNumVars+9</td></tr>
      <tr><td>\ref nSinThetaIJp1halfK</td><td>\ref nNumVars+10</td></tr>
      </table></td>
    </tr>
    <tr><th>3D (\ref nNumIntVars=15)</th></tr>
    <tr>
    <td>
      <table>
      <tr><td>Variable</td><td>Index</td></tr>
      <tr><td>\ref nP </td><td>\ref nNumVars+0</td></tr>
      <tr><td>\ref nQ0 </td><td>\ref nNumVars+1</td></tr>
      <tr><td>\ref nDenAve </td><td>\ref nNumVars+2</td></tr>
      <tr><td>\ref nDPhi </td><td>\ref nNumVars+3</td></tr>
      <tr><td>\ref nDCosThetaIJK </td><td>\ref nNumVars+4</td></tr>
      <tr><td>\ref nE </td><td>\ref nNumVars+5</td></tr>
      <tr><td>\ref nKappa </td><td>\ref nNumVars+6</td></tr>
      <tr><td>\ref nGamma </td><td>\ref nNumVars+7</td></tr>
      <tr><td>\ref nQ1</td><td>\ref nNumVars+8</td></tr>
      <tr><td>\ref nDTheta</td><td>\ref nNumVars+9</td></tr>
      <tr><td>\ref nSinThetaIJK</td><td>\ref nNumVars+10</td></tr>
      <tr><td>\ref nSinThetaIJp1halfK</td><td>\ref nNumVars+11</td></tr>
      <tr><td>\ref nCotThetaIJK</td><td>\ref nNumVars+12</td></tr>
      <tr><td>\ref nCotThetaIJp1halfK</td><td>\ref nNumVars+13</td></tr>
      <tr><td>\ref nQ2</td><td>\ref nNumVars+14</td></tr>
      </table>
    </td>
    </tr>
  </table>
  
  The variable indexes are set in \ref modelRead based on the input model.
  */
class Parameters{
  public:
    bool bEOSGammaLaw;/**<
      If true SPHERLS will use a gamma law gas instead of a 
      tabulated equation of  state. This is set in the starting model.
      */
    bool bAdiabatic;/**<
      If true SPHERLS will use adiabatic functions to calculate the energy.
      This can be used for both gamma law gas and tabulated equations of state
      (see \ref Parameters::bEOSGammaLaw).
      */
    int nTypeTurbulanceMod;/**<
      This varible indicates the type of turbulance model to be used. If 0, no turbulance model will
      be used, if 1 it will use a constant times the zoning size, and if 2 it will use the 
      Smagorinksy turbulance model which increases the value of the eddy viscosity parameter when
      there are large amounts of shear, and decrease it when there isn't.
      */
    double dPi; /**<
      The value of \f$\pi\f$.
      */
    double dSigma; /**<
      The value of \f$\sigma\f$, the Stefan-Boltzmann constant.
      */
    double dG; /**<
      The Gravitational constant \f$G\f$.
      */
    double dGamma; /**<
      The adiabatic \f$\gamma\f$, used in calculating the equation of state. If using a gamma law
      gas.
      */
    std::string sEOSFileName;/**<
      File name of equation of state table. This value is set either by the configuration file,
      SPHERLS.xml or in the model file read in. If it is specified in SPHERLS.xml it will overide
      the file name set in the model.
      */
    eos eosTable;/**<
      Holds the equation of state table. If using a tabulated equation of state.
      */
    double dA; /**<
      Artificial viscosity parameter, reasonable values range from 0 to ~3.
      */
    double dAVThreshold; /**<
      The amount of compression before AV is turned on. It is in terms of a
      velocity difference between zone sides and is in fractions of the local sound speed.*/
    double dDonorCellMultiplier; /**<
      Multiplier used to determine the faction of the sound speed at 
      which donor cell is full. e.g. a value of 1.0 means the donor cell will be full when the 
      convective velocity is equal to the sound speed. A value of 0.5 will mean that it will be full
      donor cell when the convective velocity is twice the sound speed. A value of 2.0 will mean that 
      it will use full donor cell when the convective velocity is half the sound speed.
      */
    double dDonorCellMin; /**<
      The minimum amount of donor cell allowed. Set in constructor, \ref Parameters::Parameters
      */
    double dAlpha; /**<
      This parameter controls the amount of extra mass above the outter interface.
      it is read in from the starting model, so that it will be consistent with the value used
      in calculating the starting model.
      */
    double dAlphaExtra;/**<
      */
    double dTolerance; /**<
      Amount of error to tolerate when calculating temperature from the equation of state.
      */
    int nMaxIterations; /**<
      The maximum number of iterations to try to get the the relative error in the temperture below
      \ref parameters::dTolerance.*/
    double dEddyViscosity;/**<
      Used in calculating the eddy viscosity, larger values will produce a larger value of the eddy
      viscosity, causing the rethermalization to happen at larger scales. This value should be 
      kept small, a good value is 0.17, which seems to correspond with experiments.
      */
    double dMaxConvectiveVelocity;/**<
      Holds the maximum convective velocity, it is set in the functions which calculate the 
      timestep (see \ref calDelt_R_GL, \ref calDelt_R_TEOS, \ref calDelt_RT_GL,
      \ref calDelt_RT_TEOS, \ref calDelt_RTP_GL, \ref calDelt_RTP_TEOS, \ref calDelt_CONST).
      */
    double dMaxConvectiveVelocity_c;/**<
      Holds the maximum of convective velocity divided by the sound speed. It is set in the 
      functions which calculate the timestep (see \ref calDelt_R_GL, \ref calDelt_R_TEOS,
      \ref calDelt_RT_GL, \ref calDelt_RT_TEOS, \ref calDelt_RTP_GL, \ref calDelt_RTP_TEOS,
      \ref calDelt_CONST).
      */
    double dPrt;/**<
      This is the value of the Prandtl number, a value of 0.7 is what is suggested by Lawrence D. 
      Cloutman in "The LUVD11 Large Eddy Simulation Model" April 15, 1991 a Lawrence Livermore 
      National Labratory report.
      */
    double dDEDMClampValue;/**<
      The value to use for DEDM in energy conservation equation when \ref Parameters::bDEDMClamp is
      true.
    */
    double dDEDMClampMr;/**<
      The mass above which the DEDM clamp is applied.
    */
    double dEDMClampTemperature;/**<
      The temperature at which to chose \ref Parameters::dDEDMClampMr from the stating model.
    */
    bool bDEDMClamp;/**<
      Specifies if a DEDM clamp should be used. This should only be used when starting from a model
      with out any sizable convection. It could give undesirable results if used when starting a
      calculation from a model with already established convection.
    */
    std::string sDebugProfileOutput;/**<
      output file name for debuging profile, only used if DEBUG_EQUATIONS is set to 1
    */
    
    #if DEBUG_EQUATIONS==1
    profileData profileDataDebug;/**<
      tracks and writes out profile data useful for debugging
      */
    bool bSetThisCall;/**<
      if true will call set function, if not it won't
      */
    bool bEveryJK;/**<
      if true every JK will have their own values set
      */
    #endif
    Parameters(); /**<
      Constructor for the class \ref Parameters
      */
};/**@class Parameters
  This class holds parameters and constants used for calculation.
  */
class Output{
  public:
    int nDumpFrequencyStep; /**<
      How ofter a the grid state should be written to a file according to time step index. If it is 
      1 the will state will be written every time step, if it equals 2 it will be written every 
      other time step etc. If it is 0 no dumps will be made according to the time step index.
      */
    double dDumpFrequencyTime; /**<
      How ofter a the grid state should be written to a file according to simulation time in 
      seconds. If it is 0 no dumps will be made according to simulation time.
      */
    double dTimeLastDump;/**<
      The simulation time at which the last dump was made using the \ref Output::dDumpFrequencyTime 
      criterion.
      */
    int nNumTimeStepsSinceLastDump;/**
      The number of time steps since the last model dump.
    */
    int nNumTimeStepsSinceLastPrint;/**
      The number of time steps since the last print.
    */
    bool bDump; /**<
      Should the grid state be written to a file at a frequency of \ref Output::nDumpFrequencyStep 
      timesteps, and/or every \ref Output::dDumpFrequencyTime seconds of simulation time. This is
      set to true by putting a "<dump>" node into the "SPHERLS.xml" configuration file.
      */
    bool bPrint;/**<
      Should status updates be printed to the screen.
    */
    int nPrintMode;/**<
      Sets the way in which information should be printed to the standard output during the run. If
      it is 0, it will print the standard information reporting on the progress of the code. If it 
      is 1 it will print out information to diagnose timestepping problems.
    */
    std::string sBaseOutputFileName; /**<
      Base filename used for output,
      default is "out". All model dumps, and output information will contain this
      file name and extend it to indicate their specific information. The value of 
      this variable is independent of processor \ref ProcTop::nRank.
      */
    std::ofstream *ofWatchZoneFiles; /**<
      An array of output streams of size 
      \ref Output::watchzoneList .size() which are used to write out the information of the watched
      zones.
      */
    std::vector<WatchZone> watchzoneList; /**<
      A vector used to keep information used to specify the zones to be watched.
      */
    int nPrintFrequencyStep;/**<
      How often the status is printed to the screen in time steps.*/
    double dPrintFrequencyTime;/**<
      How often the status is printed to the screen in simulation time.*/
    double dTimeLastPrint;/**<
      Simulation time when last status was printed.*/
    std::string sExeDir;/**<
      Directory where the executable is located.
      */
    void setExeDir(ProcTop &procTop);/**<
      Sets sExeDir to the directory where the current executable is located
      */
    Output(); /**<
      Constructor for this class.
      */
};/**@class Output
  This class manages information pertianing to the output of data to files.
  */
class Performance{
  public:
    double dStartTimer; /**<
      The time that the code timer was started.
      */
    double dEndTimer; /**<
      The time that the code timer was ended. The difference between
      \ref Performance::dStartTimer and \c dEndTimer gives the total run time
      */
    Performance(); /**<
      Constructor for the class \ref Performance.
      */
};/**@class Performance
  This class manages information pertianing to performace analysis of the code.
  */
class Implicit{
  public:
    int nNumImplicitZones; /**<
      The number of zones in the region near the surface which should used the implicit calculation
      of the energy equation. If zero no zones will use the implict calculation of energy.
      */
    Mat matCoeff;/**<
      Parallel coeffecient matrix (spread across all processors)
      */
    Vec vecTCorrections;/**<
      Temperature corrections solution vector (spread across all processors)
      */
    Vec vecRHS;/**<
      RHS vector (spread across all processors)
      */
    Vec vecTCorrectionsLocal;/**<
      Corrections to local temperatures only (on local processor only).
      */
    KSP kspContext;/**<
      PETSc solver context.
      */
    VecScatter vecscatTCorrections;/**<
      Scatter context, used to hold information about retrieving the distributed temperature 
      corrections from vecTCorrections and placing them into the local vector vecTCorrectionsLocal.
      */
    int nMaxNumIterations;/**<
      The maximum number of iterations to try to get the largest value of 
      \ref vecTCorrections relative to the temperature below \ref dTolerance. Ater which the 
      calculation continues.
      */
    double dTolerance;/**<
      The amount of relative error that is allowed in the calculation of the temperature with the
      implicit calculation.
      */
    int nNumRowsALocal;/**<
      The number of rows of the coeffecient matrix which is on the local processor.
      */
    int nNumRowsALocalSB;/**<
      The number or rows of the coeffecient matrix which is on the local processor, and that are in
      the surface boundary region.
      */
    int *nNumDerPerRow;/**<
      An array of size \ref nNumRowsALocal which contains the number of non-zero derivatives for a 
      given row of A.
      */
    int **nTypeDer;/**<
      An array of size \ref nNumRowsALocal by \ref nNumDerPerRow \c [q] , where \c q is a row 
      index. Thus each row of the array can have a different length. This gives the type of 
      derivative of row \c q for each derivative in that row. The value of this variable is set in 
      the function \ref initImplicitCalculation .
      */
    int ***nLocDer;/**<
      An array of size \ref nNumRowsALocal by 2 by \ref nNumDerPerRow \c [q] , where \c q is a row 
      index. This array holds the global position of the current row \c q for the current derivative
      e.g. the \c p th derivative in the \c q th row would be in row and column (\c nLocDer[q][0][p]
      ,\c nLocDer[q][1][p]). The value of this variable is set in the function 
      \ref initImplicitCalculation .
      */
    int **nLocFun;/**<
      An array of size \ref nNumRowsALocal by 3 \c [q] , where \c q is a row index. This array holds
      the local grid position of the current row \c q e.g. the (i,j,k) location of the the current
      row in the local grid. The value of this variable is set in the function 
      \ref initImplicitCalculation .
      */
    double dDerivativeStepFraction;/**<
      Dicates the size of the step that should be used to evaluate the numerical derivitves of the 
      energy equation, for solving for the temperature implicitily. This value multiplies the
      temperature to produce the step size. A good value is around 5e-7.
      */
    double dCurrentRelTError;/**<
      keeps track of the largest relative error in the calculation of the temperature
    */
    int nCurrentNumIterations;/**<
      keeps track of the number of iterations needed to converge to a solution
    */
    int nMaxNumSolverIterations;/**<
      If \ref TRACKMAXSOLVERERROR set to 1, then this will be the current maximum number of 
      iterations required for the linear equaiton solver to solve for the temperature correction
      over all iterations and time steps since the last model dump.
      */
    double dMaxErrorInRHS;/**<
      If \ref TRACKMAXSOLVERERROR set to 1, then this will be the current maximum absolute error 
      between the RHS as calculated from the solution and the coeffecient matrix, and the actual
      RHS. This value is the maximum from all values at each iteration of the solution, from each
      time step since the last model dump.
      */
    double dAverageRHS;/**<
      Holds the average value of the right hand side for the timestep where the error in the RHS is 
      the largest \ref dMaxErrorInRHS. Only set if \ref TRACKMAXSOLVERERROR is set to 1.
    */
    Implicit();/**<
      constructor the the class \ref Implicit.
      */
};/**@class Implicit
  This class holds data required for the implicit calculation.
  */
class Functions{
  public:
    void (*fpCalculateNewVelocities)(Grid&, Parameters&, Time&, ProcTop&); /**<
      Function pointer to the function used to calculate new velocities.
      */
    void (*fpCalculateNewGridVelocities)(Grid&, Parameters&, Time&, ProcTop&, MessPass&); /**<
      Function pointer to the function used to calculate new grid velocities.
      */
    void (*fpCalculateNewRadii)(Grid&, Time&); /**<
      Functin pointer to the function used to calculate new radii.
      */
    void (*fpCalculateNewDensities)(Grid&, Parameters&, Time&,ProcTop&); /**<
      Function pointer to the function used to calculate the new densities.
      */
    void (*fpCalculateNewEnergies)(Grid&, Parameters&, Time&, ProcTop&); /**<
      Function pointer to the function used to calculate the new energies.
      */
    void (*fpCalculateDeltat)(Grid&,Parameters&,Time&,ProcTop&); /**<
      Function pointer to the function used to calculate the new time step.
      */
    void (*fpCalculateAveDensities)(Grid&); /**<
      Function pointer to the function used to calculate the new average density.
      */
    void (*fpCalculateNewEOSVars)(Grid&, Parameters&);/**<
      Function pointer to the function used to calculate the new variables depending on the equation
      of state.
      */
    void (*fpCalculateNewAV)(Grid&, Parameters&);/**<
      Function pointer to the function used to calculate new Artificial viscosity.
      */
    void (*fpModelWrite)(std::string sFileName, ProcTop&, Grid&, Time&, Parameters&);/**<
      Function pointer to the function used to write out models.
      */
    void (*fpWriteWatchZones)(Output&, Grid&,Parameters&,Time&,ProcTop&);/**<
      Function pointer to the function that is used to write out watch zone files
      */
    void (*fpUpdateLocalBoundaryVelocitiesNewGrid)(ProcTop&, MessPass&,Grid&);/**<
      Function pointer to the fnction that is used to update velocities across boundaries.*/
    void (*fpImplicitSolve)(Grid&,Implicit&,Parameters&,Time&,ProcTop&,MessPass&,Functions&);/**<
      Funciton pointer to the function that is used to implicitly solve for the temperature,
      then uses the equation of state to solve for energy, opacity, and pressure.
      */
    void (*fpCalculateNewEddyVisc)(Grid&,Parameters&);/**<
      Function pointer to the function that is used to calculate the new eddy viscosity.
      */
    double (*fpImplicitEnergyFunction)(Grid&,Parameters&,Time&,double[],int,int,int);
    double (*fpImplicitEnergyFunction_SB)(Grid&,Parameters&,Time&,double[],int,int,int);
    Functions(); /**<
      Constructor for the class \ref Functions.
      */
};/**@class Functions
  This class holds function pointers used to indicate the functions which should be used to 
  calculate the various needed quantities. These functions can be different from processor to
  processor. For example \ref ProcTop::nRank=0 processor will have only 1D verions of the
  conservation equations, while the rest of the processors will have 3D versions. These functions
  will also change depending on what kind of model is being calculated and the number of dimensions
  the calculation includes.
  */
class Global{
  public:
    ProcTop procTop; /**<
      An instance of the \ref ProcTop class.
      */
    MessPass messPass; /**<
      An instance of the \ref MessPass class.
      */
    Grid grid; /**<
      An instance of the \ref Grid class.
      */
    Time time; /**<
      An instance of the \ref Time class.
      */
    Parameters parameters; /**<
      An instance of the \ref Parameters class.
      */
    Output output; /**<
      An instance of the \ref Output class.
      */
    Performance performance; /**<
      An instance of the \ref Performance class.
      */
    Functions functions; /**<
      An instance of the \ref Functions class.
      */
    Implicit implicit;/**<
      An instance of the \ref Implicit class.
      */
    Global(); /**<
      Constructor for the class \ref Global.
      */
};/**@class Global
  This class is simply a class that holds the other classes.
  */
#endif
