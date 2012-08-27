#ifndef PHYSEQUATIONS_H
#define PHYSEQUATIONS_H

/** @file
  
  Header file for \ref physEquations.cpp
*/

#include "global.h"

void setMainFunctions(Functions& functions,ProcTop &procTop,Parameters &parameters, Grid &grid
  , Time &time,Implicit &implicit);/**<
  Used to set the functions that \ref main() uses to evolve the input model.
  
  @param[out] functions is of class \ref Functions and is used to specify the functions
    called to calculate the evolution of the input model.
  @param[in] procTop is of type \ref ProcTop. \ref ProcTop::nRank is used to set different
    functions based on processor rank. For instance processor rank 1 requires 1D 
    versions of the equations.
  @param[in] parameters is of class \ref Parameters. It holds various constants and runtime 
    parameters.
  @param[in] grid of type \ref Grid. This function requires the number of dimensions, specified by
    \ref Grid::nNumDims.
  @param[in] time of type \ref Time. This function requires knowledge of the type of time setp being
    used, specified by \ref Time::bVariableTimeStep.
  @param[in] implicit of type \ref Implicit. This function needs to know if there is an implicit
    region, specified when Implicit::nNumImplicitZones>0.

  The functions are picked based on model geometry, and the physics requested or required by the 
  input model, and the configuration file. The specific functions pointers that are set are 
  described in the \ref Functions class.
  
  */
void setInternalVarInf(Grid& grid, Parameters &parameters);/**<
  This function sets the information for internal variables. While external verabile information is 
  derived from the starting model, internal variables infos are set in this function. In other words
  this function sets the values of \ref Grid::nVariables.
  
  @param[in,out] grid supplies the information needed to calculate the horizontal density average,
                 it also stores the calculated horizontally averaged density.
  @param[in] parameters is used when setting variable infos, since one needs to know if the code
             is calculating using a gamma law gas, or a tabulated equation of state.
  */
void initInternalVars(Grid& grid, ProcTop &procTop, Parameters &parameters);/**<
  This function function is used to set the initial values of the internal variables. While external
  variables are initialized from the starting model, internal variables are calculated at startup.
  
  @param[in,out] grid supplies information needed for initilizing internal variables as well as
                 storing the initilized internal variables
  @param[in] procTop contians information about processor topology
  @param[in] parameters contains parameters used in initializing the internal variables.
  */
void calNewVelocities_R(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function simply calls a function that calculate the radial velocity. Calls the function 
  \ref calNewU_R to calculate radial velocity, including only radial terms.
  
  @param[in,out] grid contains the local grid data and supplies the needed data to calculate the
  new velocities as well as holding the new velocities.
  
  @param[in] parameters contains parameters used in the calculation of the new velocities.
  @param[in] time contains time step information, current time step, and current time
  @param[in] procTop contains processor topology information
  */
void calNewVelocities_R_LES(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function simply calls a function that calculate the radial velocity. Calls the function 
  \ref calNewU_R to calculate radial velocity, including only radial terms.
  
  @param[in,out] grid contains the local grid data and supplies the needed data to calculate the
  new velocities as well as holding the new velocities.
  
  @param[in] parameters contains parameters used in the calculation of the new velocities.
  @param[in] time contains time step information, current time step, and current time
  @param[in] procTop contains processor topology information
  */
void calNewVelocities_RT(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function simply calls two other functions that calculate the radial and theta velocities.
  Calls the two functions \ref calNewU_RT and \ref calNewV_RT to calculate radial and 
  theta velocities, including both radial and theta terms.
  
  @param[in,out] grid contains the local grid data and supplies the needed data to calculate the
  new velocities as well as holding the new velocities.
  
  @param[in] parameters contains parameters used in the calculation of the new velocities.
  @param[in] time contains time step information, current time step, and current time
  @param[in] procTop contains processor topology information
  */
void calNewVelocities_RT_LES(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function simply calls two other functions that calculate the radial and theta velocities.
  Calls the two functions \ref calNewU_RT and \ref calNewV_RT to calculate radial and 
  theta velocities, including both radial and theta terms.
  
  @param[in,out] grid contains the local grid data and supplies the needed data to calculate the
  new velocities as well as holding the new velocities.
  
  @param[in] parameters contains parameters used in the calculation of the new velocities.
  @param[in] time contains time step information, current time step, and current time
  @param[in] procTop contains processor topology information
  */
void calNewVelocities_RTP(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function simply calls three other functions that calculate the radial, theta  and phi
  velocities. Calls the two functions \ref calNewU_RTP, \ref calNewV_RTP and
  \ref calNewW_RTP to calculate radial, theta, and phi velocities, including radial, theta,
  and phi terms.
  
  @param[in,out] grid contains the local grid data and supplies the needed data to calculate the
  new velocities as well as holding the new velocities.
  
  @param[in] parameters contains parameters used in the calculation of the new velocities.
  @param[in] time contains time step information, current time step, and current time
  @param[in] procTop contains processor topology information
  */
void calNewVelocities_RTP_LES(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function simply calls three other functions that calculate the radial, theta  and phi
  velocities. Calls the two functions \ref calNewU_RTP, \ref calNewV_RTP and
  \ref calNewW_RTP to calculate radial, theta, and phi velocities, including radial, theta,
  and phi terms.
  
  @param[in,out] grid contains the local grid data and supplies the needed data to calculate the
  new velocities as well as holding the new velocities.
  
  @param[in] parameters contains parameters used in the calculation of the new velocities.
  @param[in] time contains time step information, current time step, and current time
  @param[in] procTop contains processor topology information
  */
void calNewU_R(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function calculates the radial velocity, and does it by only considering the radial terms.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated radial velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
  */
void calNewU_R_LES(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function calculates the radial velocity, and does it by only considering the radial terms. 
  It also includes the terms for including real viscosity, used in the LES.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated radial velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
  */
void calNewU_RT(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function calculates the radial velocity, and does it by only considering the radial and
  theta terms.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated radial velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
  */
void calNewU_RT_LES(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function calculates the radial velocity, and does it by only considering the radial and
  theta terms. It also includes the terms for including real viscosity, used in the LES.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated radial velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
  */
void calNewU_RTP(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function calculates the radial velocity, and does it by including all radial, theta and
  phi terms.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated radial velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
*/
void calNewU_RTP_LES(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function calculates the radial velocity, and does it by including all radial, theta and
  phi terms. It also includes the terms for including real viscosity, used in the LES.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated radial velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
*/
void calNewV_RT(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function calculates the theta velocity, and does it by only considering the radial and
  theta terms.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated theta velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
  */
void calNewV_RT_LES(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function calculates the theta velocity, and does it by only considering the radial and
  theta terms. It also includes the terms for including real viscosity, used in the LES.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated theta velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
  */
void calNewV_RTP(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function calculates the theta velocity, and does it by considering the radial,
  theta, and phi terms.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated theta velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
  */
void calNewV_RTP_LES(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function calculates the theta velocity, and does it by considering the radial,
  theta, and phi terms. It also includes the terms for including real viscosity, used in the LES.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated theta velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
  */
void calNewW_RTP(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function calculates the phi velocity, and does it by only considering the radial, theta, and
  phi terms. 
  
  @param[in,out] grid contains the local grid, and will hold the newly updated theta velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
  */
void calNewW_RTP_LES(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop);/**<
  This function calculates the phi velocity, and does it by only considering the radial, theta, and
  phi terms. It also includes the terms for including real viscosity, used in the LES.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated theta velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
  */
void calNewU0_R(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop
  ,MessPass& messPass);/**<
  This function calculates the radial grid velocity, it does so by considering only the radial 
  terms
  
  @param[in,out] grid contains the local grid, and will hold the newly updated radial grid 
                 velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
  @param[in] messPass
  */
void calNewU0_RT(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop
  ,MessPass& messPass);/**<
  This function calculates the radial grid velocity, and does it by only considering the radial 
  and theta terms
  
  @param[in,out] grid contains the local grid, and will hold the newly updated radial grid 
                 velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
  @param[in,out] messPass handles data needed for message passing
  */
void calNewU0_RTP(Grid& grid,Parameters &parameters,Time &time,ProcTop &procTop
  ,MessPass& messPass);/**<
  This function calculates the radial grid velocity, and does it by considering all radial, theta 
  and phi terms
  
  @param[in,out] grid contains the local grid, and will hold the newly updated radial grid 
                 velocities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology
  @param[in,out] messPass handles data needed for message passing
  */
void calNewR(Grid& grid, Time &time);/**<
  This function calculates the radii, from the new radial grid velocities
  
  @param[in,out] grid contains the local grid, and will hold the newly updated radial velocities
  @param[in] time contains time information, e.g. time step, current time etc.
  */
void calNewD_R(Grid& grid, Parameters &parameters, Time &time,ProcTop &procTop);/**<
  This function calculates new densities using terms in the radial direction only
  
  @param[in,out] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology, uses \ref ProcTop::nRank 
             when reporting negative densities
  */
void calNewD_RT(Grid& grid, Parameters &parameters, Time &time,ProcTop &procTop);/**<
  This function calculates new densities using terms in the radial and theta directions
  
  @param[in,out] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology, uses \ref ProcTop::nRank 
             when reporting negative densities
  */
void calNewD_RTP(Grid& grid, Parameters &parameters, Time &time,ProcTop &procTop);/**<
  This function calculates new densities using terms in the radial, theta, and phi directions
  
  @param[in,out] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology, uses \ref ProcTop::nRank 
             when reporting negative densities
  */
void calNewE_R_AD(Grid& grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates new adiabatic energies using terms in the radial direction.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop
  */
void calNewE_R_NA(Grid& grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates new non-adiabatic energies using terms in the radial direction and
  includes radiation diffusion terms.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop
  */
void calNewE_R_NA_LES(Grid& grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates new non-adiabatic energies using terms in the radial direction and
  includes radiation diffusion terms. It also includes the terms for including real viscosity, used
  in the LES.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop
  */
void calNewE_RT_AD(Grid& grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates new adiabatic energies using terms in the radial and theta directions.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop
  */
void calNewE_RT_NA(Grid& grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates new non-adiabatic energies using terms in the radial and theta 
  directions and includes radiation diffusion terms.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop
  */
void calNewE_RT_NA_LES(Grid& grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates new non-adiabatic energies using terms in the radial and theta 
  directions and includes radiation diffusion terms. It also includes the terms for including real 
  viscosity, used in the LES.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop
  */
void calNewE_RTP_AD(Grid& grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates new adiabatic energies using terms in the radial, theta, and phi
  directions.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop
  */
void calNewE_RTP_NA(Grid& grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates new non-adiabatic energies using terms in the radial, theta, and phi 
  directions and includes radiation diffusion terms.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop
  */
void calNewE_RTP_NA_LES(Grid& grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates new non-adiabatic energies using terms in the radial, theta, and phi 
  directions and includes radiation diffusion terms. It also includes the terms for including real 
  viscosity, used in the LES.
  
  @param[in,out] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in] time contains time information, e.g. time step, current time etc.
  @param[in] procTop
  */
void calNewDenave_None(Grid &grid);/**<
  This function is a dumby funciton, and doesn't do anything. In the case of a 1D calculation
  the average density is undefined, and only the density is used. This is different from the case
  where the 1D region exsists on the rank 0 processor, but the grid as a whole is really 2D or 3D.
  In which case \ref calNewDenave_R should be used instead.
  
  @param[in,out] grid
  */
void calNewDenave_R(Grid& grid);/**<
  This function calculates the horizontal average density in a 3\\1D region. This really just copies 
  the density from the particular radial zone into the averaged density variable. This way it can be
  used exactly the same way in the 1D region as it is in the 3D region. This is done using the
  density in the new grid, and places the result into the new grid.
  
  @param[in,out] grid supplies the information needed to calculate the horizontal density average, 
                 it also stores the calculated horizontally averaged density.
  */
void calNewDenave_RT(Grid& grid);/**<
  This function calculates the horizontal average density in a 2D region from the new grid density
  and stores the result in the new grid.
  
  @param[in,out] grid supplies the information needed to calculate the horizontal density average, 
                 it also stores the calculated horizontally averaged density.
  */
void calNewDenave_RTP(Grid& grid);/**<
  This function calculates the horizontal average density in a 3D region from the new grid density
  and stores the result in the new grid.
  
  @param[in,out] grid supplies the information needed to calculate the horizontal density average, 
                 it also stores the calculated horizontally averaged density.
  */
void calNewP_GL(Grid& grid, Parameters &parameters);/**<
  This function calculates the pressure. It is calculated using the new values of quantities and 
  places the result in the new grid. It uses a gamma law gas give in \ref dEOS_GL to calculate the
  pressure.
  
  @param[in,out] grid supplies the input for calculating the pressure and also accepts the result
                 of the pressure calculations.
  @param[in] parameters contains parameters used in calculating the pressure, namely the adiabatic
             gamma that is used.
  */
void calNewTPKappaGamma_TEOS(Grid& grid, Parameters &parameters);/**<
  This function calculates the Temperature, pressure and opacity of a cell. It calculates it using
  the new vaules of quantities and places the result in the new grid.
  
  @param[in,out] grid supplies the input for calculating the pressure and also accepts the result
                 of the pressure calculation
  @param[in] parameters contains parameters used in calculating the pressure.
  */
void calNewPEKappaGamma_TEOS(Grid& grid, Parameters &parameters);/**<
  This function calculates the Energy, pressure and opacity of a cell. It calculates it using
  the new vaules of quantities and places the result in the new grid.
  
  @param[in,out] grid supplies the input for calculating the pressure and also accepts the result
                 of the pressure calculation
  @param[in] parameters contains parameters used in calculating the pressure.
  */
void calNewQ0_R_GL(Grid& grid, Parameters &parameters);/**<
  This funciton calculates the artificial viscosity of a cell. It calculates it using the new values
  of quantities and places the result in the new grid. It does this for the radial component of the
  viscosity only. It uses the sound speed derived from the adiabatic gamma given for the gamma law
  gas equation of state.
  
  @param[in,out] grid supplies the input for calculating the artificial viscosity and also accepts 
    the result of the artificial viscosity calculation.
  @param[in] parameters contains parameters used when calculating the artificial viscosity, namely
    the adiabatic gamma.
  */
void calNewQ0_R_TEOS(Grid& grid, Parameters &parameters);/**<
  This function calculates the artificial viscosity of a cell. It calculates it using the old values
  of quantities and places the result in the old grid. It does this for the radial component of the
  viscosity only. It uses a sound speed derived from a tabulated equaiton of state for the 
  calculation.
  
  @param[in,out] grid supplies the input for calculating the artificial viscosity and also accepts 
    the result of the artificial viscosity calculation
  @param[in] parameters contains parameters used in calculating the artificial viscosity.
  */
void calNewQ0Q1_RT_GL(Grid& grid, Parameters &parameters);/**<
  This function calculates the artificial viscosity of a cell. It calculates it using the new values
  of quantities and places the result in the new grid. It does this for the radial and theta
  componenets of the viscosity. It uses the sound speed derived from the adiabatic gamma given for 
  the gamma law gas equation of state.
  
  @param[in,out] grid supplies the input for calculating the artificial viscosity and also accepts
    the result of the artificial viscosity calculations.
  @param[in] parameters contains parameters used when calculating the artificial viscosity, namely
    the adiabatic gamma.
  */
void calNewQ0Q1_RT_TEOS(Grid& grid, Parameters &parameters);/**<
  This function calculates the artificial viscosity of a cell. It calculates it using the old values
  of quantities and places the result in the old grid. It does this for two component of the
  viscosity.
  
  @param[in,out] grid supplies the input for calculating the artificial viscosity and also accepts 
                 the result of the artificial viscosity calculation
  @param[in] parameters contains parameters used in calculating the artificial viscosity.
  */
void calNewQ0Q1Q2_RTP_GL(Grid& grid, Parameters &parameters);/**<
  This function calculates the artificial viscosity of a cell. It calculates it using the new values
  of quantities and places the result in the new grid. It does this for the radial, theta, and phi
  componenets of the viscosity. It uses the sound speed derived from the adiabatic gamma given for 
  the gamma law gas equation of state.
  
  @param[in,out] grid supplies the input for calculating the artificial viscosity and also accepts
    the result of the artificial viscosity calculations.
  @param[in] parameters contains parameters used when calculating the artificial viscosity, namely
    the adiabatic gamma.
  */
void calNewQ0Q1Q2_RTP_TEOS(Grid& grid, Parameters &parameters);/**<
  This function calculates the artificial viscosity of a cell. It calculates it using the old values
  of quantities and places the result in the old grid. It does this for the three component of the
  viscosity.
  
  @param[in,out] grid supplies the input for calculating the artificial viscosity and also accepts 
                 the result of the artificial viscosity calculation
  @param[in] parameters contains parameters used in calculating the artificial viscosity.
  */
void calNewEddyVisc_None(Grid &grid, Parameters &parameters);/**<
  This function is a empty function used as a place holder when no eddy viscosity model is being 
  used.
  
  @param[in,out] grid
  @param[in] parameters
  */
void calNewEddyVisc_R_CN(Grid &grid, Parameters &parameters);/**<
  This function calculates the eddy viscosity using a constant times the zoning with only the radial
  terms.
  
  @param[in,out] grid supplies the input for calculating the eddy viscosity.
  @param[in] parameters contains parameters used in calculating the eddy viscosity.
  */
void calNewEddyVisc_RT_CN(Grid &grid, Parameters &parameters);/**<
  This function calculates the eddy viscosity using a constant times the zoning with only the radial
  and theta terms.
  
  @param[in,out] grid supplies the input for calculating the eddy viscosity.
  @param[in] parameters contains parameters used in calculating the eddy viscosity.
  */
void calNewEddyVisc_RTP_CN(Grid &grid, Parameters &parameters);/**<
  This function calculates the eddy viscosity using a constant times the zoning with only the 
  radial, theta, and phi terms.
  
  @param[in,out] grid supplies the input for calculating the eddy viscosity.
  @param[in] parameters contains parameters used in calculating the eddy viscosity.
  */
void calNewEddyVisc_R_SM(Grid &grid, Parameters &parameters);/**<
  This function calculates the eddy viscosity with only the radial terms.
  
  @param[in,out] grid supplies the input for calculating the eddy viscosity.
  @param[in] parameters contains parameters used in calculating the eddy viscosity.
  */
void calNewEddyVisc_RT_SM(Grid &grid, Parameters &parameters);/**<
  This function calculates the eddy viscosity with only the radial and theta terms.
  
  @param[in,out] grid supplies the input for calculating the eddy viscosity.
  @param[in] parameters contains parameters used in calculating the eddy viscosity.
  */
void calNewEddyVisc_RTP_SM(Grid &grid, Parameters &parameters);/**<
  This function calculates the eddy viscosity with only the radial, theta, and phi terms.
  
  @param[in,out] grid supplies the input for calculating the eddy viscosity.
  @param[in] parameters contains parameters used in calculating the eddy viscosity.
  */
void calOldDenave_None(Grid &grid);/**<
  This function is a dumby funciton, and doesn't do anything. In the case of a 1D calculation
  the average density is undefined, and only the density is used. This is different from the case
  where the 1D region exsists on the rank 0 processor, but the grid as a whole is really 2D or 3D.
  In which case \ref calOldDenave_R should be used instead.
  */
void calOldDenave_R(Grid& grid);/**<
  This function does nothing as the averaged density is not needed in 1D calculations.
  
  @param[in,out] grid supplies the information needed to calculate the horizontal density average, 
                 it also stores the calculated horizontally averaged density.
  */
void calOldDenave_RT(Grid& grid);/**<
  This function calculates the horizontal average density in a 2D region. This function differs from
  \ref calNewDenave_RT in that it calculates the average density from the old grid density
  and stores the result in the old grid. While calNewDenave_RT calculates the average 
  density from the new grid density and places the result in the new grid.
  
  @param[in,out] grid supplies the information needed to calculate the horizontal density average, 
                 it also stores the calculated horizontally averaged density.
  */
void calOldDenave_RTP(Grid& grid);/**<
  This function calculates the horizontal average density in a 3D region. This function differs from
  \ref calNewDenave_RTP in that it calculates the average density from the old grid density
  and stores the result in the old grid. While calNewDenave_RTP calculates the average 
  density from the new grid density and places the result in the new grid.
  
  @param[in,out] grid supplies the information needed to calculate the horizontal density average, 
                 it also stores the calculated horizontally averaged density.
  */
void calOldP_GL(Grid& grid, Parameters &parameters);/**<
  This function calculates the pressure using a gamma law gas, calculate by \ref dEOS_GL.
  
  @param[in,out] grid supplies the input for calculating the pressure and also accepts the results
                 of the pressure calculations
  @param[in] parameters contains parameters used in calculating the pressure, namely the value of
             the adiabatic gamma
  */
void calOldPEKappaGamma_TEOS(Grid& grid, Parameters &parameters);/**<
  This function calculates the pressure, energy, opacity, and adiabatic index of a cell. It
  calculates it using the old vaules of quantities and places the result in the old grid. This
  function is used to initialize the internal variables pressure, energy and kappa, and is suitable
  for both 1D and 3D calculations.
  
  @param[in,out] grid supplies the input for calculating the pressure and also accepts the result
                 of the pressure calculation
  @param[in] parameters contains parameters used in calculating the pressure.
  */
void calOldQ0_R_GL(Grid& grid, Parameters &parameters);/**<
  This function calculates the artificial viscosity of a cell. It calculates it using the old vaules
  of quantities and places the result in the old grid. It does this for the radial component of the
  viscosity only. This function is used when using a gamma law gas equation of state.
  
  @param[in,out] grid supplies the input for calculating the artificial viscosity and also accepts 
                 the result of the artificial viscosity calculation
  @param[in] parameters contains parameters used in calculating the artificial viscosity.
  */
void calOldQ0_R_TEOS(Grid& grid, Parameters &parameters);/**<
  This function calculates the artificial viscosity of a cell. It calculates it using the old vaules
  of quantities and places the result in the old grid. It does this for 1D viscosity only.
  
  @param[in,out] grid supplies the input for calculating the pressure and also accepts the result
                 of the pressure calculation
  @param[in] parameters contains parameters used in calculating the artificial viscosity.
  */
void calOldQ0Q1_RT_GL(Grid& grid, Parameters &parameters);/**<
  This function calculates the artificial viscosity of a cell. It calculates it using the old vaules
  of quantities and places the result in the old grid. It does this for the two components of the
  viscosity. This function is used when using a gamma law gas equation of state.
  
  @param[in,out] grid supplies the input for calculating the artificial viscosity and also accepts 
                 the result of the artificial viscosity calculation
  @param[in] parameters contains parameters used in calculating the artificial viscosity.
  */
void calOldQ0Q1_RT_TEOS(Grid& grid, Parameters &parameters);/**<
  This function calculates the artificial viscosity of a cell. It calculates it using the old vaules
  of quantities and places the result in the old grid. It does this for two components of the
  viscosity. This function is used when using a tabulated equation of state.
  
  @param[in,out] grid supplies the input for calculating the artificial viscosity and also accepts 
                 the result of the artificial viscosity calculation
  @param[in] parameters contains parameters used in calculating the artificial viscosity.
  */
void calOldQ0Q1Q2_RTP_GL(Grid& grid, Parameters &parameters);/**<
  This function calculates the artificial viscosity of a cell. It calculates it using the old vaules
  of quantities and places the result in the old grid. It does this for the three components of the
  viscosity. This function is used when using a gamma law gas equation of state.
  
  @param[in,out] grid supplies the input for calculating the artificial viscosity and also accepts 
                 the result of the artificial viscosity calculation
  @param[in] parameters contains parameters used in calculating the artificial viscosity.
  */
void calOldQ0Q1Q2_RTP_TEOS(Grid& grid, Parameters &parameters);/**<
  This function calculates the artificial viscosity of a cell. It calculates it using the old vaules
  of quantities and places the result in the old grid. It does this for the three components of the
  viscosity. This function is used when using a tabulated equation of state.
  
  @param[in,out] grid supplies the input for calculating the artificial viscosity and also accepts 
                 the result of the artificial viscosity calculation
  @param[in] parameters contains parameters used in calculating the artificial viscosity.
  */
void calOldEddyVisc_R_CN(Grid &grid, Parameters &parameters);/**<
  Calculates the eddy viscosity using a constant times the zoning including only the radial terms.
  It puts the result into the old grid. This funciton is used to initalize the eddy viscosity when 
  the code begins execution.
  */
void calOldEddyVisc_RT_CN(Grid &grid, Parameters &parameters);/**<
  Calculates the eddy viscosity using a constant times the zoning including only the radial and 
  theta terms. It puts the result into the old grid. This funciton is used to initalize the eddy
  viscosity when the code begins execution.
  
  @param[in,out] grid supplies the input for calculating the eddy viscosity.
  @param[in] parameters contains parameters used in calculating the eddy viscosity.
  */
void calOldEddyVisc_RTP_CN(Grid &grid, Parameters &parameters);/**<
  Calculates the eddy viscosity using a constant times the zoning including the radial, theta, 
  and phi terms. It puts the result into the old grid. This funciton is used to initalize the eddy
  viscosity when the code begins execution.
  
  @param[in,out] grid supplies the input for calculating the eddy viscosity.
  @param[in] parameters contains parameters used in calculating the eddy viscosity.
  */
void calOldEddyVisc_R_SM(Grid &grid, Parameters &parameters);/**<
  Calculates the eddy viscosity including only the radial terms. It puts the result into the old 
  grid. This funciton is used to initalize the eddy viscosity when the code begins execution. It
  uses the Smagorinsky model for calculating the eddy viscosity.
  
  @param[in,out] grid supplies the input for calculating the eddy viscosity.
  @param[in] parameters contains parameters used in calculating the eddy viscosity.
  */
void calOldEddyVisc_RT_SM(Grid &grid, Parameters &parameters);/**<
  Calculates the eddy viscosity including only the radial and theta terms. It puts the result into
  the old grid. This funciton is used to initalize the eddy viscosity when the code begins
  execution.It uses the Smagorinsky model for calculating the eddy viscosity.
  
  @param[in,out] grid supplies the input for calculating the eddy viscosity.
  @param[in] parameters contains parameters used in calculating the eddy viscosity.
  */
void calOldEddyVisc_RTP_SM(Grid &grid, Parameters &parameters);/**<
  Calculates the eddy viscosity including the radial, theta, and phi terms. It puts the result into
  the old grid. This funciton is used to initalize the eddy viscosity when the code begins
  execution.It uses the Smagorinsky model for calculating the eddy viscosity.
  
  @param[in,out] grid supplies the input for calculating the eddy viscosity.
  @param[in] parameters contains parameters used in calculating the eddy viscosity.
  */
void calDelt_R_GL(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates the time step by considering the sound crossing time in the radial 
  direction only and is compatiable with a gamma law gass EOS.
  
  @param[in] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in,out] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology. This function uses 
             \ref ProcTop::nRank to pass messages.
  */
void calDelt_R_TEOS(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates the time step by considering the sound crossing time in the radial 
  direction only and is compatiable with a tabulated EOS.
  
  @param[in] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in,out] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology. This function uses 
             \ref ProcTop::nRank to pass messages.
  */
void calDelt_RT_GL(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates the time step by considering the sound crossing time in the radial 
  and theta directions only and is compatiable with a gamma law gass EOS.
  
  @param[in] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in,out] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology. This function uses 
             \ref ProcTop::nRank to pass messages.
  */
void calDelt_RT_TEOS(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates the time step by considering the sound crossing time in the radial 
  and theta directions and is compatiable with a tabulated EOS.
  
  @param[in] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in,out] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology. This function uses 
             \ref ProcTop::nRank to pass messages.
  */
void calDelt_RTP_GL(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates the time step by considering the sound crossing time in the radial, 
  theta and phi directions only and is compatiable with a gamma law gass EOS.
  
  @param[in] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in,out] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology. This function uses 
             \ref ProcTop::nRank to pass messages.
  */
void calDelt_RTP_TEOS(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function calculates the time step by considering the sound crossing time in the radial, 
  theta and phi directions and is compatiable with a tabulated EOS.
  
  @param[in] grid contains the local grid, and will hold the newly updated densities
  @param[in] parameters various parameters needed for the calculation
  @param[in,out] time contains time information, e.g. time step, current time etc.
  @param[in] procTop contains information about the processor topology. This function uses 
             \ref ProcTop::nRank to pass messages.
  */
void calDelt_CONST(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop);/**<
  This function is used when a constant tie step is desired.
  */
void implicitSolve_None(Grid &grid,Implicit &implicit,Parameters &parameters,Time &time
  , ProcTop &procTop,MessPass &messPass,Functions &functions);/**<
  This is an empty function, to be called when no implicit solution is needed. This allows the same
  code in the main program to be executed wheather or not an implicit solution is being preformed
  by setting the funciton pointer to this funciton if there is no implicit solution required.
  */
void implicitSolve_R(Grid &grid,Implicit &implicit,Parameters &parameters,Time &time
  , ProcTop &procTop,MessPass &messPass,Functions &functions);/**<
  This function solves for temperature corrections based on derivatives of the radial non-adiabatic 
  energy equation with respect to the new temperature. It then uses these derivatives as entries in
  the coeffecient matrix. The discrepancy in the balance of the energy equation with the new 
  temperature, energy, pressure, and opacity are included as the right hand side of the system of 
  equaitons. Solving this system of equaitons provides the corrections needed for the new 
  temperature. This processes is then repeated until the corrections are small. At this point the 
  new temperature is used to update the energy, pressure, and opacity in the new grid via the 
  equaiton of state.
  */
void implicitSolve_RT(Grid &grid,Implicit &implicit,Parameters &parameters,Time &time
  , ProcTop &procTop,MessPass &messPass,Functions &functions);/**<
  This function solves for temperature corrections based on derivatives of the radial-theta
  non-adiabatic energy equation with respect to the new temperature. It then uses these derivatives
  as entries in the coeffecient matrix. The discrepancy in the balance of the energy equation with 
  the new temperature, energy, pressure, and opacity are included as the right hand side of the 
  system of equaitons. Solving this system of equaitons provides the corrections needed for the new 
  temperature. This processes is then repeated until the corrections are small. At this point the 
  new temperature is used to update the energy, pressure, and opacity in the new grid via the 
  equaiton of state.
  */
void implicitSolve_RTP(Grid &grid,Implicit &implicit,Parameters &parameters,Time &time
  , ProcTop &procTop,MessPass &messPass,Functions &functions);/**<
  This function solves for temperature corrections based on derivatives of the radial-theta-phi
  non-adiabatic energy equation with respect to the new temperature. It then uses these derivatives
  as entries in the coeffecient matrix. The discrepancy in the balance of the energy equation with 
  the new temperature, energy, pressure, and opacity are included as the right hand side of the 
  system of equaitons. Solving this system of equaitons provides the corrections needed for the new 
  temperature. This processes is then repeated until the corrections are small. At this point the 
  new temperature is used to update the energy, pressure, and opacity in the new grid via the 
  equaiton of state.
  */
double dImplicitEnergyFunction_None(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k);/**<
  This is an empty function, that isn't even called when no implicit solution is needed. This safe
  guards against future addition which may need to call an empty function when no implicit solve is
  being done.
  */
double dImplicitEnergyFunction_R(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k);/**<
  This function is used to determine the agreement of the updated values at \f$n+1\f$, with
  each other in the non-adiabatic energy equation. The \c _R version of the funciton contains only 
  the radial terms, and should be used for purely radial calculations. This function can also be 
  used for calculating numerical deriviatives by varying the input temperatures.
  
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] dTemps dTemps[0]=dT_ijk_np1 is the temperature at radial position \f$(i,j,k)\f$ and 
             time \f$ n+1 \f$,dTemps[1]=dT_ip1jk_np1 is the temperature at radial position 
             \f$(i+1,j,k)\f$ and time \f$ n+1 \f$, dTemps[2]=dT_im1jk_np1 is the temperature at 
             radial position \f$(i-1,j,k)\f$ and time \f$ n+1 \f$.
  @param[in] i is the radial index to evaluate the function at.
  @param[in] j is the theta index to evaluate the function at.
  @param[in] k is the phi index to evaluate the function at.
  */
double dImplicitEnergyFunction_R_SB(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k);/**<
  This function is used to determine the agreement of the updated values at \f$n+1\f$, with
  each other in the non-adiabatic energy equation. The \c _R version of the funciton contains only 
  the radial terms, and should be used for purely radial calculations. This function can also be 
  used for calculating numerical deriviatives by varying the input temperatures. This funciton 
  differs from the version without the "_SB" suffix (\ref dImplicitEnergyFunction_R)in that it is
  tailored to the surface boundary region.
  
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] dTemps dTemps[0]=dT_ijk_np1 is the temperature at radial position \f$(i,j,k)\f$ and 
             time \f$ n+1 \f$ and time \f$ n+1 \f$, dTemps[1]=dT_im1jk_np1 is the temperature at
             radial position \f$(i-1,j,k)\f$ and time \f$ n+1 \f$.
  @param[in] i is the radial index to evaluate the function at.
  @param[in] j is the theta index to evaluate the function at.
  @param[in] k is the phi index to evaluate the function at.
  */
double dImplicitEnergyFunction_RT(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k);/**<
  This function is used to determine the agreement of the updated values at \f$n+1\f$, with
  each other in the non-adiabatic energy equation. The \c _RT version of the funciton contains only 
  the radial and theta terms, and should be used for radial-theta calculations. This function can
  also be used for calculating numerical deriviatives by varying the input temperatures.
  
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] dTemps dTemps[0]=dT_ijk_np1 is the temperature at radial position \f$(i,j,k)\f$ and 
             time \f$ n+1 \f$, dTemps[1]=dT_ip1jk_np1 is the temperature at radial position 
             \f$(i+1,j,k)\f$ and time \f$ n+1 \f$, dTemps[2]=dT_im1jk_np1 is the temperature at 
             radial position \f$(i-1,j,k)\f$ and time \f$ n+1 \f$, dTemps[3]=dT_ijp1k_np1 is the 
             temperature at radial position \f$(i,j+1,k)\f$ and time \f$ n+1 \f$,
             dTemps[4]=dT_ijm1k_np1 is the temperature at radial position \f$(i,j-1,k)\f$ and time
             \f$ n+1 \f$.
  @param[in] i is the radial index to evaluate the function at.
  @param[in] j is the theta index to evaluate the function at.
  @param[in] k is the phi index to evaluate the function at.
  */
double dImplicitEnergyFunction_RT_SB(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k);/**<
  This function is used to determine the agreement of the updated values at \f$n+1\f$, with
  each other in the non-adiabatic energy equation. The \c _RT version of the funciton contains only 
  the radial and theta terms, and should be used for radial-theta calculations. This function can
  also be used for calculating numerical deriviatives by varying the input temperatures. This 
  funciton differs from the version without the "_SB" suffix (\ref dImplicitEnergyFunction_RT)in 
  that it is tailored to the surface boundary region.
  
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] dTemps dTemps[0]=dT_ijk_np1 is the temperature at radial position \f$(i,j,k)\f$ and 
              time \f$ n+1 \f$, dTemps[1]=dT_ip1jk_np1 is the temperature at radial position 
              \f$(i+1,j,k)\f$ and time \f$ n+1 \f$, dTemps[2]=dT_im1jk_np1 is the temperature at
              radial position \f$(i-1,j,k)\f$ and time \f$ n+1 \f$, dTemps[3]=dT_ijp1k_np1 is the
              temperature at radial position \f$(i,j+1,k)\f$ and time \f$ n+1 \f$, 
              dTemps[4]=dT_ijm1k_np1 is the temperature at radial position \f$(i,j-1,k)\f$ and time 
              \f$ n+1 \f$.
  @param[in] i is the radial index to evaluate the function at.
  @param[in] j is the theta index to evaluate the function at.
  @param[in] k is the phi index to evaluate the function at.
  */
double dImplicitEnergyFunction_RTP(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k);/**<
  This function is used to determine the agreement of the updated values at \f$n+1\f$, with
  each other in the non-adiabatic energy equation. The \c _RTP version of the funciton contains 
  terms for all three directions, and should be used for calculations involving all three 
  directions. This function can also be used for calculating numerical deriviatives by varying the 
  input temperatures. This funciton differs from the version without the "_SB" suffix 
  (\ref dImplicitEnergyFunction_RT)in that it is tailored to the surface boundary region.
  
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] dTemps , dTemps[0]=dT_ijk_np1 is the temperature at radial position \f$(i,j,k)\f$ and 
             time \f$ n+1 \f$, dTemps[1]=dT_ip1jk_np1 is the temperature at radial position 
             \f$(i+1,j,k)\f$ and time \f$ n+1 \f$, dTemps[2]=dT_im1jk_np1 is the temperature at 
             radial position \f$(i-1,j,k)\f$ and time \f$ n+1 \f$, dTemps[3]=dT_ijp1k_np1 is the 
             temperature at radial position \f$(i,j+1,k)\f$ and time \f$ n+1 \f$, 
             dTemps[4]=dT_ijm1k_np1 is the temperature at radial position \f$(i,j-1,k)\f$ and time 
             \f$ n+1 \f$, dTemps[5]=dT_ijkp1_np1 is the temperature at radial position 
             \f$(i,j,k+1)\f$ and time \f$ n+1 \f$, dTemps[6]=dT_ijkm1_np1 is the temperature at 
             radial position \f$(i,j,k-1)\f$ and time \f$ n+1 \f$.
  @param[in] i is the radial index to evaluate the function at.
  @param[in] j is the theta index to evaluate the function at.
  @param[in] k is the phi index to evaluate the function at.
  */
double dImplicitEnergyFunction_RTP_SB(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k);/**<
  This function is used to determine the agreement of the updated values at \f$n+1\f$, with
  each other in the non-adiabatic energy equation. The \c _RTP version of the funciton contains 
  terms for all three directions, and should be used for calculations involving all three 
  directions. This function can also be used for calculating numerical deriviatives by varying the 
  input temperatures. This funciton differs from the version without the "_SB" suffix 
  (\ref dImplicitEnergyFunction_RT)in that it is tailored to the surface boundary region.
  
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] dTemps dTemps[0]=dT_ijk_np1 is the temperature at radial position \f$(i,j,k)\f$ and 
             time \f$ n+1 \f$, dTemps[1]=dT_ip1jk_np1 is the temperature at radial position
             \f$(i+1,j,k)\f$ and time \f$ n+1 \f$, dTemps[2]=dT_im1jk_np1 is the temperature at
             radial position \f$(i-1,j,k)\f$ and time \f$ n+1 \f$, dTemps[3]=dT_ijp1k_np1 is the 
             temperature at radial position \f$(i,j+1,k)\f$ and time \f$ n+1 \f$, 
             dTemps[4]=dT_ijm1k_np1 is the temperature at radial position \f$(i,j-1,k)\f$ and time 
             \f$ n+1 \f$, dTemps[5]=dT_ijkp1_np1 is the temperature at radial position
             \f$(i,j,k+1)\f$ and time \f$ n+1 \f$, dTemps[6]=dT_ijkm1_np1 is the temperature at
             radial position \f$(i,j,k-1)\f$ and time \f$ n+1 \f$.
  @param[in] i is the radial index to evaluate the function at.
  @param[in] j is the theta index to evaluate the function at.
  @param[in] k is the phi index to evaluate the function at.
  */
double dImplicitEnergyFunction_R_LES(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k);/**<
  This function is used to determine the agreement of the updated values at \f$n+1\f$, with
  each other in the non-adiabatic energy equation. The \c _R version of the funciton contains only 
  the radial terms, and should be used for purely radial calculations. This function can also be 
  used for calculating numerical deriviatives by varying the input temperatures.
  
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] dTemps dTemps[0]=dT_ijk_np1 is the temperature at radial position \f$(i,j,k)\f$ and 
             time \f$ n+1 \f$,dTemps[1]=dT_ip1jk_np1 is the temperature at radial position 
             \f$(i+1,j,k)\f$ and time \f$ n+1 \f$, dTemps[2]=dT_im1jk_np1 is the temperature at 
             radial position \f$(i-1,j,k)\f$ and time \f$ n+1 \f$.
  @param[in] i is the radial index to evaluate the function at.
  @param[in] j is the theta index to evaluate the function at.
  @param[in] k is the phi index to evaluate the function at.
  */
double dImplicitEnergyFunction_R_LES_SB(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k);/**<
  This function is used to determine the agreement of the updated values at \f$n+1\f$, with
  each other in the non-adiabatic energy equation. The \c _R version of the funciton contains only 
  the radial terms, and should be used for purely radial calculations. This function can also be 
  used for calculating numerical deriviatives by varying the input temperatures. This funciton 
  differs from the version without the "_SB" suffix (\ref dImplicitEnergyFunction_R)in that it is
  tailored to the surface boundary region.
  
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] dTemps dTemps[0]=dT_ijk_np1 is the temperature at radial position \f$(i,j,k)\f$ and 
             time \f$ n+1 \f$ and time \f$ n+1 \f$, dTemps[1]=dT_im1jk_np1 is the temperature at
             radial position \f$(i-1,j,k)\f$ and time \f$ n+1 \f$.
  @param[in] i is the radial index to evaluate the function at.
  @param[in] j is the theta index to evaluate the function at.
  @param[in] k is the phi index to evaluate the function at.
  */
double dImplicitEnergyFunction_RT_LES(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k);/**<
  This function is used to determine the agreement of the updated values at \f$n+1\f$, with
  each other in the non-adiabatic energy equation. The \c _RT version of the funciton contains only 
  the radial and theta terms, and should be used for radial-theta calculations. This function can
  also be used for calculating numerical deriviatives by varying the input temperatures.
  
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] dTemps dTemps[0]=dT_ijk_np1 is the temperature at radial position \f$(i,j,k)\f$ and 
             time \f$ n+1 \f$, dTemps[1]=dT_ip1jk_np1 is the temperature at radial position 
             \f$(i+1,j,k)\f$ and time \f$ n+1 \f$, dTemps[2]=dT_im1jk_np1 is the temperature at 
             radial position \f$(i-1,j,k)\f$ and time \f$ n+1 \f$, dTemps[3]=dT_ijp1k_np1 is the 
             temperature at radial position \f$(i,j+1,k)\f$ and time \f$ n+1 \f$,
             dTemps[4]=dT_ijm1k_np1 is the temperature at radial position \f$(i,j-1,k)\f$ and time
             \f$ n+1 \f$.
  @param[in] i is the radial index to evaluate the function at.
  @param[in] j is the theta index to evaluate the function at.
  @param[in] k is the phi index to evaluate the function at.
  */
double dImplicitEnergyFunction_RT_LES_SB(Grid &grid,Parameters &parameters,Time &time
  ,double dTemps[],int i,int j,int k);/**<
  This function is used to determine the agreement of the updated values at \f$n+1\f$, with
  each other in the non-adiabatic energy equation. The \c _RT version of the funciton contains only 
  the radial and theta terms, and should be used for radial-theta calculations. This function can
  also be used for calculating numerical deriviatives by varying the input temperatures. This 
  funciton differs from the version without the "_SB" suffix (\ref dImplicitEnergyFunction_RT)in 
  that it is tailored to the surface boundary region.
  
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] dTemps dTemps[0]=dT_ijk_np1 is the temperature at radial position \f$(i,j,k)\f$ and 
              time \f$ n+1 \f$, dTemps[1]=dT_ip1jk_np1 is the temperature at radial position 
              \f$(i+1,j,k)\f$ and time \f$ n+1 \f$, dTemps[2]=dT_im1jk_np1 is the temperature at
              radial position \f$(i-1,j,k)\f$ and time \f$ n+1 \f$, dTemps[3]=dT_ijp1k_np1 is the
              temperature at radial position \f$(i,j+1,k)\f$ and time \f$ n+1 \f$, 
              dTemps[4]=dT_ijm1k_np1 is the temperature at radial position \f$(i,j-1,k)\f$ and time 
              \f$ n+1 \f$.
  @param[in] i is the radial index to evaluate the function at.
  @param[in] j is the theta index to evaluate the function at.
  @param[in] k is the phi index to evaluate the function at.
  */
double dImplicitEnergyFunction_RTP_LES(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k);/**<
  This function is used to determine the agreement of the updated values at \f$n+1\f$, with
  each other in the non-adiabatic energy equation. The \c _RTP version of the funciton contains 
  terms for all three directions, and should be used for calculations involving all three 
  directions. This function can also be used for calculating numerical deriviatives by varying the 
  input temperatures. This funciton differs from the version without the "_SB" suffix 
  (\ref dImplicitEnergyFunction_RT)in that it is tailored to the surface boundary region.
  
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] dTemps , dTemps[0]=dT_ijk_np1 is the temperature at radial position \f$(i,j,k)\f$ and 
             time \f$ n+1 \f$, dTemps[1]=dT_ip1jk_np1 is the temperature at radial position 
             \f$(i+1,j,k)\f$ and time \f$ n+1 \f$, dTemps[2]=dT_im1jk_np1 is the temperature at 
             radial position \f$(i-1,j,k)\f$ and time \f$ n+1 \f$, dTemps[3]=dT_ijp1k_np1 is the 
             temperature at radial position \f$(i,j+1,k)\f$ and time \f$ n+1 \f$, 
             dTemps[4]=dT_ijm1k_np1 is the temperature at radial position \f$(i,j-1,k)\f$ and time 
             \f$ n+1 \f$, dTemps[5]=dT_ijkp1_np1 is the temperature at radial position 
             \f$(i,j,k+1)\f$ and time \f$ n+1 \f$, dTemps[6]=dT_ijkm1_np1 is the temperature at 
             radial position \f$(i,j,k-1)\f$ and time \f$ n+1 \f$.
  @param[in] i is the radial index to evaluate the function at.
  @param[in] j is the theta index to evaluate the function at.
  @param[in] k is the phi index to evaluate the function at.
  */
double dImplicitEnergyFunction_RTP_LES_SB(Grid &grid,Parameters &parameters,Time &time
  ,double dTemps[],int i,int j,int k);/**<
  This function is used to determine the agreement of the updated values at \f$n+1\f$, with
  each other in the non-adiabatic energy equation. The \c _RTP version of the funciton contains 
  terms for all three directions, and should be used for calculations involving all three 
  directions. This function can also be used for calculating numerical deriviatives by varying the 
  input temperatures. This funciton differs from the version without the "_SB" suffix 
  (\ref dImplicitEnergyFunction_RT)in that it is tailored to the surface boundary region.
  
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] dTemps dTemps[0]=dT_ijk_np1 is the temperature at radial position \f$(i,j,k)\f$ and 
             time \f$ n+1 \f$, dTemps[1]=dT_ip1jk_np1 is the temperature at radial position
             \f$(i+1,j,k)\f$ and time \f$ n+1 \f$, dTemps[2]=dT_im1jk_np1 is the temperature at
             radial position \f$(i-1,j,k)\f$ and time \f$ n+1 \f$, dTemps[3]=dT_ijp1k_np1 is the 
             temperature at radial position \f$(i,j+1,k)\f$ and time \f$ n+1 \f$, 
             dTemps[4]=dT_ijm1k_np1 is the temperature at radial position \f$(i,j-1,k)\f$ and time 
             \f$ n+1 \f$, dTemps[5]=dT_ijkp1_np1 is the temperature at radial position
             \f$(i,j,k+1)\f$ and time \f$ n+1 \f$, dTemps[6]=dT_ijkm1_np1 is the temperature at
             radial position \f$(i,j,k-1)\f$ and time \f$ n+1 \f$.
  @param[in] i is the radial index to evaluate the function at.
  @param[in] j is the theta index to evaluate the function at.
  @param[in] k is the phi index to evaluate the function at.
  */
double dEOS_GL(double dRho, double dE, Parameters parameters);/**<
  Calculates the pressure from the energy and density using a \f$\gamma\f$-law gas.
  
  @param[in] dRho the density of a cell
  @param[in] dE the energy of a cell
  @param[in] parameters contians various parameters, including \f$\gamma\f$ needed to calculate the
             pressure.
  @return the pressure
  
  This version of \ref dEOS_GL uses the same value of \f$\gamma\f$ through out the model. The
  equation of state is given by \f$\rho(\gamma-1)E\f$.
  */
void initDonorFracAndMaxConVel_R_GL(Grid &grid, Parameters &parameters);/**<
  Initializes the donor fraction, and the maximum convective velocity when starting a calculation.
  The donor fraction is used to determine the amount of upwinded donor cell to use in advection 
  terms. The maximum convective velocity is used for calculation of constant eddy viscosity
  parameter. This version of the fuction is for 1D, gamma law calculations.
  */
void initDonorFracAndMaxConVel_R_TEOS(Grid &grid, Parameters &parameters);/**<
  Initializes the donor fraction, and the maximum convective velocity when starting a calculation.
  The donor fraction is used to determine the amount of upwinded donor cell to use in advection 
  terms. The maximum convective velocity is used for calculation of constant eddy viscosity
  parameter. This version of the fuction is for 1D, tabulated equation of state calculations.
  */
void initDonorFracAndMaxConVel_RT_GL(Grid &grid, Parameters &parameters);/**<
  Initializes the donor fraction, and the maximum convective velocity when starting a calculation.
  The donor fraction is used to determine the amount of upwinded donor cell to use in advection 
  terms. The maximum convective velocity is used for calculation of constant eddy viscosity
  parameter. This version of the fuction is for 2D, gamma law calculations.
  */
void initDonorFracAndMaxConVel_RT_TEOS(Grid &grid, Parameters &parameters);/**<
  Initializes the donor fraction, and the maximum convective velocity when starting a calculation.
  The donor fraction is used to determine the amount of upwinded donor cell to use in advection 
  terms. The maximum convective velocity is used for calculation of constant eddy viscosity
  parameter. This version of the fuction is for 2D, tabulated equation of state calculations.
  */
void initDonorFracAndMaxConVel_RTP_GL(Grid &grid, Parameters &parameters);/**<
  Initializes the donor fraction, and the maximum convective velocity when starting a calculation.
  The donor fraction is used to determine the amount of upwinded donor cell to use in advection 
  terms. The maximum convective velocity is used for calculation of constant eddy viscosity
  parameter. This version of the fuction is for 3D, gamma law calculations.
  */
void initDonorFracAndMaxConVel_RTP_TEOS(Grid &grid, Parameters &parameters);/**<
  Initializes the donor fraction, and the maximum convective velocity when starting a calculation.
  The donor fraction is used to determine the amount of upwinded donor cell to use in advection 
  terms. The maximum convective velocity is used for calculation of constant eddy viscosity
  parameter. This version of the fuction is for 3D, tabulated equation of state calculations.
  */
inline double dET4(Parameters &parameters,double dEddyVisc_ijk_np1half, double dRho_ijk_np1half
  ,double dLengthScale4_ijk_np1half);/**<
  This is an additional turbulance term to be added to the energy equation.
  */
#endif
