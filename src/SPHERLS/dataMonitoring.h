#ifndef DATAMONITORING_H
#define DATAMONITORING_H

/** @file
  
  Header file for \ref dataMonitoring.cpp
*/


#include <string>
#include "xmlParser.h"
#include "global.h"

void initWatchZones(XMLNode xParent,ProcTop &procTop, Grid &grid, Output &output
  , Parameters &parameters, Time &time);/**<
  Reads in watchzones set in configuration file "SPHERLS.xml". A list is created on each processor 
  containing the watchzones on that processor's local grid. It also opens file streams for each 
  watchzone and writes out a header.
  
  @param[in] xParent
  @param[in] procTop
  @param[in] grid
  @param[in,out] output
  @param[in] parameters
  @param[in] time
  */
void writeWatchZones_R_GL(Output &output, Grid &grid, Parameters &parameters, Time &time
  , ProcTop &procTop);/**<
  Writes out the information for each watchzone specified in "SPHERLS.xml" in the case of a 1D
  gamma-law gas.
  
  @param[in,out] output
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] procTop
  */
void writeWatchZones_R_TEOS(Output &output, Grid &grid, Parameters &parameters, Time &time
  , ProcTop &procTop);/**<
  Writes out the information for each watchzone specified in "SPHERLS.xml" in the case of a 1D
  tabulated equation of state.
  
  @param[in,out] output
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] procTop
  */
void writeWatchZones_RT_GL(Output &output, Grid &grid, Parameters &parameters, Time &time
  , ProcTop &procTop);/**<
  Writes out the information for each watchzone specified in "SPHERLS.xml" in the case of a 2D
  gamma-law gas.
  
  @param[in,out] output
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] procTop
  */
void writeWatchZones_RT_TEOS(Output &output, Grid &grid, Parameters &parameters, Time &time
  , ProcTop &procTop);/**<
  Writes out the information for each watchzone specified in "SPHERLS.xml" in the case of a 2D
  tabulated equation of state.
  
  @param[in,out] output
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] procTop
  */
void writeWatchZones_RTP_GL(Output &output, Grid &grid, Parameters &parameters, Time &time
  , ProcTop &procTop);/**<
  Writes out the information for each watchzone specified in "SPHERLS.xml" in the case of a 3D
  gamma-law gas.
  
  @param[in,out] output
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] procTop
  */
void writeWatchZones_RTP_TEOS(Output &output, Grid &grid, Parameters &parameters, Time &time
  , ProcTop &procTop);/**<
  Writes out the information for each watchzone specified in "SPHERLS.xml" in the case of a 3D
  tabulated equation of state.
  
  @param[in,out] output
  @param[in] grid
  @param[in] parameters
  @param[in] time
  @param[in] procTop
  */
void finWatchZones(Output &output);/**<
  Closes the files opened for writting out the watchzones
  
  @param[in] output
  */
#endif
