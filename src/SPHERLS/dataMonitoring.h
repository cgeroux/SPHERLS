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
void initPeakKE(ProcTop &procTop, Output &output, PeakKETracking &peakKETracking, Time &time);/**<
  Opens the output file for tracking the peak kinetic energy, and writes out a header.
  
  @param[in] procTop
  @param[in] output
  @param[in,out] peakKETracking
  @param[in] time
  */
void writePeakKE_R(Grid &grid, Parameters &parameters, PeakKETracking &peakKETracking
  , ProcTop &procTop,Time &time);/**<
  This function calculates, and tracks the peak kinetic energy of the model, and averages it over
  three full periods. It does this by watch for 6 peaks in the kinetic energy and then averaging 
  them. It uses 6 because there are 2 peaks every period, when during contraction, and one during
  expansion. And writes the peak kinetic energy, the averaged peak kinetic energy. It does this in 
  the case of a 1D model.
  
  @param[in] grid
  @param[in] parameters
  @param[in,out] peakKETracking
  @param[in] procTop
  @param[in] time
  */
void writePeakKE_RTP(Grid &grid, Parameters &parameters, PeakKETracking &peakKETracking
  , ProcTop &procTop,Time &time);/**<
  This function calculates, and tracks the peak kinetic energy of the model, and averages it over
  three full periods. It does this by watch for 6 peaks in the kinetic energy and then averaging 
  them. It uses 6 because there are 2 peaks every period, when during contraction, and one during
  expansion. And writes the peak kinetic energy, the averaged peak kinetic energy. It does this in 
  the case of a 2D and 3D model.
  
  @param[in] grid
  @param[in] parameters
  @param[in,out] peakKETracking
  @param[in] procTop
  @param[in] time
  */
void finPeakKE(PeakKETracking &peakKETracking);/**<
  Finishes tracking the peak kinetic energy by flushing and closing the output files.
  
  @param[in,out] peakKETracking
  */
bool bFileExists(std::string sFilename);/**<
  Tests if the file exists by attempting to open the file for reading, if
  it fails it returns false, if it succeeds it returns true. This does not
  take into consideration permissions but that is ok for this project.
  
  @param[in] sFilename file name of the file to check if it exists or not
  @return returns true of false depending on weather the file exsists
*/
#endif
