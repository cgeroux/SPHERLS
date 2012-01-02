#ifndef WATCHZONE_H
#define WATCHZONE_H

/** @file
  
  This file holds the definition of the watchzone class.
  
*/

#include <string>
#include <fstream>

class WatchZone{
  public:
  int i;
  int j;
  int k;
  std::string sOutFileName;
  
  WatchZone();
  WatchZone(int nIIn, int nJIn, int nKIn);
  WatchZone(int nIIn, int nJIn, int nKIn, std::string sFileName);
};/**@class WatchZone
  This class contains information used to monitor a particular zone of the grid.
*/

#endif
