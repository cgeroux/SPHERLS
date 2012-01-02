/**
  @file
  
  This file holds the implementation of the watchzone class.
*/
#include "watchzone.h"
#include "exception2.h"
#include <sstream>

WatchZone::WatchZone(){
  i=0;
  j=0;
  k=0;
  sOutFileName="watchzone.txt";
}
WatchZone::WatchZone(int nIIn, int nJIn, int nKIn){
  i=nIIn;
  j=nJIn;
  k=nKIn;
}
WatchZone::WatchZone(int nIIn, int nJIn, int nKIn, std::string sFileName){
  i=nIIn;
  j=nJIn;
  k=nKIn;
  sOutFileName=sFileName;
}
