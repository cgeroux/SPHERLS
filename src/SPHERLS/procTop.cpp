/**
  @file
  
  Implementation file for the ProcTop class
  
*/

#include "procTop.h"
#include <cstring>
ProcTop::ProcTop(){
  //initialize
  nProcDims=NULL;
  nPeriodic=NULL;
  nCoords=NULL;
  nNumNeighbors=0;
  nNeighborRanks=NULL;
  nNumRadialNeighbors=0;
  nRadialNeighborRanks=NULL;
  nRadialNeighborNeighborIDs=NULL;
}
