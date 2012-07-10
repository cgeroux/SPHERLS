#include "fileExists.h"
#include <fstream>

bool bFileExists(std::string sFilename){
  std::ifstream ifTest;
  ifTest.open(sFilename.c_str(),std::ios::in);
  if(!ifTest){
    return false;//doesn't exsist
  }
  else{
    ifTest.close();
    return true;//does exsist
  }
}
