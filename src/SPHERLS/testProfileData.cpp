#include "profileData.h"
#include <iostream>
#include <string>

int main(){
  profileData myKeepData;
  
  myKeepData.test();
  myKeepData.clear();
  myKeepData.toFile("test_ProfileData_cleared.txt",0.0);
  return 0;
}
