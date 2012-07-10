#ifndef FILEEXISTS_H
#define FILEEXISTS_H
#include <string>

bool bFileExists(std::string sFilename);/**<
  Tests if the file exists by attempting to open the file for reading, if
  it fails it returns false, if it succeeds it returns true. This does not
  take into consideration permissions but that is ok for this project.
  
  @param[in] sFilename file name of the file to check if it exists or not
  @return returns true of false depending on weather the file exsists
*/
#endif
