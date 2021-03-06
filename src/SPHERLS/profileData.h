#ifndef PROFILEDATA_H
#define PROFILEDATA_H

/** @file
  
  Header file for \ref keepMax.cpp
*/

#include <string>
#include <vector>
#include <map>
#include <limits>
#include "time.h"
#include "procTop.h"
#include <fstream>

class profileData{
  public:
    void set(std::string sName,unsigned int nZone,double dValue);/**<
      Sets a new bit of data to dValue, identified by sName in radial zone nZone.
      */
    void set(std::string sName,unsigned int nZone,int nValue);/**<
      Sets a new bit of data to nValue, identified by sName in radial zone nZone.
      */
    void setSum(std::string sName,unsigned int nZone,double dValue);/**<
      If the value is already set it will add to it
      */
    void setSum(std::string sName,unsigned int nZone,int nValue);/**<
      If the value is already set it will add to it
      */
    void setAve(std::string sName,unsigned int nZone,double dValue);/**<
      If the value is already set it will add to it keeping track of the times
      it was added, to later compute an average when writing to a file.
      */
    void setAve(std::string sName,unsigned int nZone,int nValue);/**<
      If the value is already set it will add to it keeping track of the times
      it was added, to later compute an average when writing to a file.
      */
    void setMax(std::string sName,unsigned int nZone,double dValue);/**<
      If the value is already set it will set it to which ever is largest, the current value or the
      new value I am trying to set it to
      */
    void setMax(std::string sName,unsigned int nZone,int nValue);/**<
      If the value is already set it will set it to which ever is largest, the current value or the
      new value.
      */
    void setMin(std::string sName,unsigned int nZone,double dValue);/**<
      If the value is already set it will set it to which ever is smallest, the
      current value or the new value.
    */
    void setMin(std::string sName,unsigned int nZone,int nValue);/**<
      If the value is already set it will set it to which ever is smallest, the
      current value or the new value.
      */
    void setMaxAbs(std::string sName,unsigned int nZone,double dValue);/**<
      If the value is already set it will set it to which ever has the largest absolute value, the
      current value or the new value.
      */
    void setMaxAbs(std::string sName,unsigned int nZone,int nValue);/**<
      If the value is already set it will set it to which ever has the largest absolute value, the
      current value or the new value.
      */
    void toFile(std::string sFileName,Time time,ProcTop procTop);/**<
      Prints the data to a file in the same format as the radial profiles generated by SPHERLSanal
      */
    void clear();/**<
      Resets values to their initial values. It doesn't free any memory.
    */
    unsigned int nMaxNumZones();/**<
      Returns the maximum number of zones found under a key
      */
    bool test();/**<
      Runs a series of tests to insure that the functions are doing what they should be. Returns
      true if all tests passed, returns false other wise.
      */
    profileData();/**<
      Constructor for class
      */
  private:
    double dInitValue;
    int nInitValue;
    int nWidthIntField;
    int nWidthDoubleField;
    int nPrecision;
    std::map<std::string,std::vector<double> > dProfileData;
    std::map<std::string,std::vector<int> > nDoubleProfileDataCount;
    std::map<std::string,std::vector<int> > nProfileData;
    std::map<std::string,std::vector<int> > nIntegerProfileDataCount;
    int mergeOverLap(std::fstream& ifIn,int nFirstZone
      ,std::vector<std::string> sIntColumnNames, std::vector<std::string> sDoubleColumnNames
      ,ProcTop &procTop,Time& time);/**<
      returns the last line of the file pointed to by sFileName
    */
    int getFirstZoneInTable();/**<
      returns the first radial zone with a set value in the table, if no set values found it returns
      -1
    */
    bool bAlreadySetInt(std::string sName,unsigned int nZone);/**<
      Checks to see if an integer with name sName in zone nZone has already been set.
      */
    bool bAlreadySetDouble(std::string sName, unsigned int nZone);/**<
      Checks to see if a double with name sName in zone nZone has already been set.
      */
    bool inKeysInt(std::string sName);/**<
      Checks to see if sName is already in the set of integer keys, if so returns true else returns
      false.
      */
    bool inKeysDouble(std::string sName);/**<
      Checks to see if sName is already in the set of integer keys, if so returns true else returns
      false.
      */
    std::vector<std::string> getIntKeys();/**<
      returns a vector of strings containing the key names for the integer array
    */
    std::vector<std::string> getDoubleKeys();/**<
      returns a vector of strings containing the key names for the double array
    */
};/**@class profileData
  Class for tracking data on a radial profile, useful for debugging
*/

#endif
