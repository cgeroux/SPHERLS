#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include "exception2.h"
#include "procTop.h"
#include "profileData.h"
#include "fileExists.h"
#include <cstring>
#include <stdlib.h>

profileData::profileData(){
  dInitValue=NAN;
  nInitValue=std::numeric_limits<int>::min();
  nWidthIntField=12;
  nWidthDoubleField=25;
  nPrecision=16;
}
void profileData::set(std::string sName,unsigned int nZone,double dValue){
  
  unsigned int nCurrentSize=dProfileData[sName].size();
  if (nZone>=nCurrentSize){
    
    unsigned int nCurrent=nCurrentSize;
    while(nCurrent<nZone+1){//add until we get to that zone
      dProfileData[sName].push_back(dInitValue);
      nCurrent++;
    }
    dProfileData[sName][nZone]=dValue;
  }
  else{
    dProfileData[sName][nZone]=dValue;
  }
}
void profileData::set(std::string sName,unsigned int nZone,int nValue){
  
  unsigned int nCurrentSize=nProfileData[sName].size();
  if (nZone>=nCurrentSize){
    
    unsigned int nCurrent=nCurrentSize;
    while(nCurrent<nZone+1){//add until we get to that zone
      nProfileData[sName].push_back(nInitValue);
      nCurrent++;
    }
    nProfileData[sName][nZone]=nValue;
  }
  else{
    nProfileData[sName][nZone]=nValue;
  }
}
void profileData::setSum(std::string sName,unsigned int nZone,double dValue){
  
  //check to see if the value is already set
  if(bAlreadySetDouble(sName,nZone)){//if already set just add to it
    dProfileData[sName][nZone]=dProfileData[sName][nZone]+dValue;
  }
  else{//if not already set just set it
    set(sName,nZone,dValue);
  }
}
void profileData::setSum(std::string sName,unsigned int nZone,int nValue){
  
  //check to see if the value is already set
  if(bAlreadySetInt(sName,nZone)){//if already set just add to it
    nProfileData[sName][nZone]=nProfileData[sName][nZone]+nValue;
  }
  else{//if not already set just set it
    set(sName,nZone,nValue);
  }
}
void profileData::setMax(std::string sName,unsigned int nZone,double dValue){
  
  //check to see if the value is already set
  if(bAlreadySetDouble(sName,nZone)){//if already set just add to it
    if(dValue>dProfileData[sName][nZone]){
      dProfileData[sName][nZone]=dValue;
    }
  }
  else{//if not already set just set it
    set(sName,nZone,dValue);
  }
}
void profileData::setMax(std::string sName,unsigned int nZone,int nValue){
  
  //check to see if the value is already set
  if(bAlreadySetInt(sName,nZone)){//if already set just add to it
    if(nValue>nProfileData[sName][nZone]){
      nProfileData[sName][nZone]=nValue;
    }
  }
  else{//if not already set just set it
    set(sName,nZone,nValue);
  }
}
void profileData::setMaxAbs(std::string sName,unsigned int nZone,double dValue){
  
  //check to see if the value is already set
  if(bAlreadySetDouble(sName,nZone)){//if already set just add to it
    if(fabs(dValue)>fabs(dProfileData[sName][nZone])){
      dProfileData[sName][nZone]=dValue;
    }
  }
  else{//if not already set just set it
    set(sName,nZone,dValue);
  }
}
void profileData::setMaxAbs(std::string sName,unsigned int nZone,int nValue){
  
  //check to see if the value is already set
  if(bAlreadySetInt(sName,nZone)){//if already set just add to it
    if(abs(nValue)>abs(nProfileData[sName][nZone])){
      nProfileData[sName][nZone]=nProfileData[sName][nZone]+nValue;
    }
  }
  else{//if not already set just set it
    set(sName,nZone,nValue);
  }
}
void profileData::toFile(std::string sFileName,Time time,ProcTop procTop){
  
  //package up local header
  
  //FOR INTEGER HEADERS
  
  //get local integer header string
  std::string sHeaderInt;
  sHeaderInt="Zone#";
  for(std::map<std::string, std::vector<int> >::iterator it = nProfileData.begin();
    it !=nProfileData.end(); ++it) {
    sHeaderInt+=" ";
    sHeaderInt+=it->first;
  }
  
  //get length of each header string on each processor
  int* nIntHeaderLengths=new int[procTop.nNumProcs];
  int nLocalIntHeaderLength=sHeaderInt.size();
  nIntHeaderLengths[procTop.nRank]=nLocalIntHeaderLength;
  int i;
  for(i=0;i<procTop.nNumProcs;i++){
    MPI::COMM_WORLD.Bcast(&nIntHeaderLengths[i],1,MPI::INT,i);
  }
  
  //get headers from each processor
  char** cHeadersInt=new char*[procTop.nNumProcs];
  for(i=0;i<procTop.nNumProcs;i++){
    cHeadersInt[i]=new char[nIntHeaderLengths[i]+1];
  }
  strcpy(cHeadersInt[procTop.nRank],sHeaderInt.c_str());
  for(i=0;i<procTop.nNumProcs;i++){
    MPI::COMM_WORLD.Bcast(cHeadersInt[i],nIntHeaderLengths[i]+1,MPI::CHAR,i);
  }
  
  
  //FOR DOUBLE HEADERS
  
  //get local double header string
  std::string sHeaderDouble;
  for(std::map<std::string, std::vector<double> >::iterator it = dProfileData.begin();
    it !=dProfileData.end(); ++it) {
    sHeaderDouble+=" ";
    sHeaderDouble+=it->first;
  }
  
  //get length of each header string on each processor
  int* nDoubleHeaderLengths=new int[procTop.nNumProcs];
  int nLocalDoubleHeaderLength=sHeaderDouble.size();
  nDoubleHeaderLengths[procTop.nRank]=nLocalDoubleHeaderLength;
  for(i=0;i<procTop.nNumProcs;i++){
    MPI::COMM_WORLD.Bcast(&nDoubleHeaderLengths[i],1,MPI::INT,i);
  }
  
  //get headers from each processor
  char** cHeadersDouble=new char*[procTop.nNumProcs];
  for(i=0;i<procTop.nNumProcs;i++){
    cHeadersDouble[i]=new char[nDoubleHeaderLengths[i]+1];
  }
  strcpy(cHeadersDouble[procTop.nRank],sHeaderDouble.c_str());
  for(i=0;i<procTop.nNumProcs;i++){
    MPI::COMM_WORLD.Bcast(cHeadersDouble[i],nDoubleHeaderLengths[i]+1,MPI::CHAR,i);
  }
  
  
  //get a list of all columns for integers
  std::vector<std::string> sIntColumnNames;
  int j;
  bool bAddToList;
  for(i=0;i<procTop.nNumProcs;i++){
    
    //convert character array in to a vector of strings
    std::stringstream ssIntHeader;
    ssIntHeader<<cHeadersInt[i];
    std::string sTemp;
    while(!ssIntHeader.eof()){//get entries from list
      ssIntHeader>>sTemp;
      
      //check to see if we have it already
      bAddToList=true;
      for (j=0;j<sIntColumnNames.size();j++){
        if(sIntColumnNames[j].compare(sTemp)==0){//found it in the list already
          bAddToList=false;
          break;
        }
      }
      if(bAddToList){
        sIntColumnNames.push_back(sTemp);
      }
    }
  }
  
  //get a list of all columns for doubles
  std::vector<std::string> sDoubleColumnNames;
  for(i=0;i<procTop.nNumProcs;i++){
    
    //convert character array in to a vector of strings
    std::stringstream ssDoubleHeader;
    ssDoubleHeader<<cHeadersDouble[i];
    std::string sTemp;
    while(!ssDoubleHeader.eof()){//get entries from list
      ssDoubleHeader>>sTemp;
      
      //check to see if we have it already
      bAddToList=true;
      for (j=0;j<sDoubleColumnNames.size();j++){
        if(sDoubleColumnNames[j].compare(sTemp)==0){//found it in the list already
          bAddToList=false;
          break;
        }
      }
      if(bAddToList&&sTemp.compare("")!=0){
        sDoubleColumnNames.push_back(sTemp);
      }
    }
  }
  
  //post a blocking receive from inner radial neighbour
  char cDummy;
  for(i=0;i<procTop.nNumRadialNeighbors;i++){
    if(procTop.nCoords[procTop.nRank][0]>procTop.nCoords[procTop.nRadialNeighborRanks[i]][0]){/*
      If the neighbour has a radial coordinate smaller than the current processor wait to receive
      a message from it*/
      MPI::COMM_WORLD.Recv(&cDummy,1,MPI::CHAR,procTop.nRadialNeighborRanks[i],2);
    }
  }
  
  //open file
  std::fstream ofOut;
  
  //set scientific format
  ofOut.unsetf(std::ios::fixed);
  ofOut.setf(std::ios::scientific);
  ofOut.precision(nPrecision);
  bool bAppend=false;
  int nLastZone=0;
  
  //check if the file already exists or not
  if(bFileExists(sFileName)){//if the file exists
    if(procTop.nRank!=0){//and we are not the first processor append to the end
      bAppend=true;
      
      //open for reading and writing (i.e. appending)
      ofOut.open(sFileName.c_str(),std::ios::out|std::ios::in);
      if(!ofOut.good()){
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": error opening the file \""
          <<sFileName<<"\"\n";
        throw exception2(ssTemp.str(),OUTPUT);
      }
      
      //figure out what is the first zone in table
      int nFirstZone=getFirstZoneInTable();
      
      //merge two tables in over lapping region
      nLastZone=mergeOverLap(ofOut,nFirstZone,sIntColumnNames, sDoubleColumnNames,procTop,time);
    }
  }
  
  if(!bAppend){//open a new file and add the headers
    
    //close in case it is already open
    ofOut.close();
    
    //open for writing only
    ofOut.open(sFileName.c_str(),std::ios::out);
    if(!ofOut.good()){
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": error opening the file \""
        <<sFileName<<"\"\n";
      throw exception2(ssTemp.str(),OUTPUT);
    }
    
    //print header
    ofOut<<"time= "<<time.dt<<" [s]"<<std::endl;
    
    //print all integer names
    for(i=0;i<sIntColumnNames.size();i++){
      ofOut<<std::setw(nWidthIntField)<<sIntColumnNames[i];
    }
    
    for(i=0;i<sDoubleColumnNames.size();i++){
      ofOut<<std::setw(nWidthDoubleField)<<sDoubleColumnNames[i];
    }
    ofOut<<std::endl;
  }
  
  //print out all data to a file, integers then doubles in the same order as the headers
  //starting with radial zone nLastZone+1
  unsigned int nNumZones=nMaxNumZones();
  if(nLastZone!=0){
    nLastZone++;
  }
  int nLastPrintedZone=0;
  for(unsigned int i=nLastZone;i<nNumZones;i++){
    
    //write zone number
    ofOut<<std::setw(nWidthIntField)<<i;
    
    //write all integer values
    for(int j=1;j<sIntColumnNames.size();j++){
      if(inKeysInt(sIntColumnNames[j])){
        if(i<nProfileData[sIntColumnNames[j]].size()){
          if(nProfileData[sIntColumnNames[j]][i]!=nInitValue){
            ofOut<<std::setw(nWidthIntField)<<nProfileData[sIntColumnNames[j]][i];
          }
          else{
            ofOut<<std::setw(nWidthIntField)<<"-";
          }
        }
        else{
          ofOut<<std::setw(nWidthIntField)<<"-";
        }
      }
      else{
        ofOut<<std::setw(nWidthIntField)<<"-";
      }
    }
    
    //write all double values
    for(int j=0;j<sDoubleColumnNames.size();j++){
      if(inKeysDouble(sDoubleColumnNames[j])){
        if(i<dProfileData[sDoubleColumnNames[j]].size()){
          if(!isnan(dProfileData[sDoubleColumnNames[j]][i])){
            ofOut<<std::setw(nWidthDoubleField)<<dProfileData[sDoubleColumnNames[j]][i];
          }
          else{
            ofOut<<std::setw(nWidthDoubleField)<<"-";
          }
        }
        else{
          ofOut<<std::setw(nWidthDoubleField)<<"-";
        }
      }
      else{
        ofOut<<std::setw(nWidthDoubleField)<<"-";
      }
    }
    ofOut<<std::endl;
    nLastPrintedZone=i;
  }
  ofOut.close();
  
  //post a blocking send to outer radial neighbour to let it know writing is done
  for(i=0;i<procTop.nNumRadialNeighbors;i++){
    if(procTop.nCoords[procTop.nRank][0]<procTop.nCoords[procTop.nRadialNeighborRanks[i]][0]){/*
      if current processor has a radial neighbour at inside post a receive*/
      MPI::COMM_WORLD.Send(&cDummy,1,MPI::CHAR,procTop.nRadialNeighborRanks[i],2);
    }
  }
}
void profileData::clear(){
  int i;
  int nSize;
  
  //reset integer values
  for(std::map<std::string, std::vector<int> >::iterator it = nProfileData.begin();
    it !=nProfileData.end(); ++it) {
    nSize=it->second.size();
    for(i=0;i<nSize;i++){
      it->second[i]=nInitValue;
    }
  }
  
  //reset double values
  for(std::map<std::string, std::vector<double> >::iterator it = dProfileData.begin();
    it !=dProfileData.end(); ++it) {
    nSize=it->second.size();
    for(i=0;i<nSize;i++){
      it->second[i]=dInitValue;
    }
  }
}
unsigned int profileData::nMaxNumZones(){
  unsigned int nMaxNumZones=0;
  for(std::map<std::string, std::vector<int> >::iterator it = nProfileData.begin();
    it !=nProfileData.end(); ++it) {
    if (it->second.size()>nMaxNumZones){
      nMaxNumZones=it->second.size();
    }
  }
  for(std::map<std::string, std::vector<double> >::iterator it = dProfileData.begin();
    it !=dProfileData.end(); ++it) {
    if (it->second.size()>nMaxNumZones){
      nMaxNumZones=it->second.size();
    }
  }
  return nMaxNumZones;
}
bool profileData::bAlreadySetInt(std::string sName, unsigned int nZone){
  if(inKeysInt(sName)){
    if(nZone<nProfileData[sName].size()){//if inside the array
      if (nProfileData[sName][nZone]!=nInitValue){
          return true;
      }
    }
  }
  return false;
}
bool profileData::bAlreadySetDouble(std::string sName, unsigned int nZone){
  if(inKeysDouble(sName)){
    if(nZone<dProfileData[sName].size()){//if inside the array
      if (!isnan(dProfileData[sName][nZone])){
        return true;
      }
    }
  }
  return false;
}
bool profileData::inKeysInt(std::string sName){
  for(std::map<std::string, std::vector<int> >::iterator it = nProfileData.begin();
    it !=nProfileData.end(); ++it) {
    if (sName.compare(it->first)==0){
      return true;
    }
  }
  return false;
}
bool profileData::inKeysDouble(std::string sName){
  for(std::map<std::string, std::vector<double> >::iterator it = dProfileData.begin();
    it !=dProfileData.end(); ++it) {
    if (sName.compare(it->first)==0){
      return true;
    }
  }
  return false;
}
std::vector<std::string> profileData::getIntKeys(){
  std::vector<std::string> sKeys;
  for(std::map<std::string, std::vector<int> >::iterator it = nProfileData.begin();
    it !=nProfileData.end(); ++it) {
    sKeys.push_back(it->first);
  }
  return sKeys;
}
std::vector<std::string> profileData::getDoubleKeys(){
  std::vector<std::string> sKeys;
  for(std::map<std::string, std::vector<double> >::iterator it = dProfileData.begin();
    it !=dProfileData.end(); ++it) {
    sKeys.push_back(it->first);
  }
  return sKeys;
}
bool profileData::test(){
  bool bAllTestsPassed=true;
  
  //test 0: test a simple set of a double
  set("stuff1",10,2.3);
  if (dProfileData["stuff1"][10]!=2.3){
    std::cout<<"test 0: couldn't set \"stuff1\" in zone 10 equal to 2.3"<<std::endl;
    bAllTestsPassed=false;
  }
  
  //test 1: test a simple set of a double in a different variable at a larger zone
  set("stuff2",12,2.5);
  if (dProfileData["stuff2"][12]!=2.5){
    std::cout<<"test 1: couldn't set \"stuff2\" in zone 12 equal to 2.5"<<std::endl;
    bAllTestsPassed=false;
  }
  
  //test 2: test a simple set of a double in the same variable in an earlier zone
  set("stuff1",8,2.7);
  if (dProfileData["stuff1"][8]!=2.7){
    std::cout<<"test 2: couldn't set \"stuff1\" in zone 8 equal to 2.7"<<std::endl;
    bAllTestsPassed=false;
  }
  
  //test 3: setting an integer
  set("stuff3",12,2);
  if (nProfileData["stuff3"][12]!=2){
    std::cout<<"test 3: couldn't set \"stuff3\" in zone 12 equal to 2"<<std::endl;
    bAllTestsPassed=false;
  }
  
  //test 4: setting an integer
  set("stuff3",14,4);
  if (nProfileData["stuff3"][14]!=4){
    std::cout<<"test 4: couldn't set \"stuff3\" in zone 14 equal to 4"<<std::endl;
    bAllTestsPassed=false;
  }
  
  //test 5: try testing to see if a variable has been set
  if(bAlreadySetInt("stuff3",14)!=true){
    std::cout<<"test 5: finding that \"stuff3\" in zone 14 hasn't been set when it has."<<std::endl;
    bAllTestsPassed=false;
  }
  
  //test 6: 
  if(bAlreadySetInt("stuff3",2)!=false){
    std::cout<<"test 6: finding that \"stuff3\" in zone 2 has been set when it hasn't."<<std::endl;
    bAllTestsPassed=false;
  }
  
  //test 7: test that I can do sums on integers already set
  setSum("stuff3",14,8);
  if(nProfileData["stuff3"][14]!=12){
    std::cout<<"test 7: Adding 8 to \"stuff3\" in zone 14 already set to 4 resulted in"
      <<dProfileData["stuff1"][14]<<" expecting 12"<<std::endl;
    bAllTestsPassed=false;
  }
  
  //test 8: test that I can do sums on doubles already set
  setSum("stuff1",8,6.2);
  if(dProfileData["stuff1"][8]!=8.9){//this is a little bit risky because of rounding
    std::cout<<"test 8: Adding 6.2 to \"stuff1\" in zone 8 already set to 2.7 resulted in"
      <<dProfileData["stuff1"][8]<<" expecting 8.9"<<std::endl;
    bAllTestsPassed=false;
  }
  
  //test 9: test that it finds that names are in keys properly
  if(!inKeysInt("stuff3")){
    std::cout<<"test 9: didn't find \"stuff3\" in integer keys while it has been set"<<std::endl;
  }
  
  //test 10: test that it finds that names are in keys properly
  if(inKeysInt("blob")){
    std::cout<<"test 10: found \"blob\" in integer keys while it has not been set"<<std::endl;
  }
  
  //test 11: test that it finds that names are in keys properly
  if(!inKeysDouble("stuff1")){
    std::cout<<"test 11: didn't find \"stuff1\" in double keys while it has been set"<<std::endl;
  }
  
  //test 12: test that it finds that names are in keys properly
  if(inKeysInt("blob2")){
    std::cout<<"test 12: found \"blob\" in double keys while it has not been set"<<std::endl;
  }
  
  //test 13: test setting max
  setMax("stuff1",8,9.9);
  if(dProfileData["stuff1"][8]!=9.9){
    std::cout<<"test 13: found \"stuff1\" in zone 8 had a value of "<<dProfileData["stuff1"][8]
      <<" expecting 9.9"<<std::endl;
  }
  
  //test 14: 
  setMax("stuff1",8,8.9);
  if(dProfileData["stuff1"][8]!=9.9){
    std::cout<<"test 14: found \"stuff1\" in zone 8 had a value of "<<dProfileData["stuff1"][8]
      <<" expecting 9.9"<<std::endl;
  }
  
  //at some point should test the contents of the file
  //toFile("test_profileData.txt",0.0);
  
  if(bAllTestsPassed){
    std::cout<<"all tests passed"<<std::endl;
    return true;
  }
  return false;
}
int profileData::getFirstZoneInTable(){
  int i;
  
  //get first double zone
  int nFirstZone=-1;
  for(std::map<std::string, std::vector<double> >::iterator it = dProfileData.begin();
    it !=dProfileData.end(); ++it) {
    for(i=0;i<it->second.size();i++){
      if(!isnan(it->second[i])){
        if(nFirstZone==-1||nFirstZone>i){
          nFirstZone=i;
          break;//break out of loop for this column
        }
      }
    }
  }
  
  //see if integer has a value earlier
  for(std::map<std::string, std::vector<int> >::iterator it = nProfileData.begin();
    it !=nProfileData.end(); ++it) {
    for(i=0;i<it->second.size();i++){
      if((it->second[i])!=nInitValue){
        if(nFirstZone==-1||nFirstZone>i){
          nFirstZone=i;
        }
        break;//break out of loop for this column
      }
    }
  }
  
  return nFirstZone;
}
int profileData::mergeOverLap(std::fstream& ifIn,int nFirstZone
  ,std::vector<std::string> sIntColumnNames, std::vector<std::string> sDoubleColumnNames
  ,ProcTop &procTop,Time& time){
  
  std::string sLine0="";
  std::stringstream ssTemp;
  ssTemp.unsetf(std::ios::fixed);
  ssTemp.setf(std::ios::scientific);
  ssTemp.precision(nPrecision);
  int nZoneNum=0;
  int nColumn;
  int nNumIntColumns=sIntColumnNames.size();
  int nNumDoubleColumns=sDoubleColumnNames.size();
  int nIntValue;
  int nTotalcolumns=nNumIntColumns+nNumDoubleColumns;
  char* cIntString=new char[nWidthIntField+1];
  char* cDoubleString=new char[nWidthDoubleField+1];
  char cTemp[]="c";//usef for reading new line
  bool bMoveToNextColumn;
  
  //read first couple lines, they are the time and headers
  std::getline(ifIn,sLine0);
  std::getline(ifIn,sLine0);
  
  while(!ifIn.eof()){
    
    if(nZoneNum>=nFirstZone-1){//they overlap
      
      ifIn.tellp();//seems to be need for some unclear reason, it seems to allow
      //switching form reading to writting.
      ifIn.read(cIntString,nWidthIntField);
      
      cIntString[nWidthIntField]='\0';
      ssTemp.str(cIntString);
      ssTemp>>nZoneNum;
      
      //read in column by column
      nColumn=1;
      if(ifIn.eof()){//got end of file read, first read after end of line character read, shoudl test for eof
        break;
      }
      
      while(nColumn<nTotalcolumns){
        
        //do integer columns
        if(nColumn<nNumIntColumns){
          
          //if local table has a set integer replace existing integer
          bMoveToNextColumn=true;
          std::vector<std::string> sIntKeys=getIntKeys();
          if(inKeysInt(sIntColumnNames[nColumn])){
            if(nProfileData[sIntColumnNames[nColumn]][nZoneNum]!=nInitValue){
              
              //replace value
              ssTemp.str("");
              ssTemp.clear();
              ssTemp<<std::setw(nWidthIntField)<<nProfileData[sIntColumnNames[nColumn]][nZoneNum];
              ifIn.tellp();/*seems to be need for some unclear reason, it seems to allow
                switching form reading to writting.*/
              ifIn.write(ssTemp.str().c_str(),nWidthIntField);
              bMoveToNextColumn=false;
            }
          }
          if(bMoveToNextColumn){
            
            //move file points to next column
            ifIn.tellp();/*seems to be need for some unclear reason, it seems to allow
              switching form reading to writting.*/
            ifIn.read(cIntString,nWidthIntField);
          }
        }
        else if(nColumn>=nNumIntColumns){
          
          //do double columns
          bMoveToNextColumn=true;
          if(inKeysDouble(sDoubleColumnNames[nColumn-nNumIntColumns])){
            if(!isnan(dProfileData[sDoubleColumnNames[nColumn-nNumIntColumns]][nZoneNum])){
              
              //replace value
              ssTemp.str("");
              ssTemp.clear();
              ssTemp<<std::setw(nWidthDoubleField)
                <<dProfileData[sDoubleColumnNames[nColumn-nNumIntColumns]][nZoneNum];
              ifIn.tellp();/*seems to be need for some unclear reason, it seems to allow
                switching form reading to writting.*/
              ifIn.write(ssTemp.str().c_str(),nWidthDoubleField);
              bMoveToNextColumn=false;
            }
          }
          if(bMoveToNextColumn){
            
            //move to next column
            ifIn.tellp();/*seems to be need for some unclear reason, it seems to allow
                switching form reading to writting.*/
            ifIn.read(cDoubleString,nWidthDoubleField);
            cDoubleString[nWidthDoubleField]='\0';
          }
        }
        nColumn++;
      }
      ifIn.tellp();/*seems to be need for some unclear reason, it seems to allow
        switching form reading to writting.*/
      ifIn.read(cTemp,1);//read in end line characater
    }
    else{//if before overlap region read lines and get zone number
      
      //get line
      std::getline(ifIn,sLine0);
      
      //get zone number
      ssTemp.str(sLine0);//set string
      ssTemp>>nZoneNum;
    }
  }
  ifIn.clear();//clear eof bit
  return nZoneNum;
}
