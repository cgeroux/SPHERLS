/** @file
  
  Implements the eos (equation of state) class defined in \ref eos.h
*/
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <unistd.h>

#include "eos.h"
#include "exception2.h"

eos::eos(){//empty constructor
  nNumT=0;
  nNumRho=0;
  dLogP=NULL;
  dLogE=NULL;
  dLogKappa=NULL;
  setExePath();
}
eos& eos::operator=(const eos & rhs){//assignment operator
  if (this !=&rhs){
    //deallocate old memory
    for(int i=0;i<nNumRho;i++){
      delete [] dLogP[i];
      delete [] dLogE[i];
      delete [] dLogKappa[i];
    }
    delete [] dLogP;
    delete [] dLogE;
    delete [] dLogKappa;
    
    nNumT=rhs.nNumT;
    nNumRho=rhs.nNumRho;
    dLogRhoMin=rhs.dLogRhoMin;
    dLogTMin=rhs.dLogTMin;
    dLogRhoDelta=rhs.dLogRhoDelta;
    dLogTDelta=rhs.dLogTDelta;
    
    //allocate new memory and set it's values
    dLogP=new double*[nNumRho];
    dLogE=new double*[nNumRho];
    dLogKappa=new double*[nNumRho];
    for(int i=0;i<nNumRho;i++){
      dLogP[i]=new double[nNumT];
      dLogE[i]=new double[nNumT];
      dLogKappa[i]=new double[nNumT];
      for(int j=0;j<nNumT;j++){
        dLogP[i][j]=rhs.dLogP[i][j];
        dLogE[i][j]=rhs.dLogE[i][j];
        dLogKappa[i][j]=rhs.dLogKappa[i][j];
      }
    }
  }
  return *this;
}
eos::eos(const eos &ref){//copy constructor
  nNumRho=ref.nNumRho;
  nNumT=ref.nNumT;
  dXMassFrac=ref.dXMassFrac;
  dYMassFrac=ref.dYMassFrac;
  dLogRhoMin=ref.dLogRhoMin;
  dLogRhoDelta=ref.dLogRhoDelta;
  dLogTMin=ref.dLogTMin;
  dLogTDelta=ref.dLogTDelta;
  dLogP=new double*[nNumRho];
  dLogE=new double*[nNumRho];
  dLogKappa=new double*[nNumRho];
  for(int i=0;i<nNumRho;i++){
    dLogP[i]=new double[nNumT];
    dLogE[i]=new double[nNumT];
    dLogKappa[i]=new double[nNumT];
    for(int j=0;j<nNumT;j++){
      dLogP[i][j]=ref.dLogP[i][j];
      dLogE[i][j]=ref.dLogE[i][j];
      dLogKappa[i][j]=ref.dLogKappa[i][j];
    }
  }
}
eos::~eos(){//destructor
  for(int i=0;i<nNumRho;i++){
    delete [] dLogP[i];
    delete [] dLogE[i];
    delete [] dLogKappa[i];
  }
  delete [] dLogP;
  delete [] dLogE;
  delete [] dLogKappa;
}
void eos::readAscii(std::string sFileName)throw(exception2){
  
  //open file
  std::ifstream ifIn;
  ifIn.open(sFileName.c_str());
  if(!ifIn.is_open()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": error opening the file \""<<sFileName.c_str()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //delete memory
  for(int i=0;i<nNumRho;i++){
    delete [] dLogP[i];
    delete [] dLogE[i];
    delete [] dLogKappa[i];
  }
  delete [] dLogP;
  delete [] dLogE;
  delete [] dLogKappa;
  
  //read in sizes
  ifIn>>nNumRho>>nNumT;
  
  //read in file
  ifIn>>dXMassFrac>>dYMassFrac>>dLogRhoMin>>dLogRhoDelta>>dLogTMin>>dLogTDelta;
  
  //allocate space/read in file
  dLogP=new double*[nNumRho];
  dLogE=new double*[nNumRho];
  dLogKappa=new double*[nNumRho];
  std::string sLogP;
  std::string sLogE;
  std::string sLogKappa;
  char* cEnd;
  for(int i=0;i<nNumRho;i++){
    dLogP[i]=new double[nNumT];
    dLogE[i]=new double[nNumT];
    dLogKappa[i]=new double[nNumT];
    for(int j=0;j<nNumT;j++){
      ifIn>>sLogP>>sLogE>>sLogKappa;
      dLogP[i][j]=strtod(sLogP.c_str(),&cEnd);
      dLogE[i][j]=strtod(sLogE.c_str(),&cEnd);
      dLogKappa[i][j]=strtod(sLogKappa.c_str(),&cEnd);
      
      //check that reading went ok
      if(!ifIn.good()){
        std::cout<<"line="<<(i*j+2)<<std::endl;
        std::cout<<dLogP[i][j]<<" "<<dLogE[i][j]<<" "<<dLogKappa[i][j]<<std::endl;
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
          <<": error reading from file \""<<sFileName.c_str()<<"\"\n";
        throw exception2(ssTemp.str(),INPUT);
      }
    }
  }
  
  ifIn.close();
}
void eos::readBobsAscii(std::string sFileName)throw(exception2){
  
  //open file
  std::ifstream ifIn;
  ifIn.open(sFileName.c_str());
  if(!ifIn.is_open()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": error opening the file \""<<sFileName.c_str()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //delete memory
  for(int i=0;i<nNumRho;i++){
    delete [] dLogP[i];
    delete [] dLogE[i];
    delete [] dLogKappa[i];
  }
  delete [] dLogP;
  delete [] dLogE;
  delete [] dLogKappa;
  
  //Bob's format
  
  //read in composition
  ifIn>>dXMassFrac>>dYMassFrac;
  
  //read in size
  ifIn>>dLogTMin>>dLogTDelta>>dLogRhoMin>>dLogRhoDelta>>nNumT>>nNumRho;
  
  
  //allocate space/read in file
  dLogP=new double*[nNumRho];
  dLogE=new double*[nNumRho];
  dLogKappa=new double*[nNumRho];
  for(int i=0;i<nNumRho;i++){
    dLogP[i]=new double[nNumT];
    dLogE[i]=new double[nNumT];
    dLogKappa[i]=new double[nNumT];
    for(int j=0;j<nNumT;j++){
      ifIn>>dLogP[i][j]>>dLogE[i][j]>>dLogKappa[i][j];
    }
  }
  
  //check that reading went ok
  if(!ifIn.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": error reading from file \""<<sFileName.c_str()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  ifIn.close();
}
void eos::writeAscii(std::string sFileName)throw(exception2){
  
  //open file
  std::ofstream ofOut;
  ofOut.open(sFileName.c_str());
  if(!ofOut.is_open()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": error opening the file \""<<sFileName.c_str()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //write out ascii file
  ofOut.precision(16);
  ofOut.unsetf(std::ios::fixed);
  ofOut.setf(std::ios::scientific);
  ofOut<<nNumRho<<" "<<nNumT<<std::endl;
  ofOut<<dXMassFrac<<" "<<dYMassFrac<<" "<<dLogRhoMin<<" "<<dLogRhoDelta<<" "<<dLogTMin<<" "<<dLogTDelta<<" "<<std::endl;
  for(int i=0;i<nNumRho;i++){
    for(int j=0;j<nNumT;j++){
      ofOut<<dLogP[i][j]<<" "<<dLogE[i][j]<<" "<<dLogKappa[i][j]<<std::endl;
    }
  }
  ofOut.close();
}
void eos::readBin(std::string sFileName)throw(exception2){
  
  //test to see if it is relative to the executable directory
  std::string sTemp;
  std::stringstream ssTemp;
  if (sFileName.substr(0,1)!="/" 
    && sFileName.substr(0,2)!="./"){
    
    //if relative to executable directory add executable directory
    sTemp=sExePath+"/"+sFileName;
  }
  else{
    sTemp=sFileName;
  }
  sFileName=sTemp;
  
  //open file
  std::ifstream ifIn;
  ifIn.open(sFileName.c_str(),std::ios::binary);
  if(!ifIn.is_open()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": error opening the file \""<<sFileName.c_str()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //delete memory
  for(int i=0;i<nNumRho;i++){
    delete [] dLogP[i];
    delete [] dLogE[i];
    delete [] dLogKappa[i];
  }
  delete [] dLogP;
  delete [] dLogE;
  delete [] dLogKappa;
  
  //read in file
  ifIn.read((char*)(&nNumRho),sizeof(int));
  ifIn.read((char*)(&nNumT),sizeof(int));
  ifIn.read((char*)(&dXMassFrac),sizeof(double));
  ifIn.read((char*)(&dYMassFrac),sizeof(double));
  ifIn.read((char*)(&dLogRhoMin),sizeof(double));
  ifIn.read((char*)(&dLogRhoDelta),sizeof(double));
  ifIn.read((char*)(&dLogTMin),sizeof(double));
  ifIn.read((char*)(&dLogTDelta),sizeof(double));
  
  //allocate space/read in file
  dLogP=new double*[nNumRho];
  dLogE=new double*[nNumRho];
  dLogKappa=new double*[nNumRho];
  for(int i=0;i<nNumRho;i++){
    dLogP[i]=new double[nNumT];
    dLogE[i]=new double[nNumT];
    dLogKappa[i]=new double[nNumT];
    ifIn.read((char*)(dLogP[i]),nNumT*sizeof(double));
    ifIn.read((char*)(dLogE[i]),nNumT*sizeof(double));
    ifIn.read((char*)(dLogKappa[i]),nNumT*sizeof(double));
  }
  
  //check that reading went ok
  if(!ifIn.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": error reading from the file \""<<sFileName.c_str()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  ifIn.close();
}
void eos::writeBin(std::string sFileName)throw(exception2){
  
  //open file
  std::ofstream ofOut;
  ofOut.open(sFileName.c_str(),std::ios::binary);
  if(!ofOut.is_open()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": error opening the file \""<<sFileName.c_str()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  ofOut.write((char*)(&nNumRho),sizeof(int));
  ofOut.write((char*)(&nNumT),sizeof(int));
  ofOut.write((char*)(&dXMassFrac),sizeof(double));
  ofOut.write((char*)(&dYMassFrac),sizeof(double));
  ofOut.write((char*)(&dLogRhoMin),sizeof(double));
  ofOut.write((char*)(&dLogRhoDelta),sizeof(double));
  ofOut.write((char*)(&dLogTMin),sizeof(double));
  ofOut.write((char*)(&dLogTDelta),sizeof(double));
  for(int i=0;i<nNumRho;i++){
    ofOut.write((char*)(dLogP[i]),nNumT*sizeof(double));
    ofOut.write((char*)(dLogE[i]),nNumT*sizeof(double));
    ofOut.write((char*)(dLogKappa[i]),nNumT*sizeof(double));
  }
  ofOut.close();
}
double eos::dGetPressure(double dT, double dRho)throw(exception2){
  
  //calculate logs of dT and dRho
  if(dRho<0){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": density ("<<dRho<<") is negative.\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho-1)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT-1)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nILower and nIUpper
  double dRhoFrac=(dLogRho-dLogRhoLower)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated pressure
  double dP_j  =(dLogP[nIUpper][nJLower]-dLogP[nILower][nJLower])*dRhoFrac+dLogP[nILower][nJLower];
  double dP_jp1=(dLogP[nIUpper][nJUpper]-dLogP[nILower][nJUpper])*dRhoFrac+dLogP[nILower][nJUpper];
  double dP=pow(10.0,((dP_jp1-dP_j)*dTFrac+dP_j));
  if (std::isnan(dP)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the pressure at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  return dP;
}
double eos::dGetEnergy(double dT, double dRho)throw(exception2){
  
  //calculate logs of dT and dRho
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nILower and nIUpper
  double dRhoFrac=(dLogRho-dLogRhoLower)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated energy
  double dE_j  =(dLogE[nIUpper][nJLower]-dLogE[nILower][nJLower])*dRhoFrac+dLogE[nILower][nJLower];
  double dE_jp1=(dLogE[nIUpper][nJUpper]-dLogE[nILower][nJUpper])*dRhoFrac+dLogE[nILower][nJUpper];
  double dE=pow(10.0,((dE_jp1-dE_j)*dTFrac+dE_j));
  if (std::isnan(dE)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the energy at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  return dE;
}
double eos::dGetOpacity(double dT, double dRho)throw(exception2){
  
  //calculate logs of dT and dRho
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nILower and nIUpper
  double dRhoFrac=(dLogRho-dLogRhoLower)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated opacity
  double dKappa_j  =(dLogKappa[nIUpper][nJLower]-dLogKappa[nILower][nJLower])*dRhoFrac
    +dLogKappa[nILower][nJLower];
  double dKappa_jp1=(dLogKappa[nIUpper][nJUpper]-dLogKappa[nILower][nJUpper])*dRhoFrac
    +dLogKappa[nILower][nJUpper];
  double dKappa=pow(10.0,((dKappa_jp1-dKappa_j)*dTFrac+dKappa_j));
  if (std::isnan(dKappa)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the opacity at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more "
      <<"values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  return dKappa;
}
double eos::dDRhoDP(double dT,double dRho)throw(exception2){
  
  //calculate logs of dT and dRho
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated partial derivative of density with respect to pressure holding
  //temperature constant
  double dP_i  =(dLogP[nILower][nJUpper]-dLogP[nILower][nJLower])*dTFrac+dLogP[nILower][nJLower];
  double dP_ip1=(dLogP[nIUpper][nJUpper]-dLogP[nIUpper][nJLower])*dTFrac+dLogP[nIUpper][nJLower];
  double dDRhoDP=(pow(10.0,dLogRhoUpper)-pow(10.0,dLogRhoLower))/(pow(10.0,dP_ip1)-pow(10.0,dP_i));
  if (std::isnan(dDRhoDP)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for DRhoDP at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more values"
      <<" used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  return dDRhoDP;
}
double eos::dSoundSpeed(double dT,double dRho)throw(exception2){
  
  //calculate logs of dT and dRho
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's, density varys with i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's, temperature varys with j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nILower and nIUpper
  double dRhoFrac=(dLogRho-dLogRhoLower)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated pressures at upper and lower temperatures
  double dP_j  =(dLogP[nIUpper][nJLower]-dLogP[nILower][nJLower])*dRhoFrac+dLogP[nILower][nJLower];
  double dP_jp1=(dLogP[nIUpper][nJUpper]-dLogP[nILower][nJUpper])*dRhoFrac+dLogP[nILower][nJUpper];
  
  //calculate interpolated pressures at upper and lower densities
  double dP_i  =(dLogP[nILower][nJUpper]-dLogP[nILower][nJLower])*dTFrac+dLogP[nILower][nJLower];
  double dP_ip1=(dLogP[nIUpper][nJUpper]-dLogP[nIUpper][nJLower])*dTFrac+dLogP[nIUpper][nJLower];
  
  //calculate interpolated energy at upper and lower temperatures
  double dE_j  =(dLogE[nIUpper][nJLower]-dLogE[nILower][nJLower])*dRhoFrac+dLogE[nILower][nJLower];
  double dE_jp1=(dLogE[nIUpper][nJUpper]-dLogE[nILower][nJUpper])*dRhoFrac+dLogE[nILower][nJUpper];
  
  //calculate dlnP/dlnT at constant density
  double dDlnPDlnT=(dP_jp1-dP_j)/(dLogTUpper-dLogTLower);
  
  //calculate dlnP/dlnRho at constant temperature
  double dDlnPDlnRho=(dP_ip1-dP_i)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate dE/dT at constant density, equal to C_v (specific heat at constant volume)
  double dDEDT=(pow(10.0,dE_jp1)-pow(10.0,dE_j))/(pow(10.0,dLogTUpper)-pow(10.0,dLogTLower));
  
  //calculate interpolated pressure
  double dP=pow(10.0,((dP_jp1-dP_j)*dTFrac+dP_j));
  
  //calculate Gamma3 - 1
  double dGamma3m1=dP/(dRho*dT*dDEDT)*dDlnPDlnT;
  
  //calculate Gamma1
  double dGamma1=dDlnPDlnT*dGamma3m1+dDlnPDlnRho;
  
  //calculate speed of sound
  double dC=sqrt(dGamma1*dP/dRho);
  if (std::isnan(dC)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the sound speed at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or"
      <<" more values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  return dC;
}
void eos::getEKappa(double dT, double dRho, double &dE, double &dKappa)throw(exception2){
  
  //calculate logs of dT and dRho
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nILower and nIUpper
  double dRhoFrac=(dLogRho-dLogRhoLower)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated energy
  double dE_j  =(dLogE[nIUpper][nJLower]-dLogE[nILower][nJLower])*dRhoFrac+dLogE[nILower][nJLower];
  double dE_jp1=(dLogE[nIUpper][nJUpper]-dLogE[nILower][nJUpper])*dRhoFrac+dLogE[nILower][nJUpper];
  dE=pow(10.0,((dE_jp1-dE_j)*dTFrac+dE_j));
  if (std::isnan(dE)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the energy at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate interpolated opacity
  double dKappa_j  =(dLogKappa[nIUpper][nJLower]-dLogKappa[nILower][nJLower])*dRhoFrac
    +dLogKappa[nILower][nJLower];
  double dKappa_jp1=(dLogKappa[nIUpper][nJUpper]-dLogKappa[nILower][nJUpper])*dRhoFrac
    +dLogKappa[nILower][nJUpper];
  dKappa=pow(10.0,((dKappa_jp1-dKappa_j)*dTFrac+dKappa_j));
  if (std::isnan(dKappa)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the opacity at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
}
void eos::getPEKappa(double dT, double dRho, double &dP, double &dE, double &dKappa)throw(exception2){
  
  //calculate logs of dT and dRho
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nILower and nIUpper
  double dRhoFrac=(dLogRho-dLogRhoLower)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated pressure
  double dP_j  =(dLogP[nIUpper][nJLower]-dLogP[nILower][nJLower])*dRhoFrac+dLogP[nILower][nJLower];
  double dP_jp1=(dLogP[nIUpper][nJUpper]-dLogP[nILower][nJUpper])*dRhoFrac+dLogP[nILower][nJUpper];
  dP=pow(10.0,((dP_jp1-dP_j)*dTFrac+dP_j));
  if (std::isnan(dP)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the pressure at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate interpolated energy
  double dE_j  =(dLogE[nIUpper][nJLower]-dLogE[nILower][nJLower])*dRhoFrac+dLogE[nILower][nJLower];
  double dE_jp1=(dLogE[nIUpper][nJUpper]-dLogE[nILower][nJUpper])*dRhoFrac+dLogE[nILower][nJUpper];
  dE=pow(10.0,((dE_jp1-dE_j)*dTFrac+dE_j));
  if (std::isnan(dE)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the energy at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more "
      <<"values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate interpolated opacity
  double dKappa_j  =(dLogKappa[nIUpper][nJLower]-dLogKappa[nILower][nJLower])*dRhoFrac
    +dLogKappa[nILower][nJLower];
  double dKappa_jp1=(dLogKappa[nIUpper][nJUpper]-dLogKappa[nILower][nJUpper])*dRhoFrac
    +dLogKappa[nILower][nJUpper];
  dKappa=pow(10.0,((dKappa_jp1-dKappa_j)*dTFrac+dKappa_j));
  if (std::isnan(dKappa)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the opacity at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
}
void eos::getPEKappaGamma(double dT, double dRho, double &dP, double &dE, double &dKappa
  ,double &dGamma)throw(exception2){
  
  //calculate logs of dT and dRho
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nILower and nIUpper
  double dRhoFrac=(dLogRho-dLogRhoLower)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated log10 pressure at upper and lower temperatures
  double dP_j  =(dLogP[nIUpper][nJLower]-dLogP[nILower][nJLower])*dRhoFrac+dLogP[nILower][nJLower];
  double dP_jp1=(dLogP[nIUpper][nJUpper]-dLogP[nILower][nJUpper])*dRhoFrac+dLogP[nILower][nJUpper];
  
  //calculate interpolated log10 energy at upper and lower temperatures
  double dE_j  =(dLogE[nIUpper][nJLower]-dLogE[nILower][nJLower])*dRhoFrac+dLogE[nILower][nJLower];
  double dE_jp1=(dLogE[nIUpper][nJUpper]-dLogE[nILower][nJUpper])*dRhoFrac+dLogE[nILower][nJUpper];
  
  //calculate interpolated log10 opacity at upper and lower temperatures
  double dKappa_j  =(dLogKappa[nIUpper][nJLower]-dLogKappa[nILower][nJLower])*dRhoFrac
    +dLogKappa[nILower][nJLower];
  double dKappa_jp1=(dLogKappa[nIUpper][nJUpper]-dLogKappa[nILower][nJUpper])*dRhoFrac
    +dLogKappa[nILower][nJUpper];
  
  //calculate interpolated log pressures at upper and lower densities
  double dP_i  =(dLogP[nILower][nJUpper]-dLogP[nILower][nJLower])*dTFrac+dLogP[nILower][nJLower];
  double dP_ip1=(dLogP[nIUpper][nJUpper]-dLogP[nIUpper][nJLower])*dTFrac+dLogP[nIUpper][nJLower];
  
  //calculate dlnP/dlnT at constant density
  double dDlnPDlnT=(dP_jp1-dP_j)/(dLogTUpper-dLogTLower);
  
  //calculate dlnP/dlnRho at constant temperature
  double dDlnPDlnRho=(dP_ip1-dP_i)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate dE/dT at constant density, equal to C_v (specific heat at constant volume)
  double dDEDT=(pow(10.0,dE_jp1)-pow(10.0,dE_j))/(pow(10.0,dLogTUpper)-pow(10.0,dLogTLower));
  
  //calculate interpolated energy
  dE=pow(10.0,((dE_jp1-dE_j)*dTFrac+dE_j));
  if (std::isnan(dE)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the energy at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more values used in the"
      <<"interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate interpolated pressure
  dP=pow(10.0,((dP_jp1-dP_j)*dTFrac+dP_j));
  if (std::isnan(dP)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the pressure at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate interpolated opacity
  dKappa=pow(10.0,((dKappa_jp1-dKappa_j)*dTFrac+dKappa_j));
  if (std::isnan(dKappa)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the opacity at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate Gamma3-1
  double dGamma3m1=dP/(dRho*dT*dDEDT)*dDlnPDlnT;
  
  //calculate Gamma1
  dGamma=dDlnPDlnT*dGamma3m1+dDlnPDlnRho;
  if (std::isnan(dGamma)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the gamma at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
}
void eos::getPEKappaGammaCp(double dT, double dRho, double &dP, double &dE, double &dKappa
  ,double &dGamma, double &dC_p)throw(exception2){
  
  //calculate logs of dT and dRho
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nILower and nIUpper
  double dRhoFrac=(dLogRho-dLogRhoLower)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated log10 pressure at upper and lower temperatures
  double dP_j  =(dLogP[nIUpper][nJLower]-dLogP[nILower][nJLower])*dRhoFrac+dLogP[nILower][nJLower];
  double dP_jp1=(dLogP[nIUpper][nJUpper]-dLogP[nILower][nJUpper])*dRhoFrac+dLogP[nILower][nJUpper];
  
  //calculate interpolated log10 energy at upper and lower temperatures
  double dE_j  =(dLogE[nIUpper][nJLower]-dLogE[nILower][nJLower])*dRhoFrac+dLogE[nILower][nJLower];
  double dE_jp1=(dLogE[nIUpper][nJUpper]-dLogE[nILower][nJUpper])*dRhoFrac+dLogE[nILower][nJUpper];
  
  //calculate interpolated log10 opacity at upper and lower temperatures
  double dKappa_j  =(dLogKappa[nIUpper][nJLower]-dLogKappa[nILower][nJLower])*dRhoFrac
    +dLogKappa[nILower][nJLower];
  double dKappa_jp1=(dLogKappa[nIUpper][nJUpper]-dLogKappa[nILower][nJUpper])*dRhoFrac
    +dLogKappa[nILower][nJUpper];
  
  //calculate interpolated log pressures at upper and lower densities
  double dP_i  =(dLogP[nILower][nJUpper]-dLogP[nILower][nJLower])*dTFrac+dLogP[nILower][nJLower];
  double dP_ip1=(dLogP[nIUpper][nJUpper]-dLogP[nIUpper][nJLower])*dTFrac+dLogP[nIUpper][nJLower];
  
  //calculate dlnP/dlnT at constant density
  double dDlnPDlnT=(dP_jp1-dP_j)/(dLogTUpper-dLogTLower);
  
  //calculate dlnP/dlnRho at constant temperature
  double dDlnPDlnRho=(dP_ip1-dP_i)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate dE/dT at constant density, equal to C_v (specific heat at constant volume)
  double dDEDT=(pow(10.0,dE_jp1)-pow(10.0,dE_j))/(pow(10.0,dLogTUpper)-pow(10.0,dLogTLower));
  
  //calculate interpolated energy
  dE=pow(10.0,((dE_jp1-dE_j)*dTFrac+dE_j));
  if (std::isnan(dE)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the energy at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate interpolated pressure
  dP=pow(10.0,((dP_jp1-dP_j)*dTFrac+dP_j));
  if (std::isnan(dP)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the pressure at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate interpolated opacity
  dKappa=pow(10.0,((dKappa_jp1-dKappa_j)*dTFrac+dKappa_j));
  if (std::isnan(dKappa)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the opacity at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate Gamma3 - 1
  double dGamma3m1=dP/(dRho*dT*dDEDT)*dDlnPDlnT;
  
  //calculate Gamma1
  dGamma=dDlnPDlnT*dGamma3m1+dDlnPDlnRho;
  if (std::isnan(dGamma)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the gamma at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate dE/dT at constant density, equal to C_v (specific heat at constant volume)
  double dC_v=dE/dT*(dE_jp1-dE_j)/(dLogTUpper-dLogTLower);
  
  //calculate dE/dT at constant pressure, equal to C_p (specific heat at constant pressure)
  dC_p=dGamma*dC_v/dDlnPDlnRho;
  if (std::isnan(dC_p)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the specific heat at constant pressure at (rho,T)=("<<dRho<<","<<dT
      <<"), indicating that one or more values used in the interpolation are outside the calculated"
      <<" grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
}
void eos::getPKappaGamma(double dT, double dRho, double &dP, double &dKappa,double &dGamma)throw(exception2){
  
  //calculate logs of dT and dRho
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nILower and nIUpper
  double dRhoFrac=(dLogRho-dLogRhoLower)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated log10 pressure at upper and lower temperatures
  double dP_j  =(dLogP[nIUpper][nJLower]-dLogP[nILower][nJLower])*dRhoFrac+dLogP[nILower][nJLower];
  double dP_jp1=(dLogP[nIUpper][nJUpper]-dLogP[nILower][nJUpper])*dRhoFrac+dLogP[nILower][nJUpper];
  
  //calculate interpolated log10 energy at upper and lower temperatures
  double dE_j  =(dLogE[nIUpper][nJLower]-dLogE[nILower][nJLower])*dRhoFrac+dLogE[nILower][nJLower];
  double dE_jp1=(dLogE[nIUpper][nJUpper]-dLogE[nILower][nJUpper])*dRhoFrac+dLogE[nILower][nJUpper];
  
  //calculate interpolated log10 opacity at upper and lower temperatures
  double dKappa_j  =(dLogKappa[nIUpper][nJLower]-dLogKappa[nILower][nJLower])*dRhoFrac
    +dLogKappa[nILower][nJLower];
  double dKappa_jp1=(dLogKappa[nIUpper][nJUpper]-dLogKappa[nILower][nJUpper])*dRhoFrac
    +dLogKappa[nILower][nJUpper];
  
  //calculate interpolated log pressures at upper and lower densities
  double dP_i  =(dLogP[nILower][nJUpper]-dLogP[nILower][nJLower])*dTFrac+dLogP[nILower][nJLower];
  double dP_ip1=(dLogP[nIUpper][nJUpper]-dLogP[nIUpper][nJLower])*dTFrac+dLogP[nIUpper][nJLower];
  
  //calculate dlnP/dlnT at constant density
  double dDlnPDlnT=(dP_jp1-dP_j)/(dLogTUpper-dLogTLower);
  
  //calculate dlnP/dlnRho at constant temperature
  double dDlnPDlnRho=(dP_ip1-dP_i)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate dE/dT at constant density, equal to C_v (specific heat at constant volume)
  double dDEDT=(pow(10.0,dE_jp1)-pow(10.0,dE_j))/(pow(10.0,dLogTUpper)-pow(10.0,dLogTLower));
  
  //calculate interpolated energy
  //dE=pow(10.0,((dE_jp1-dE_j)*dTFrac+dE_j));
  
  //calculate interpolated pressure
  dP=pow(10.0,((dP_jp1-dP_j)*dTFrac+dP_j));
  if (std::isnan(dP)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the pressure at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate interpolated opacity
  dKappa=pow(10.0,((dKappa_jp1-dKappa_j)*dTFrac+dKappa_j));
  if (std::isnan(dKappa)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the opacity at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate Gamma3 - 1
  double dGamma3m1=dP/(dRho*dT*dDEDT)*dDlnPDlnT;
  
  //calculate Gamma1
  dGamma=dDlnPDlnT*dGamma3m1+dDlnPDlnRho;
  if (std::isnan(dGamma)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the gamma at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
}
void eos::gamma1DelAdC_v(double dT,double dRho,double &dGamma1, double &dDelAd,double &dC_v)throw(exception2){
  
  //calculate logs of dT and dRho
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's, density varys with i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's, temperature varys with j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nILower and nIUpper
  double dRhoFrac=(dLogRho-dLogRhoLower)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated pressures at upper and lower temperatures
  double dP_j  =(dLogP[nIUpper][nJLower]-dLogP[nILower][nJLower])*dRhoFrac+dLogP[nILower][nJLower];
  double dP_jp1=(dLogP[nIUpper][nJUpper]-dLogP[nILower][nJUpper])*dRhoFrac+dLogP[nILower][nJUpper];
  
  //calculate interpolated pressures at upper and lower densities
  double dP_i  =(dLogP[nILower][nJUpper]-dLogP[nILower][nJLower])*dTFrac+dLogP[nILower][nJLower];
  double dP_ip1=(dLogP[nIUpper][nJUpper]-dLogP[nIUpper][nJLower])*dTFrac+dLogP[nIUpper][nJLower];
  
  //calculate interpolated energy at upper and lower temperatures
  double dE_j  =(dLogE[nIUpper][nJLower]-dLogE[nILower][nJLower])*dRhoFrac+dLogE[nILower][nJLower];
  double dE_jp1=(dLogE[nIUpper][nJUpper]-dLogE[nILower][nJUpper])*dRhoFrac+dLogE[nILower][nJUpper];
  
  //calculate dlnP/dlnT at constant density
  double dDlnPDlnT=(dP_jp1-dP_j)/(dLogTUpper-dLogTLower);
  
  //calculate dlnP/dlnRho at constant temperature
  double dDlnPDlnRho=(dP_ip1-dP_i)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate interpolated pressure and energy
  double dP=pow(10.0,((dP_jp1-dP_j)*dTFrac+dP_j));
  double dE=pow(10.0,((dE_jp1-dE_j)*dTFrac+dE_j));
  
  //calculate dE/dT at constant density, equal to C_v (specific heat at constant volume)
  dC_v=dE/dT*(dE_jp1-dE_j)/(dLogTUpper-dLogTLower);
  if (std::isnan(dC_v)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the specific heat at constant volume at (rho,T)=("<<dRho<<","<<dT
      <<"), indicating that one or more values used in the interpolation are outside the "
      <<"calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate Gamma3 - 1
  double dGamma3m1=dP/(dRho*dT*dC_v)*dDlnPDlnT;
  
  //calculate Gamma1
  dGamma1=dDlnPDlnT*dGamma3m1+dDlnPDlnRho;
  if (std::isnan(dGamma1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the gamma at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate DelAd
  dDelAd=dGamma3m1/dGamma1;
  if (std::isnan(dDelAd)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the adiatabtic temperature gradient at (rho,T)=("<<dRho<<","<<dT
      <<"), indicating that one or more values used in the interpolation are outside the "
      <<"calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
}
void eos::getPAndDRhoDP(double dT,double dRho,double &dP, double &dDRhoDP)throw(exception2){
  
  //calculate logs of dT and dRho
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nILower and nIUpper
  double dRhoFrac=(dLogRho-dLogRhoLower)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated pressure
  double dLogP_i  =(dLogP[nILower][nJUpper]-dLogP[nILower][nJLower])*dTFrac+dLogP[nILower][nJLower];
  double dLogP_ip1=(dLogP[nIUpper][nJUpper]-dLogP[nIUpper][nJLower])*dTFrac+dLogP[nIUpper][nJLower];
  dDRhoDP=(pow(10.0,dLogRhoUpper)-pow(10.0,dLogRhoLower))/(pow(10.0,dLogP_ip1)-pow(10.0,dLogP_i));
  if (std::isnan(dDRhoDP)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the DRhoDP at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate interpolated pressure
  dP=pow(10.0,((dLogP_ip1-dLogP_i)*dRhoFrac+dLogP_i));
  if (std::isnan(dP)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the pressure at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
}
void eos::getEAndDTDE(double dT,double dRho,double &dE, double &dDTDE)throw(exception2){
  
  //calculate logs of dT and dRho
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nILower and nIUpper
  double dRhoFrac=(dLogRho-dLogRhoLower)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated energy
  double dLogE_j  =(dLogE[nIUpper][nJLower]-dLogE[nILower][nJLower])*dRhoFrac+dLogE[nILower][nJLower];
  double dLogE_jp1=(dLogE[nIUpper][nJUpper]-dLogE[nILower][nJUpper])*dRhoFrac+dLogE[nILower][nJUpper];
  dDTDE=(pow(10.0,dLogTUpper)-pow(10.0,dLogTLower))/(pow(10.0,dLogE_jp1)-pow(10.0,dLogE_j));
  if (std::isnan(dDTDE)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the DTDE at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate interpolated energy
  dE=pow(10.0,((dLogE_jp1-dLogE_j)*dTFrac+dLogE_j));
  if (std::isnan(dE)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": got nan for the energy at (rho,T)=("<<dRho<<","<<dT<<"), indicating that one or more"
      <<" values used in the interpolation are outside the calculated grid points.\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
}
void eos::getDlnPDlnTDlnPDlnPDEDT(double dT, double dRho, double &dDlnPDlnT,
  double &dDlnPDlnRho, double &dDEDT)throw(exception2){
  
  //calculate logs of dT and dRho
  double dLogRho=log10(dRho);
  double dLogT=log10(dT);
  
  //if density too low
  if(dLogRho<dLogRhoMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\" is lower than the minimum density in the table, \""<<dLogRhoMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //if temperature too low
  if(dLogT<dLogTMin){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\" is lower than the minimum log temperature in the table, \""<<dLogTMin<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate maximum values of grid
  double dLogRhoMax=dLogRhoMin+double(nNumRho)*dLogRhoDelta;
  double dLogTMax=dLogTMin+double(nNumT)*dLogTDelta;
  
  //calculate independent quantities at bracketing i's
  int nILower=int((dLogRho-dLogRhoMin)/dLogRhoDelta);
  int nIUpper=nILower+1;
  double dLogRhoLower=dLogRhoMin+double(nILower)*dLogRhoDelta;
  double dLogRhoUpper=dLogRhoMin+double(nIUpper)*dLogRhoDelta;
  
  //if density too high
  if(dLogRho>dLogRhoMax||nIUpper>(nNumRho-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log density to interpolate to, \""<<dLogRho
      <<"\"("<<nIUpper<<") is higher than the maximum density in the table, \""<<dLogRhoMax
      <<"\"("<<nNumRho-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate independent quantities at bracketing j's
  int nJLower=int((dLogT-dLogTMin)/dLogTDelta);
  int nJUpper=nJLower+1;
  double dLogTLower=dLogTMin+double(nJLower)*dLogTDelta;
  double dLogTUpper=dLogTMin+double(nJUpper)*dLogTDelta;
  
  //if temperature too high
  if(dLogT>dLogTMax||nJUpper>(nNumT-1)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the log temperature to interpolate to, \""<<dLogT
      <<"\"("<<nJUpper<<") is higher than the maximum temperature in the table, \""<<dLogTMax
      <<"\"("<<nNumT-1<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //calculate fractional distance between nILower and nIUpper
  double dRhoFrac=(dLogRho-dLogRhoLower)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate fractional distance between nJLower and nJUpper
  double dTFrac=(dLogT-dLogTLower)/(dLogTUpper-dLogTLower);
  
  //calculate interpolated log10 pressure at upper and lower temperatures
  double dP_j  =(dLogP[nIUpper][nJLower]-dLogP[nILower][nJLower])*dRhoFrac+dLogP[nILower][nJLower];
  double dP_jp1=(dLogP[nIUpper][nJUpper]-dLogP[nILower][nJUpper])*dRhoFrac+dLogP[nILower][nJUpper];
  
  //calculate interpolated log10 energy at upper and lower temperatures
  double dE_j  =(dLogE[nIUpper][nJLower]-dLogE[nILower][nJLower])*dRhoFrac+dLogE[nILower][nJLower];
  double dE_jp1=(dLogE[nIUpper][nJUpper]-dLogE[nILower][nJUpper])*dRhoFrac+dLogE[nILower][nJUpper];
  
  //calculate interpolated log10 opacity at upper and lower temperatures
  double dKappa_j  =(dLogKappa[nIUpper][nJLower]-dLogKappa[nILower][nJLower])*dRhoFrac
    +dLogKappa[nILower][nJLower];
  double dKappa_jp1=(dLogKappa[nIUpper][nJUpper]-dLogKappa[nILower][nJUpper])*dRhoFrac
    +dLogKappa[nILower][nJUpper];
  
  //calculate interpolated log pressures at upper and lower densities
  double dP_i  =(dLogP[nILower][nJUpper]-dLogP[nILower][nJLower])*dTFrac+dLogP[nILower][nJLower];
  double dP_ip1=(dLogP[nIUpper][nJUpper]-dLogP[nIUpper][nJLower])*dTFrac+dLogP[nIUpper][nJLower];
  
  //calculate dlnP/dlnT at constant density
  dDlnPDlnT=(dP_jp1-dP_j)/(dLogTUpper-dLogTLower);
  
  //calculate dlnP/dlnRho at constant temperature
  dDlnPDlnRho=(dP_ip1-dP_i)/(dLogRhoUpper-dLogRhoLower);
  
  //calculate dE/dT at constant density, equal to C_v (specific heat at constant volume)
  dDEDT=(pow(10.0,dE_jp1)-pow(10.0,dE_j))/(pow(10.0,dLogTUpper)-pow(10.0,dLogTLower));
}
void eos::setExePath(){
  /*This method might not be 100% portable, may need to look into other 
  solutions if problems arise with this not being reliable*/
  
  char buff[1024];
  ssize_t len = readlink("/proc/self/exe", buff, sizeof(buff)-1);
  if (len != -1) {
    buff[len] = '\0';
    sExePath=std::string(buff);
    
    //find the first "/" from the end
    unsigned pos=sExePath.find_last_of("/");
    
    //keep from the begging to the location of the last "/" to remove the name
    //of the executable
    sExePath=sExePath.substr(0,pos);
    
    //check to see if the last directory is "bin" if so remove that also
    //as installed versions put the exe's into the bin directory and sExePath
    //should point the top level directory.
    pos=sExePath.find_last_of("/");
    std::string sBin=sExePath.substr(pos+1,3);
    
    //if installed remove bin directory
    if(sBin.compare("bin")==0){
      sExePath=sExePath.substr(0,pos);
    }
    
    
  } else {
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": error determining executable path"<<std::endl;
    throw exception2(ssTemp.str(),OUTPUT);
  }
}
