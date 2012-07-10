/**
  @file
  
  This file holds functions used for examining the grid data during execution.
  This includes initializing structures, handlng watching zones during the execution of
  the program, opening files to write out the peak kinetic energy, etc.
*/

#include <mpi.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "watchzone.h"
#include "exception2.h"
#include "xmlFunctions.h"
#include "dataMonitoring.h"
#include "global.h"
#include "fileExists.h"

using namespace std;

void initWatchZones(XMLNode xParent,ProcTop &procTop, Grid &grid, Output &output
  , Parameters &parameters, Time &time){
  
  XMLNode xWatchZones=getXMLNodeNoThrow(xParent,"watchZones",0);
  
  if(!xWatchZones.isEmpty()){//if there are watch zones set
    int w=0;
    XMLNode xTemp=getXMLNodeNoThrow(xWatchZones,"watchZone",w);
    while(!xTemp.isEmpty()){
      
      //get zone location
      int nITemp;
      int nJTemp;
      int nKTemp;
      getXMLAttribute(xTemp,"x0",nITemp);
      getXMLAttribute(xTemp,"x1",nJTemp);
      getXMLAttribute(xTemp,"x2",nKTemp);
      
      //check that it is within global dimensions
      bool bOnGlobalGrid=true;
      int nGhostCells0=1;
      int nGhostCells1=0;
      int nGhostCells2=0;
      if(grid.nNumDims>1){
        nGhostCells1=1;
      }
      if(grid.nNumDims>2){
        nGhostCells2=1;
      }
      if(!(nITemp<(grid.nGlobalGridDims[0]+2*nGhostCells0*grid.nNumGhostCells))){//is it outside model?
        bOnGlobalGrid=false;
      }
      getXMLAttribute(xTemp,"x1",nJTemp);
      if(!(nJTemp<(grid.nGlobalGridDims[1]+2*nGhostCells1*grid.nNumGhostCells))){//is it outside model?
        bOnGlobalGrid=false;
      }
      getXMLAttribute(xTemp,"x2",nKTemp);
      if(!(nKTemp<(grid.nGlobalGridDims[2]+2*nGhostCells2*grid.nNumGhostCells))){//is it outside model?
        bOnGlobalGrid=false;
      }
      if(nITemp<grid.nNum1DZones+grid.nNumGhostCells){
        if(nKTemp!=0||nJTemp!=0){
          bOnGlobalGrid=false;
        }
      }
      
      //check that it is on current processor
      bool bOnLocalGrid=true;
      
      //check x0-direction
      int nStartX0=0;
      int nEndX0=0;
      for(int p=0;p<procTop.nNumProcs;p++){
        if( (procTop.nCoords[p][0]<procTop.nCoords[procTop.nRank][0])
          &&(procTop.nCoords[p][1]==procTop.nCoords[procTop.nRank][1]||procTop.nCoords[p][1]==-1)
          &&(procTop.nCoords[p][2]==procTop.nCoords[procTop.nRank][2]||procTop.nCoords[p][2]==-1)){
          nStartX0+=grid.nLocalGridDims[p][grid.nD][0];
        }
      }
      if(nStartX0!=0){//if not first proc in line, add ghost cells of inner boundary
        nStartX0+=grid.nNumGhostCells;
      }
      nEndX0=nStartX0+grid.nLocalGridDims[procTop.nRank][grid.nD][0];
      if(nStartX0==0){//if first proc in line, add ghost cells to end
        nEndX0+=grid.nNumGhostCells;
      }
      if(procTop.nCoords[procTop.nRank][0]==procTop.nProcDims[0]-1){//if last processor add outer ghost cells
        nEndX0+=grid.nNumGhostCells;
      }
      if(nITemp<nStartX0||nITemp>=nEndX0){
        bOnLocalGrid=false;
      }
      
      //check x1-direction
      int nStartX1=0;
      int nEndX1=0;
      for(int p=0;p<procTop.nNumProcs;p++){
        if( (procTop.nCoords[p][0]==procTop.nCoords[procTop.nRank][0]||procTop.nCoords[p][0]==-1)
          &&(procTop.nCoords[p][1]<procTop.nCoords[procTop.nRank][1])
          &&(procTop.nCoords[p][2]==procTop.nCoords[procTop.nRank][2]||procTop.nCoords[p][2]==-1)){
          nStartX1+=grid.nLocalGridDims[p][grid.nD][1];
        }
      }
      if(nStartX1!=0){//if not first proc in line, add ghost cells of inner boundary
        nStartX1+=grid.nNumGhostCells;
      }
      nEndX1=nStartX1+grid.nLocalGridDims[procTop.nRank][grid.nD][1];
      if(nStartX1==0){//if first proc in line, add ghost cells to end
        nEndX1+=grid.nNumGhostCells;
      }
      if(procTop.nCoords[procTop.nRank][1]==procTop.nProcDims[1]-1){//if last processor add outer ghost cells
        nEndX1+=grid.nNumGhostCells;
      }
      if(nJTemp<nStartX1||nJTemp>=nEndX1){
        bOnLocalGrid=false;
      }
      
      //check x2-direction
      int nStartX2=0;
      int nEndX2=0;
      for(int p=0;p<procTop.nNumProcs;p++){
        if( (procTop.nCoords[p][0]==procTop.nCoords[procTop.nRank][0]||procTop.nCoords[p][0]==-1)
          &&(procTop.nCoords[p][1]==procTop.nCoords[procTop.nRank][1]||procTop.nCoords[p][0]==-1)
          &&(procTop.nCoords[p][2]<procTop.nCoords[procTop.nRank][2])){
          nStartX2+=grid.nLocalGridDims[p][grid.nD][2];
        }
      }
      if(nStartX2!=0){//if not first proc in line, add ghost cells of inner boundary
        nStartX2+=grid.nNumGhostCells;
      }
      nEndX2=nStartX2+grid.nLocalGridDims[procTop.nRank][grid.nD][2];
      if(nStartX2==0){//if first proc in line, add ghost cells to end
        nEndX2+=grid.nNumGhostCells;
      }
      if(procTop.nCoords[procTop.nRank][2]==procTop.nProcDims[2]-1){//if last processor add outer ghost cells
        nEndX2+=grid.nNumGhostCells;
      }
      if(nKTemp<nStartX2||nKTemp>=nEndX2){
        bOnLocalGrid=false;
      }
      
      /*Note: the above if's don't work if the model is 1D*/
      if(bOnGlobalGrid && bOnLocalGrid){//if on local and global grid add to list
        
        //set output file name
        std::stringstream ssTemp;
        ssTemp<<output.sBaseOutputFileName<<"_watchZone_r"<<nITemp<<"_t"<<nJTemp<<"_p"<<nKTemp<<".txt";
        
        //adjust position to fall on local grid
        if(nStartX0!=0){
          nStartX0-=grid.nNumGhostCells;
        }
        if(nStartX1!=0){
          nStartX1-=grid.nNumGhostCells;
        }
        if(nStartX2!=0){
          nStartX2-=grid.nNumGhostCells;
        }
        nITemp-=nStartX0;
        nJTemp-=nStartX1;
        nKTemp-=nStartX2;
        
        //save watch zone
        output.watchzoneList.push_back(WatchZone (nITemp,nJTemp,nKTemp,ssTemp.str()));
      }
      if(!bOnGlobalGrid){
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
        <<": WARNING the watch zone ("<<nITemp<<","<<nJTemp<<","<<nKTemp<<") is not being included "
        "since it isn't on the global grid.\n";
      }
      w++;
      
      //get next watch zone
      xTemp=getXMLNodeNoThrow(xWatchZones,"watchZone",w);
    }
  }
  
  //open files for watch zones
  output.ofWatchZoneFiles=new std::ofstream[output.watchzoneList.size()];
  for(unsigned int i=0;i<output.watchzoneList.size();i++){
    bool bAppend=bFileExists(output.watchzoneList[i].sOutFileName);
    if(time.nTimeStepIndex!=0&&bAppend){//append to end of file
      //open file and go to the start of the line for the current time step
      output.ofWatchZoneFiles[i].open(output.watchzoneList[i].sOutFileName.c_str(),ios::in|ios::out);
      output.ofWatchZoneFiles[i].seekp((time.nTimeStepIndex+2)*(9+23*23)*sizeof(char));
    }
    else{//open a new file
      output.ofWatchZoneFiles[i].open(output.watchzoneList[i].sOutFileName.c_str());
    }
    if(!output.ofWatchZoneFiles[i].good()){//didn't open properly
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
        <<": error opening the file \""<<output.watchzoneList[i].sOutFileName<<"\"\n";
      throw exception2(ssTemp.str(),OUTPUT);
    }
    
    output.ofWatchZoneFiles[i].precision(14);
    output.ofWatchZoneFiles[i].setf(std::ios::scientific);
    int nWidthOutputField=23;
    if(time.nTimeStepIndex==0||!bAppend){//write out file header
      stringstream ssHeader;
      ssHeader<<"zone= ("<<output.watchzoneList[i].i<<","<<output.watchzoneList[i].j<<","
      <<output.watchzoneList[i].k<<")";
      if(parameters.bEOSGammaLaw){
        ssHeader<<" gamma= "<<parameters.dGamma;
      }
      output.ofWatchZoneFiles[i]<<std::setw(8+23*23)<<std::left<<ssHeader.str()<<std::endl;
      output.ofWatchZoneFiles[i]<<std::right;
      output.ofWatchZoneFiles[i]
        <<std::setw(8)<<"j(1)"
        <<std::setw(nWidthOutputField)<<"t[s](2)"
        <<std::setw(nWidthOutputField)<<"u_ip1half[cm/s](3)"
        <<std::setw(nWidthOutputField)<<"u_im1half[cm/s](4)"
        <<std::setw(nWidthOutputField)<<"u_0_ip1half[cm/s](5)"
        <<std::setw(nWidthOutputField)<<"u_0_im1half[cm/s](6)"
        <<std::setw(nWidthOutputField)<<"q0[dyne/cm^2](7)"
        <<std::setw(nWidthOutputField)<<"v_jp1half[cm/s](8)"
        <<std::setw(nWidthOutputField)<<"v_jm1half[cm/s](9)"
        <<std::setw(nWidthOutputField)<<"q1[dyne/cm^2](10)"
        <<std::setw(nWidthOutputField)<<"w_kp1half[cm/s](11)"
        <<std::setw(nWidthOutputField)<<"w_km1half[cm/s](12)"
        <<std::setw(nWidthOutputField)<<"q2[dyne/cm^2](13)"
        <<std::setw(nWidthOutputField)<<"R_ip1half[cm](14)"
        <<std::setw(nWidthOutputField)<<"R_im1half[cm](15)"
        <<std::setw(nWidthOutputField)<<"Density[g/cm^3](16)"
        <<std::setw(nWidthOutputField)<<"Den_ave[g/cm^3](17)"
        <<std::setw(nWidthOutputField)<<"Eddy_Visc(18)"
        <<std::setw(nWidthOutputField)<<"E[erg/g](19)"
        <<std::setw(nWidthOutputField)<<"P[dyne/cm^2](20)"
        <<std::setw(nWidthOutputField)<<"T[K](21)"
        <<std::setw(nWidthOutputField)<<"DMr(t=0)[g](22)"
        <<std::setw(nWidthOutputField)<<"Del_MCalc[g](23)"
        <<std::setw(nWidthOutputField)<<"Rel_Error_Del_M(24)"
        <<std::endl;
    }
  }
}
void writeWatchZones_R_GL(Output &output, Grid &grid, Parameters &parameters, Time &time
  , ProcTop &procTop){
  for(unsigned int i=0;i<output.watchzoneList.size();i++){
    
    //calculate zone mass
    double dM=grid.dLocalGridOld[grid.nDM][output.watchzoneList[i].i][0][0];
    double dMCalc=4.0/3.0*parameters.dPi*grid.dLocalGridOld[grid.nD][output.watchzoneList[i].i][0][0]
      *(pow(grid.dLocalGridOld[grid.nR][output.watchzoneList[i].i+grid.nCenIntOffset[0]][0][0],3.0)
      -pow(grid.dLocalGridOld[grid.nR][output.watchzoneList[i].i-1+grid.nCenIntOffset[0]][0][0],3.0));
    
    int nWidthOutputField=23;
    int nIInt=output.watchzoneList[i].i+grid.nCenIntOffset[0];
    output.ofWatchZoneFiles[i]
      <<std::setw(8)<<time.nTimeStepIndex
      <<std::setw(nWidthOutputField)<<time.dt
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU][nIInt][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU][nIInt-1][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nQ0][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<"-"//V_ijp1halfk
      <<std::setw(nWidthOutputField)<<"-"//V_ijm1halfk
      <<std::setw(nWidthOutputField)<<"-"//Q1
      <<std::setw(nWidthOutputField)<<"-"//W_ijkp1half
      <<std::setw(nWidthOutputField)<<"-"//W_ijkm1half
      <<std::setw(nWidthOutputField)<<"-"//Q2
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nR][nIInt][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nD][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<"-"
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nD][output.watchzoneList[i].i][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nE][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nP][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<"-"//T
      <<std::setw(nWidthOutputField)<<dM
      <<std::setw(nWidthOutputField)<<dMCalc
      <<std::setw(nWidthOutputField)<<(dM-dMCalc)/dM
      <<std::endl;
  }
}
void writeWatchZones_R_TEOS(Output &output, Grid &grid, Parameters &parameters, Time &time
  , ProcTop &procTop){
  for(unsigned int i=0;i<output.watchzoneList.size();i++){
    
    //calculate zone mass
    double dM=grid.dLocalGridOld[grid.nDM][output.watchzoneList[i].i][0][0];
    double dMCalc=4.0/3.0*parameters.dPi*grid.dLocalGridOld[grid.nD][output.watchzoneList[i].i][0][0]
      *(pow(grid.dLocalGridOld[grid.nR][output.watchzoneList[i].i+grid.nCenIntOffset[0]][0][0],3.0)
      -pow(grid.dLocalGridOld[grid.nR][output.watchzoneList[i].i-1+grid.nCenIntOffset[0]][0][0],3.0));
    
    int nWidthOutputField=23;
    int nIInt=output.watchzoneList[i].i+grid.nCenIntOffset[0];
    output.ofWatchZoneFiles[i]
      <<std::setw(8)<<time.nTimeStepIndex
      <<std::setw(nWidthOutputField)<<time.dt
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU][nIInt][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU][nIInt-1][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nQ0][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<"-"//V_ijp1halfk
      <<std::setw(nWidthOutputField)<<"-"//V_ijm1halfk
      <<std::setw(nWidthOutputField)<<"-"//Q1
      <<std::setw(nWidthOutputField)<<"-"//W_ijkp1half
      <<std::setw(nWidthOutputField)<<"-"//W_ijkm1half
      <<std::setw(nWidthOutputField)<<"-"//Q2
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nR][nIInt][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nD][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nD][output.watchzoneList[i].i][0][0]
      <<std::setw(nWidthOutputField)<<"-"
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nE][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nP][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)
        <<grid.dLocalGridOld[grid.nT][output.watchzoneList[i].i][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<dM
      <<std::setw(nWidthOutputField)<<dMCalc
      <<std::setw(nWidthOutputField)<<(dM-dMCalc)/dM
      <<std::endl;
  }
}
void writeWatchZones_RT_GL(Output &output, Grid &grid, Parameters &parameters, Time &time
  , ProcTop &procTop){
  for(unsigned int i=0;i<output.watchzoneList.size();i++){
    
    //calculate zone mass
    double dM=grid.dLocalGridOld[grid.nDM][output.watchzoneList[i].i][0][0];
    double dMCalc=4.0/3.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][output.watchzoneList[i].i][0][0]
      *(pow(grid.dLocalGridOld[grid.nR][output.watchzoneList[i].i+grid.nCenIntOffset[0]][0][0],3.0)
      -pow(grid.dLocalGridOld[grid.nR][output.watchzoneList[i].i-1+grid.nCenIntOffset[0]][0][0],3.0));
    
    int nWidthOutputField=23;
    int nIInt=output.watchzoneList[i].i+grid.nCenIntOffset[0];
    int nJInt=output.watchzoneList[i].j+grid.nCenIntOffset[1];
    output.ofWatchZoneFiles[i]
      <<std::setw(8)<<time.nTimeStepIndex
      <<std::setw(nWidthOutputField)<<time.dt
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU][nIInt][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU][nIInt-1][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nQ0][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nV][output.watchzoneList[i].i][nJInt][output.watchzoneList[i].k];
    if(nJInt-1<0){
      output.ofWatchZoneFiles[i]
        <<std::setw(nWidthOutputField)<<"-";//doesn't work when at inner edge
    }
    else{
      output.ofWatchZoneFiles[i]
        <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nV][output.watchzoneList[i].i][nJInt-1][output.watchzoneList[i].k];
    }
    output.ofWatchZoneFiles[i]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nQ1][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<"-"//W_ijkp1half
      <<std::setw(nWidthOutputField)<<"-"//W_ijkm1half
      <<std::setw(nWidthOutputField)<<"-"//Q2
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nR][nIInt][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nD][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k];
    if(parameters.nTypeTurbulanceMod>0){
      output.ofWatchZoneFiles[i]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nEddyVisc][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k];
    }
    else{
      output.ofWatchZoneFiles[i]<<std::setw(nWidthOutputField)<<"-";
    }
    output.ofWatchZoneFiles[i]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nDenAve][output.watchzoneList[i].i][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nE][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nP][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<"-"//T
      <<std::setw(nWidthOutputField)<<dM
      <<std::setw(nWidthOutputField)<<dMCalc
      <<std::setw(nWidthOutputField)<<(dM-dMCalc)/dM
      <<std::endl;
  }
}
void writeWatchZones_RT_TEOS(Output &output, Grid &grid, Parameters &parameters, Time &time
  , ProcTop &procTop){
  for(unsigned int i=0;i<output.watchzoneList.size();i++){
    
    //calculate zone mass
    double dM=grid.dLocalGridOld[grid.nDM][output.watchzoneList[i].i][0][0];
    double dMCalc=4.0/3.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][output.watchzoneList[i].i][0][0]
      *(pow(grid.dLocalGridOld[grid.nR][output.watchzoneList[i].i+grid.nCenIntOffset[0]][0][0],3.0)
      -pow(grid.dLocalGridOld[grid.nR][output.watchzoneList[i].i-1+grid.nCenIntOffset[0]][0][0],3.0));
    
    int nWidthOutputField=23;
    int nIInt=output.watchzoneList[i].i+grid.nCenIntOffset[0];
    int nJInt=output.watchzoneList[i].j+grid.nCenIntOffset[1];
    output.ofWatchZoneFiles[i]
      <<std::setw(8)<<time.nTimeStepIndex
      <<std::setw(nWidthOutputField)<<time.dt
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nU][nIInt][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nU][nIInt-1][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nQ0][output.watchzoneList[i].i][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nV][output.watchzoneList[i].i][nJInt][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField);
    if(nJInt-1<0){
      output.ofWatchZoneFiles[i]
        <<"-";//doesn't work when at inner edge
    }
    else{
      output.ofWatchZoneFiles[i]
        <<grid.dLocalGridOld[grid.nV][output.watchzoneList[i].i][nJInt-1][output.watchzoneList[i].k];
    }
    output.ofWatchZoneFiles[i]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nQ1][output.watchzoneList[i].i][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<"-"//W_ijkp1half
      <<std::setw(nWidthOutputField)<<"-"//W_ijkm1half
      <<std::setw(nWidthOutputField)<<"-"//Q2
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nR][nIInt][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nD][output.watchzoneList[i].i][output.watchzoneList[i].j]
      [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nDenAve][output.watchzoneList[i].i][0][0];
    if(parameters.nTypeTurbulanceMod>0){
      output.ofWatchZoneFiles[i]
        <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nEddyVisc]
        [output.watchzoneList[i].i][output.watchzoneList[i].j][output.watchzoneList[i].k];
    }
    else{
      output.ofWatchZoneFiles[i]
        <<std::setw(nWidthOutputField)<<"-";
    }
    output.ofWatchZoneFiles[i]<<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nE][output.watchzoneList[i].i][output.watchzoneList[i].j]
      [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nP][output.watchzoneList[i].i][output.watchzoneList[i].j]
      [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nT][output.watchzoneList[i].i][output.watchzoneList[i].j]
      [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<dM
      <<std::setw(nWidthOutputField)<<dMCalc
      <<std::setw(nWidthOutputField)<<(dM-dMCalc)/dM
      <<std::endl;
  }
}
void writeWatchZones_RTP_GL(Output &output, Grid &grid, Parameters &parameters, Time &time
  , ProcTop &procTop){
  for(unsigned int i=0;i<output.watchzoneList.size();i++){
    
    //calculate zone mass
    double dM=grid.dLocalGridOld[grid.nDM][output.watchzoneList[i].i][0][0];
    double dMCalc=4.0/3.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][output.watchzoneList[i].i][0][0]
      *(pow(grid.dLocalGridOld[grid.nR][output.watchzoneList[i].i+grid.nCenIntOffset[0]][0][0],3.0)
      -pow(grid.dLocalGridOld[grid.nR][output.watchzoneList[i].i-1+grid.nCenIntOffset[0]][0][0],3.0));
    
    int nWidthOutputField=23;
    int nIInt=output.watchzoneList[i].i+grid.nCenIntOffset[0];
    int nJInt=output.watchzoneList[i].j+grid.nCenIntOffset[1];
    int nKInt=output.watchzoneList[i].k+grid.nCenIntOffset[2];
    output.ofWatchZoneFiles[i]
      <<std::setw(8)<<time.nTimeStepIndex
      <<std::setw(nWidthOutputField)<<time.dt
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU][nIInt][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU][nIInt-1][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nQ0][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nV][output.watchzoneList[i].i][nJInt][output.watchzoneList[i].k];
    if(nJInt-1==-1){
      output.ofWatchZoneFiles[i]
        <<std::setw(nWidthOutputField)<<"-";//not defined when on edge
    }
    else{
      output.ofWatchZoneFiles[i]
        <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nV][output.watchzoneList[i].i][nJInt-1][output.watchzoneList[i].k];
    }
    output.ofWatchZoneFiles[i]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nQ1][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nW][output.watchzoneList[i].i][output.watchzoneList[i].j][nKInt];
    if(nKInt-1==-1){
      output.ofWatchZoneFiles[i]
        <<std::setw(nWidthOutputField)<<"-";//not defined when on edge
    }
    else{
      output.ofWatchZoneFiles[i]
        <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nW][output.watchzoneList[i].i][output.watchzoneList[i].j][nKInt-1];
    }
    output.ofWatchZoneFiles[i]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nQ2][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nR][nIInt][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nD][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k];
    if(parameters.nTypeTurbulanceMod>0){
      output.ofWatchZoneFiles[i]
        <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nEddyVisc]
        [output.watchzoneList[i].i][output.watchzoneList[i].j][output.watchzoneList[i].k];
    }
    else{
      output.ofWatchZoneFiles[i]
        <<std::setw(nWidthOutputField)<<"-";
    }
    output.ofWatchZoneFiles[i]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nDenAve][output.watchzoneList[i].i][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nE][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nP][output.watchzoneList[i].i][output.watchzoneList[i].j]
        [output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<"-"//T
      <<std::setw(nWidthOutputField)<<dM
      <<std::setw(nWidthOutputField)<<dMCalc
      <<std::setw(nWidthOutputField)<<(dM-dMCalc)/dM
      <<std::endl;
  }
}
void writeWatchZones_RTP_TEOS(Output &output, Grid &grid, Parameters &parameters, Time &time
  , ProcTop &procTop){
  for(unsigned int i=0;i<output.watchzoneList.size();i++){
    
    //calculate zone mass
    double dM=grid.dLocalGridOld[grid.nDM][output.watchzoneList[i].i][0][0];
    double dMCalc=4.0/3.0*parameters.dPi
      *grid.dLocalGridOld[grid.nDenAve][output.watchzoneList[i].i][0][0]
      *(pow(grid.dLocalGridOld[grid.nR][output.watchzoneList[i].i+grid.nCenIntOffset[0]][0][0],3.0)
      -pow(grid.dLocalGridOld[grid.nR][output.watchzoneList[i].i-1
      +grid.nCenIntOffset[0]][0][0],3.0));
    
    int nWidthOutputField=23;
    int nIInt=output.watchzoneList[i].i+grid.nCenIntOffset[0];
    int nJInt=output.watchzoneList[i].j+grid.nCenIntOffset[1];
    int nKInt=output.watchzoneList[i].k+grid.nCenIntOffset[2];
    output.ofWatchZoneFiles[i]
      <<std::setw(8)<<time.nTimeStepIndex
      <<std::setw(nWidthOutputField)<<time.dt
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nU][nIInt][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nU][nIInt-1][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nQ0][output.watchzoneList[i].i][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nV][output.watchzoneList[i].i][nJInt][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField);
    if(nJInt-1==-1){
      output.ofWatchZoneFiles[i]<<"-";//not defined when on edge
    }
    else{
      output.ofWatchZoneFiles[i]
        <<grid.dLocalGridOld[grid.nV][output.watchzoneList[i].i][nJInt-1][output.watchzoneList[i].k];
    }
    output.ofWatchZoneFiles[i]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nQ1][output.watchzoneList[i].i][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nW][output.watchzoneList[i].i][output.watchzoneList[i].j][nKInt]
      <<std::setw(nWidthOutputField);
    if(nKInt-1==-1){
      output.ofWatchZoneFiles[i]<<"-";//not defined when on edge
    }
    else{
      output.ofWatchZoneFiles[i]
        <<grid.dLocalGridOld[grid.nW][output.watchzoneList[i].i][output.watchzoneList[i].j][nKInt-1];
    }
    output.ofWatchZoneFiles[i]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nQ2][output.watchzoneList[i].i][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nR][nIInt][0][0]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nD][output.watchzoneList[i].i][output.watchzoneList[i].j]
      [output.watchzoneList[i].k];
    if(parameters.nTypeTurbulanceMod>0){
      output.ofWatchZoneFiles[i]
        <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nEddyVisc]
        [output.watchzoneList[i].i][output.watchzoneList[i].j][output.watchzoneList[i].k];
    }
    else{
      output.ofWatchZoneFiles[i]<<std::setw(nWidthOutputField)<<"-";
    }
    output.ofWatchZoneFiles[i]
      <<std::setw(nWidthOutputField)<<grid.dLocalGridOld[grid.nDenAve][output.watchzoneList[i].i][0][0]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nE][output.watchzoneList[i].i][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nP][output.watchzoneList[i].i][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)
      <<grid.dLocalGridOld[grid.nT][output.watchzoneList[i].i][output.watchzoneList[i].j][output.watchzoneList[i].k]
      <<std::setw(nWidthOutputField)<<dM
      <<std::setw(nWidthOutputField)<<dMCalc
      <<std::setw(nWidthOutputField)<<(dM-dMCalc)/dM
      <<std::endl;
  }
}
void finWatchZones(Output &output){
  for(unsigned int i=0;i<output.watchzoneList.size();i++){
    output.ofWatchZoneFiles[i].flush();//write out buffer
    output.ofWatchZoneFiles[i].close();
  }
}

