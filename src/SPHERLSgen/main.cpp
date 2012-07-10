/** @file
  
  This code is used to generate the starting model for SPHERLS
  
  \todo Want to make printing model to the screen an option in the configuration file, and the 
  default should be not print the model.
  
  \todo It would also make sense to print to a binary model rather than an ascii model. It could 
  however be an option with the default being to print to a binary file.
*/

#include "exception2.h"
#include "xmlParser.h"
#include "xmlFunctions.h"
#include "eos.h"
#include "main.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <vector>
#include <limits>
#include <string>
#include <algorithm>
#include <iomanip>

int main(){
  try{
    
    //read in settings of config file "config.xml"
    readConfig("SPHERLSgen.xml","data");
    
  }
  catch(exception2& eTemp){
    std::cout<<eTemp.getMsg();
  }
  return 0;
}
void readConfig(std::string sConfigFileName,std::string sStartNode){
  
  //open file
  XMLNode xData=openXMLFile(sConfigFileName,sStartNode);
  
  //GET CONSTANTS
  
  //get the gravitational constant
  if(!getXMLValueNoThrow(xData,"G",0,dG)){
    dG=6.67259e-08;//default
  }
  
  //get mass of the sun
  if(!getXMLValueNoThrow(xData,"M-sun",0,dMSun)){
    dMSun=1.9891e+33;//default
  }
  
  //get raduis of the sun
  if(!getXMLValueNoThrow(xData,"R-sun",0,dRSun)){
    dRSun=6.958e+10;//default
  }
  
  //get luminosity of the sun
  if(!getXMLValueNoThrow(xData,"L-sun",0,dLSun)){
    dLSun=3.839e33;
  }
  
  //get the Stefan-Boltzman constant
  if(!getXMLValueNoThrow(xData,"sigma",0,dSigma)){
    dSigma=5.6704e-5;
  }
  
  //get first model
  XMLNode xModel=getXMLNodeNoThrow(xData, "model",0);
  
  if(xModel.isEmpty()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": no \"model\" node found under \"data\" node. Need at least one model node.\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  int nCount=0;
  int nCount2=0;
  XMLNode xMDeltaDelta;
  XMLNode xMDeltaPicking;
  MDeltaDelta mDeltaDeltaTemp;
  std::string sModelType;
  std::string sEOSType;
  while(!xModel.isEmpty()){
    
    //get model type
    getXMLAttribute(xModel,"type",sModelType);
    if(sModelType.compare("sedov")==0){
      
      //no alpha for sedov model
      dAlpha=0.0;
      
      /////////////////////////////////////////////
      //OUTPUT OPTIONS
      
      XMLNode xOutput=getXMLNode(xModel,"output",0);
      
      //get name of file to write model to
      getXMLValue(xOutput,"fileName",0,sOutPutfile);
      std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": generating sedov model \""<<sOutPutfile<<"\" ...\n";
      
      //get output file type, binary or ascii
      if(!getXMLValueNoThrow(xOutput,"binary",0,bBinaryOutput)){
        bBinaryOutput=true;//default
      }
      
      //get wheather to write to screen
      if(!getXMLValueNoThrow(xOutput,"writeToScreen",0,bWriteToScreen)){
        bWriteToScreen=false;//default value
      }
      
      //get precision
      if(!getXMLValueNoThrow(xOutput,"precision",0,nPrecision)){
        nPrecision=16;
      }
      
      //get dTimeStepFactor
      getXMLValue(xOutput,"timeStepFactor",0,dTimeStepFactor);
      
      
      ////////////////////////////////////////
      //GET DIMENSIONS OF MODEL
      
      //switch to dimensions node
      XMLNode xDims=getXMLNode(xModel,"dimensions",0);
      
      //get theta dimension
      getXMLValue(xDims,"num-theta",0,nNumTheta);
      if(nNumTheta<1){
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
          <<": number of theta zones \"num-theta\" must be 1 or greater\n";
        throw exception2(ssTemp.str(),INPUT);
      }
      
      //get phi dimension
      getXMLValue(xDims,"num-phi",0,nNumPhi);
      if(nNumPhi<1){
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
          <<": number of phi zones \"num-phi\" must be 1 or greater\n";
        throw exception2(ssTemp.str(),INPUT);
      }
      
      //get delta theta
      getXMLValue(xDims,"delta-theta",0,dDeltaTheta);
      
      //get delta phi
      getXMLValue(xDims,"delta-phi",0,dDeltaPhi);
      
      //get number of ghost cells
      getXMLValue(xDims,"num-ghost-cells",0,nNumGhostCells);
      
      //set number of dimensions
      if(nNumTheta==1&&nNumPhi==1){
        nNumDims=1;
      }
      else if(nNumTheta>1&&nNumPhi==1){
        nNumDims=2;
      }
      else if(nNumTheta>1&&nNumPhi>1){
        nNumDims=3;
      }
      else if(nNumTheta==1&&nNumPhi>1){//only support 2D in theta direction
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
          <<": 2D simulations use only radial and theta directions, try switching number of "
          <<"theta and phi zones.\n";
        throw exception2(ssTemp.str(),INPUT);
      }
      
      
      //RADIAL INDEPENDENT VARIABLE
      
      //get deleta r
      XMLNode xIndepVar=getXMLNode(xDims,"radIndepVar",0);
      getXMLValue(xIndepVar,"r-delta",0,dRDelta);
      
      //get minimum radius
      getXMLValue(xIndepVar,"r-min",0,dRMin);
      
      //get number of raidial zones
      getXMLValue(xIndepVar,"num-r",0,nNumR);
      
      //get number of 1D zones at center
      getXMLValue(xIndepVar,"num-1D",0,nNumZones1D);
      
      
      ////////////////////////////////////////
      //GET PERIODICITY
      
      //switch to periodic node
      XMLNode xPeriodic=getXMLNode(xModel,"periodic",0);
      
      //get x-periodicity
      getXMLValue(xPeriodic,"x0",0,nPeriodic[0]);
      if(nPeriodic[0]!=0){//x0 periodicity not yet supported
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
          <<": periodicity in x-direction not yet implemented\n";
        throw exception2(ssTemp.str(),INPUT);
      }
      
      //get y-periodicity
      getXMLValue(xPeriodic,"x1",0,nPeriodic[1]);
      if(nPeriodic[1]!=0&&nNumTheta<nNumGhostCells){//must have enough zones for periodic BC
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
          <<": WARNING number of theta zones is "<<nNumTheta
          <<" which is too small to support periodicity with "<<nNumGhostCells
          <<" ghost cells at boundary. Unsetting periodic boundary condition.\n";
        nPeriodic[1]=0;
      }
      
      //get z-periodicity
      getXMLValue(xPeriodic,"x2",0,nPeriodic[2]);
      if(nPeriodic[2]!=0&&nNumPhi<nNumGhostCells){//must have enough zones for periodic BC
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
          <<": WARNING number of phi zones is "<<nNumPhi
          <<" which is too small to support periodicity with "<<nNumGhostCells
          <<" ghost cells at boundary. Unsetting periodic boundary condition.\n";
        nPeriodic[2]=0;
      }
      
      
      //GET STATE
      XMLNode xState=getXMLNode(xModel,"state",0);
      
      //get gamma
      getXMLValue(xState,"gamma",0,dGamma);
      bGammaLawEOS=true;
      
      //get number of central cells
      getXMLValue(xState,"num-shells-center",0,nNumCellsCent);
      
      //get central energy
      getXMLValue(xState,"eng-cent",0,dEngCent);
      
      //get energy in for rest of grid
      getXMLValue(xState,"eng",0,dEng);
      
      //get density of grid
      getXMLValue(xState,"rho",0,dRho);
      
      //generate sedov model
      generateModel_SEDOV();
      
      //write to screen
      if(bWriteToScreen){
        writeModelToScreen_GL();
      }
      
      //write to file
      if(bBinaryOutput){
        if(nNumDims==1){
          writeModel_Bin_R_GL();
        }
        else if(nNumDims==2){
          writeModel_Bin_RT_GL();
        }
        else if(nNumDims==3){
          writeModel_Bin_RTP_GL();
        }
      }
      else{
        if(nNumDims==1){
          writeModel_R_GL();
        }
        else if(nNumDims==2){
          writeModel_RT_GL();
        }
        else if(nNumDims==3){
          writeModel_RTP_GL();
        }
      }

    }
    else if(sModelType.compare("stellar")==0){
      
      
      /////////////////////////////////////////////
      //OUTPUT OPTIONS
      
      XMLNode xOutput=getXMLNode(xModel,"output",0);
      
      //get name of file to write model to
      getXMLValue(xOutput,"fileName",0,sOutPutfile);
      std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": generating stellar model \""
        <<sOutPutfile<<"\" ...\n";
      
      //get output file type, binary or ascii
      if(!getXMLValueNoThrow(xOutput,"binary",0,bBinaryOutput)){
        bBinaryOutput=true;//default
      }
      
      //get wheather to write to screen
      if(!getXMLValueNoThrow(xOutput,"writeToScreen",0,bWriteToScreen)){
        bWriteToScreen=false;//default value
      }
      
      //get precision
      if(!getXMLValueNoThrow(xOutput,"precision",0,nPrecision)){
        nPrecision=16;
      }
      
      //get dTimeStepFactor
      getXMLValue(xOutput,"timeStepFactor",0,dTimeStepFactor);
      
      
      ////////////////////////////////////////
      //GET EQUATION OF STATE
      XMLNode xEOS=getXMLNode(xModel,"EOS",0);
      getXMLAttribute(xEOS,"type",sEOSType);
      if(sEOSType.compare("gammaLaw")==0){//eosTable overides gamma
        
        //get gamma
        getXMLValue(xEOS,"gamma",0,dGamma);
        
        //get name of file containing internal energy profile
        std::string sEProFileName;
        getXMLValue(xEOS,"E-pro",0,sEProFileName);
        
        //read in energy profile
        readEnergyProfile_GL(sEProFileName);
        
        //get dRSurf
        getXMLValue(xEOS,"R-surf",0,dRSurf);
      }
      else if(sEOSType.compare("table")==0){//if no gamma and no eosTable specified 
        
        //get file name for eos table
        getXMLValue(xEOS,"eosTable",0,sEOSFile);
        
        //read in equation of state
        eosTable.readBin(sEOSFile);
        
        //get dTeff
        getXMLValue(xEOS,"T-eff",0,dTeff);
        
        //get dL
        getXMLValue(xEOS,"L",0,dL);
        
        //get allowed tolerance
        getXMLValue(xEOS,"tolerance",0,dTolerance);
        
        //don't use gamma-law EOS
        bGammaLawEOS=false;
      }
      else{
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
            <<": EOS node has an unknown type in model "<<nCount<<".\n";
          throw exception2(ssTemp.str(),INPUT);
      }
      
      
      ////////////////////////////////////////
      //GET DIMENSIONS OF MODEL
      
      //switch to dimensions node
      XMLNode xDims=getXMLNode(xModel,"dimensions",0);
      
      //get theta dimension
      getXMLValue(xDims,"num-theta",0,nNumTheta);
      if(nNumTheta<1){
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
          <<": number of theta zones \"num-theta\" must be 1 or greater in model"<<nCount<<"\n";
        throw exception2(ssTemp.str(),INPUT);
      }
      
      //get phi dimension
      getXMLValue(xDims,"num-phi",0,nNumPhi);
      if(nNumPhi<1){
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
          <<": number of phi zones \"num-phi\" must be 1 or greater in model "<<nCount<<"\n";
        throw exception2(ssTemp.str(),INPUT);
      }
      
      //get delta theta
      getXMLValue(xDims,"delta-theta",0,dDeltaTheta);
      
      //get delta phi
      getXMLValue(xDims,"delta-phi",0,dDeltaPhi);
      
      //get number of ghost cells
      getXMLValue(xDims,"num-ghost-cells",0,nNumGhostCells);
      
      //set number of dimensions
      if(nNumTheta==1&&nNumPhi==1){
        nNumDims=1;
      }
      else if(nNumTheta>1&&nNumPhi==1){
        nNumDims=2;
      }
      else if(nNumTheta>1&&nNumPhi>1){
        nNumDims=3;
      }
      else if(nNumTheta==1&&nNumPhi>1){//only support 2D in theta direction
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
          <<": 2D simulations use only radial and theta directions, try switching number of "
          <<"theta and phi zones in model "<<nCount<<".\n";
        throw exception2(ssTemp.str(),INPUT);
      }
      
      
      //RADIAL INDEPENDENT VARIABLE
      
      //get total mass
      XMLNode xIndepVar=getXMLNode(xDims,"radIndepVar",0);
      getXMLValue(xIndepVar,"M-total",0,dMTotal);
      
      //get delta mass init
      getXMLValue(xIndepVar,"M-delta-init",0,dMDelta);
      
      //get method for determining % change in delta mass
      xMDeltaPicking=getXMLNodeNoThrow(xIndepVar,"M-delta-picking",0);
      if(!xMDeltaPicking.isEmpty()){
        
        std::string sDeltaMPicking;
        getXMLAttribute(xMDeltaPicking,"type",sDeltaMPicking);
        if(sDeltaMPicking.compare("auto")==0){
          bAutoDeltaM=true;
        }
        else if (sDeltaMPicking.compare("manual")==0){
          bAutoDeltaM=false;
        }
        else{
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
            <<": \"M-delta-picking\" node attribute \"type\" must be either  \"auto\" or \"manual\""
            <<" node under the "<<nCount<<"th model.\n";
          throw exception2(ssTemp.str(),INPUT);
        }
      }
      if(!bAutoDeltaM){
        
        //get % change in delta mass
        xMDeltaDelta=getXMLNodeNoThrow(xMDeltaPicking, "M-delta-delta",0);
        if(xMDeltaDelta.isEmpty()){
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
            <<": no \"M-delta-delta\" node found under \"radIndepVar\" node, under \"dimensions\""
            <<" node under the "<<nCount<<"th model.\n";
          throw exception2(ssTemp.str(),INPUT);
        }
        nCount2=0;
        while(!xMDeltaDelta.isEmpty()){
          
          //get type of stop
          getXMLAttribute(xMDeltaDelta,"stopType",mDeltaDeltaTemp.sStopType);
          
          //if type is T and model is adiabatic, halt
          if(mDeltaDeltaTemp.sStopType.compare("T")==0&&sEOSType.compare("gammaLaw")==0){
            std::stringstream ssTemp;
            ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
              <<":"<<nCount2<<"th \"M-delta-delta\" node of "<<nCount<<"th model has a \"stopType\""
              <<" of \"T\" but model uses a gamma-law gas, must have a tabulated equation of state to"
              <<" use a \"stopType\" of \"T\". Perhpas use a \"stoptype\" of \"R\" or use a different"
              <<" equation of state.\n";
            throw exception2(ssTemp.str(),INPUT);
          }
          
          //get value of stop
          getXMLAttribute(xMDeltaDelta,"stopValue",mDeltaDeltaTemp.dStopValue);
          
          /**\todo need to check that T is increasing, and R is decreasing. This will get tricky if
          R and T types are mixed.*/
          
          //get value of MDeltaDelta
          getXMLValue(xMDeltaPicking,"M-delta-delta",nCount2,mDeltaDeltaTemp.dMDeltaDelta);
          
          //add to vector
          vecMDeltaDeltaList.push_back(mDeltaDeltaTemp);
          
          //get next node
          nCount2++;
          xMDeltaDelta=getXMLNodeNoThrow(xMDeltaPicking, "M-delta-delta",nCount2);
        }
      }
      //get dAlpha
      getXMLValue(xIndepVar,"alpha",0,dAlpha);
      
      //get number of 1D zones at center
      getXMLValue(xIndepVar,"num-1D",0,nNumZones1D);
      
      ////////////////////////////////////////
      //GET PERIODICITY
      
      //set defaults
      nPeriodic[0]=0;
      nPeriodic[1]=1;
      nPeriodic[2]=1;
      
      //switch to periodic node if there is one
      XMLNode xPeriodic=getXMLNodeNoThrow(xModel,"periodic",0);
      
      //if there is a periodic node get values
      if(!xPeriodic.isEmpty()){
        getXMLValue(xPeriodic,"x0",0,nPeriodic[0]);
        getXMLValue(xPeriodic,"x1",0,nPeriodic[1]);
        getXMLValue(xPeriodic,"x2",0,nPeriodic[2]);
      }
      
      //check periodicity against grid sizes
      if(nPeriodic[0]!=0){//x0 periodicity not yet supported
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
          <<": periodicity in x-direction not yet implemented\n";
        throw exception2(ssTemp.str(),INPUT);
      }
      if(nPeriodic[1]!=0&&nNumTheta<nNumGhostCells){//must have enough zones for periodic BC
        if(!xPeriodic.isEmpty()){//don't need to tell user about this unless they tried to set it
          std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
            <<": WARNING number of theta zones is "<<nNumTheta
            <<" which is too small to support periodicity with "<<nNumGhostCells
            <<" ghost cells at boundary. Unsetting periodic boundary condition.\n";
        }
        nPeriodic[1]=0;
      }
      if(nPeriodic[2]!=0&&nNumPhi<nNumGhostCells){//must have enough zones for periodic BC
        if(!xPeriodic.isEmpty()){//don't need to tell user about this unless they tried to set it
          std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
            <<": WARNING number of phi zones is "<<nNumPhi
            <<" which is too small to support periodicity with "<<nNumGhostCells
            <<" ghost cells at boundary. Unsetting periodic boundary condition.\n";
        }
        nPeriodic[2]=0;
      }
      
      //GET VELOCITY DISTRIBUTION
      XMLNode xVelDist=getXMLNode(xModel,"velocityDist",0);
      std::string sType;
      getXMLAttribute(xVelDist,"type",sUDistType);
      if(sUDistType=="POLY"){
        int nIndex=0;
        XMLNode xTemp=getXMLNodeNoThrow(xVelDist,"term",nIndex);
        while(!xTemp.isEmpty()){
          
          //READ IN TERM
          term tTemp;
          
          //get coeffecient
          getXMLValue(xTemp,"c",0,tTemp.dCoeff);
          
          //get power
          getXMLValueNoThrow(xTemp,"p",0,tTemp.dPower);
          
          //add term
          vectVelDist.push_back(tTemp);
          
          //get next term
          nIndex++;
          xTemp=getXMLNodeNoThrow(xVelDist,"term",nIndex);
        }
        if(vectVelDist.size()==0){//no terms found
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
            <<": need at least one \"term\" node in the \"velocityDist\" node!\n";
          throw exception2(ssTemp.str(),INPUT);
        }
      }
      else if(sUDistType=="PRO"){
        
        //get file name 
        std::string sProfileFileName;
        getXMLValue(xVelDist,"fileName",0,sProfileFileName);
        
        //read u profile
        readUProfile(sProfileFileName);
        
        //get surface velocity value
        getXMLValue(xVelDist,"uSurf",0,dUSurf);
        
      }
      else{
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
            <<": velocityDist node has an unknown type.\n";
          throw exception2(ssTemp.str(),INPUT);
      }
      
      //generate stellar model
      if(bGammaLawEOS){
        
        //generate model
        generateModel_GL();
        
        //check for velocity perturbations and apply them
        XMLNode xPerturb=getXMLNodeNoThrow(xVelDist,"perturb",0);
        int nPerturbation=0;
        while(!xPerturb.isEmpty()){
          
          
          //figure out type of preturbaiton
          std::string sPerturbType;
          getXMLAttribute(xPerturb,"type",sPerturbType);
          if(sPerturbType=="torus"){
            
            //ready perturbation info, and apply a torus velocity perturbation
            pretubeVelocityTorus(xPerturb);
          }
          else{
            std::stringstream ssTemp;
            ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
              <<": unknown perturbation type, \""<<sPerturbType<<"\", in "<<nPerturbation
              <<"th preturbation in velocityDist node with type \""<<sUDistType<<"\"\n";
            throw exception2(ssTemp.str(),INPUT);
          }
          
          //get next perturbation
          nPerturbation++;
          xPerturb=getXMLNodeNoThrow(xVelDist,"pretrub",nPerturbation);
        }
        
        //write to screen
        if(bWriteToScreen){
          writeModelToScreen_GL();
        }
        
        //write to file
        if(bBinaryOutput){
          if(nNumDims==1){
            writeModel_Bin_R_GL();
          }
          else if(nNumDims==2){
            writeModel_Bin_RT_GL();
          }
          else if(nNumDims==3){
            writeModel_Bin_RTP_GL();
          }
        }
        else{
          if(nNumDims==1){
            writeModel_R_GL();
          }
          else if(nNumDims==2){
            writeModel_RT_GL();
          }
          else if(nNumDims==3){
            writeModel_RTP_GL();
          }
        }
      }
      else{
        
        //generate model
        generateModel_TEOS();
        
        //check for velocity perturbations and apply them
        XMLNode xPerturb=getXMLNodeNoThrow(xVelDist,"perturb",0);
        int nPerturbation=0;
        while(!xPerturb.isEmpty()){
          
          
          //figure out type of preturbaiton
          std::string sPerturbType;
          getXMLAttribute(xPerturb,"type",sPerturbType);
          if(sPerturbType=="torus"){
            
            //ready perturbation info, and apply a torus velocity perturbation
            pretubeVelocityTorus(xPerturb);
          }
          else{
            std::stringstream ssTemp;
            ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
              <<": unknown perturbation type, \""<<sPerturbType<<"\", in "<<nPerturbation
              <<"th preturbation in velocityDist node with type \""<<sUDistType<<"\"\n";
            throw exception2(ssTemp.str(),INPUT);
          }
          
          //get next perturbation
          nPerturbation++;
          xPerturb=getXMLNodeNoThrow(xVelDist,"pretrub",nPerturbation);
        }
        
        //write to screen
        if(bWriteToScreen){
          writeModelToScreen_TEOS();
        }
        
        //write to file
        if(bBinaryOutput){
          if(nNumDims==1){
            writeModel_Bin_R_TEOS();
          }
          else if(nNumDims==2){
            writeModel_Bin_RT_TEOS();
          }
          else if(nNumDims==3){
            writeModel_Bin_RTP_TEOS();
          }
        }
        else{
          if(nNumDims==1){
            writeModel_R_TEOS();
          }
          else if(nNumDims==2){
            writeModel_RT_TEOS();
          }
          else if(nNumDims==3){
            writeModel_RTP_TEOS();
          }
        }
      }
      
    }
    else{
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
        <<": unknown model type \""<<sModelType<<"\"\n";
      throw exception2(ssTemp.str(),INPUT);
    }
    //get next model
    nCount++;
    xModel=getXMLNodeNoThrow(xData,"model",nCount);
    
    //erase vectors
    vecdM.erase(vecdM.begin(),vecdM.end());
    vecdMDel.erase(vecdMDel.begin(),vecdMDel.end());
    vecdP.erase(vecdP.begin(),vecdP.end());
    vecdE.erase(vecdE.begin(),vecdE.end());
    vecdRho.erase(vecdRho.begin(),vecdRho.end());
    vecdR.erase(vecdR.begin(),vecdR.end());
    vecdT.erase(vecdT.begin(),vecdT.end());
    vecdKappa.erase(vecdKappa.begin(),vecdKappa.end());
    vectVelDist.erase(vectVelDist.begin(),vectVelDist.end());
  }
}
void readUProfile(std::string sProfileFileName){
  
  //attempt to open the file
  std::ifstream ifIn;
  ifIn.open(sProfileFileName.c_str());
  if(!ifIn.good()){//file not ready for reading
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sProfileFileName<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //get number of points
  ifIn>>nNumUProPoints;
  dUPro =new double[nNumUProPoints];
  dUProR=new double[nNumUProPoints];
  for(unsigned int i=0;i<nNumUProPoints;i++){
    ifIn>>dUProR[i];
    ifIn>>dUPro[i];
  }
  ifIn.close();
}
void generateModel_SEDOV(){
  
  //calculate first shell
  calculateFirstShell_SEDOV();
  
  //calculate the rest of the shells
  for(int i=nNumR-2;i>=0;i--){
    calculateShell_SEDOV(i);
  }
  
  makeVelocityDist_SEDOV();
}
void generateModel_TEOS(){

  //calculate first shell
  calculateFirstShell_TEOS();
  
  //calculate the rest of the shells
  int nShell=1;
  bool bContinue=true;
  while(bContinue){//while not deep enough
    calculateShell_TEOS(nShell);
    if(vecMDeltaDeltaList[vecMDeltaDeltaList.size()-1].sStopType.compare("R")==0){
      if(vecdR[nShell]<=vecMDeltaDeltaList[vecMDeltaDeltaList.size()-1].dStopValue*dRSun){//radius is deep enough
        bContinue=false;
      }
    }
    else if(vecMDeltaDeltaList[vecMDeltaDeltaList.size()-1].sStopType.compare("T")==0){
      if(vecdT[nShell]>=vecMDeltaDeltaList[vecMDeltaDeltaList.size()-1].dStopValue){//temperature is deep enough
        bContinue=false;
      }
    }
    nShell++;
  }
  
  makeVelocityDist();
  
}
void generateModel_GL(){
  
  //calculate first shell
  calculateFirstShell_GL();
  
  //calculate the rest of the shells
  int nShell=1;
  bool bContinue=true;
  while(bContinue){//while not deep enough
    calculateShell_GL(nShell);
    if(vecMDeltaDeltaList[vecMDeltaDeltaList.size()-1].sStopType.compare("R")==0){
      if(vecdR[nShell]<=vecMDeltaDeltaList[vecMDeltaDeltaList.size()-1].dStopValue*dRSun){//radius is deep enough
        bContinue=false;
      }
    }
    else if(vecMDeltaDeltaList[vecMDeltaDeltaList.size()-1].sStopType.compare("T")==0){
      if(vecdT[nShell]>=vecMDeltaDeltaList[vecMDeltaDeltaList.size()-1].dStopValue){//temperature is deep enough
        bContinue=false;
      }
    }
    nShell++;
  }
  
  makeVelocityDist();
}
void calculateFirstShell_SEDOV(){
  
  //mass (-1/2)
  vecdM.push_back(0.0);
  
  //get Mass at (1/2)
  vecdM.push_back(0.0);
  
  //Radius (-1/2)
  vecdR.push_back(dRMin+dRDelta*double(nNumR));
  
  //Radius (+1/2)
  vecdR.push_back(dRMin+dRDelta*double(nNumR-1));
  
  //delta M (0)
  double dV=1.333333333333333*dPi*(pow(vecdR[0],3)-pow(vecdR[1],3));
  vecdMDel.push_back(-1.0*dV*dRho);
  
  //density
  vecdRho.push_back(dRho);
  
  //set energy at (0)
  vecdE.push_back(dEng);
  
  //set pressure
  vecdP.push_back((dGamma-1.0)*dRho*vecdE[0]);
  
}
void calculateFirstShell_TEOS(){
  
  //mass (-1/2)
  vecdM.push_back(dMTotal*dMSun);
  
  //temeperature (0)
  double dTSurf=pow(2.0,-0.25)*dTeff;
  vecdT.push_back(dTSurf);
  
  //radius (-1/2), slight inconsistancey here with where the temperature and radius are defined
  double dRSurf=sqrt(dL*dLSun/(4.0*dPi*dSigma*pow(dTeff,4.0)));
  vecdR.push_back(dRSurf);
  
  //delta M (0)
  vecdMDel.push_back(-1.0*dMDelta*dMSun);
  
  //pressure at (0)
  vecdP.push_back(dG*vecdM[0]/(-8.0*dPi*pow(vecdR[0],4.0))*(0.5+dAlpha)*vecdMDel[0]);
  
  //calculate density at (0)
  //density will be low, start between first two grid points of table
  double dRho=pow(10.0,(eosTable.dLogRhoMin+0.5*eosTable.dLogRhoDelta));
  double dError=std::numeric_limits<double>::max();
  double dP;
  double dDRhoDP;
  
  int nIteration=0;
  while(dError>dTolerance){
    
    //get pressure, to try to match to dP_N, and derivative of density w.r.t. pressure which used to
    //calculate correction to density
    eosTable.getPAndDRhoDP(vecdT[0],dRho,dP,dDRhoDP);
    
    //corrected density
    dRho=(vecdP[0]-dP)*dDRhoDP+dRho;
    
    //how far off was the pressure
    dError=fabs((vecdP[0]-dP)/vecdP[0]);
    nIteration++;
  }
  vecdRho.push_back(dRho);
  
  double dE;
  double dKappa;
  eosTable.getEKappa(vecdT[0],vecdRho[0],dE,dKappa);
  
  //get energy at (0)
  vecdE.push_back(dE);
  
  //get opcacity at (0)
  vecdKappa.push_back(dKappa);
  
  //get Mass at (1/2)
  vecdM.push_back(vecdM[0]+vecdMDel[0]);
  
  //get radius at (1/2)
  vecdR.push_back(pow(3.0/(4.0*dPi)*vecdMDel[0]/vecdRho[0]+pow(vecdR[0],3.0)
    ,0.33333333333333333333333333333333));
}
void calculateFirstShell_GL(){
  //calculated differently from other cells because it needs to have the inner 
  //interfaces calculated awell, and in adition initializes the vectors
  //using the data provided in the config file
  
  //inner interface quantities
  vecdM.   push_back(dMTotal*dMSun);
  vecdR.   push_back(dRSurf*dRSun);
  
  //centered quantities
  vecdMDel.push_back(-1.0*dMDelta*dMSun);
  double dM_i=vecdM[0]+vecdMDel[0]*0.5;
  //double dM_im1=vecdM[0];
  double dPTemp=-0.125*dG*vecdM[0]/(dPi*pow(vecdR[0],4.0))*vecdMDel[0]*(0.5+dAlpha);
  
  vecdP.   push_back(dPTemp);
  double dIntVar=log10(1.0-dM_i/dMTotal/dMSun);
  vecdE.   push_back(interpolateE_GL(dIntVar));//last point in profile is at the surface
  vecdRho. push_back(EOS_GL(vecdP[0],vecdE[0]));//use gamma law
  
  //outer interface quantities
  vecdM.   push_back(vecdM[0]+vecdMDel[0]);
  vecdR.   push_back(pow(3.0*vecdMDel[0]*0.25/vecdRho[0]/dPi+pow(vecdR[0],3.0)
    ,0.33333333333333333333333333333333));//first outer radius
}
void calculateShell_SEDOV(int nShell){
  
  //mass (-1/2)
  vecdM.push_back(0.0);
  
  //Radius (+1/2)
  vecdR.push_back(dRMin+double(nShell)*dRDelta);
  
  //delta M (0)
  double dV=1.333333333333333*dPi*(pow(vecdR[nNumR-1-nShell],3)-pow(vecdR[nNumR-nShell],3));
  vecdMDel.push_back(-1.0*dV*dRho);
  
  //density
  vecdRho.push_back(dRho);
  
  //set energy at (0)
  if(nShell<nNumCellsCent){//if in central region
    double dVCent=1.3333333333333*dPi*pow(double(nNumCellsCent)*dRDelta+dRMin,3);
    double dMCent=dVCent*dRho;
    vecdE.push_back(dEngCent/dMCent);
  }
  else{//if outside central region
    vecdE.push_back(dEng);
  }
  
  //get Mass at (1/2)
  vecdM.push_back(0.0);
  
  //set pressure
  vecdP.push_back((dGamma-1.0)*dRho*vecdE[nNumR-1-nShell]);
}
void calculateShell_TEOS(unsigned int nShell){
  
  //calculate delta M for current shell (nShell)
  dMDeltaDelta=vecMDeltaDeltaList[0].dMDeltaDelta;
  for(unsigned int i=1;i<vecMDeltaDeltaList.size();i++){
    
    //test temperature, test radius
    if(vecMDeltaDeltaList[i].sStopType.compare("T")==0){// i is T
      if(vecMDeltaDeltaList[i-1].sStopType.compare("T")==0){//both i and i-1 are T
        if(vecMDeltaDeltaList[i].dStopValue>=vecdT[nShell-1]
          && vecMDeltaDeltaList[i-1].dStopValue<vecdT[nShell-1]){//was last zone in range
          dMDeltaDelta=vecMDeltaDeltaList[i].dMDeltaDelta;
        }
      }
      else if(vecMDeltaDeltaList[i-1].sStopType.compare("R")==0){//i is T and i-1 is R
        if(vecMDeltaDeltaList[i].dStopValue>=vecdT[nShell-1]
          && vecMDeltaDeltaList[i-1].dStopValue*dRSun<vecdR[nShell]){//was last zone in range
          dMDeltaDelta=vecMDeltaDeltaList[i].dMDeltaDelta;
        }
      }
    }
    else if(vecMDeltaDeltaList[i].sStopType.compare("R")==0){//i is R
      if(vecMDeltaDeltaList[i-1].sStopType.compare("T")==0){//i is R and i-1 are T
        if(vecMDeltaDeltaList[i].dStopValue*dRSun>=vecdR[nShell]
          && vecMDeltaDeltaList[i-1].dStopValue<vecdT[nShell-1]){//was last zone in range
          dMDeltaDelta=vecMDeltaDeltaList[i].dMDeltaDelta;
        }
      }
      else if(vecMDeltaDeltaList[i-1].sStopType.compare("R")==0){//both i and i-1 are R
        if(vecMDeltaDeltaList[i].dStopValue*dRSun>=vecdR[nShell]
          && vecMDeltaDeltaList[i-1].dStopValue*dRSun<vecdR[nShell]){//was last zone in range
          dMDeltaDelta=vecMDeltaDeltaList[i].dMDeltaDelta;
        }
      }
    }
  }
  vecdMDel.push_back(vecdMDel[nShell-1]*(1.0+dMDeltaDelta));
  if(fabs(vecdMDel[nShell-1]/vecdM[nShell-1])<=std::numeric_limits<double>::epsilon()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": delta M over M ("<<vecdMDel[nShell-1]/vecdM[nShell-1]
      <<") for shell "<<nShell
      <<" is smaller than is representable by double precision on current machine. Check "
      <<"\"M-delta-delta\" nodes.\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //initial guesses at T and rho
  double dT=vecdT[nShell-1];
  double dRho=vecdRho[nShell-1];
  double dRhoCorrection=dRho;
  double dTCorrection=dT;
  double dRhoError=dRhoCorrection/dRho;
  double dTError  =dTCorrection  /dT;
  double dCorrectionFrac=0.1;
  
  //keep going if error temperature or error in density is too big
  while(fabs(dRhoError)>dTolerance||fabs(dTError)>dTolerance){
    
    //calculate coeffecients of the equations to solve for T and rho corrections
    double dError=10.0;
    double dB1=-1.0*dMomentumCons(dT,dRho,nShell);
    double dB2=-1.0*dEnergyCons(dT,dRho,nShell);
    double dH=1.0e-10*dT;
    double dA11=dPartialDerivativeVar1(dT,dRho,nShell,dMomentumCons,dH,dError);
    dH=1.0e-10*dRho;
    double dA12=dPartialDerivativeVar2(dT,dRho,nShell,dMomentumCons,dH,dError);
    dH=1.0e-10*dT;
    double dA21=dPartialDerivativeVar1(dT,dRho,nShell,dEnergyCons,dH,dError);
    dH=1.0e-10*dRho;
    double dA22=dPartialDerivativeVar2(dT,dRho,nShell,dEnergyCons,dH,dError);
    
    //solve for the corrections
    double dFracA21A11=dA21/dA11;
    dRhoCorrection=(dB2-dFracA21A11*dB1)/(dA22-dFracA21A11*dA12);
    dTCorrection  =dB1/dA11-dA12/dA11*dRhoCorrection;
    
    while(dRho+dCorrectionFrac*dRhoCorrection<0.0){//keep it from decreasing too fast and becoming negative
      dCorrectionFrac=dCorrectionFrac*0.5;
    }
    while(dT  +dCorrectionFrac*dTCorrection<0.0){//keep it from decreasing too fast and becoming negative
      dCorrectionFrac=dCorrectionFrac*0.5;
    }
    
    //apply corrections
    dT  =dT  +dCorrectionFrac*dTCorrection;
    dRho=dRho+dCorrectionFrac*dRhoCorrection;
    
    //calculate relative error
    dTError  =dTCorrection  /dT;
    dRhoError=dRhoCorrection/dRho;
  }
  
  //temperature at (nShell)
  vecdT.push_back(dT);
  
  //density at (nShell)
  vecdRho.push_back(dRho);
  
  //get P, E, and Kappa
  double dP;
  double dE;
  double dKappa;
  eosTable.getPEKappa(dT,dRho,dP,dE,dKappa);
  
  //pressure at (nShell)
  vecdP.push_back(dP);
  
  //energy at (nShell)
  vecdE.push_back(dE);
  
  //opacity at (nShell)
  vecdKappa.push_back(dKappa);
  
  //mass at (nShell+1/2)
  vecdM.push_back(vecdM[nShell]+vecdMDel[nShell]);
  if(vecdM[nShell+1]<0.0){//gone past center of star
    //print out what we have so far
    vecdR.push_back(std::numeric_limits<float>::quiet_NaN());//set radius to a nan
    if(nNumDims==1){
      writeModel_R_TEOS();
    }
    else if(nNumDims==2){
      writeModel_RT_TEOS();
    }
    else if(nNumDims==3){
      writeModel_RTP_TEOS();
    }
    
    //throw exception
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": vecdM["<<nShell+1<<"]="<<vecdM[nShell+1]<<", we have gone past the center of the star!\n"
      <<"\tTry adjusting \"R-stop\", \"M-delta-delta\", and or \"M-delta-init\". See output model "
      <<"for details.\n";
    throw exception2(ssTemp.str(),CALCULATION);
  }
  
  //radius at (nShell+1/2)
  vecdR.push_back(pow(3.0/(4.0*dPi)*vecdMDel[nShell]/vecdRho[nShell]+pow(vecdR[nShell],3.0)
    ,0.33333333333333333333333333333333));
}
void calculateShell_GL(unsigned int nShell){
  
  //calculate mass, note: mass is one zone a head of the rest
  vecdMDel.push_back(vecdMDel[nShell-1]*(1.0+dMDeltaDelta));
  vecdM.   push_back(vecdM[nShell]+vecdMDel[nShell]);
  if(vecdM[nShell+1]<0.0){//gone past center of star
    //print out what we have so far
    vecdR.push_back(std::numeric_limits<float>::quiet_NaN());//set radius to a nan
    if(nNumDims==1){
      writeModel_R_GL();
    }
    else if(nNumDims==2){
      writeModel_RT_GL();
    }
    else if(nNumDims==3){
      writeModel_RTP_GL();
    }
    
    //throw exception
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": vecdM["<<nShell+1<<"]="<<vecdM[nShell+1]<<", we have gone past the center of the star!\n"
      <<"\tTry adjusting \"R-stop\", \"M-delta-delta\", and or \"M-delta-init\". See output model "
      <<"for details.\n";
    throw exception2(ssTemp.str(),CALCULATION);
  }
  
  //calculate centered quantities
  double dM_i    =(vecdM[nShell+1]+vecdM[nShell])/2.0;
  //double dM_im1  =(vecdM[nShell]+vecdM[nShell-1])/2.0;
  double dR_4=pow(vecdR[nShell],4.0);
  vecdP.   push_back(vecdP[nShell-1]-dG*vecdM[nShell]*0.25/dPi/dR_4*(vecdMDel[nShell]
    +vecdMDel[nShell-1])*0.5);
  double dIntVar=log10(1.0-dM_i/dMTotal/dMSun);
  vecdE.   push_back(interpolateE_GL(dIntVar));//last point in profile is at the surface
  vecdRho. push_back(EOS_GL(vecdP[nShell],vecdE[nShell]));//use gamma law
  
  //outer interface quantities
  vecdR.   push_back(pow(3.0*vecdMDel[nShell]*0.25/vecdRho[nShell]/dPi+pow(vecdR[nShell],3)
    ,0.33333333333333333333333333333333));//first outer radius
  
  //double dTerm1=(vecdP[nShell]-vecdP[nShell-1])/(dM_i-dM_im1);
  //double dTerm2=-1.0*dG*vecdM[nShell]/4.0/dPi/pow(vecdR[nShell],4);
  //std::cout<<dTerm1<<" "<<dTerm2<<std::endl;
}
double dTimeStep_TEOS(){
  
  double dTimeStep=1.0e200;
  
  //set up theta ranges
  double dStartTheta=(90.0-dDeltaTheta*double((nNumTheta+2*nNumGhostCells)/2))*dPi/180.0;  
  if(2*(nNumTheta/2)!=nNumTheta){
    dStartTheta-=dDeltaTheta*dPi/360.0;
  }
  dStartTheta+=dDeltaTheta*dPi/360.0;
  unsigned int nNumTheta2=nNumTheta+2*nNumGhostCells+1;

  for(unsigned int i=0;i<vecdRho.size()-nNumGhostCells;i++){
    
    //calculate sound speed
    double dC=eosTable.dSoundSpeed(vecdT[i],vecdRho[i]);
    
    //calculate time step in r-direction
    double dTimeStepR=(vecdR[i]-vecdR[i+1])/dC;
    
    //constant delta theta
    double dTimeStepTheta=dDeltaTheta*vecdR[i+1]/dC;
    
    //constant delta phi
    double dTimeStepPhi=1.0e200;
    for(unsigned int j=0;j<nNumTheta2;j++){
      double dTheta=dStartTheta+(double(j)*dDeltaTheta)*dPi/180.0;
      double dTemp=vecdR[j+1]*sin(dTheta)*dDeltaPhi/dC;
      if(dTemp<dTimeStepPhi){
        dTimeStepPhi=dTemp;
      }
    }
    
    if(dTimeStepR<dTimeStep){
      dTimeStep=dTimeStepR;
    }
    if(dTimeStepTheta<dTimeStep){
      dTimeStep=dTimeStepTheta;
    }
    if(dTimeStepPhi<dTimeStep){
      dTimeStep=dTimeStepPhi;
    }
  }
  return dTimeStep;
}
double dTimeStep_GL(){
  
  double dTimeStep=1.0e200;
  
  //set up theta ranges
  double dStartTheta=(90.0-dDeltaTheta*double((nNumTheta+2*nNumGhostCells)/2))*dPi/180.0;  
  if(2*(nNumTheta/2)!=nNumTheta){
    dStartTheta-=dDeltaTheta*dPi/360.0;
  }
  dStartTheta+=dDeltaTheta*dPi/360.0;
  unsigned int nNumTheta2=nNumTheta+2*nNumGhostCells+1;

  for(unsigned int i=0;i<vecdRho.size()-nNumGhostCells;i++){
    
    //calculate sound speed
    double dC=sqrt(dGamma*vecdP[i]/vecdRho[i]);
    
    //calculate time step in r-direction
    double dTimeStepR=(vecdR[i]-vecdR[i+1])/dC;
    
    //constant delta theta
    double dTimeStepTheta=dDeltaTheta*vecdR[i+1]/dC;
    
    //constant delta phi
    double dTimeStepPhi=1.0e200;
    for(unsigned int j=0;j<nNumTheta2;j++){
      double dTheta=dStartTheta+(double(j)*dDeltaTheta)*dPi/180.0;
      double dTemp=vecdR[j+1]*sin(dTheta)*dDeltaPhi/dC;
      if(dTemp<dTimeStepPhi){
        dTimeStepPhi=dTemp;
      }
    }
    
    if(dTimeStepR<dTimeStep){
      dTimeStep=dTimeStepR;
    }
    if(dTimeStepTheta<dTimeStep){
      dTimeStep=dTimeStepTheta;
    }
    if(dTimeStepPhi<dTimeStep){
      dTimeStep=dTimeStepPhi;
    }
  }
  return dTimeStep;
}
void writeModel_R_TEOS(){
  
  std::ofstream ofOut;
  ofOut.open(sOutPutfile.c_str());
  if(!ofOut.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sOutPutfile<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //file type and version
  ofOut<<"a 1"<<std::endl;
  
  //set double output precision
  ofOut.precision(nPrecision);
  ofOut.unsetf(std::ios::fixed);
  ofOut.setf(std::ios::scientific);
  
  //write out start time
  ofOut<<0.0<<" "<<0<<std::endl;
  
  //write out time step
  ofOut<<dTimeStep_TEOS()*dTimeStepFactor<<std::endl;
  
  //write out time step again for new timestep
  ofOut<<dTimeStep_TEOS()*dTimeStepFactor<<std::endl;
  
  //write out alpha
  ofOut<<dAlpha<<std::endl;
  
  //write out the equation of state file name and it's length
  ofOut<<sEOSFile.size()<<" "<<sEOSFile<<std::endl;
  
  //write place holder for artificial viscosity
  double dTemp=0.0;
  ofOut<<dTemp<<std::endl;
  
  //write place holder for artificial viscosity threshold
  ofOut<<dTemp<<std::endl;
  
  //output dimensions
  ofOut<<vecdP.size()-2*nNumGhostCells<<" "<<nNumTheta<<" "<<nNumPhi<<std::endl;
  
  //write out periodicity
  ofOut<<nPeriodic[0]<<" "<<nPeriodic[1]<<" "<<nPeriodic[2]<<std::endl;
  
  //write out number of 1D zones
  ofOut<<nNumZones1D<<std::endl;
  
  //write out number of ghost cells
  ofOut<<nNumGhostCells<<std::endl;
  
  //write out number of variables
  ofOut<<7<<std::endl;
  
  //write out variable info
  //output interior mass info (M,1)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<0<<"  ";//radial centered, theta and phi not defined, time independent
  
  //output delta M info (DM,2)
  ofOut<<0<<" "<<-1<<" "<<-1<<" "<<0<<"  ";//radial centered, theta and phi not defined, time independent
  
  //output radius info (R,3)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi not defined, time dependent
  
  //output density info (D,4)
  ofOut<<0<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//r,theta, and phi centered, time dependent
  
  //output radial velocity (U,5)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi centered, time dependent
  
  //output radial grid velocity (U0,6)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi centered, time dependent
  
  //output temperature (T,7)
  ofOut<<0<<" "<<-1<<" "<<-1<<" "<<1<<"  "<<std::endl;//r, theta, and phi centered, time dependent
  
  //write out variables
  ofOut.precision(nPrecision);
  ofOut.unsetf(std::ios::fixed);
  ofOut.setf(std::ios::scientific);
  
  //write out M_r - 1D
  int nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdM[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out Del M
  for(int i=vecdMDel.size()-1;i>=0;i--){//start at the center and work outward
    ofOut<<-1.0*vecdMDel[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out r
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdR[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out rho
  //write out 1D region
  for(int i=vecdRho.size()-1;i>=0;i--){//start at the center and work outward
    ofOut<<vecdRho[i]<<" "<<std::endl<<std::endl;
  }
  
  //write out u
  //write out 1D region
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart-=1;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<dU[i][0][0]<<" "<<std::endl<<std::endl;
  }

  //write out u_0
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<dU0[i]<<" "<<std::endl<<std::endl;//same as radial velocity
  }
  ofOut<<std::endl;
  
  //write out T - 3D
  nStart=vecdT.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdT.size()-2;
  }
  for(int i=vecdT.size()-1;i>=0;i--){//start at the center and work outward
    ofOut<<vecdT[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  ofOut.close();
}
void writeModel_R_GL(){
  
  std::ofstream ofOut;
  ofOut.open(sOutPutfile.c_str());
  if(!ofOut.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sOutPutfile<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //file type and version
  ofOut<<"a 1"<<std::endl;
  
  //set double output precision
  ofOut.precision(nPrecision);
  ofOut.unsetf(std::ios::fixed);
  ofOut.setf(std::ios::scientific);
  
  //write out start time
  ofOut<<0.0<<" "<<0<<std::endl;
  
  //write out time step for old step
  ofOut<<dTimeStep_GL()*dTimeStepFactor<<std::endl;
  
  //write out time step for next step
  ofOut<<dTimeStep_GL()*dTimeStepFactor<<std::endl;
  
  //write out alpha
  ofOut<<dAlpha<<std::endl;
  
  ofOut<<0<<" "<<dGamma<<std::endl;
  
  //write place holder for artificial viscosity
  double dTemp=0.0;
  ofOut<<dTemp<<std::endl;
  
  //write place holder for artificial viscosity threshold
  ofOut<<dTemp<<std::endl;
  
  //output dimensions
  ofOut<<vecdP.size()-2*nNumGhostCells<<" "<<nNumTheta<<" "<<nNumPhi<<std::endl;
  
  //write out periodicity
  ofOut<<nPeriodic[0]<<" "<<nPeriodic[1]<<" "<<nPeriodic[2]<<std::endl;
  
  //write out number of 1D zones
  ofOut<<nNumZones1D<<std::endl;
  
  //write out number of ghost cells
  ofOut<<nNumGhostCells<<std::endl;
  
  //write out number of variables
  ofOut<<7<<std::endl;
  
  //write out variable info
  //output interior mass info (M_r,1)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<0<<"  ";//radial centered, theta and phi not defined, time independent
  
  //output delta M info (DM,2)
  ofOut<<0<<" "<<-1<<" "<<-1<<" "<<0<<"  ";//radial centered, theta and phi not defined, time independent
  
  //output radius info (r,3)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi not defined, time dependent
  
  //output density info (rho,4)
  ofOut<<0<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//r,theta, and phi centered, time dependent
  
  //output radial velocity (u,5)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi centered, time dependent
  
  //output radial grid velocity (u_0,6)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi centered, time dependent

  //output internal energy info (E,7)
  ofOut<<0<<" "<<-1<<" "<<-1<<" "<<1<<"  "<<std::endl;//r, theta, and phi centered, time dependent
  
  //write out variables
  
  //write out M_r - 1D
  int nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdM[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out Del M
  for(int i=vecdMDel.size()-1;i>=0;i--){//start at the center and work outward
    ofOut<<-1.0*vecdMDel[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out r
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdR[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out rho
  nStart=vecdRho.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdRho.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdRho[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out u
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<dU[i][0][0]<<" "<<std::endl<<std::endl;//same as radial velocity
  }
  ofOut<<std::endl;
  
  //write out u_0 - 1D
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<dU0[i]<<" "<<std::endl<<std::endl;//same as radial velocity
  }
  ofOut<<std::endl;
  
  //write out E
  nStart=vecdE.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdE.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdE[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  ofOut.close();
}
void writeModel_RT_TEOS(){
  
  std::ofstream ofOut;
  ofOut.open(sOutPutfile.c_str());
  if(!ofOut.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sOutPutfile<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //file type and version
  ofOut<<"a 1"<<std::endl;
  
  //set double output precision
  ofOut.precision(nPrecision);
  ofOut.unsetf(std::ios::fixed);
  ofOut.setf(std::ios::scientific);
  
  //write out start time
  ofOut<<0.0<<" "<<0<<std::endl;
  
  //write out time step for old step
  ofOut<<dTimeStep_TEOS()*dTimeStepFactor<<std::endl;
  
  //write out time step for new step
  ofOut<<dTimeStep_TEOS()*dTimeStepFactor<<std::endl;
  
  //write out alpha
  ofOut<<dAlpha<<std::endl;
  
  //write out the equation of state file name and it's length
  ofOut<<sEOSFile.size()<<" "<<sEOSFile<<std::endl;
  
  //write place holder for artificial viscosity
  double dTemp=0.0;
  ofOut<<dTemp<<std::endl;
  
  //write place holder for artificial viscosity threshold
  ofOut<<dTemp<<std::endl;
  
  //output dimensions
  ofOut<<vecdP.size()-2*nNumGhostCells<<" "<<nNumTheta<<" "<<nNumPhi<<std::endl;
  
  //write out periodicity
  ofOut<<nPeriodic[0]<<" "<<nPeriodic[1]<<" "<<nPeriodic[2]<<std::endl;
  
  //write out number of 1D zones
  ofOut<<nNumZones1D<<std::endl;
  
  //write out number of ghost cells
  ofOut<<nNumGhostCells<<std::endl;
  
  //write out number of variables
  ofOut<<9<<std::endl;
  
  //write out variable info
  //output interior mass info (M)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<0<<"  ";//radial centered, theta and phi not defined, time independent
  
  //output theta (THETA)
  ofOut<<-1<<" "<<1<<" "<<-1<<" "<<0<<"  ";//radial undefined, theta interface, phi undefined, time independent
  
  //output delta M info (DM)
  ofOut<<0<<" "<<-1<<" "<<-1<<" "<<0<<"  ";//radial centered, theta and phi not defined, time independent
  
  //output radius info (R)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi not defined, time dependent
  
  //output density info (D)
  ofOut<<0<<" "<<0<<" "<<-1<<" "<<1<<"  ";//r,theta, and phi centered, time dependent
  
  //output radial velocity (U)
  ofOut<<1<<" "<<0<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi centered, time dependent
  
  //output radial grid velocity (U0)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi centered, time dependent
  
  //output theta velocity (V)
  ofOut<<0<<" "<<1<<" "<<-1<<" "<<1<<"  ";//theta interface, r and phi centered, time dependent
  
  //output temperature (T)
  ofOut<<0<<" "<<0<<" "<<-1<<" "<<1<<"  "<<std::endl;//r, theta, and phi centered, time dependent
  
  //write out variables
  
  //write out M_r - 1D
  int nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdM[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out theta - 1D
  double dStartTheta=(90.0-dDeltaTheta*double((nNumTheta+2*nNumGhostCells)/2))*dPi/180.0;
  if(2*(nNumTheta/2)!=nNumTheta){
    dStartTheta-=dDeltaTheta*dPi/360.0;
  }
  unsigned int nNumThetaInt=nNumTheta+2*nNumGhostCells+1;
  if(nPeriodic[1]==1){
    dStartTheta+=dDeltaTheta*dPi/180.0;
    nNumThetaInt-=1;
  }
  for(unsigned int i=0;i<nNumThetaInt;i++){
    ofOut<<dStartTheta+(double(i)*dDeltaTheta)*dPi/180.0<<" "<<std::endl;
  }
  ofOut<<std::endl;
  ofOut<<std::endl;
  
  //write out Del M - 1D
  for(int i=vecdMDel.size()-1;i>=0;i--){//start at the center and work outward
    ofOut<<-1.0*vecdMDel[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out r - 1D
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdR[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out rho - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<vecdRho[i]<<" "<<std::endl<<std::endl;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      ofOut<<vecdRho[i]<<" ";//spherically symetric
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out u - 3D
  //write out 1D region
  nStart=vecdR.size()-1;
  unsigned int nEnd=vecdR.size()-nNumZones1D-1-nNumGhostCells;//needs an extra zone since it is interface centered
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart-=1;
  }
  for(unsigned int i=nStart;i>=nEnd;i--){//start at the center and work outward
    ofOut<<dU[i][0][0]<<" "<<std::endl<<std::endl;
  }
  
  //write out 3D region
  for(int i=nEnd-1;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      ofOut<<dU[i][j][0]<<" ";//spherically symetric
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out u_0 - 1D
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<dU0[i]<<" "<<std::endl<<std::endl;//same as radial velocity
  }
  ofOut<<std::endl;
  
  //write out v - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<dV[i][0][0]<<" "<<std::endl<<std::endl;
  }
  
  int nInt=1;//if not periodic add inner interface
  if(nPeriodic[1]==1){//if periodic don't add inner interface
    nInt=0;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells+nInt;j++){
      ofOut<<dV[i][j][0]<<" ";//spherically symetric
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;

  //write out T - 3D
  //write out 1D region
  for(unsigned int i=vecdT.size()-1;i>=vecdT.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<vecdT[i]<<" "<<std::endl<<std::endl;
  }
  
  //write out 3D region
  for(int i=vecdT.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      ofOut<<vecdT[i]<<" ";//spherically symetric
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  ofOut.close();
}
void writeModel_RT_GL(){
  
  std::ofstream ofOut;
  ofOut.open(sOutPutfile.c_str());
  if(!ofOut.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sOutPutfile<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //file type and version
  ofOut<<"a 1"<<std::endl;
  
  //set double output precision
  ofOut.precision(nPrecision);
  ofOut.unsetf(std::ios::fixed);
  ofOut.setf(std::ios::scientific);
  
  //write out start time
  ofOut<<0.0<<" "<<0<<std::endl;
  
  //write out time step for old step
  ofOut<<dTimeStep_GL()*dTimeStepFactor<<std::endl;
  
  //write out time step for new step
  ofOut<<dTimeStep_GL()*dTimeStepFactor<<std::endl;
  
  //write out alpha
  ofOut<<dAlpha<<std::endl;
  
  ofOut<<0<<" "<<dGamma<<std::endl;
  
  //write place holder for artificial viscosity
  double dTemp=0.0;
  ofOut<<dTemp<<std::endl;
  
  //write place holder for artificial viscosity threshold
  ofOut<<dTemp<<std::endl;
  
  //output dimensions
  ofOut<<vecdP.size()-2*nNumGhostCells<<" "<<nNumTheta<<" "<<nNumPhi<<std::endl;
  
  //write out periodicity
  ofOut<<nPeriodic[0]<<" "<<nPeriodic[1]<<" "<<nPeriodic[2]<<std::endl;
  
  //write out number of 1D zones
  ofOut<<nNumZones1D<<std::endl;
  
  //write out number of ghost cells
  ofOut<<nNumGhostCells<<std::endl;
  
  //write out number of variables
  ofOut<<9<<std::endl;
  
  //write out variable info
  //output interior mass info (M_r,0)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<0<<"  ";//radial centered, theta and phi not defined, time independent
  
  //output theta (theta,1)
  ofOut<<-1<<" "<<1<<" "<<-1<<" "<<0<<"  ";//radial undefined, theta interface, phi undefined, time independent
  
  //output delta M info (DM,2)
  ofOut<<0<<" "<<-1<<" "<<-1<<" "<<0<<"  ";//radial centered, theta and phi not defined, time independent
  
  //output radius info (r,3)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi not defined, time dependent
  
  //output density info (rho,4)
  ofOut<<0<<" "<<0<<" "<<-1<<" "<<1<<"  ";//r,theta, and phi centered, time dependent
  
  //output radial velocity (u,5)
  ofOut<<1<<" "<<0<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi centered, time dependent
  
  //output radial grid velocity (u_0,6)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi centered, time dependent
  
  //output theta velocity (v,7)
  ofOut<<0<<" "<<1<<" "<<-1<<" "<<1<<"  ";//radial center, theta interface and phi centered, time dependent
  
  //output internal energy info (E,8)
  ofOut<<0<<" "<<0<<" "<<-1<<" "<<1<<"  "<<std::endl;//r, theta, and phi centered, time dependent
  
  //write out variables
  
  //write out M_r
  int nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdM[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out theta - 1D
  double dStartTheta=(90.0-dDeltaTheta*double((nNumTheta+2*nNumGhostCells)/2))*dPi/180.0;  
  if(2*(nNumTheta/2)!=nNumTheta){
    dStartTheta-=dDeltaTheta*dPi/360.0;
  }
  unsigned int nNumThetaInt=nNumTheta+2*nNumGhostCells+1;
  if(nPeriodic[1]==1){
    dStartTheta+=dDeltaTheta*dPi/180.0;
    nNumThetaInt-=1;
  }
  for(unsigned int i=0;i<nNumThetaInt;i++){
    ofOut<<dStartTheta+(double(i)*dDeltaTheta)*dPi/180.0<<" "<<std::endl;
  }
  ofOut<<std::endl;
  ofOut<<std::endl;
  
  //write out Del M - 1D
  for(int i=vecdMDel.size()-1;i>=0;i--){//start at the center and work outward
    ofOut<<-1.0*vecdMDel[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out r
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdR[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out rho
  
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<vecdRho[i]<<" "<<std::endl<<std::endl;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      ofOut<<vecdRho[i]<<" ";//spherically symetric
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out u
  //write out 1D region
  nStart=vecdR.size()-1;
  unsigned int nEnd=vecdR.size()-nNumZones1D-1-nNumGhostCells;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart-=1;
  }
  for(unsigned int i=nStart;i>=nEnd;i--){//start at the center and work outward
    ofOut<<dU[i][0][0]<<" "<<std::endl<<std::endl;
  }
  
  //write out 3D region
  for(int i=nEnd-1;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      ofOut<<dU[i][j][0]<<" ";
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out u_0
  //write out 1D region
  nStart=vecdR.size()-1;
  nEnd=vecdR.size()-nNumZones1D-1-nNumGhostCells;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart-=1;
  }
  for(unsigned int i=nStart;i>=nEnd;i--){//start at the center and work outward
    ofOut<<dU0[i]<<" "<<std::endl<<std::endl;
  }
  
  //write out 3D region
  for(int i=nEnd-1;i>=0;i--){//start at the center and work outward
    ofOut<<dU0[i]<<" ";
    ofOut<<std::endl;//next theta
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out v - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<dV[i][0][0]<<" "<<std::endl<<std::endl;
  }
  
  int nInt=1;//if not periodic add inner interface
  if(nPeriodic[1]==1){//if periodic don't add inner interface
    nInt=0;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells+nInt;j++){
      ofOut<<dV[i][j][0]<<" ";//spherically symetric
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out E
  
  //write out 1D region
  for(unsigned int i=vecdE.size()-1;i>=vecdE.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<vecdE[i]<<" "<<std::endl<<std::endl;
  }
  
  //write out 3D region
  for(int i=vecdE.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      ofOut<<vecdE[i]<<" ";//spherically symetric
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  ofOut.close();
}
void writeModel_RTP_TEOS(){
  
  std::ofstream ofOut;
  ofOut.open(sOutPutfile.c_str());
  if(!ofOut.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sOutPutfile<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //file type and version
  ofOut<<"a 1"<<std::endl;
  
  //set double output precision
  ofOut.precision(nPrecision);
  ofOut.unsetf(std::ios::fixed);
  ofOut.setf(std::ios::scientific);
  
  //write out start time
  ofOut<<0.0<<" "<<0<<std::endl;
  
  //write out time step for old step
  ofOut<<dTimeStep_TEOS()*dTimeStepFactor<<std::endl;
  
  //write out time step for new step
  ofOut<<dTimeStep_TEOS()*dTimeStepFactor<<std::endl;
  
  //write out alpha
  ofOut<<dAlpha<<std::endl;
  
  //write out the equation of state file name and it's length
  ofOut<<sEOSFile.size()<<" "<<sEOSFile<<std::endl;
  
  //write place holder for artificial viscosity
  double dTemp=0.0;
  ofOut<<dTemp<<std::endl;
  
  //write place holder for artificial viscosity threshold
  ofOut<<dTemp<<std::endl;
  
  //output dimensions
  ofOut<<vecdP.size()-2*nNumGhostCells<<" "<<nNumTheta<<" "<<nNumPhi<<std::endl;
  
  //write out periodicity
  ofOut<<nPeriodic[0]<<" "<<nPeriodic[1]<<" "<<nPeriodic[2]<<std::endl;
  
  //write out number of 1D zones
  ofOut<<nNumZones1D<<std::endl;
  
  //write out number of ghost cells
  ofOut<<nNumGhostCells<<std::endl;
  
  //write out number of variables
  ofOut<<11<<std::endl;
  
  //write out variable info
  //output interior mass info (M)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<0<<"  ";//radial centered, theta and phi not defined, time independent
  
  //output theta (THETA)
  ofOut<<-1<<" "<<1<<" "<<-1<<" "<<0<<"  ";//radial undefined, theta interface, phi undefined, time independent
  
  //output phi (PHI)
  ofOut<<-1<<" "<<-1<<" "<<1<<" "<<0<<"  ";//radial and theta undefined, phi interface, time independent
  
  //output delta M info (DM)
  ofOut<<0<<" "<<-1<<" "<<-1<<" "<<0<<"  ";//radial centered, theta and phi not defined, time independent
  
  //output radius info (R)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi not defined, time dependent
  
  //output density info (D)
  ofOut<<0<<" "<<0<<" "<<0<<" "<<1<<"  ";//r,theta, and phi centered, time dependent
  
  //output radial velocity (U)
  ofOut<<1<<" "<<0<<" "<<0<<" "<<1<<"  ";//radial interface, theta and phi centered, time dependent
  
  //output radial grid velocity (U0)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi centered, time dependent
  
  //output theta velocity (V)
  ofOut<<0<<" "<<1<<" "<<0<<" "<<1<<"  ";//theta interface, r and phi centered, time dependent
  
  //output phi velocity info (W)
  ofOut<<0<<" "<<0<<" "<<1<<" "<<1<<"  ";//r and theta centered, phi interface, time dependent
  
  //output temperature (T)
  ofOut<<0<<" "<<0<<" "<<0<<" "<<1<<"  "<<std::endl;//r, theta, and phi centered, time dependent
  
  //write out variables
  
  //write out M_r - 1D
  int nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdM[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out theta - 1D
  double dStartTheta=(90.0-dDeltaTheta*double((nNumTheta+2*nNumGhostCells)/2))*dPi/180.0;
  if(2*(nNumTheta/2)!=nNumTheta){
    dStartTheta-=dDeltaTheta*dPi/360.0;
  }
  unsigned int nNumThetaInt=nNumTheta+2*nNumGhostCells+1;
  if(nPeriodic[1]==1){
    dStartTheta+=dDeltaTheta*dPi/180.0;
    nNumThetaInt-=1;
  }
  for(unsigned int i=0;i<nNumThetaInt;i++){
    ofOut<<dStartTheta+(double(i)*dDeltaTheta)*dPi/180.0<<" "<<std::endl;
  }
  ofOut<<std::endl;
  ofOut<<std::endl;
  
  //write out phi's
  
  //write out phi's
  double dStartPhi=(0.0-dDeltaPhi*double((nNumPhi+2*nNumGhostCells)/2))*dPi/180.0;
  if(2*(nNumPhi/2)!=nNumPhi){
    dStartPhi-=dDeltaPhi*dPi/360.0;
  }
  unsigned int nNumPhiInt=nNumPhi+2*nNumGhostCells+1;
  if(nPeriodic[2]==1){
    dStartPhi+=dDeltaPhi*dPi/180.0;
    nNumPhiInt-=1;
  }
  for(unsigned int i=0;i<nNumPhiInt;i++){
    ofOut<<dStartPhi+double(i)*dDeltaPhi*dPi/180.0<<" ";
  }
  ofOut<<std::endl<<std::endl<<std::endl;
  
  //write out Del M - 1D
  for(int i=vecdMDel.size()-1;i>=0;i--){//start at the center and work outward
    ofOut<<-1.0*vecdMDel[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out r - 1D
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdR[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out rho - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<vecdRho[i]<<" "<<std::endl<<std::endl;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut<<vecdRho[i]<<" ";//spherically symetric
      }
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out u - 3D
  //write out 1D region
  nStart=vecdR.size()-1;
  unsigned int nEnd=vecdR.size()-nNumZones1D-1-nNumGhostCells;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart-=1;
  }
  for(unsigned int i=nStart;i>=nEnd;i--){//start at the center and work outward
    ofOut<<dU[i][0][0]<<" "<<std::endl<<std::endl;
  }
  
  //write out 3D region
  for(int i=nEnd-1;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut<<dU[i][j][k]<<" ";//spherically symetric
      }
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out u_0 - 1D
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<dU0[i]<<" "<<std::endl<<std::endl;//same as radial velocity
  }
  ofOut<<std::endl;
  
  //write out v - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<dV[i][0][0]<<" "<<std::endl<<std::endl;
  }
  
  int nInt=1;//if not periodic add inner interface
  if(nPeriodic[1]==1){//if periodic don't add inner interface
    nInt=0;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells+nInt;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut<<dV[i][j][k]<<" ";//spherically symetric
      }
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out w - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<dW[i][0][0]<<" "<<std::endl<<std::endl;
  }
  
  nInt=1;//if not periodic add inner interface
  if(nPeriodic[2]==1){//if periodic don't add inner interface
    nInt=0;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells+nInt;k++){
        ofOut<<dW[i][j][k]<<" ";//spherically symetric
      }
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out T - 3D
  //write out 1D region
  for(unsigned int i=vecdT.size()-1;i>=vecdT.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<vecdT[i]<<" "<<std::endl<<std::endl;
  }
  
  //write out 3D region
  for(int i=vecdT.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut<<vecdT[i]<<" ";//spherically symetric
      }
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  ofOut.close();
}
void writeModel_RTP_GL(){
  
  std::ofstream ofOut;
  ofOut.open(sOutPutfile.c_str());
  if(!ofOut.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sOutPutfile<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //file type and version
  ofOut<<"a 1"<<std::endl;
  
  //set double output precision
  ofOut.precision(nPrecision);
  ofOut.unsetf(std::ios::fixed);
  ofOut.setf(std::ios::scientific);
  
  //write out start time
  ofOut<<0.0<<" "<<0<<std::endl;
  
  //write out time step for old step
  ofOut<<dTimeStep_GL()*dTimeStepFactor<<std::endl;
  
  //write out time step for new step
  ofOut<<dTimeStep_GL()*dTimeStepFactor<<std::endl;
  
  //write out alpha
  ofOut<<dAlpha<<std::endl;
  
  ofOut<<0<<" "<<dGamma<<std::endl;
  
  //write place holder for artificial viscosity
  double dTemp=0.0;
  ofOut<<dTemp<<std::endl;
  
  //write place holder for artificial viscosity threshold
  ofOut<<dTemp<<std::endl;
  
  //output dimensions
  ofOut<<vecdP.size()-2*nNumGhostCells<<" "<<nNumTheta<<" "<<nNumPhi<<std::endl;
  
  //write out periodicity
  ofOut<<nPeriodic[0]<<" "<<nPeriodic[1]<<" "<<nPeriodic[2]<<std::endl;
  
  //write out number of 1D zones
  ofOut<<nNumZones1D<<std::endl;
  
  //write out number of ghost cells
  ofOut<<nNumGhostCells<<std::endl;
  
  //write out variable info
  
  //write out number of variables
  ofOut<<11<<std::endl;
  
  //output interior mass info (M_r,1)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<0<<"  ";//radial centered, theta and phi not defined, time independent
  
  //output theta (theta,2)
  ofOut<<-1<<" "<<1<<" "<<-1<<" "<<0<<"  ";//radial undefined, theta interface, phi undefined, time independent
  
  //output phi (phi,3)
  ofOut<<-1<<" "<<-1<<" "<<1<<" "<<0<<"  ";//radial and theta undefined, phi interface, time independent
  
  //output delta M info (DM,4)
  ofOut<<0<<" "<<-1<<" "<<-1<<" "<<0<<"  ";//radial centered, theta and phi not defined, time independent
  
  //output radius info (r,5)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi not defined, time dependent
  
  //output density info (rho,6)
  ofOut<<0<<" "<<0<<" "<<0<<" "<<1<<"  ";//r,theta, and phi centered, time dependent
  
  //output radial velocity (u,7)
  ofOut<<1<<" "<<0<<" "<<0<<" "<<1<<"  ";//radial interface, theta and phi centered, time dependent
  
  //output radial grid velocity (u_0,8)
  ofOut<<1<<" "<<-1<<" "<<-1<<" "<<1<<"  ";//radial interface, theta and phi centered, time dependent

  //output theta velocity (v,9)
  ofOut<<0<<" "<<1<<" "<<0<<" "<<1<<"  ";//theta interface, r and phi centered, time dependent
  
  //output phi velocity info (w,10)
  ofOut<<0<<" "<<0<<" "<<1<<" "<<1<<"  ";//r and theta centered, phi interface, time dependent

  //output internal energy info (E,11)
  ofOut<<0<<" "<<0<<" "<<0<<" "<<1<<"  "<<std::endl;//r, theta, and phi centered, time dependent
  
  //write out variables
  
  //write out M_r - 1D
  int nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdM[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out theta - 1D
  double dStartTheta=(90.0-dDeltaTheta*double((nNumTheta+2*nNumGhostCells)/2))*dPi/180.0;  
  if(2*(nNumTheta/2)!=nNumTheta){
    dStartTheta-=dDeltaTheta*dPi/360.0;
  }
  unsigned int nNumThetaInt=nNumTheta+2*nNumGhostCells+1;
  if(nPeriodic[1]==1){
    dStartTheta+=dDeltaTheta*dPi/180.0;
    nNumThetaInt-=1;
  }
  for(unsigned int i=0;i<nNumThetaInt;i++){
    ofOut<<dStartTheta+(double(i)*dDeltaTheta)*dPi/180.0<<" "<<std::endl;
  }
  ofOut<<std::endl;
  ofOut<<std::endl;
  
  //write out phi's
  double dStartPhi=(0.0-dDeltaPhi*double((nNumPhi+2*nNumGhostCells)/2))*dPi/180.0;
  if(2*(nNumPhi/2)!=nNumPhi){
    dStartPhi-=dDeltaPhi*dPi/360.0;
  }
  unsigned int nNumPhiInt=nNumPhi+2*nNumGhostCells+1;
  if(nPeriodic[2]==1){
    dStartPhi+=dDeltaPhi*dPi/180.0;
    nNumPhiInt-=1;
  }
  
  for(unsigned int i=0;i<nNumPhiInt;i++){
    ofOut<<dStartPhi+double(i)*dDeltaPhi*dPi/180.0<<" ";
  }
  ofOut<<std::endl<<std::endl<<std::endl;
  
  //write out Del M - 1D
  for(int i=vecdMDel.size()-1;i>=0;i--){//start at the center and work outward
    ofOut<<-1.0*vecdMDel[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out r - 1D
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<vecdR[i]<<" "<<std::endl<<std::endl;
  }
  ofOut<<std::endl;
  
  //write out rho - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<vecdRho[i]<<" "<<std::endl<<std::endl;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut<<vecdRho[i]<<" ";//spherically symetric
      }
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out u - 3D
  //write out 1D region
  nStart=vecdR.size()-1;
  unsigned int nEnd=vecdR.size()-nNumZones1D-1-nNumGhostCells;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart-=1;
  }
  for(unsigned int i=nStart;i>=nEnd;i--){//start at the center and work outward
    ofOut<<dU[i][0][0]<<" "<<std::endl<<std::endl;
  }
  
  //write out 3D region
  for(int i=nEnd-1;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut<<dU[i][j][k]<<" ";//spherically symetric
      }
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out u_0 - 1D
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut<<dU0[i]<<" "<<std::endl<<std::endl;//same as radial velocity
  }
  ofOut<<std::endl;
  
  //write out v - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<dV[i][0][0]<<" "<<std::endl<<std::endl;
  }
  
  int nInt=1;//if not periodic add inner interface
  if(nPeriodic[1]==1){//if periodic don't add inner interface
    nInt=0;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells+nInt;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut<<dV[i][j][k]<<" ";//spherically symetric
      }
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out w - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<dW[i][0][0]<<" "<<std::endl<<std::endl;
  }
  
  nInt=1;//if not periodic add inner interface
  if(nPeriodic[2]==1){//if periodic don't add inner interface
    nInt=0;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells+nInt;k++){
        ofOut<<dW[i][j][k]<<" ";//spherically symetric
      }
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;
  
  //write out E - 3D
  //write out 1D region
  for(unsigned int i=vecdE.size()-1;i>=vecdE.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut<<vecdE[i]<<" "<<std::endl<<std::endl;
  }
  
  //write out 3D region
  for(int i=vecdE.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut<<vecdE[i]<<" ";//spherically symetric
      }
      ofOut<<std::endl;//next theta
    }
    ofOut<<std::endl;//next r
  }
  ofOut<<std::endl;

  ofOut.close();
}
void writeModel_Bin_R_TEOS(){
  
  std::ofstream ofOut;
  ofOut.open(sOutPutfile.c_str(),std::ios::binary);
  if(!ofOut.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sOutPutfile<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //file type
  char cTemp='b';
  ofOut.write((char*)(&cTemp),sizeof(char));
  
  //file version
  int nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out start time
  double dTemp=0.0;
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //write out timestep index
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out time step
  dTemp=dTimeStep_TEOS()*dTimeStepFactor;
  ofOut.write((char*)(&dTemp),sizeof(double));//for old step
  
  ofOut.write((char*)(&dTemp),sizeof(double));//for new step
  
  //write out alpha
  ofOut.write((char*)(&dAlpha),sizeof(double));
  
  //write out gamma index
  nTemp=sEOSFile.size();
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out gamma
  ofOut.write(sEOSFile.c_str(),sEOSFile.size()*sizeof(char));
  
  //write place holder for artificial viscosity
  dTemp=0.0;
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //write place holder for artificial viscosity threshold
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //output dimensions
  nTemp=vecdP.size()-2*nNumGhostCells;
  ofOut.write((char*)(&nTemp),sizeof(int));
  ofOut.write((char*)(&nNumTheta),sizeof(int));
  ofOut.write((char*)(&nNumPhi),sizeof(int));
  
  //write out periodicity
  ofOut.write((char*)(&nPeriodic[0]),sizeof(int));
  ofOut.write((char*)(&nPeriodic[1]),sizeof(int));
  ofOut.write((char*)(&nPeriodic[2]),sizeof(int));
  
  //write out number of 1D zones
  ofOut.write((char*)(&nNumZones1D),sizeof(int));
  
  //write out number of ghost cells
  ofOut.write((char*)(&nNumGhostCells),sizeof(int));
  
  //write out number of variables
  nTemp=7;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out variable info
  //output interior mass info (M_r,1)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output delta M info (DM,2)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radius info (r,3)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output density info (rho,4)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial velocity (u,5)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial grid velocity (u_0,6)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));

  //output temperature energy info (T,7)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out variables
  
  //write out M_r - 1D
  int nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdM[i]),sizeof(double));
  }
  
  //write out Del M
  for(int i=vecdMDel.size()-1;i>=0;i--){//start at the center and work outward
    dTemp=-1.0*vecdMDel[i];
    ofOut.write((char*)(&dTemp),sizeof(double));
  }
  
  //write out r
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdR[i]),sizeof(double));
  }
  
  //write out rho
  nStart=vecdRho.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdRho.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdRho[i]),sizeof(double));
  }
  
  //write out u
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&dU[i][0][0]),sizeof(double));
  }
  
  //write out u_0 - 1D
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&dU0[i]),sizeof(double));
  }
  
  //write out T
  nStart=vecdT.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdT.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdT[i]),sizeof(double));
  }
  ofOut.close();
}
void writeModel_Bin_R_GL(){
  
  std::ofstream ofOut;
  ofOut.open(sOutPutfile.c_str(),std::ios::binary);
  if(!ofOut.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sOutPutfile<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //file type
  char cTemp='b';
  ofOut.write((char*)(&cTemp),sizeof(char));
  
  //file version
  int nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out start time
  double dTemp=0.0;
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //write out timestep index
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out time step
  dTemp=dTimeStepFactor*dTimeStep_GL();
  ofOut.write((char*)(&dTemp),sizeof(double));//for old step
  ofOut.write((char*)(&dTemp),sizeof(double));//for new step
  
  //write out alpha
  ofOut.write((char*)(&dAlpha),sizeof(double));
  
  //write out gamma index
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out gamma
  ofOut.write((char*)(&dGamma),sizeof(double));
  
  //write place holder for artificial viscosity
  dTemp=0.0;
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //write place holder for artificial viscosity threshold
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //output dimensions
  nTemp=vecdP.size()-2*nNumGhostCells;
  ofOut.write((char*)(&nTemp),sizeof(int));
  ofOut.write((char*)(&nNumTheta),sizeof(int));
  ofOut.write((char*)(&nNumPhi),sizeof(int));
  
  //write out periodicity
  ofOut.write((char*)(&nPeriodic[0]),sizeof(int));
  ofOut.write((char*)(&nPeriodic[1]),sizeof(int));
  ofOut.write((char*)(&nPeriodic[2]),sizeof(int));
  
  //write out number of 1D zones
  ofOut.write((char*)(&nNumZones1D),sizeof(int));
  
  //write out number of ghost cells
  ofOut.write((char*)(&nNumGhostCells),sizeof(int));
  
  //write out number of variables
  nTemp=7;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out variable info
  //output interior mass info (M_r,1)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output delta M info (DM,2)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radius info (r,3)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output density info (rho,4)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial velocity (u,5)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial grid velocity (u_0,6)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));

  //output internal energy info (E,7)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out variables
  
  //write out M_r - 1D
  int nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdM[i]),sizeof(double));
  }
  
  //write out Del M
  for(int i=vecdMDel.size()-1;i>=0;i--){//start at the center and work outward
    dTemp=-1.0*vecdMDel[i];
    ofOut.write((char*)(&dTemp),sizeof(double));
  }
  
  //write out r
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdR[i]),sizeof(double));
  }
  
  //write out rho
  nStart=vecdRho.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdRho.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdRho[i]),sizeof(double));
  }
  
  //write out u
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&dU[i][0][0]),sizeof(double));
  }
  
  //write out u_0 - 1D
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&dU[i][0][0]),sizeof(double));
  }
  
  //write out E
  nStart=vecdE.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdE.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdE[i]),sizeof(double));
  }
  ofOut.close();
}
void writeModel_Bin_RT_TEOS(){
  
  std::ofstream ofOut;
  ofOut.open(sOutPutfile.c_str(),std::ios::binary);
  if(!ofOut.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sOutPutfile<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //file type
  char cTemp='b';
  ofOut.write((char*)(&cTemp),sizeof(char));
  
  //file version
  int nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out start time
  double dTemp=0.0;
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //write out timestep index
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out time step
  dTemp=dTimeStepFactor*dTimeStep_TEOS();
  ofOut.write((char*)(&dTemp),sizeof(double));//for old step
  ofOut.write((char*)(&dTemp),sizeof(double));//for new step
  
  //write out alpha
  ofOut.write((char*)(&dAlpha),sizeof(double));
  
  //write out gamma index
  nTemp=sEOSFile.size();
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out gamma
  ofOut.write((char*)sEOSFile.c_str(),sEOSFile.size()*sizeof(char));
  
  //write place holder for artificial viscosity
  dTemp=0.0;
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //write place holder for artificial viscosity threshold
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //output dimensions
  nTemp=vecdP.size()-2*nNumGhostCells;
  ofOut.write((char*)(&nTemp),sizeof(int));
  ofOut.write((char*)(&nNumTheta),sizeof(int));
  ofOut.write((char*)(&nNumPhi),sizeof(int));
  
  //write out periodicity
  ofOut.write((char*)(&nPeriodic[0]),sizeof(int));
  ofOut.write((char*)(&nPeriodic[1]),sizeof(int));
  ofOut.write((char*)(&nPeriodic[2]),sizeof(int));
  
  //write out number of 1D zones
  ofOut.write((char*)(&nNumZones1D),sizeof(int));
  
  //write out number of ghost cells
  ofOut.write((char*)(&nNumGhostCells),sizeof(int));
  
  //write out number of variables
  nTemp=9;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out variable info
  //output interior mass info (M_r,0)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output interior mass info (theta,1)
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output delta M info (DM,2)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radius info (r,3)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output density info (rho,4)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial velocity (u,5)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial grid velocity (u_0,6)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  
  //output radial grid velocity (u_0,7)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));

  //output internal energy info (T,8)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out variables
  
  //write out M_r - 1D
  int nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdM[i]),sizeof(double));
  }
  
  //write out theta - 1D
  double dStartTheta=(90.0-dDeltaTheta*double((nNumTheta+2*nNumGhostCells)/2))*dPi/180.0;  
  if(2*(nNumTheta/2)!=nNumTheta){
    dStartTheta-=dDeltaTheta*dPi/360.0;
  }
  unsigned int nNumThetaInt=nNumTheta+2*nNumGhostCells+1;
  if(nPeriodic[1]==1){
    dStartTheta+=dDeltaTheta*dPi/180.0;
    nNumThetaInt-=1;
  }
  for(unsigned int i=0;i<nNumThetaInt;i++){
    dTemp=dStartTheta+(double(i)*dDeltaTheta)*dPi/180.0;
    ofOut.write((char*)(&dTemp),sizeof(double));
  }
  
  //write out Del M
  for(int i=vecdMDel.size()-1;i>=0;i--){//start at the center and work outward
    dTemp=-1.0*vecdMDel[i];
    ofOut.write((char*)(&dTemp),sizeof(double));
  }
  
  //write out r
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdR[i]),sizeof(double));
  }
  
  //write out rho
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdRho[i]),sizeof(double));
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      ofOut.write((char*)(&vecdRho[i]),sizeof(double));//spherically symetric
    }
  }
  
  //write out u
  //write out 1D region
  nStart=vecdR.size()-1;
  unsigned int nEnd=vecdR.size()-nNumZones1D-1-nNumGhostCells;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart-=1;
  }
  for(unsigned int i=nStart;i>=nEnd;i--){//start at the center and work outward
    ofOut.write((char*)(&dU[i][0][0]),sizeof(double));
  }
  
  //write out 3D region
  for(int i=nEnd-1;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      ofOut.write((char*)(&dU[i][j][0]),sizeof(double));//spherically symetric
    }
  }
  
  //write out u_0 - 1D
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&dU0[i]),sizeof(double));
  }
  
 //write out v - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&dV[i][0][0]),sizeof(double));
  }
  
  int nInt=1;//if not periodic add inner interface
  if(nPeriodic[1]==1){//if periodic don't add inner interface
    nInt=0;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells+nInt;j++){
      ofOut.write((char*)(&dV[i][j][0]),sizeof(double));//spherically symetric
    }
  }
  
  //write out T
  
  //write out 1D region
  for(unsigned int i=vecdT.size()-1;i>=vecdT.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdT[i]),sizeof(double));
  }
  
  //write out 3D region
  for(int i=vecdT.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      ofOut.write((char*)(&vecdT[i]),sizeof(double));//spherically symetric
    }
  }
  ofOut.close();
}
void writeModel_Bin_RT_GL(){
  
  std::ofstream ofOut;
  ofOut.open(sOutPutfile.c_str(),std::ios::binary);
  if(!ofOut.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sOutPutfile<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //file type
  char cTemp='b';
  ofOut.write((char*)(&cTemp),sizeof(char));
  
  //file version
  int nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out start time
  double dTemp=0.0;
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //write out timestep index
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out time step
  dTemp=dTimeStepFactor*dTimeStep_GL();
  ofOut.write((char*)(&dTemp),sizeof(double));//for old step
  ofOut.write((char*)(&dTemp),sizeof(double));//for new step
  
  //write out alpha
  ofOut.write((char*)(&dAlpha),sizeof(double));
  
  //write out gamma index
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out gamma
  ofOut.write((char*)(&dGamma),sizeof(double));
  
  //write place holder for artificial viscosity
  dTemp=0.0;
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //write place holder for artificial viscosity threshold
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //output dimensions
  nTemp=vecdP.size()-2*nNumGhostCells;
  ofOut.write((char*)(&nTemp),sizeof(int));
  ofOut.write((char*)(&nNumTheta),sizeof(int));
  ofOut.write((char*)(&nNumPhi),sizeof(int));
  
  //write out periodicity
  ofOut.write((char*)(&nPeriodic[0]),sizeof(int));
  ofOut.write((char*)(&nPeriodic[1]),sizeof(int));
  ofOut.write((char*)(&nPeriodic[2]),sizeof(int));
  
  //write out number of 1D zones
  ofOut.write((char*)(&nNumZones1D),sizeof(int));
  
  //write out number of ghost cells
  ofOut.write((char*)(&nNumGhostCells),sizeof(int));
  
  //write out number of variables
  nTemp=9;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out variable info
  //output interior mass info (M_r,0)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output interior mass info (theta,1)
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output delta M info (DM,2)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radius info (r,3)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output density info (rho,4)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial velocity (u,5)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial grid velocity (u_0,6)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  
  //output radial grid velocity (v,7)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));

  //output internal energy info (E,8)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out variables
  
  //write out M_r - 1D
  int nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdM[i]),sizeof(double));
  }
  
  //write out theta - 1D
  double dStartTheta=(90.0-dDeltaTheta*double((nNumTheta+2*nNumGhostCells)/2))*dPi/180.0;  
  if(2*(nNumTheta/2)!=nNumTheta){
    dStartTheta-=dDeltaTheta*dPi/360.0;
  }
  unsigned int nNumThetaInt=nNumTheta+2*nNumGhostCells+1;
  if(nPeriodic[1]==1){
    dStartTheta+=dDeltaTheta*dPi/180.0;
    nNumThetaInt-=1;
  }
  for(unsigned int i=0;i<nNumThetaInt;i++){
    dTemp=dStartTheta+(double(i)*dDeltaTheta)*dPi/180.0;
    ofOut.write((char*)(&dTemp),sizeof(double));
  }
  
  //write out Del M
  for(int i=vecdMDel.size()-1;i>=0;i--){//start at the center and work outward
    dTemp=-1.0*vecdMDel[i];
    ofOut.write((char*)(&dTemp),sizeof(double));
  }
  
  //write out r
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdR[i]),sizeof(double));
  }
  
  //write out rho
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdRho[i]),sizeof(double));
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      ofOut.write((char*)(&vecdRho[i]),sizeof(double));//spherically symetric
    }
  }
  
  //write out u
  //write out 1D region
  nStart=vecdR.size()-1;
  unsigned int nEnd=vecdR.size()-nNumZones1D-1-nNumGhostCells;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart-=1;
  }
  for(unsigned int i=nStart;i>=nEnd;i--){//start at the center and work outward
    ofOut.write((char*)(&dU[i][0][0]),sizeof(double));
  }
  
  //write out 3D region
  for(int i=nEnd-1;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      ofOut.write((char*)(&dU[i][j][0]),sizeof(double));//spherically symetric
    }
  }
  
  //write out u_0 - 1D
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&dU0[i]),sizeof(double));
  }
  
 //write out v - 3D
  //write out 1D region
  dTemp=0.0;
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&dV[i][0][0]),sizeof(double));
  }
  
  int nInt=1;//if not periodic add inner interface
  if(nPeriodic[1]==1){//if periodic don't add inner interface
    nInt=0;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells+nInt;j++){
      ofOut.write((char*)(&dV[i][j][0]),sizeof(double));//spherically symetric
    }
  }
  
  //write out E
  
  //write out 1D region
  for(unsigned int i=vecdE.size()-1;i>=vecdE.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdE[i]),sizeof(double));
  }
  
  //write out 3D region
  for(int i=vecdE.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      ofOut.write((char*)(&vecdE[i]),sizeof(double));//spherically symetric
    }
  }
  ofOut.close();
}
void writeModel_Bin_RTP_TEOS(){
  
  std::ofstream ofOut;
  ofOut.open(sOutPutfile.c_str(),std::ios::binary);
  if(!ofOut.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sOutPutfile<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //file type
  char cTemp='b';
  ofOut.write((char*)(&cTemp),sizeof(char));
  
  //file version
  int nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out start time
  double dTemp=0.0;
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //write out timestep index
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out time step
  dTemp=dTimeStepFactor*dTimeStep_TEOS();
  ofOut.write((char*)(&dTemp),sizeof(double));//for old step
  ofOut.write((char*)(&dTemp),sizeof(double));//for new step
  
  //write out alpha
  ofOut.write((char*)(&dAlpha),sizeof(double));
  
  //write out gamma index
  nTemp=sEOSFile.size();
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out gamma
  ofOut.write((char*)sEOSFile.c_str(),sEOSFile.size()*sizeof(char));
  
  //write place holder for artificial viscosity
  dTemp=0.0;
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //write place holder for artificial viscosity threshold
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //output dimensions
  nTemp=vecdP.size()-2*nNumGhostCells;
  ofOut.write((char*)(&nTemp),sizeof(int));
  ofOut.write((char*)(&nNumTheta),sizeof(int));
  ofOut.write((char*)(&nNumPhi),sizeof(int));
  
  //write out periodicity
  ofOut.write((char*)(&nPeriodic[0]),sizeof(int));
  ofOut.write((char*)(&nPeriodic[1]),sizeof(int));
  ofOut.write((char*)(&nPeriodic[2]),sizeof(int));
  
  //write out number of 1D zones
  ofOut.write((char*)(&nNumZones1D),sizeof(int));
  
  //write out number of ghost cells
  ofOut.write((char*)(&nNumGhostCells),sizeof(int));
  
  //write out number of variables
  nTemp=11;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out variable info
  //output interior mass info (M_r,0)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output interior mass info (theta,1)
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output interior mass info (phi,2)
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output delta M info (DM,3)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radius info (r,4)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output density info (rho,5)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial velocity (u,6)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial grid velocity (u_0,7)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial grid velocity (v,8)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial grid velocity (w,9)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output internal energy info (T,10)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out variables
  
  //write out M_r - 1D
  int nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdM[i]),sizeof(double));
  }
  
  //write out theta - 1D
  double dStartTheta=(90.0-dDeltaTheta*double((nNumTheta+2*nNumGhostCells)/2))*dPi/180.0;  
  if(2*(nNumTheta/2)!=nNumTheta){
    dStartTheta-=dDeltaTheta*dPi/360.0;
  }
  unsigned int nNumThetaInt=nNumTheta+2*nNumGhostCells+1;
  if(nPeriodic[1]==1){
    dStartTheta+=dDeltaTheta*dPi/180.0;
    nNumThetaInt-=1;
  }
  for(unsigned int i=0;i<nNumThetaInt;i++){
    dTemp=dStartTheta+(double(i)*dDeltaTheta)*dPi/180.0;
    ofOut.write((char*)(&dTemp),sizeof(double));
  }
  
  //write out phi's
  
  
  //write out phi's
  double dStartPhi=(0.0-dDeltaPhi*double((nNumPhi+2*nNumGhostCells)/2))*dPi/180.0;
  if(2*(nNumPhi/2)!=nNumPhi){
    dStartPhi-=dDeltaPhi*dPi/360.0;
  }
  unsigned int nNumPhiInt=nNumPhi+2*nNumGhostCells+1;
  if(nPeriodic[1]==1){
    dStartPhi+=dDeltaPhi*dPi/180.0;
    nNumPhiInt-=1;
  }
  for(unsigned int i=0;i<nNumPhiInt;i++){
    dTemp=dStartPhi+double(i)*dDeltaPhi*dPi/180.0;
    ofOut.write((char*)(&dTemp),sizeof(double));
  }
  
  //write out Del M
  for(int i=vecdMDel.size()-1;i>=0;i--){//start at the center and work outward
    dTemp=-1.0*vecdMDel[i];
    ofOut.write((char*)(&dTemp),sizeof(double));
  }
  
  //write out r
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdR[i]),sizeof(double));
  }
  
  //write out rho
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdRho[i]),sizeof(double));
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut.write((char*)(&vecdRho[i]),sizeof(double));//spherically symetric
      }
    }
  }
  
  //write out u
  //write out 1D region
  nStart=vecdR.size()-1;
  unsigned int nEnd=vecdR.size()-nNumZones1D-1-nNumGhostCells;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart-=1;
  }
  for(unsigned int i=nStart;i>=nEnd;i--){//start at the center and work outward
    ofOut.write((char*)(&dU[i][0][0]),sizeof(double));
  }
  
  //write out 3D region
  for(int i=nEnd-1;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut.write((char*)(&dU[i][j][k]),sizeof(double));//spherically symetric
      }
    }
  }
  
  //write out u_0 - 1D
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&dU0[i]),sizeof(double));
  }
  
 //write out v - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&dV[i][0][0]),sizeof(double));
  }
  
  int nInt=1;//if not periodic add inner interface
  if(nPeriodic[1]==1){//if periodic don't add inner interface
    nInt=0;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells+nInt;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut.write((char*)(&dV[i][j][k]),sizeof(double));//spherically symetric
      }
    }
  }
  
 //write out w - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&dW[i][0][0]),sizeof(double));
  }
  
  nInt=1;//if not periodic add inner interface
  if(nPeriodic[2]==1){//if periodic don't add inner interface
    nInt=0;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells+nInt;k++){
        ofOut.write((char*)(&dW[i][j][k]),sizeof(double));//spherically symetric
      }
    }
  }
  
  //write out T
  
  //write out 1D region
  for(unsigned int i=vecdT.size()-1;i>=vecdT.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdT[i]),sizeof(double));
  }
  
  //write out 3D region
  for(int i=vecdT.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut.write((char*)(&vecdT[i]),sizeof(double));//spherically symetric
      }
    }
  }
  ofOut.close();
}
void writeModel_Bin_RTP_GL(){
  
  std::ofstream ofOut;
  ofOut.open(sOutPutfile.c_str(),std::ios::binary);
  if(!ofOut.good()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sOutPutfile<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //file type
  char cTemp='b';
  ofOut.write((char*)(&cTemp),sizeof(char));
  
  //file version
  int nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out start time
  double dTemp=0.0;
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //write out timestep index
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out time step
  dTemp=dTimeStepFactor*dTimeStep_GL();
  ofOut.write((char*)(&dTemp),sizeof(double));//for old step
  ofOut.write((char*)(&dTemp),sizeof(double));//for new step
  
  //write out alpha
  ofOut.write((char*)(&dAlpha),sizeof(double));
  
  //write out gamma index
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out gamma
  ofOut.write((char*)(&dGamma),sizeof(double));
  
  //write place holder for artificial viscosity
  dTemp=0.0;
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //write place holder for artificial viscosity threshold
  ofOut.write((char*)(&dTemp),sizeof(double));
  
  //output dimensions
  nTemp=vecdP.size()-2*nNumGhostCells;
  ofOut.write((char*)(&nTemp),sizeof(int));
  ofOut.write((char*)(&nNumTheta),sizeof(int));
  ofOut.write((char*)(&nNumPhi),sizeof(int));
  
  //write out periodicity
  ofOut.write((char*)(&nPeriodic[0]),sizeof(int));
  ofOut.write((char*)(&nPeriodic[1]),sizeof(int));
  ofOut.write((char*)(&nPeriodic[2]),sizeof(int));
  
  //write out number of 1D zones
  ofOut.write((char*)(&nNumZones1D),sizeof(int));
  
  //write out number of ghost cells
  ofOut.write((char*)(&nNumGhostCells),sizeof(int));
  
  //write out number of variables
  nTemp=11;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out variable info
  //output interior mass info (M_r,0)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output interior mass info (theta,1)
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output interior mass info (phi,2)
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output delta M info (DM,3)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radius info (r,4)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output density info (rho,5)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial velocity (u,6)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial grid velocity (u_0,7)
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=-1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial grid velocity (v,8)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output radial grid velocity (w,9)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //output internal energy info (E,10)
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  nTemp=1;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out variables
  
  //write out M_r - 1D
  int nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdM[i]),sizeof(double));
  }
  
  //write out theta - 1D
  double dStartTheta=(90.0-dDeltaTheta*double((nNumTheta+2*nNumGhostCells)/2))*dPi/180.0;  
  if(2*(nNumTheta/2)!=nNumTheta){
    dStartTheta-=dDeltaTheta*dPi/360.0;
  }
  unsigned int nNumThetaInt=nNumTheta+2*nNumGhostCells+1;
  if(nPeriodic[1]==1){
    dStartTheta+=dDeltaTheta*dPi/180.0;
    nNumThetaInt-=1;
  }
  for(unsigned int i=0;i<nNumThetaInt;i++){
    dTemp=dStartTheta+(double(i)*dDeltaTheta)*dPi/180.0;
    ofOut.write((char*)(&dTemp),sizeof(double));
  }
  
  //write out phi's
  //write out phi's
  double dStartPhi=(0.0-dDeltaPhi*double((nNumPhi+2*nNumGhostCells)/2))*dPi/180.0;
  if(2*(nNumPhi/2)!=nNumPhi){
    dStartPhi-=dDeltaPhi*dPi/360.0;
  }
  unsigned int nNumPhiInt=nNumPhi+2*nNumGhostCells+1;
  if(nPeriodic[1]==1){
    dStartPhi+=dDeltaPhi*dPi/180.0;
    nNumPhiInt-=1;
  }
  for(unsigned int i=0;i<nNumPhiInt;i++){
    dTemp=dStartPhi+double(i)*dDeltaPhi*dPi/180.0;
    ofOut.write((char*)(&dTemp),sizeof(double));
  }
  
  //write out Del M
  for(int i=vecdMDel.size()-1;i>=0;i--){//start at the center and work outward
    dTemp=-1.0*vecdMDel[i];
    ofOut.write((char*)(&dTemp),sizeof(double));
  }
  
  //write out r
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdR[i]),sizeof(double));
  }
  
  //write out rho
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdRho[i]),sizeof(double));
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut.write((char*)(&vecdRho[i]),sizeof(double));//spherically symetric
      }
    }
  }
  
  //write out u
  //write out 1D region
  nStart=vecdR.size()-1;
  unsigned int nEnd=vecdR.size()-nNumZones1D-1-nNumGhostCells;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart-=1;
  }
  for(unsigned int i=nStart;i>=nEnd;i--){//start at the center and work outward
    ofOut.write((char*)(&dU[i][0][0]),sizeof(double));
  }
  
  //write out 3D region
  for(int i=nEnd-1;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut.write((char*)(&dU[i][j][k]),sizeof(double));//spherically symetric
      }
    }
  }
  
  //write out u_0 - 1D
  nStart=vecdR.size()-1;
  if(nPeriodic[0]==1){//don't need inner interface if periodic
    nStart=vecdR.size()-2;
  }
  for(int i=nStart;i>=0;i--){//start at the center and work outward
    ofOut.write((char*)(&dU0[i]),sizeof(double));
  }
  
 //write out v - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&dV[i][0][0]),sizeof(double));
  }
  
  int nInt=1;//if not periodic add inner interface
  if(nPeriodic[1]==1){//if periodic don't add inner interface
    nInt=0;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells+nInt;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut.write((char*)(&dV[i][j][k]),sizeof(double));//spherically symetric
      }
    }
  }
  
 //write out w - 3D
  //write out 1D region
  for(unsigned int i=vecdRho.size()-1;i>=vecdRho.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&dW[i][0][0]),sizeof(double));
  }
  
  nInt=1;//if not periodic add inner interface
  if(nPeriodic[2]==1){//if periodic don't add inner interface
    nInt=0;
  }
  
  //write out 3D region
  for(int i=vecdRho.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells+nInt;k++){
        ofOut.write((char*)(&dW[i][j][k]),sizeof(double));//spherically symetric
      }
    }
  }
  
  //write out E
  
  //write out 1D region
  for(unsigned int i=vecdE.size()-1;i>=vecdE.size()-nNumZones1D-nNumGhostCells;i--){//start at the center and work outward
    ofOut.write((char*)(&vecdE[i]),sizeof(double));
  }
  
  //write out 3D region
  for(int i=vecdE.size()-nNumZones1D-1-nNumGhostCells;i>=0;i--){//start at the center and work outward
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        ofOut.write((char*)(&vecdE[i]),sizeof(double));//spherically symetric
      }
    }
  }
  ofOut.close();
}
void writeModelToScreen_TEOS(){
  int nWidth=nPrecision+9;
  std::cout<<std::setw(3)<<"#"
    <<std::setw(nWidth)<<"M_r"
    <<std::setw(nWidth)<<"Del_M"
    <<std::setw(nWidth)<<"r"
    <<std::setw(nWidth)<<"rho"
    <<std::setw(nWidth)<<"E"
    <<std::setw(nWidth)<<"T"
    <<std::setw(nWidth)<<"P"
    <<std::setw(nWidth)<<"u"<<std::endl;
  std::cout.setf(std::ios::scientific);
  std::cout.precision(nPrecision);
  for(unsigned int i=0;i<vecdP.size();i++){
    std::cout<<std::setw(3)<<i
      <<std::setw(nWidth)<<vecdM[i]
      <<std::setw(nWidth)<<vecdMDel[i]
      <<std::setw(nWidth)<<vecdR[i]
      <<std::setw(nWidth)<<vecdRho[i]
      <<std::setw(nWidth)<<vecdE[i]
      <<std::setw(nWidth)<<vecdT[i]
      <<std::setw(nWidth)<<vecdP[i]
      <<std::setw(nWidth)<<dU[i][0][0]<<std::endl;
  }
  std::cout<<std::setw(3)<<vecdP.size()
    <<std::setw(nWidth)<<vecdM[vecdP.size()]
    <<std::setw(nWidth)<<"-"
    <<std::setw(nWidth)<<vecdR[vecdP.size()]
    <<std::setw(nWidth)<<"-"
    <<std::setw(nWidth)<<"-"
    <<std::setw(nWidth)<<"-"
    <<std::setw(nWidth)<<"-"
    <<std::setw(nWidth)<<dU[vecdP.size()][0][0]<<std::endl;
}
void writeModelToScreen_GL(){
  int nWidth=nPrecision+9;
  std::cout<<std::setw(3)<<"#"
    <<std::setw(nWidth)<<"M_r"
    <<std::setw(nWidth)<<"Del_M"
    <<std::setw(nWidth)<<"r"
    <<std::setw(nWidth)<<"rho"
    <<std::setw(nWidth)<<"E"
    <<std::setw(nWidth)<<"P"
    <<std::setw(nWidth)<<"u"<<std::endl;
  std::cout.setf(std::ios::scientific);
  std::cout.precision(nPrecision);
  for(unsigned int i=0;i<vecdP.size();i++){
    std::cout<<std::setw(3)<<i
      <<std::setw(nWidth)<<vecdM[i]
      <<std::setw(nWidth)<<vecdMDel[i]
      <<std::setw(nWidth)<<vecdR[i]
      <<std::setw(nWidth)<<vecdRho[i]
      <<std::setw(nWidth)<<vecdE[i]
      <<std::setw(nWidth)<<vecdP[i]
      <<std::setw(nWidth)<<dU[i][0][0]<<std::endl;
  }
  std::cout<<std::setw(3)<<vecdP.size()
    <<std::setw(nWidth)<<vecdM[vecdP.size()]
    <<std::setw(nWidth)<<"-"
    <<std::setw(nWidth)<<vecdR[vecdP.size()]
    <<std::setw(nWidth)<<"-"
    <<std::setw(nWidth)<<"-"
    <<std::setw(nWidth)<<"-"
    <<std::setw(nWidth)<<dU[vecdP.size()][0][0]<<std::endl;
}
double interpolateU(double dIntVar){
  
  //find bracketing points
  unsigned int i=0;
  while(dUProR[i]<dIntVar){//bracketing values will be at i-1 and i
    i++;
    if(i>=nNumUProPoints){
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
        <<": the value for dIntVar="<<dIntVar<<" is out side the velocity profile,\n"
        <<"  try a different velocity profile\n";
      throw exception2(ssTemp.str(),CALCULATION);
    }
  }
  if(i==0){
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
        <<": the value for dIntVar="<<dIntVar<<" is out side the velocity profile,\n"
        <<"  perhaps try a larger \"R-stop\", or a different velocity profile\n";
      throw exception2(ssTemp.str(),CALCULATION);
  }
  //do linear interpolation
  return (dUPro[i]-dUPro[i-1])/(dUProR[i]-dUProR[i-1])*(dIntVar-dUProR[i-1])+dUPro[i-1];
}
void makeVelocityDist(){
  
  //RADIAL VELOCITY, U
  
  //allocate space
  dU=new double**[vecdR.size()];
  dU0=new double[vecdR.size()];
  
  if(sUDistType=="POLY"){
    
    //set radial velocity
    for(int i=vecdR.size()-1;i>=0;i--){//should go same way as R
      
      //calculate new velocity at new radius
      double dVelocity=0.0;
      for(unsigned int n=0;n<vectVelDist.size();n++){
        double dRFrac=(vecdR[i]-vecdR[nNumGhostCells])/(vecdR[vecdR.size()
          -nNumGhostCells-1]-vecdR[nNumGhostCells]);
        if(dRFrac<0.0){//set to zero when at inner ghost cells
          dRFrac=0.0;
        }
        dVelocity+=vectVelDist[n].dCoeff*pow(dRFrac,vectVelDist[n].dPower);
      }
      
      //allocate more memory and assign values
      dU[i]=new double* [nNumTheta+2*nNumGhostCells];
      for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
        dU[i][j]=new double[nNumPhi+2*nNumGhostCells];
        for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
          dU[i][j][k]=dVelocity;
        }
      }
      dU0[i]=dVelocity;
    }
    
    //set theta velocity
    dV=new double**[vecdRho.size()];
    int nInt=1;
    if(nPeriodic[1]==1){
      nInt=0;
    }
    for(unsigned int i=0;i<vecdRho.size();i++){
      dV[i]=new double*[nNumTheta+2*nNumGhostCells+nInt];
      for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells+nInt;j++){
        dV[i][j]=new double[nNumTheta+2*nNumGhostCells];
        for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
          dV[i][j][k]=0.0;
        }
      }
    }
    
    //set phi velocity
    dW=new double**[vecdRho.size()];
    nInt=1;
    if(nPeriodic[2]==1){
      nInt=0;
    }
    for(unsigned int i=0;i<vecdRho.size();i++){
      dW[i]=new double*[nNumTheta+2*nNumGhostCells];
      for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
        dW[i][j]=new double[nNumTheta+2*nNumGhostCells+nInt];
        for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells+nInt;k++){
          dW[i][j][k]=0.0;
        }
      }
    }
  }
  else if (sUDistType=="PRO"){
    
    //set radial velocity
    for(unsigned int i=0;i<vecdR.size()-3;i++){//should go same way as R
      
      double dIntVar=vecdR[i]/vecdR[0];
      double dVelocity=dUSurf*interpolateU(dIntVar);
      
      //allocate more memory and assign values
      dU[i]=new double* [nNumTheta+2*nNumGhostCells];
      for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
        dU[i][j]=new double[nNumPhi+2*nNumGhostCells];
        for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
          dU[i][j][k]=dVelocity;
        }
      }
      dU0[i]=dVelocity;
    }
    for(unsigned int i=vecdR.size()-3;i<vecdR.size();i++){//should go same way as R
      //allocate more memory and assign values
      dU[i]=new double* [nNumTheta+2*nNumGhostCells];
      for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
        dU[i][j]=new double[nNumPhi+2*nNumGhostCells];
        for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
          dU[i][j][k]=0.0;
        }
      }
      dU0[i]=0.0;
    }
    
    //set theta velocity
    dV=new double**[vecdRho.size()];
    int nInt=1;
    if(nPeriodic[1]==1){
      nInt=0;
    }
    for(unsigned int i=0;i<vecdRho.size();i++){
      dV[i]=new double*[nNumTheta+2*nNumGhostCells+nInt];
      for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells+nInt;j++){
        dV[i][j]=new double[nNumTheta+2*nNumGhostCells];
        for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
          dV[i][j][k]=0.0;
        }
      }
    }
    
    //set phi velocity
    dW=new double**[vecdRho.size()];
    nInt=1;
    if(nPeriodic[2]==1){
      nInt=0;
    }
    for(unsigned int i=0;i<vecdRho.size();i++){
      dW[i]=new double*[nNumTheta+2*nNumGhostCells];
      for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
        dW[i][j]=new double[nNumTheta+2*nNumGhostCells+nInt];
        for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells+nInt;k++){
          dW[i][j][k]=0.0;
        }
      }
    }
  }
}
void pretubeVelocityTorus(XMLNode xPerturb){
  
  if(nNumDims!=3){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": torus perturbation can only be applied to 3D models, and the current model is a "
      <<nNumDims<<"D model\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //GET TORUS PERTURBATION INFO
  
  //offset of torus from center of grid location of center of torus
  double dR_torus;
  getXMLValue(xPerturb,"r_cen_off",0,dR_torus);  
  double dTheta_torus;
  getXMLValue(xPerturb,"theta_cen_off",0,dTheta_torus);  
  double dPhi_torus;
  getXMLValue(xPerturb,"phi_cen_off",0,dPhi_torus);
  
  //get inner and outer radii of torus
  double dCen_torus_radius;
  getXMLValue(xPerturb,"radius_cen",0,dCen_torus_radius);
  double dOutter_torus_radius;
  getXMLValue(xPerturb,"radius_outter",0,dOutter_torus_radius);
  
  //get width of guassian and amplitude of preturbation
  double dWidthGuassian;
  getXMLValue(xPerturb,"width_guassian",0,dWidthGuassian);
  double dAmplitude;
  getXMLValue(xPerturb,"amplitude",0,dAmplitude);
  
  
  //PERTURB VELOCITY
  
  //set starting theta of model
  double dStartTheta=(90.0-dDeltaTheta*double((nNumTheta+2*nNumGhostCells)/2))*dPi/180.0;  
  unsigned int nTest=2*int(nNumTheta/2.0);
  if(nTest!=nNumTheta){//if an uneven number of theta
    dStartTheta-=dDeltaTheta*dPi/360.0;
  }
  
  //set start of phi's
  double dStartPhi=(0.0-dDeltaPhi*double((nNumPhi+2*nNumGhostCells)/2))*dPi/180.0;
  if(2*(nNumPhi/2)!=nNumPhi){
    dStartPhi-=dDeltaPhi*dPi/360.0;
  }
  
  //set a couple conversion factors
  double dRadPerDegree=dPi/180.0;
  double dHWHMconversion=1.0/sqrt(2.0*log(2.0));
  
  //set torus absolute position
  double dR_t=dR_torus;
  double dTheta_t=(90.0)+dTheta_torus;
  double dPhi_t=(0.0)+dPhi_torus;
  double dX_t=dR_t*sin(dTheta_t*dRadPerDegree)*cos(dPhi_t*dRadPerDegree);
  double dY_t=dR_t*sin(dTheta_t*dRadPerDegree)*sin(dPhi_t*dRadPerDegree);
  double dZ_t=dR_t*cos(dTheta_t*dRadPerDegree);
  
  //RADIAL VELOCITY, U
  for(unsigned int i=0;i<vecdR.size();i++){//should go same way as R
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      double dThetaCenter=(2.0*dStartTheta+(double(2*j+1)*dDeltaTheta*dRadPerDegree))*0.5;
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        double dPhiCenter=(2.0*dStartPhi+double(2*k+1)*dDeltaPhi*dRadPerDegree)*0.5;
        
        //get location x,y,z of current velocity component
        double dX=vecdR[i]*sin(dThetaCenter)*cos(dPhiCenter);
        double dY=vecdR[i]*sin(dThetaCenter)*sin(dPhiCenter);
        double dZ=vecdR[i]*cos(dThetaCenter);
        
        //get torus angle1
        double dAngle1=-10.0;
        if(dY>dY_t&&dZ>dZ_t){//(1)
          dAngle1=atan( (dZ-dZ_t)/(dY-dY_t) );
        }
        else if(dY<dY_t&&dZ>dZ_t){//(2)
          dAngle1=atan( (dY_t-dY)/(dZ-dZ_t) )+dPi*0.5;
        }
        else if(dY<dY_t&&dZ<dZ_t){//(3)
          dAngle1=atan( (dZ_t-dZ)/(dY_t-dY) )+dPi;
        }
        else if(dY>dY_t&&dZ<dZ_t){//(4)
          dAngle1=atan( (dY-dY_t)/(dZ_t-dZ) )+1.5*dPi;
        }
        else if(dY>dY_t&&dZ==dZ_t){//(5)
          dAngle1=0.0;
        }
        else if(dY==dY_t&&dZ>dZ_t){//(6)
          dAngle1=dPi*0.5;
        }
        else if(dY<dY_t&&dZ==dZ_t){//(7)
          dAngle1=dPi;
        }
        else if(dY==dY_t&&dZ<dZ_t){//(8)
          dAngle1=1.5*dPi;
        }
        
        //get torus angle2
        double dX_tt=dX_t;
        double dY_tt=dY_t+dCen_torus_radius*cos(dAngle1);
        double dZ_tt=dZ_t+dCen_torus_radius*sin(dAngle1);
        double dB=sqrt( pow( (dY_tt-dY),2) + pow( (dZ_tt-dZ),2));
        double dR=sqrt( pow( (dY_t-dY),2) + pow( (dZ_t-dZ),2));
        double dAngle2=-10.0;
        if(dR>dCen_torus_radius&&dX>dX_tt){//1
          dAngle2=atan( (dX-dX_tt)/dB );
        }
        else if(dR<dCen_torus_radius&&dX>dX_tt){//2
          dAngle2=atan( dB/(dX-dX_tt) )+dPi*0.5;
        }
        else if(dR<dCen_torus_radius&&dX<dX_tt){//3
          dAngle2=atan( (dX_tt-dX)/dB )+dPi;
        }
        else if(dR>dCen_torus_radius&&dX<dX_tt){//4
          dAngle2=atan( dB/(dX_tt-dX) )+dPi*1.5;
        }
        else if(dR>dCen_torus_radius&&dX==dX_tt){//5
          dAngle2=0.0;
        }
        else if(dR==dCen_torus_radius&&dX>dX_tt){//6
          dAngle2=dPi*0.5;
        }
        else if(dR<dCen_torus_radius&&dX==dX_tt){//7
          dAngle2=dPi;
        }
        else if(dR==dCen_torus_radius&&dX<dX_tt){//7
          dAngle2=dPi*1.5;
        }
        
        //get closest point on surface of torus
        double dY_st=(dCen_torus_radius+dOutter_torus_radius*cos(dAngle2))*cos(dAngle1)+dY_t;
        double dZ_st=(dCen_torus_radius+dOutter_torus_radius*cos(dAngle2))*sin(dAngle1)+dZ_t;
        double dX_st=dOutter_torus_radius*sin(dAngle2)+dX_t;
        
        //get distance from torus using parametric equation for the torus
        double dD_sq=pow((dX-dX_st),2)+pow((dY-dY_st),2)+pow((dZ-dZ_st),2);
        
        //calculate velocity perturbation from gaussian
        double dC_sq=pow(dWidthGuassian*dHWHMconversion,2);
        double dVAmplitude=dAmplitude*exp(-1.0*dD_sq/(2.0*dC_sq));
        double dTheta_prime=atan(dR/dX);
        double dBeta=dPi*0.5-dTheta_prime-dAngle2;
        double dV_r=-1.0*dVAmplitude*sin(dBeta);
        //double dV_r=-1.0*1e5*sin(dBeta);
        
        //preturb the velocity
        dU[i][j][k]=dU[i][j][k]+dV_r;
      }
    }
  }
    
  //THETA VELOCITY, V
  int nInt=1;
  if(nPeriodic[1]==1){
    nInt=0;
  }
  for(unsigned int i=0;i<vecdRho.size();i++){
    double dRCen=(vecdR[i]+vecdR[i+1])*0.5;
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells+nInt;j++){
      double dTheta=(dStartTheta+double(j+(1-nInt))*dDeltaTheta*dRadPerDegree);
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        double dPhiCenter=(2.0*dStartPhi+double(2*k+1)*dDeltaPhi*dRadPerDegree)*0.5;
        //get location x,y,z
        double dX=dRCen*sin(dTheta)*cos(dPhiCenter);
        double dY=dRCen*sin(dTheta)*sin(dPhiCenter);
        double dZ=dRCen*cos(dTheta);
        
        //get torus angle1
        double dAngle1=-10.0;
        if(dY>dY_t&&dZ>dZ_t){//(1)
          dAngle1=atan( (dZ-dZ_t)/(dY-dY_t) );
        }
        else if(dY<dY_t&&dZ>dZ_t){//(2)
          dAngle1=atan( (dY_t-dY)/(dZ-dZ_t) )+dPi*0.5;
        }
        else if(dY<dY_t&&dZ<dZ_t){//(3)
          dAngle1=atan( (dZ_t-dZ)/(dY_t-dY) )+dPi;
        }
        else if(dY>dY_t&&dZ<dZ_t){//(4)
          dAngle1=atan( (dY-dY_t)/(dZ_t-dZ) )+1.5*dPi;
        }
        else if(dY>dY_t&&dZ==dZ_t){//(5)
          dAngle1=0.0;
        }
        else if(dY==dY_t&&dZ>dZ_t){//(6)
          dAngle1=dPi*0.5;
        }
        else if(dY<dY_t&&dZ==dZ_t){//(7)
          dAngle1=dPi;
        }
        else if(dY==dY_t&&dZ<dZ_t){//(8)
          dAngle1=1.5*dPi;
        }
        
        //get torus angle2
        double dX_tt=dX_t;
        double dY_tt=dY_t+dCen_torus_radius*cos(dAngle1);
        double dZ_tt=dZ_t+dCen_torus_radius*sin(dAngle1);
        double dB=sqrt( pow( (dY_tt-dY),2) + pow( (dZ_tt-dZ),2));
        double dR=sqrt( pow( (dY_t-dY),2) + pow( (dZ_t-dZ),2));
        double dAngle2=-10.0;
        if(dR>dCen_torus_radius&&dX>dX_tt){//1
          dAngle2=atan( (dX-dX_tt)/dB );
        }
        else if(dR<dCen_torus_radius&&dX>dX_tt){//2
          dAngle2=atan( dB/(dX-dX_tt) )+dPi*0.5;
        }
        else if(dR<dCen_torus_radius&&dX<dX_tt){//3
          dAngle2=atan( (dX_tt-dX)/dB )+dPi;
        }
        else if(dR>dCen_torus_radius&&dX<dX_tt){//4
          dAngle2=atan( dB/(dX_tt-dX) )+dPi*1.5;
        }
        else if(dR>dCen_torus_radius&&dX==dX_tt){//5
          dAngle2=0.0;
        }
        else if(dR==dCen_torus_radius&&dX>dX_tt){//6
          dAngle2=dPi*0.5;
        }
        else if(dR<dCen_torus_radius&&dX==dX_tt){//7
          dAngle2=dPi;
        }
        else if(dR==dCen_torus_radius&&dX<dX_tt){//7
          dAngle2=dPi*1.5;
        }
        
        //get closest point on surface of torus
        double dY_st=(dCen_torus_radius+dOutter_torus_radius*cos(dAngle2))*cos(dAngle1)+dY_t;
        double dZ_st=(dCen_torus_radius+dOutter_torus_radius*cos(dAngle2))*sin(dAngle1)+dZ_t;
        double dX_st=dOutter_torus_radius*sin(dAngle2)+dX_t;
        
        //get distance from torus using parametric equation for the torus
        double dD_sq=pow((dX-dX_st),2)+pow((dY-dY_st),2)+pow((dZ-dZ_st),2);
        
        //calculate velocity perturbation from gaussian
        double dC_sq=pow(dWidthGuassian*dHWHMconversion,2);
        double dVAmplitude=dAmplitude*exp(-1.0*dD_sq/(2.0*dC_sq));
        double dTheta_prime=atan(dR/dX);
        double dBeta=dPi*0.5-dTheta_prime-dAngle2;
        double dV_theta=-1.0*dVAmplitude*cos(dBeta)*sin(dAngle1);
        
        //preturb the velocity
        dV[i][j][k]+=dV_theta;
      }
    }
  }
    
  //PHI VELOCITY, W
  nInt=1;
  if(nPeriodic[2]==1){
    nInt=0;
  }
  for(unsigned int i=0;i<vecdRho.size();i++){
    double dRCen=(vecdR[i]+vecdR[i+1])*0.5;
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      double dThetaCenter=(2.0*dStartTheta+(double(2*j+1)*dDeltaTheta*dRadPerDegree))*0.5;
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells+nInt;k++){
       double dPhi=dStartPhi+double(k+(1-nInt))*dDeltaPhi*dRadPerDegree;
        
        //get location x,y,z
        double dX=dRCen*sin(dThetaCenter)*cos(dPhi);
        double dY=dRCen*sin(dThetaCenter)*sin(dPhi);
        double dZ=dRCen*cos(dThetaCenter);
        
        //get torus angle1
        double dAngle1=-10.0;
        if(dY>dY_t&&dZ>dZ_t){//(1)
          dAngle1=atan( (dZ-dZ_t)/(dY-dY_t) );
        }
        else if(dY<dY_t&&dZ>dZ_t){//(2)
          dAngle1=atan( (dY_t-dY)/(dZ-dZ_t) )+dPi*0.5;
        }
        else if(dY<dY_t&&dZ<dZ_t){//(3)
          dAngle1=atan( (dZ_t-dZ)/(dY_t-dY) )+dPi;
        }
        else if(dY>dY_t&&dZ<dZ_t){//(4)
          dAngle1=atan( (dY-dY_t)/(dZ_t-dZ) )+1.5*dPi;
        }
        else if(dY>dY_t&&dZ==dZ_t){//(5)
          dAngle1=0.0;
        }
        else if(dY==dY_t&&dZ>dZ_t){//(6)
          dAngle1=dPi*0.5;
        }
        else if(dY<dY_t&&dZ==dZ_t){//(7)
          dAngle1=dPi;
        }
        else if(dY==dY_t&&dZ<dZ_t){//(8)
          dAngle1=1.5*dPi;
        }
        
        //get torus angle2
        double dX_tt=dX_t;
        double dY_tt=dY_t+dCen_torus_radius*cos(dAngle1);
        double dZ_tt=dZ_t+dCen_torus_radius*sin(dAngle1);
        double dB=sqrt( pow( (dY_tt-dY),2) + pow( (dZ_tt-dZ),2));
        double dR=sqrt( pow( (dY_t-dY),2) + pow( (dZ_t-dZ),2));
        double dAngle2=-10.0;
        if(dR>dCen_torus_radius&&dX>dX_tt){//1
          dAngle2=atan( (dX-dX_tt)/dB );
        }
        else if(dR<dCen_torus_radius&&dX>dX_tt){//2
          dAngle2=atan( dB/(dX-dX_tt) )+dPi*0.5;
        }
        else if(dR<dCen_torus_radius&&dX<dX_tt){//3
          dAngle2=atan( (dX_tt-dX)/dB )+dPi;
        }
        else if(dR>dCen_torus_radius&&dX<dX_tt){//4
          dAngle2=atan( dB/(dX_tt-dX) )+dPi*1.5;
        }
        else if(dR>dCen_torus_radius&&dX==dX_tt){//5
          dAngle2=0.0;
        }
        else if(dR==dCen_torus_radius&&dX>dX_tt){//6
          dAngle2=dPi*0.5;
        }
        else if(dR<dCen_torus_radius&&dX==dX_tt){//7
          dAngle2=dPi;
        }
        else if(dR==dCen_torus_radius&&dX<dX_tt){//7
          dAngle2=dPi*1.5;
        }
        
        //get closest point on surface of torus
        double dY_st=(dCen_torus_radius+dOutter_torus_radius*cos(dAngle2))*cos(dAngle1)+dY_t;
        double dZ_st=(dCen_torus_radius+dOutter_torus_radius*cos(dAngle2))*sin(dAngle1)+dZ_t;
        double dX_st=dOutter_torus_radius*sin(dAngle2)+dX_t;
        
        //get distance from torus using parametric equation for the torus
        double dD_sq=pow((dX-dX_st),2)+pow((dY-dY_st),2)+pow((dZ-dZ_st),2);
        
        //calculate velocity perturbation from gaussian
        double dC_sq=pow(dWidthGuassian*dHWHMconversion,2);
        double dVAmplitude=dAmplitude*exp(-1.0*dD_sq/(2.0*dC_sq));
        double dTheta_prime=atan(dR/dX);
        double dBeta=dPi*0.5-dTheta_prime-dAngle2;
        double dV_phi=dVAmplitude*cos(dBeta)*cos(dAngle1);
        
        //preturb the velocity
        dW[i][j][k]+=dV_phi;
        
      }
    }
  }
}
void makeVelocityDist_SEDOV(){
  
  //allocate space
  dU=new double**[vecdR.size()];
  dU0=new double[vecdR.size()];
  
  //set radial velocity
  for(unsigned int i=vecdR.size()-1;i>=0;i--){//should go same way as R
    
    //calculate new velocity at new radius
    double dVelocity=0.0;
    if(i>vecdR.size()-1-nNumCellsCent&&i!=vecdR.size()-1){
      dVelocity=1.0e4;
    }
    //allocate more memory and assign values
    dU[i]=new double* [nNumTheta+2*nNumGhostCells];
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      dU[i][j]=new double[nNumPhi+2*nNumGhostCells];
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        dU[i][j][k]=dVelocity;
      }
    }
    dU0[i]=dVelocity;
  }
  
  
  
  
  //RADIAL VELOCITY, U
  
  //allocate space
  dU=new double**[vecdR.size()];
  dU0=new double[vecdR.size()];
  
  //set radial velocity
  for(unsigned int i=vecdR.size()-1;i>=0;i--){//should go same way as R
    
    //calculate new velocity at new radius
    double dVelocity=0.0;
    if(i>vecdR.size()-1-nNumCellsCent&&i!=vecdR.size()-1){
      dVelocity=1.0e4;
    }
    
    //allocate more memory and assign values
    dU[i]=new double* [nNumTheta+2*nNumGhostCells];
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      dU[i][j]=new double[nNumPhi+2*nNumGhostCells];
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        dU[i][j][k]=dVelocity;
      }
    }
    dU0[i]=dVelocity;
  }
    
  //set theta velocity
  dV=new double**[vecdRho.size()];
  int nInt=1;
  if(nPeriodic[1]==1){
    nInt=0;
  }
  for(unsigned int i=0;i<vecdRho.size();i++){
    dV[i]=new double*[nNumTheta+2*nNumGhostCells+nInt];
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells+nInt;j++){
      dV[i][j]=new double[nNumTheta+2*nNumGhostCells];
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells;k++){
        dV[i][j][k]=0.0;
      }
    }
  }
    
  //set phi velocity
  dW=new double**[vecdRho.size()];
  nInt=1;
  if(nPeriodic[2]==1){
    nInt=0;
  }
  for(unsigned int i=0;i<vecdRho.size();i++){
    dW[i]=new double*[nNumTheta+2*nNumGhostCells];
    for(unsigned int j=0;j<nNumTheta+2*nNumGhostCells;j++){
      dW[i][j]=new double[nNumTheta+2*nNumGhostCells+nInt];
      for(unsigned int k=0;k<nNumPhi+2*nNumGhostCells+nInt;k++){
        dW[i][j][k]=0.0;
      }
    }
  }
}
double dMomentumCons(double dT, double dRho,int nShell){
  double dP_nShellm1=vecdP[nShell-1];
  double dP_nShell  =eosTable.dGetPressure(dT,dRho);
  double dGravTerm=dG*vecdM[nShell]/(4.0*dPi*pow(vecdR[nShell],4));
  double dPressureTerm=(dP_nShell-dP_nShellm1)/(vecdMDel[nShell]+vecdMDel[nShell-1])*2.0;
  return dGravTerm+dPressureTerm;
}
double dEnergyCons(double dT, double dRho,int nShell){
  double dKappa_nShellm1=vecdKappa[nShell-1];
  double dKappa_nShell=eosTable.dGetOpacity(dT,dRho);
  double dT_nShellP4=pow(dT,4);
  double dT_nShellm1P4=pow(vecdT[nShell-1],4);
  double dKappa_nShellm1half=(dT_nShellP4/dKappa_nShell+dT_nShellm1P4/dKappa_nShellm1)
    /(dT_nShellP4+dT_nShellm1P4);
  double dFluxTerm=-64.0*dPi*dPi*dSigma*pow(vecdR[nShell],4)/3.0*dKappa_nShellm1half
    *(dT_nShellP4-dT_nShellm1P4)/(vecdMDel[nShell]+vecdMDel[nShell-1])*2.0;
  return dFluxTerm/dLSun-dL;
}
double dPartialDerivativeVar1(double dVar1, double dVar2,int nShell,
  double (*dFunction)(double dVar1, double dVar2, int nShell),double dH, double &dError){
  return (dFunction(dVar1+dH,dVar2,nShell)-dFunction(dVar1-dH,dVar2,nShell))/(2.0*dH);
}
double dPartialDerivativeVar2(double dVar1, double dVar2,int nShell,
  double (*dFunction)(double dVar1, double dVar2, int nShell),double dH, double &dError){
  return (dFunction(dVar1,dVar2+dH,nShell)-dFunction(dVar1,dVar2-dH,nShell))/(2.0*dH);
}
void readEnergyProfile_GL(std::string sProfileFileName){
  
  //attempt to open the file
  std::ifstream ifIn;
  ifIn.open(sProfileFileName.c_str());
  if(!ifIn.good()){//file not ready for reading
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": the file \""<<sProfileFileName<<"\" didn't open properly.\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //get number of points
  ifIn>>nNumEProPoints;
  dEPro =new double[nNumEProPoints];
  dEProM=new double[nNumEProPoints];
  for(int i=0;i<nNumEProPoints;i++){
    ifIn>>dEProM[i];
    ifIn>>dEPro[i];
  }
  ifIn.close();
}
void deleteEnergyProfile_GL(){
  delete [] dEPro;
  delete [] dEProM;
}
void cleanUp(){
  deleteEnergyProfile_GL();
}
double EOS_GL(double dP, double dE,int i){
  if(i>53&&i<59){
    return dP/(1.2-1.0)/dE;
  }
  else{
    return dP/(dGamma-1.0)/dE;
  }
}
double EOS_GL(double dP, double dE){
  return dP/(dGamma-1.0)/dE;
}
double interpolateE_GL(double dIntVar){
  
  //find bracketing points
  int i=1;
  while(dEProM[i]>dIntVar){//bracketing values will be at i-1 and i
    i++;
    if(i>=nNumEProPoints){
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
        <<": the value for dIntVar="<<dIntVar
        <<" is out side the energy profile, perhaps try a larger \"M-delta-init\"\n";
      throw exception2(ssTemp.str(),CALCULATION);
    }
  }
  //do linear interpolation
  return (dEPro[i]-dEPro[i-1])/(dEProM[i]-dEProM[i-1])*(dIntVar-dEProM[i-1])+dEPro[i-1];
}
