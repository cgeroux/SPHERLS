/** 
  @file
  
  This file holds functions for manipulating data. This includes initializing
  the program, parsing the configuration file "config.xml", allocating memory
  for the model to be read in, reading in the input model, etc.
*/

#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <fenv.h>//linux
#include "dataManipulation.h"
#include "global.h"
#include "xmlFunctions.h"
#include "exception2.h"
#include "dataMonitoring.h"
#include "physEquations.h"
#include <string>
#include "fileExists.h"

void init(ProcTop &procTop,Grid &grid,Output &output,Time &time,Parameters &parameters
  ,MessPass &messPass,Performance &performance,Implicit &implicit
  ,int nNumArgs,char* cArgs[]){
  
  //find out number of processes
  procTop.nNumProcs=MPI::COMM_WORLD.Get_size();
  
  //find out process rank
  procTop.nRank=MPI::COMM_WORLD.Get_rank();
  
  //start timer
  if(procTop.nRank==0){
    performance.dStartTimer=MPI::Wtime();
  }
  
  //turn on floating point exceptions
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  
  
  //SET INITIAL VALUES FROM "SPHERLS.xml"
  
  //open configuration file using the data node
  XMLNode xData=openXMLFile("SPHERLS.xml","data");
  
  
  //READ IN PROCESSOR TOPOLOGY
  
  //setup array to hold demensions of processors
  procTop.nProcDims=new int[3];
  
  //move into dimSize node
  XMLNode xProcDims=getXMLNode(xData,"procDims",0);
  
  //get size of processor dims in x0 direction
  getXMLValue(xProcDims,"x0",0,procTop.nProcDims[0]);
  
  //get size of processor dims in x1 direction
  getXMLValue(xProcDims,"x1",0,procTop.nProcDims[1]);
  
  //get size of processor dims in x2 direction
  getXMLValue(xProcDims,"x2",0,procTop.nProcDims[2]);
  
  //make sure we have the right number of prcessors for this setup
  int nTotalProcs=(procTop.nProcDims[0]-1)*procTop.nProcDims[1]*procTop.nProcDims[2]+1;
  if(nTotalProcs!=procTop.nNumProcs){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": (procTop.nProcDims[0]-1)*procTop.nProcDims[1]*procTop.nProcDims[2]+1="<<nTotalProcs
      <<", not equal to the number of processors="<<procTop.nNumProcs<<std::endl;
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //get output file name
  getXMLValue(xData,"outputName",0,output.sBaseOutputFileName);
  
  //get debug outputfile name if there is one set
  if(!getXMLValueNoThrow(xData,"debugProfileOutput",0,parameters.sDebugProfileOutput)){
    parameters.sDebugProfileOutput=output.sBaseOutputFileName+"_debug";
  };
  
  //get starting model file name
  std::string sStartModel;
  getXMLValue(xData,"startModel",0,sStartModel);
  
  //get eos node
  XMLNode xEOS=getXMLNode(xData,"eos",0);
  
  //read in file name for equation of state, to override starting model's eos file
  getXMLValueNoThrow(xEOS,"eosFile",0,parameters.sEOSFileName);
  
  //get if using the turbulance model or not
  XMLNode xTurbModel=getXMLNode(xData,"turbMod",0);
  if(!xTurbModel.isEmpty()){
    std::string sTemp;
    getXMLValue(xTurbModel,"type",0,sTemp);
    if(sTemp.compare("smagorinsky")==0){//use smagorinsky eddy viscosity
      parameters.nTypeTurbulanceMod=2;
    }
    else if(sTemp.compare("constant")==0){//using a constant eddy viscosity
      parameters.nTypeTurbulanceMod=1;
    }
    else{//unknown eddy viscosity type
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
        <<": Unknown type, \""<<sTemp
        <<"\", of turbulance model, try either \"smagorinsky\" or \"constant\"."<<std::endl;
      throw exception2(ssTemp.str(),INPUT);
    }
    getXMLValue(xTurbModel,"eddyVisc",0,parameters.dEddyViscosity);
  }
  else{//if no node found
    parameters.nTypeTurbulanceMod=0;//not using a turbulance model
  }
  
  //switch to dedm node if there is one
  XMLNode xDEDM=getXMLNodeNoThrow(xData,"dedm",0);
  if(!xDEDM.isEmpty()){
    parameters.bDEDMClamp=true;
    
    //check for a ./DEDMClamp.dat file
    bool bIsDEDMClampfile=bFileExists("./DEDMClamp.dat");
    if(!bIsDEDMClampfile){
      
      //get temperature to set DEDM clamp
      getXMLAttribute(xDEDM,"temperature",parameters.dEDMClampTemperature);
    }
    else{
      
      //open file
      std::ifstream ifDEDMClampFile;
      ifDEDMClampFile.open("./DEDMClamp.dat");
      
      //get M_r to set the clamp at
      ifDEDMClampFile>>parameters.dDEDMClampMr;
      
      //get DEDM of clamp
      ifDEDMClampFile>>parameters.dDEDMClampValue;
    }
  }
  else{
    parameters.bDEDMClamp=false;
  }
  
  //read in model
  modelRead(sStartModel,procTop,grid,time,parameters);
  
  //switch to model dump node
  XMLNode xDump=getXMLNodeNoThrow(xData,"dumps",0);
  
  //set time of last model dump equal to the current simulation time
  output.dTimeLastDump=time.dt;
  
  //get dump fequencies
  output.nDumpFrequencyStep=0;//not dumping according to number of time steps
  output.dDumpFrequencyTime=0.0;//not dumping according to simulation time
  
  if(!xDump.isEmpty()){
    output.bDump=true;
    
    //get dump frequencies
    XMLNode xFrequency1=getXMLNodeNoThrow(xDump,"frequency",0);
    if(!xFrequency1.isEmpty()){//no frequency node found
      std::string sType;
      getXMLAttribute(xFrequency1,"type",sType);
      if(sType.compare("timeSteps")==0){
        getXMLValue(xDump,"frequency",0,output.nDumpFrequencyStep);
      }
      else if(sType.compare("seconds")==0){
        getXMLValue(xDump,"frequency",0,output.dDumpFrequencyTime);
      }
      else{
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
          <<": unknown attribute, \""<<sType<<"\" in first frequency node."<<std::endl;
        throw exception2(ssTemp.str(),INPUT);
      }
      
      //get second dump frequency
      XMLNode xFrequency2=getXMLNodeNoThrow(xDump,"frequency",1);
      if(!xFrequency2.isEmpty()){
        getXMLAttribute(xFrequency2,"type",sType);
        if(sType.compare("timeSteps")==0){
          getXMLValue(xDump,"frequency",1,output.nDumpFrequencyStep);
        }
        else if(sType.compare("seconds")==0){
          getXMLValue(xDump,"frequency",1,output.dDumpFrequencyTime);
        }
        else{
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
            <<": unknown attribute, \""<<sType<<"\" in second frequency node."<<std::endl;
          throw exception2(ssTemp.str(),INPUT);
        }
      }
    }
    else{
      std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
        <<": WARNING no \"frequency\" node found under \"dump\" node, no model dumps will be made!"
        <<std::endl;
      output.bDump=false;//no dumps
    }
  }
  else{
    output.bDump=false;
  }
  
  //switch to status print node
  XMLNode xPrint=getXMLNodeNoThrow(xData,"prints",0);
  
  //set time of last model dump equal to the current simulation time
  output.dTimeLastPrint=time.dt;
  
  //get print fequencies
  output.nPrintFrequencyStep=0;//not printing according to number of time steps
  output.dPrintFrequencyTime=0.0;//not printing according to simulation time
  
  if(!xPrint.isEmpty()){
    output.bPrint=true;
    
    //get print mode
    std::string sType;
    if(getXMLAttributeNoThrow(xPrint,"type",sType)){//if attribute is set
      if(sType.compare("normal")==0){
        output.nPrintMode=0;
      }
      else if(sType.compare("timeStepInfo")==0){
        output.nPrintMode=1;
      }
      else{
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
          <<": unknown attribute, \""<<sType<<"\" in first print node under \"data\" node."
          <<std::endl;
        throw exception2(ssTemp.str(),INPUT);
      }
    }
    else{//if attribute isn't set, assume normal output
      output.nPrintMode=0;//normal print mode
      std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
        <<": WARNING: type attribute not given under \"print\" node, assuming normal print type"
        <<std::endl;
    }
    
    //get print frequencies
    XMLNode xFrequency1=getXMLNodeNoThrow(xPrint,"frequency",0);
    if(xFrequency1.isEmpty()){//no frequency node found
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
        <<": no \"frequency\" node found under \"prints\" node."
        <<std::endl;
      throw exception2(ssTemp.str(),INPUT);
    }
    getXMLAttribute(xFrequency1,"type",sType);
    if(sType.compare("timeSteps")==0){
      getXMLValue(xPrint,"frequency",0,output.nPrintFrequencyStep);
    }
    else if(sType.compare("seconds")==0){
      getXMLValue(xPrint,"frequency",0,output.dPrintFrequencyTime);
    }
    else{
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
        <<": unknown attribute, \""<<sType<<"\" in first frequency node under \"prints\" node."
        <<std::endl;
      throw exception2(ssTemp.str(),INPUT);
    }
    
    //get second print frequency
    XMLNode xFrequency2=getXMLNodeNoThrow(xPrint,"frequency",1);
    if(!xFrequency2.isEmpty()){
      getXMLAttribute(xFrequency2,"type",sType);
      if(sType.compare("timeSteps")==0){
        getXMLValue(xPrint,"frequency",1,output.nPrintFrequencyStep);
      }
      else if(sType.compare("seconds")==0){
        getXMLValue(xPrint,"frequency",1,output.dPrintFrequencyTime);
      }
      else{
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
          <<": unknown attribute, \""<<sType<<"\" in second frequency node under prints node."
          <<std::endl;
        throw exception2(ssTemp.str(),INPUT);
      }
    }
  }
  else{
    output.bPrint=false;
  }
  
  //get time node
  XMLNode xTime=getXMLNode(xData,"time",0);
  
  //get end time
  bool bGotEnd=true;
  if(!getXMLValueNoThrow(xTime,"endTime",0,time.dEndTime)){
    time.dEndTime=std::numeric_limits<double>::max();
    bGotEnd=false;
  }
  
  //get end timestep
  if(!getXMLValueNoThrow(xTime,"endTimeStep",0,time.nEndTimeStep)){
    time.nEndTimeStep=std::numeric_limits<int>::max();/*set large so it won't be triggered, if it is
      still triggered at this large value we will so have problems anyhow and need to stop.*/
    if(!bGotEnd){//should have some way of knowing when to stop
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
        <<": no \"endTime\" node or \"endTimeStep\" node found under \"time\" node. Must specify at "
        <<"least one of these (need to know when to stop)."<<std::endl;
      throw exception2(ssTemp.str(),INPUT);
    }
  }
  
  //get time factor
  bool bNoTimeStepFactor=false;
  if(!getXMLValueNoThrow(xTime,"timeStepFactor",0,time.dTimeStepFactor)){
    bNoTimeStepFactor=true;
  }
  
  //get percent chnage allowed per time step
  time.dPerChange=1.0e-1;//default is 10%, probably need something an order or two smaller
  getXMLValueNoThrow(xTime,"percentChangePerTimeStep",0,time.dPerChange);
  
  //get constant time step
  if(getXMLValueNoThrow(xTime,"timeStep",0,time.dConstTimeStep)){
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": WARNING: using constant time step of "<<time.dConstTimeStep<<"s "<<std::endl;
    time.bVariableTimeStep=false;
  }
  else{
    time.dConstTimeStep=0;
    time.bVariableTimeStep=true;
  }
  
  //check that timestep parameters have been chosen
  if(time.bVariableTimeStep&&bNoTimeStepFactor){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": no \"timeStep\" node or \"timeStepFactor\" found under \"time\" node. Must specify at "
      <<"least one of these."<<std::endl;
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //get artificial viscosity paramter
  getXMLValue(xData,"av",0,parameters.dA);
  
  //get extra alpha
  parameters.dAlphaExtra=0.0;
  getXMLValueNoThrow(xData,"extraAlpha",0,parameters.dAlphaExtra);
  
  //get extra alpha
  parameters.dDonorCellMultiplier=1.0;
  getXMLValueNoThrow(xData,"donorMult",0,parameters.dDonorCellMultiplier);
  
  //get A.V. threshold value
  getXMLValue(xData,"av-threshold",0,parameters.dAVThreshold);
  
  //read in equation of state if using a tabulated equation of state
  if(!parameters.bEOSGammaLaw){
    parameters.eosTable.readBin(parameters.sEOSFileName);
    
    //get tolerance for interated quantities
    getXMLValue(xEOS,"tolerance",0,parameters.dTolerance);
    
    //get maximum number of iternations for achieving allowed tolerance
    getXMLValue(xEOS,"max-iterations",0,parameters.nMaxIterations);
  }
  
  //get if the calculation is to be adiabatic or non-adiabatic
  if(!getXMLValueNoThrow(xData,"adiabatic",0,parameters.bAdiabatic)){
    parameters.bAdiabatic=false;
  }
  
  //get implicit info
  XMLNode xImplicit=getXMLNodeNoThrow(xData,"implicit",0);
  if(!xImplicit.isEmpty()){
    
    //get number of zones at the surface to treat implicitly
    getXMLValue(xImplicit,"numImplicitZones",0,implicit.nNumImplicitZones);
    if(implicit.nNumImplicitZones>grid.nGlobalGridDims[0]-grid.nNum1DZones){
      std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
        <<": WARNING: number of implicit zones ("<<implicit.nNumImplicitZones
        <<") is larger than the total number of radial multi-D zones in the model ("
        <<grid.nGlobalGridDims[0]-grid.nNum1DZones
        <<") setting number of implicit zones equal to the number of radial zones in multi-D region.\n";
      implicit.nNumImplicitZones=grid.nGlobalGridDims[0]-grid.nNum1DZones;
    }
    if(implicit.nNumImplicitZones<0){
      std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
        <<": WARNING number of implicit zones ("<<implicit.nNumImplicitZones
        <<") is negative setting to zero.\n";
      implicit.nNumImplicitZones=0;
    }
    
    //get fraciton of temperature to use for step size in numerical derivatives
    getXMLValue(xImplicit,"derivativeStepFraction",0,implicit.dDerivativeStepFraction);
    
    //get maximum number of iterations to use for the implicit calculation
    getXMLValue(xImplicit,"max-iterations",0,implicit.nMaxNumIterations);
    
    //get tolerance to use when calculating the temperature using the implicit calculation method.
    getXMLValue(xImplicit,"tolerance",0,implicit.dTolerance);
  }
  else{
    implicit.nNumImplicitZones=0;
  }
  
  //initialize boundary updates
  initUpdateLocalBoundaries(procTop, grid, messPass,implicit);
  
  //initialize internal variables
  initInternalVars(grid,procTop,parameters);
  
  //initilize implicit calculation
  if(implicit.nNumImplicitZones>0){
    initImplicitCalculation(implicit, grid, procTop,nNumArgs,cArgs);
  }
  
  //parse, and initialize watch zones
  initWatchZones(xData, procTop,grid,output,parameters,time);
}
void setupLocalGrid(ProcTop &procTop, Grid &grid){
  
  //set coordinates for all processors
  int nRankCur=1;
  procTop.nCoords=new int*[procTop.nNumProcs];
  for(int i=1;i<procTop.nProcDims[0];i++){
    for(int j=0;j<procTop.nProcDims[1];j++){
      for(int k=0;k<procTop.nProcDims[2];k++){
        procTop.nCoords[nRankCur]=new int[3];
        procTop.nCoords[nRankCur][0]=i;
        procTop.nCoords[nRankCur][1]=j;
        procTop.nCoords[nRankCur][2]=k;
        nRankCur++;
      }
    }
  }
  procTop.nCoords[0]=new int[3];
  procTop.nCoords[0][0]=0;
  procTop.nCoords[0][1]=-1;//matches all y
  procTop.nCoords[0][2]=-1;//matches all z
  
  //calculate grid sizes for all processors
  grid.nLocalGridDims=new int**[procTop.nNumProcs];
  
  //for processors other than 0
  for(int p=1;p<procTop.nNumProcs;p++){
    grid.nLocalGridDims[p]=new int*[grid.nNumVars+grid.nNumIntVars];
    
    //calculate grid size and remainder on all processors in 3D region
    for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
      grid.nLocalGridDims[p][n]=new int[3];
      int *nRemainder=new int[3];
      for(int l=0;l<3;l++){
        
        if(l==0){/*remove processor 0 from radial grid size, and remove number of 1D zones form 
          global dimension*/
          grid.nLocalGridDims[p][n][l]=int((grid.nGlobalGridDims[l]-grid.nNum1DZones)
            /(procTop.nProcDims[l]-1));
          nRemainder[l]=(grid.nGlobalGridDims[l]-grid.nNum1DZones)%(procTop.nProcDims[l]-1);
        }
        else{
          grid.nLocalGridDims[p][n][l]=int((grid.nGlobalGridDims[l])/procTop.nProcDims[l]);
          nRemainder[l]=(grid.nGlobalGridDims[l])%procTop.nProcDims[l];
        }
        if(grid.nVariables[n][l]==-1){//if not defined in direction l
          grid.nLocalGridDims[p][n][l]=1;
          nRemainder[l]=0;
        }
        //if interface quantity, first processor in direction l, and not periodic
        if(grid.nVariables[n][l]==1&&procTop.nCoords[p][l]==0&&procTop.nPeriodic[l]==0){
          grid.nLocalGridDims[p][n][l]++;//add space for innner interface
          
        }
      }
    
      //allocate remainder of grid evenly to other processors
      for(int l=0;l<3;l++){
        if(procTop.nCoords[p][l]>=procTop.nProcDims[l]-nRemainder[l]){
          grid.nLocalGridDims[p][n][l]++;
        }
      }
      delete [] nRemainder;
    }
  }
  
  //for processsor 0
  grid.nLocalGridDims[0]=new int*[grid.nNumVars+grid.nNumIntVars];
  for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
    grid.nLocalGridDims[0][n]=new int[3];
    
    grid.nLocalGridDims[0][n][0]=grid.nNum1DZones;
    grid.nLocalGridDims[0][n][1]=1;//only one theta
    grid.nLocalGridDims[0][n][2]=1;//only one phi
    for(int l=0;l<3;l++){
      //if interface quantity, first processor in direction l, and not periodic
      if(grid.nVariables[n][l]==1&&(procTop.nCoords[0][l]==-1||procTop.nCoords[0][l]==0)
        &&procTop.nPeriodic[l]==0){
        grid.nLocalGridDims[0][n][l]++;//add space for innner interface
      }
      if(grid.nVariables[n][l]==-1){
        grid.nLocalGridDims[0][n][l]=1;
      }
    }
    if(grid.nVariables[n][0]==-1){
      grid.nLocalGridDims[0][n][0]=0;
      grid.nLocalGridDims[0][n][1]=0;
      grid.nLocalGridDims[0][n][2]=0;/*if variable not defined in radial direciton, processor 0 
        dosn't care about it*/
    }
  }
  
  //allocate memory for local grids
  if(procTop.nRank==0){// 1D region doesn't need ghost cells in theta and phi directions
    grid.dLocalGridOld=new double***[grid.nNumVars+grid.nNumIntVars];
    grid.dLocalGridNew=new double***[grid.nNumVars+grid.nNumIntVars];
    for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
      
      //allocate radial grid memory for old and new grid
      int nGhostCellsX=1;
      if(grid.nVariables[n][0]==-1){
        nGhostCellsX=0;
      }
      grid.dLocalGridOld[n]=new double**[grid.nLocalGridDims[procTop.nRank][n][0]
        +2*nGhostCellsX*grid.nNumGhostCells];
      grid.dLocalGridNew[n]=new double**[grid.nLocalGridDims[procTop.nRank][n][0]
        +2*nGhostCellsX*grid.nNumGhostCells];
      
      //do inner part of old and new grid
      for(int i=0;i<grid.nLocalGridDims[procTop.nRank][n][0]+nGhostCellsX*grid.nNumGhostCells;i++){
        grid.dLocalGridOld[n][i]=new double*[grid.nLocalGridDims[procTop.nRank][n][1]];
        grid.dLocalGridNew[n][i]=new double*[grid.nLocalGridDims[procTop.nRank][n][1]];
        for(int j=0;j<grid.nLocalGridDims[procTop.nRank][n][1];j++){
          grid.dLocalGridOld[n][i][j]=new double[grid.nLocalGridDims[procTop.nRank][n][2]];
          grid.dLocalGridNew[n][i][j]=new double[grid.nLocalGridDims[procTop.nRank][n][2]];
        }
      }
      
      //expand out last grid.nNumGhostCells to hold data from adjacent 3D grid, to later be averaged
      int nStartX=0;
      int nEndX=0;
      int nSizeY=1;
      int nSizeZ=1;
      if(grid.nVariables[n][0]!=-1){
        nStartX=grid.nLocalGridDims[procTop.nRank][n][0]+grid.nNumGhostCells;
        nEndX=grid.nLocalGridDims[procTop.nRank][n][0]+2*grid.nNumGhostCells;
      }
      if(grid.nVariables[n][1]!=-1){
        nSizeY=grid.nGlobalGridDims[1]+grid.nVariables[n][1];
      }
      else if(grid.nVariables[n][1]==-1){
        nSizeY=procTop.nProcDims[1];//if not defined in that y-direction
                            //allow space for each neighboring processor to send data
      }
      if(grid.nVariables[n][2]!=-1){
        nSizeZ=grid.nGlobalGridDims[2]+grid.nVariables[n][2];
      }
      else if(grid.nVariables[n][2]==-1){
        nSizeZ=procTop.nProcDims[2];//if not defined in that z-direction
                            //allow space for each neighboring processor to send data
      }
      for(int i=nStartX;i<nEndX;i++){
        grid.dLocalGridOld[n][i]=new double*[nSizeY];
        grid.dLocalGridNew[n][i]=new double*[nSizeY];
        for(int j=0;j<nSizeY;j++){
          grid.dLocalGridOld[n][i][j]=new double[nSizeZ];
          grid.dLocalGridNew[n][i][j]=new double[nSizeZ];
        }
      }
    }
  }
  else{// 3D region
    grid.dLocalGridOld=new double***[grid.nNumVars+grid.nNumIntVars];
    grid.dLocalGridNew=new double***[grid.nNumVars+grid.nNumIntVars];
    for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
      int nSizeX=1;
      int nSizeY=1;
      int nSizeZ=1;
      if(grid.nVariables[n][0]!=-1){
        nSizeX=grid.nLocalGridDims[procTop.nRank][n][0]+2*grid.nNumGhostCells;
      }
      if(grid.nVariables[n][1]!=-1){
        if(grid.nNumDims>1){//only need ghost cells if greater than 1D
          nSizeY=grid.nLocalGridDims[procTop.nRank][n][1]+2*grid.nNumGhostCells;
        }
        else{
          nSizeY=grid.nLocalGridDims[procTop.nRank][n][1];
        }
      }
      if(grid.nVariables[n][2]!=-1){
        if(grid.nNumDims>2){
          nSizeZ=grid.nLocalGridDims[procTop.nRank][n][2]+2*grid.nNumGhostCells;
        }
        else{
          nSizeZ=grid.nLocalGridDims[procTop.nRank][n][2];
        }
      }
      grid.dLocalGridOld[n]=new double**[nSizeX];
      grid.dLocalGridNew[n]=new double**[nSizeX];
      for(int i=0;i<nSizeX;i++){
        grid.dLocalGridOld[n][i]=new double*[nSizeY];
        grid.dLocalGridNew[n][i]=new double*[nSizeY];
        for(int j=0;j<nSizeY;j++){
          grid.dLocalGridOld[n][i][j]=new double[nSizeZ];
          grid.dLocalGridNew[n][i][j]=new double[nSizeZ];
        }
      }
    }
  }
  
  //set offset for interface centered quantities
  grid.nCenIntOffset=new int[3];
  for(int l=0;l<3;l++){
    if(procTop.nCoords[procTop.nRank][l]==0&&procTop.nPeriodic[l]==0){
      grid.nCenIntOffset[l]=1;
    }
    else{
      grid.nCenIntOffset[l]=0;
    }
  }
  //SET LOCAL GRID START POSITIONS
  //x0-direction
  grid.nGlobalGridPositionLocalGrid[0]=0;
  for(int p=0;p<procTop.nNumProcs;p++){
    if( (procTop.nCoords[p][0]<procTop.nCoords[procTop.nRank][0])
      &&(procTop.nCoords[p][1]==procTop.nCoords[procTop.nRank][1]||procTop.nCoords[p][1]==-1)
      &&(procTop.nCoords[p][2]==procTop.nCoords[procTop.nRank][2]||procTop.nCoords[p][2]==-1)){
      grid.nGlobalGridPositionLocalGrid[0]+=grid.nLocalGridDims[p][grid.nD][0];
    }
  }
  if(grid.nGlobalGridPositionLocalGrid[0]!=0){//if not first proc in line, add ghost cells of inner boundary
    grid.nGlobalGridPositionLocalGrid[0]+=grid.nNumGhostCells;
  }
  
  //x1-direction
  grid.nGlobalGridPositionLocalGrid[1]=0;
  for(int p=0;p<procTop.nNumProcs;p++){
    if( (procTop.nCoords[p][0]==procTop.nCoords[procTop.nRank][0]||procTop.nCoords[p][0]==-1)
      &&(procTop.nCoords[p][1]<procTop.nCoords[procTop.nRank][1])
      &&(procTop.nCoords[p][2]==procTop.nCoords[procTop.nRank][2]||procTop.nCoords[p][2]==-1)){
      grid.nGlobalGridPositionLocalGrid[1]+=grid.nLocalGridDims[p][grid.nD][1];
    }
  }
  if(grid.nGlobalGridPositionLocalGrid[1]!=0){//if not first proc in line, add ghost cells of inner boundary
    grid.nGlobalGridPositionLocalGrid[1]+=grid.nNumGhostCells;
  }
  
  //x2-direction
  grid.nGlobalGridPositionLocalGrid[2]=0;
  for(int p=0;p<procTop.nNumProcs;p++){
    if( (procTop.nCoords[p][0]==procTop.nCoords[procTop.nRank][0]||procTop.nCoords[p][0]==-1)
      &&(procTop.nCoords[p][1]==procTop.nCoords[procTop.nRank][1]||procTop.nCoords[p][0]==-1)
      &&(procTop.nCoords[p][2]<procTop.nCoords[procTop.nRank][2])){
      grid.nGlobalGridPositionLocalGrid[2]+=grid.nLocalGridDims[p][grid.nD][2];
    }
  }
  if(grid.nGlobalGridPositionLocalGrid[2]!=0){//if not first proc in line, add ghost cells of inner boundary
    grid.nGlobalGridPositionLocalGrid[2]+=grid.nNumGhostCells;
  }
}
void fin(bool bWriteCurrentStateToFile, Time &time, Output &output,ProcTop
  &procTop,Grid& grid,Parameters &parameters,Functions &functions
  ,Performance& performance,Implicit& implicit){
  
  
  //wait for all processors to finish before quiting
  MPI::COMM_WORLD.Barrier();
  
  if(bWriteCurrentStateToFile){
    
    //write out last model
    std::stringstream ssNumber;
    ssNumber.precision(10);
    ssNumber.unsetf(std::ios::fixed);
    ssNumber.setf(std::ios::scientific);
    ssNumber<<time.dt;
    ssNumber.str(ssNumber.str().erase(1,1));
    std::stringstream ssFileNameOut;
    ssFileNameOut<<output.sBaseOutputFileName<<"_t"<<std::setfill('0')<<std::setw(8)
      <<time.nTimeStepIndex;
          if(procTop.nRank==0){
            std::cout.setf(std::ios::scientific);
            std::cout.precision(14);
            if(output.nPrintMode==0){
              std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank<<":"
                <<std::endl<<"  Dumping model to file: "<<ssFileNameOut.str()<<std::endl;
              std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank<<":"
                <<std::endl
                <<"  For current time step "<<time.nTimeStepIndex<<" at time t="
                <<time.dt<<" [s] with time step dt=";
              std::cout.precision(5);
              std::cout<<time.dDeltat_np1half<<" [s]"<<std::endl;
              std::cout
                <<"    * max( (Del Rho)_t/Rho )   ="<<time.dDelRho_t_Rho_max<<std::endl
                <<"    * max( (Del T)_t/T )       ="<<time.dDelT_t_T_max<<std::endl
                <<"    * max( (Del UmU0)_t/UmU0 ) ="<<time.dDelUmU0_t_UmU0_max<<std::endl;
              if(grid.nNumDims>1){//if 2D or more
                std::cout<<"    * max( (Del V)_t/V )       ="<<time.dDelV_t_V_max<<std::endl;
              }
              if(grid.nNumDims>2){//if 3D or more
                std::cout<<"    * max( (Del W)_t/W )       ="<<time.dDelW_t_W_max<<std::endl;
              }
              if(grid.nNumDims>1){
                std::cout<<"    * max(convective velocity) ="
                  <<parameters.dMaxConvectiveVelocity/1.0e5<<" [km/s]"<<std::endl;
              }
              if(implicit.nNumImplicitZones>0){
                std::cout
                  <<"  : Over last "<<output.nNumTimeStepsSinceLastPrint<<" time steps :\n"
                  <<"    * Largest number of iterations for temperature calculation was "
                  <<implicit.nCurrentNumIterations<<std::endl
                  <<"    * Largest relative error in calculation of temperature was     "
                  <<implicit.dCurrentRelTError<<std::endl;
                #if TRACKMAXSOLVERERROR==1
                  std::cout
                    <<"    * Largest number of iterations for solver was                  "
                    <<implicit.nMaxNumSolverIterations<<std::endl
                    <<"    * Largest absolute error in RHS of system was                  "
                    <<implicit.dMaxErrorInRHS<<std::endl
                    <<"    * with average RHS magnitude of                                "
                    <<implicit.dAverageRHS<<std::endl;
                #endif
              }
            }
            else if(output.nPrintMode==1){//print out time step diagnostic info
              std::cout<<time.nTimeStepIndex
                <<" "<<time.dt
                <<" "<<time.dDeltat_np1half
                <<" "<<time.dDelRho_t_Rho_max
                <<" "<<time.dDelT_t_T_max
                <<" "<<time.dDelUmU0_t_UmU0_max
                <<" "<<time.dDelV_t_V_max<<std::endl;
            }
          }
    functions.fpModelWrite(ssFileNameOut.str(),procTop,grid,time,parameters);
    
    #if DEBUG_EQUATIONS==1
    std::stringstream ssFileNameProOut;
    ssFileNameProOut<<parameters.sDebugProfileOutput<<"_t"<<std::setfill('0')
      <<std::setw(8)<<time.nTimeStepIndex<<"_pro.txt";
    parameters.profileDataDebug.toFile(ssFileNameProOut.str(),time,procTop);
    parameters.profileDataDebug.clear();
    #endif
  }
  
  //finish other tasks
  finWatchZones(output);
  
  //report on performance
  if(procTop.nRank==0){
    
    //get end time
    performance.dEndTimer=MPI::Wtime();
    
    //set floating point format
    std::cout.precision(10);
    std::cout.unsetf(std::ios::fixed);
    std::cout.setf(std::ios::scientific);
    
    //write out run time
    std::cout<<"Run time for proc "<<procTop.nRank<<" is "
      <<(performance.dEndTimer-performance.dStartTimer)<<" [s]"<<std::endl;
  }
}
void modelWrite_GL(std::string sFileName,ProcTop &procTop, Grid &grid, Time &time
  , Parameters &parameters){
  
  //set file name to be the sFilename-procTop.nRank, where sFileName should be the same
  // for each processor
  std::ostringstream ossFileName;
  ossFileName<<sFileName<<"-"<<procTop.nRank;
  
  //open file
  std::ofstream ofOut;
  ofOut.open(ossFileName.str().c_str(),std::ios::binary);
  if(!ofOut.is_open()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": error opening the file "<<ossFileName.str().c_str()<<std::endl;
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //write out file type as binary
  char cTemp='b';
  ofOut.write((char*)(&cTemp),sizeof(char));
  
  //write out file version
  int nTemp=DUMP_VERSION;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out time
  ofOut.write((char*)(&time.dt),sizeof(double));
  
  //write out time step index
  ofOut.write((char*)(&time.nTimeStepIndex),sizeof(int));
  
  //write out last time step
  ofOut.write((char*)(&time.dDeltat_nm1half),sizeof(double));
  
  //write out alpha
  ofOut.write((char*)(&parameters.dAlpha),sizeof(double));
  
  //write out using a gamma law
  nTemp=0;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out gamma
  ofOut.write((char*)(&parameters.dGamma),sizeof(double));
  
  //write out artificial viscosity
  ofOut.write((char*)(&parameters.dA),sizeof(double));
  
  //write out artificial viscosity threshold
  ofOut.write((char*)(&parameters.dAVThreshold),sizeof(double));
  
  if(procTop.nRank==0){
    
    //write out processor dimensions
    ofOut.write((char*)(procTop.nProcDims),3*sizeof(int));
    
    //write out processor coordinates
    ofOut.write((char*)(procTop.nCoords[procTop.nRank]),3*sizeof(int));
    
    //write out periodicity
    ofOut.write((char*)(procTop.nPeriodic),3*sizeof(int));
    
    //write out number of variables
    ofOut.write((char*)(&grid.nNumVars),sizeof(int));
    
    //write out variable info
    for(int i=0;i<grid.nNumVars;i++){
      ofOut.write((char*)(grid.nVariables[i]),(4)*sizeof(int));
    }
    
    //write out number of 1D zones
    ofOut.write((char*)(&grid.nNum1DZones),sizeof(int));
    
    //write out global grid size
    ofOut.write((char*)(grid.nGlobalGridDims),3*sizeof(int));
    
    //write out localgrid size
    for(int i=0;i<grid.nNumVars;i++){
      ofOut.write((char*)(grid.nLocalGridDims[procTop.nRank][i]),3*sizeof(int));
    }
    
    //write number of ghostcells
    ofOut.write((char*)(&grid.nNumGhostCells),sizeof(int));
    
    //write out processor local grid
    for(int n=0;n<grid.nNumVars;n++){
      
      int nGhostCellsX=1;
      if(grid.nVariables[n][0]==-1){
        nGhostCellsX=0;
      }
      
      //write out inner grid
      for(int i=0;i<grid.nLocalGridDims[procTop.nRank][n][0]+nGhostCellsX*grid.nNumGhostCells;i++){
        for(int j=0;j<grid.nLocalGridDims[procTop.nRank][n][1];j++){
          //have to write out the multidimensional array one row at a time since
          // there is no garantee that from one row to another the memory is
          // contiguous
          ofOut.write((char*)(grid.dLocalGridOld[n][i][j])
            ,(grid.nLocalGridDims[procTop.nRank][n][2])*sizeof(double));
        }
      }
      
      //write out outter grid
      int nSizeY=grid.nGlobalGridDims[1];
      int nSizeZ=grid.nGlobalGridDims[2];
      if(grid.nVariables[n][1]==-1){
        nSizeY=procTop.nProcDims[1];//if not defined in that y-direction
                            //allow space for each neighboring processor to send data
      }
      if(grid.nVariables[n][2]==-1){
        nSizeZ=procTop.nProcDims[2];//if not defined in that z-direction
                            //allow space for each neighboring processor to send data
      }
      for(int i=grid.nLocalGridDims[procTop.nRank][n][0]+nGhostCellsX*grid.nNumGhostCells;
        i<grid.nLocalGridDims[procTop.nRank][n][0]+nGhostCellsX*2*grid.nNumGhostCells;i++){
        for(int j=0;j<nSizeY;j++){
          //have to write out the multidimensional array one row at a time since there is no 
          //garantee that from one row to another the memory is contiguous
          ofOut.write((char*)(grid.dLocalGridOld[n][i][j]),nSizeZ*sizeof(double));
        }
      }
    }
    ofOut.flush();
    ofOut.close();
  }
  else{
    
    //write out processor coordinates
    ofOut.write((char*)(procTop.nCoords[procTop.nRank]),3*sizeof(int));
    
    //write out number of variables
    ofOut.write((char*)(&grid.nNumVars),sizeof(int));
    
    //write out variable info
    for(int i=0;i<grid.nNumVars;i++){
      ofOut.write((char*)(grid.nVariables[i]),(4)*sizeof(int));
    }
    
    //write out localgrid size
    for(int i=0;i<grid.nNumVars;i++){
      ofOut.write((char*)(grid.nLocalGridDims[procTop.nRank][i]),3*sizeof(int));
    }
    
    //write number of ghostcells
    ofOut.write((char*)(&grid.nNumGhostCells),sizeof(int));
    
    //write out processor local grid
    for(int n=0;n<grid.nNumVars;n++){
      int nGhostCellsX=1;
      int nGhostCellsY=1;
      int nGhostCellsZ=1;
      if(grid.nVariables[n][0]==-1){
        nGhostCellsX=0;
      }
      if(grid.nVariables[n][1]==-1){
        nGhostCellsY=0;
      }
      if(grid.nVariables[n][2]==-1){
        nGhostCellsZ=0;
      }
      for(int i=0;i<grid.nLocalGridDims[procTop.nRank][n][0]
        +nGhostCellsX*2*grid.nNumGhostCells;i++){
        for(int j=0;j<grid.nLocalGridDims[procTop.nRank][n][1]
          +nGhostCellsY*2*grid.nNumGhostCells;j++){
          //have to write out the multidimensional array one row at a time since there is no
          //garantee that from one row to another the memory is contiguous
          ofOut.write((char*)(grid.dLocalGridOld[n][i][j])
            ,(grid.nLocalGridDims[procTop.nRank][n][2]
            +nGhostCellsZ*2*grid.nNumGhostCells)*sizeof(double));
        }
      }
    }
    ofOut.flush();
    ofOut.close();
  }
}
void modelWrite_TEOS(std::string sFileName,ProcTop &procTop, Grid &grid, Time &time
  , Parameters &parameters){
  
  //set file name to be the sFilename-procTop.nRank, where sFileName should be the same
  // for each processor
  std::ostringstream ossFileName;
  ossFileName<<sFileName<<"-"<<procTop.nRank;
  
  //open file
  std::ofstream ofOut;
  ofOut.open(ossFileName.str().c_str(),std::ios::binary);
  if(!ofOut.is_open()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": error opening the file "<<ossFileName.str().c_str()<<std::endl;
    throw exception2(ssTemp.str(),OUTPUT);
  }
  
  //write out file type
  char cTemp='b';
  ofOut.write((char*)(&cTemp),sizeof(char));
  
  //write out file version
  int nTemp=DUMP_VERSION;
  ofOut.write((char*)(&nTemp),sizeof(int));
  
  //write out time
  ofOut.write((char*)(&time.dt),sizeof(double));
  
  //write out time step index
  ofOut.write((char*)(&time.nTimeStepIndex),sizeof(int));
  
  //write out last time step
  ofOut.write((char*)(&time.dDeltat_nm1half),sizeof(double));
  
  //write out last time step
  ofOut.write((char*)(&time.dDeltat_np1half),sizeof(double));
  
  //write out alpha
  ofOut.write((char*)(&parameters.dAlpha),sizeof(double));
  
  //write out size of equation of state string
  int nEOSFileNameLen=parameters.sEOSFileName.length();
  ofOut.write((char*)(&nEOSFileNameLen),sizeof(int));
  
  //write out equation of state string
  ofOut.write((char*)(parameters.sEOSFileName.c_str())
    ,nEOSFileNameLen*sizeof(char));
  
  //write out artificial viscosity
  ofOut.write((char*)(&parameters.dA),sizeof(double));
  
  //write out artificial viscosity threshold
  ofOut.write((char*)(&parameters.dAVThreshold),sizeof(double));
  
  if(procTop.nRank==0){
    
    //write out processor dimensions
    ofOut.write((char*)(procTop.nProcDims),3*sizeof(int));
    
    //write out processor coordinates
    ofOut.write((char*)(procTop.nCoords[procTop.nRank]),3*sizeof(int));
    
    //write out preiodicity
    ofOut.write((char*)(procTop.nPeriodic),3*sizeof(int));
    
    //write out number of variables
    ofOut.write((char*)(&grid.nNumVars),sizeof(int));
    
    //write out variable info
    for(int i=0;i<grid.nNumVars;i++){
      ofOut.write((char*)(grid.nVariables[i]),(4)*sizeof(int));
    }
    
    //write out number of 1D zones
    ofOut.write((char*)(&grid.nNum1DZones),sizeof(int));
    
    //write out global grid size
    ofOut.write((char*)(grid.nGlobalGridDims),3*sizeof(int));
    
    //write out localgrid size
    for(int i=0;i<grid.nNumVars;i++){
      ofOut.write((char*)(grid.nLocalGridDims[procTop.nRank][i]),3*sizeof(int));
    }
    
    //write number of ghostcells
    ofOut.write((char*)(&grid.nNumGhostCells),sizeof(int));
    
    //write out processor local grid
    for(int n=0;n<grid.nNumVars;n++){
      
      int nGhostCellsX=1;
      if(grid.nVariables[n][0]==-1){
        nGhostCellsX=0;
      }
      
      //write out inner grid
      for(int i=0;i<grid.nLocalGridDims[procTop.nRank][n][0]+nGhostCellsX*grid.nNumGhostCells;i++){
        for(int j=0;j<grid.nLocalGridDims[procTop.nRank][n][1];j++){
          //have to write out the multidimensional array one row at a time since
          // there is no garantee that from one row to another the memory is
          // contiguous
          ofOut.write((char*)(grid.dLocalGridOld[n][i][j])
            ,(grid.nLocalGridDims[procTop.nRank][n][2])*sizeof(double));
        }
      }
      
      //write out outter grid
      int nSizeY=grid.nGlobalGridDims[1];
      int nSizeZ=grid.nGlobalGridDims[2];
      if(grid.nVariables[n][1]==-1){
        nSizeY=procTop.nProcDims[1];//if not defined in that y-direction
                            //allow space for each neighboring processor to send data
      }
      if(grid.nVariables[n][2]==-1){
        nSizeZ=procTop.nProcDims[2];//if not defined in that z-direction
                            //allow space for each neighboring processor to send data
      }
      for(int i=grid.nLocalGridDims[procTop.nRank][n][0]+nGhostCellsX*grid.nNumGhostCells;
        i<grid.nLocalGridDims[procTop.nRank][n][0]+nGhostCellsX*2*grid.nNumGhostCells;i++){
        for(int j=0;j<nSizeY;j++){
          //have to write out the multidimensional array one row at a time since there is no 
          //garantee that from one row to another the memory is contiguous
          ofOut.write((char*)(grid.dLocalGridOld[n][i][j]),nSizeZ*sizeof(double));
        }
      }
    }
    ofOut.flush();
    ofOut.close();
  }
  else{
    
    //write out processor coordinates
    ofOut.write((char*)(procTop.nCoords[procTop.nRank]),3*sizeof(int));
    
    //write out number of variables
    ofOut.write((char*)(&grid.nNumVars),sizeof(int));
    
    //write out variable info
    for(int i=0;i<grid.nNumVars;i++){
      ofOut.write((char*)(grid.nVariables[i]),(4)*sizeof(int));
    }
    
    //write out localgrid size
    for(int i=0;i<grid.nNumVars;i++){
      ofOut.write((char*)(grid.nLocalGridDims[procTop.nRank][i]),3*sizeof(int));
    }
    
    //write number of ghostcells
    ofOut.write((char*)(&grid.nNumGhostCells),sizeof(int));
    
    //write out processor local grid
    for(int n=0;n<grid.nNumVars;n++){
      int nGhostCellsX=1;
      int nGhostCellsY=1;
      int nGhostCellsZ=1;
      if(grid.nVariables[n][0]==-1){
        nGhostCellsX=0;
      }
      if(grid.nVariables[n][1]==-1){
        nGhostCellsY=0;
      }
      if(grid.nVariables[n][2]==-1){
        nGhostCellsZ=0;
      }
      for(int i=0;i<grid.nLocalGridDims[procTop.nRank][n][0]
        +nGhostCellsX*2*grid.nNumGhostCells;i++){
        for(int j=0;j<grid.nLocalGridDims[procTop.nRank][n][1]
          +nGhostCellsY*2*grid.nNumGhostCells;j++){
          //have to write out the multidimensional array one row at a time since there is no
          //garantee that from one row to another the memory is contiguous
          ofOut.write((char*)(grid.dLocalGridOld[n][i][j])
            ,(grid.nLocalGridDims[procTop.nRank][n][2]
            +nGhostCellsZ*2*grid.nNumGhostCells)*sizeof(double));
        }
      }
    }
    ofOut.flush();
    ofOut.close();
  }
}
void modelRead(std::string sFileName,ProcTop &procTop, Grid &grid, Time &time
  , Parameters &parameters){
  
  //open file
  std::ifstream ifIn;
  ifIn.open(sFileName.c_str(),std::ios::binary);
  if(!ifIn.is_open()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": error opening the file \""<<sFileName.c_str()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //set up array to hold size of dimensions
  grid.nGlobalGridDims=new int[3];
  
  //check file type
  char cTemp;
  ifIn.read((char*)(&cTemp),sizeof(char));
  if(cTemp!='b'){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": file \""<<sFileName<<"\" is not a binary file.\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //check file version
  int nTemp;
  ifIn.read((char*)(&nTemp),sizeof(int));
  if(nTemp!=DUMP_VERSION){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": file \""<<sFileName<<"\" has version \""<<nTemp
      <<"\" which is not the same as the supported version \"DUMP_VERSION\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //read in time
  ifIn.read((char*)(&time.dt),sizeof(double));
  
  //read in time step index
  ifIn.read((char*)(&time.nTimeStepIndex),sizeof(int));
  
  //read in last time step
  ifIn.read((char*)(&time.dDeltat_nm1half),sizeof(double));
  
  //read in last time step
  ifIn.read((char*)(&time.dDeltat_np1half),sizeof(double));
  
  //set other time values to reasonable initial values
  time.dDeltat_n=(time.dDeltat_nm1half+time.dDeltat_np1half)*0.5;
  
  //read in dAlpha
  ifIn.read((char*)(&parameters.dAlpha),sizeof(double));
  
  //read in gammalaw gas
  int nGammaLaw;
  ifIn.read((char*)(&nGammaLaw),sizeof(int));
  
  if(nGammaLaw==0){//if zero use a gamma law gas
    ifIn.read((char*)(&parameters.dGamma),sizeof(double));
    parameters.bEOSGammaLaw=true;
    if(parameters.sEOSFileName.size()!=0){//if we have an EOS file name from configuration file
      std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": WARNING the equation of state file specified in the configuration file, \""
      <<parameters.sEOSFileName<<"\" is being ignored since the model specifies a gamma law gas.\n";
    }
  }
  else{
    char *cBuffer=new char[nGammaLaw+1];//if nGammaLaw not zero then it is the size of the string following it
    ifIn.read(cBuffer,(nGammaLaw)*sizeof(char));
    cBuffer[nGammaLaw]='\0';
    std::string sTemp=cBuffer;
    delete [] cBuffer;
    if(parameters.sEOSFileName.size()==0){//sEOSFileName already read form configuration file if specified
      parameters.sEOSFileName=sTemp;
    }
    if(sTemp!=parameters.sEOSFileName){//using a different EOS file warn about this
      std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": WARNING the equation of state file specified in the model, \""<<sTemp
      <<"\" is not the same as that specified in the configuration file, \""
      <<parameters.sEOSFileName<<"\". Using file given in configuration file for EOS.\n";
    }
    parameters.bEOSGammaLaw=false;
  }
  
  //write out artificial viscosity
  ifIn.read((char*)(&parameters.dA),sizeof(double));
  
  //write out artificial viscosity threshold
  ifIn.read((char*)(&parameters.dAVThreshold),sizeof(double));
  
  //read in grid dimensions
  ifIn.read((char*)(grid.nGlobalGridDims),3*sizeof(int));
  
  //allocate memory to hold periodicity, and read it in
  procTop.nPeriodic=new int[3];
  ifIn.read((char*)(procTop.nPeriodic),3*sizeof(int));
  
  //read in number of 1D zones
  ifIn.read((char*)(&grid.nNum1DZones),sizeof(int));
  if(grid.nGlobalGridDims[0]>grid.nNum1DZones&&procTop.nNumProcs<=1){//need at least 2 processors
    /**\todo At some point should get it working with only 1 processor*/
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": file \""<<sFileName<<"\"\n"
      <<"  has "<<grid.nNum1DZones<<" 1D zones and a greater number of global grid zones ("
      <<grid.nGlobalGridDims[0]
      <<") in direction 0.\n  The number of processors ("<<procTop.nNumProcs
      <<") is not greater than 1. At least 2 processors are requiredfor for 3D grids.\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  ifIn.read((char*)(&grid.nNumGhostCells),sizeof(int));
  ifIn.read((char*)(&grid.nNumVars),sizeof(int));
  
  //set number of dimensions
  if(grid.nGlobalGridDims[0]>grid.nNumGhostCells&&grid.nGlobalGridDims[1]==1
    &&grid.nGlobalGridDims[2]==1){//1D radial
    grid.nNumDims=1;
  }
  else if(grid.nGlobalGridDims[0]>grid.nNumGhostCells&&grid.nGlobalGridDims[1]>1
    &&grid.nGlobalGridDims[2]==1){//2D radial-theta
    grid.nNumDims=2;
  }
  else if(grid.nGlobalGridDims[0]>grid.nNumGhostCells&&grid.nGlobalGridDims[1]>1
    &&grid.nGlobalGridDims[2]>1){//3D radial-theta-phi
    grid.nNumDims=3;
  }
  else{
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": the specified global grid dims ("<<grid.nGlobalGridDims[0]<<","<<grid.nGlobalGridDims[1]
      <<","<<grid.nGlobalGridDims[2]<<") do not define a 1D radial, 2D radial-theta, "
      <<"or 3D radial-theta-phi grid.\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //dimension checking
  if(grid.nGlobalGridDims[0]<grid.nNumGhostCells){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": the size of the radial dimension ("<<grid.nGlobalGridDims[0]
      <<" zones) is less than the number of ghost cells ("<<grid.nNumGhostCells<<")\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  if(grid.nGlobalGridDims[1]<grid.nNumGhostCells&&grid.nNumDims>1){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": the size of the theta dimension ("<<grid.nGlobalGridDims[1]
      <<" zones) is less than the number of ghost cells ("<<grid.nNumGhostCells
      <<") and the number of dimensions is 2. Try changing the number of theta to above "
      <<grid.nNumGhostCells<<" or to 1 for a 1D calculation.\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  if(grid.nGlobalGridDims[2]<grid.nNumGhostCells&&grid.nNumDims>2){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": the size of the phi dimension ("<<grid.nGlobalGridDims[2]
      <<" zones) is less than the number of ghost cells ("<<grid.nNumGhostCells
      <<") and the number of dimensions is 3. Try changing the number of phi to above "
      <<grid.nNumGhostCells<<" or to 1 for a 2D or 1D calculation.\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //set number of internal variables, and indexes of external and internal variables
  if(parameters.bEOSGammaLaw){
    if(grid.nNumDims==1){
      grid.nNumIntVars=4;
      
      grid.nM                 = 0;
      grid.nDM                = 1;
      grid.nR                 = 2;
      grid.nD                 = 3;
      grid.nU                 = 4;
      grid.nU0                = 5;
      grid.nE                 = 6;
      grid.nP                 = grid.nNumVars+0;
      grid.nQ0                = grid.nNumVars+1;
      grid.nDenAve            = grid.nNumVars+2;
      grid.nDonorCellFrac     = grid.nNumVars+3;
    }
    else if(grid.nNumDims==2){
      grid.nNumIntVars=9;
      
      grid.nM                 = 0;
      grid.nTheta             = 1;
      grid.nDM                = 2;
      grid.nR                 = 3;
      grid.nD                 = 4;
      grid.nU                 = 5;
      grid.nU0                = 6;
      grid.nV                 = 7;
      grid.nE                 = 8;
      grid.nP                 = grid.nNumVars+0;
      grid.nQ0                = grid.nNumVars+1;
      grid.nDenAve            = grid.nNumVars+2;
      grid.nDCosThetaIJK      = grid.nNumVars+3;
      grid.nQ1                = grid.nNumVars+4;
      grid.nDTheta            = grid.nNumVars+5;
      grid.nSinThetaIJK       = grid.nNumVars+6;
      grid.nSinThetaIJp1halfK = grid.nNumVars+7;
      grid.nDonorCellFrac     = grid.nNumVars+8;
    }
    else if(grid.nNumDims==3){
      grid.nNumIntVars=13;
      
      grid.nM                 = 0;
      grid.nTheta             = 1;
      grid.nPhi               = 2;
      grid.nDM                = 3;
      grid.nR                 = 4;
      grid.nD                 = 5;
      grid.nU                 = 6;
      grid.nU0                = 7;
      grid.nV                 = 8;
      grid.nW                 = 9;
      grid.nE                 = 10;
      grid.nP                 = grid.nNumVars+0;
      grid.nQ0                = grid.nNumVars+1;
      grid.nDenAve            = grid.nNumVars+2;
      grid.nDPhi              = grid.nNumVars+3;
      grid.nDCosThetaIJK      = grid.nNumVars+4;
      grid.nQ1                = grid.nNumVars+5;
      grid.nDTheta            = grid.nNumVars+6;
      grid.nSinThetaIJK       = grid.nNumVars+7;
      grid.nSinThetaIJp1halfK = grid.nNumVars+8;
      grid.nCotThetaIJK       = grid.nNumVars+9;
      grid.nCotThetaIJp1halfK = grid.nNumVars+10;
      grid.nQ2                = grid.nNumVars+11;
      grid.nDonorCellFrac    = grid.nNumVars+12;
    }
  }
  else{
    if(parameters.nTypeTurbulanceMod>0){//uses a turbulance model
      if(grid.nNumDims==1){
        grid.nNumIntVars=8;
        
        grid.nM                 = 0;
        grid.nDM                = 1;
        grid.nR                 = 2;
        grid.nD                 = 3;
        grid.nU                 = 4;
        grid.nU0                = 5;
        grid.nT                 = 6;
        grid.nP                 = grid.nNumVars+0;
        grid.nQ0                = grid.nNumVars+1;
        grid.nE                 = grid.nNumVars+2;
        grid.nKappa             = grid.nNumVars+3;
        grid.nGamma             = grid.nNumVars+4;
        grid.nDenAve            = grid.nNumVars+5;
        grid.nEddyVisc          = grid.nNumVars+6;
        grid.nDonorCellFrac     = grid.nNumVars+7;
      }
      else if(grid.nNumDims==2){
        grid.nNumIntVars=15;
        
        grid.nM                 = 0;
        grid.nTheta             = 1;
        grid.nDM                = 2;
        grid.nR                 = 3;
        grid.nD                 = 4;
        grid.nU                 = 5;
        grid.nU0                = 6;
        grid.nV                 = 7;
        grid.nT                 = 8;
        grid.nP                 = grid.nNumVars+0;
        grid.nQ0                = grid.nNumVars+1;
        grid.nDenAve            = grid.nNumVars+2;
        grid.nDCosThetaIJK      = grid.nNumVars+3;
        grid.nE                 = grid.nNumVars+4;
        grid.nKappa             = grid.nNumVars+5;
        grid.nGamma             = grid.nNumVars+6;
        grid.nQ1                = grid.nNumVars+7;
        grid.nDTheta            = grid.nNumVars+8;
        grid.nSinThetaIJK       = grid.nNumVars+9;
        grid.nSinThetaIJp1halfK = grid.nNumVars+10;
        grid.nCotThetaIJK       = grid.nNumVars+11;
        grid.nCotThetaIJp1halfK = grid.nNumVars+12;
        grid.nEddyVisc          = grid.nNumVars+13;
        grid.nDonorCellFrac     = grid.nNumVars+14;
      }
      else if(grid.nNumDims==3){
        grid.nNumIntVars=17;
        
        grid.nM                 = 0;
        grid.nTheta             = 1;
        grid.nPhi               = 2;
        grid.nDM                = 3;
        grid.nR                 = 4;
        grid.nD                 = 5;
        grid.nU                 = 6;
        grid.nU0                = 7;
        grid.nV                 = 8;
        grid.nW                 = 9;
        grid.nT                 = 10;
        grid.nP                 = grid.nNumVars+0;
        grid.nQ0                = grid.nNumVars+1;
        grid.nDenAve            = grid.nNumVars+2;
        grid.nDPhi              = grid.nNumVars+3;
        grid.nDCosThetaIJK      = grid.nNumVars+4;
        grid.nE                 = grid.nNumVars+5;
        grid.nKappa             = grid.nNumVars+6;
        grid.nGamma             = grid.nNumVars+7;
        grid.nQ1                = grid.nNumVars+8;
        grid.nDTheta            = grid.nNumVars+9;
        grid.nSinThetaIJK       = grid.nNumVars+10;
        grid.nSinThetaIJp1halfK = grid.nNumVars+11;
        grid.nCotThetaIJK       = grid.nNumVars+12;
        grid.nCotThetaIJp1halfK = grid.nNumVars+13;
        grid.nQ2                = grid.nNumVars+14;
        grid.nEddyVisc          = grid.nNumVars+15;
        grid.nDonorCellFrac     = grid.nNumVars+16;
      }
    }
    else{//no turulance model
      if(grid.nNumDims==1){
        grid.nNumIntVars=6;
        
        grid.nM                 = 0;
        grid.nDM                = 1;
        grid.nR                 = 2;
        grid.nD                 = 3;
        grid.nU                 = 4;
        grid.nU0                = 5;
        grid.nT                 = 6;
        grid.nP                 = grid.nNumVars+0;
        grid.nQ0                = grid.nNumVars+1;
        grid.nE                 = grid.nNumVars+2;
        grid.nKappa             = grid.nNumVars+3;
        grid.nGamma             = grid.nNumVars+4;
        grid.nDonorCellFrac     = grid.nNumVars+5;
      }
      else if(grid.nNumDims==2){
        grid.nNumIntVars=12;
        
        grid.nM                 = 0;
        grid.nTheta             = 1;
        grid.nDM                = 2;
        grid.nR                 = 3;
        grid.nD                 = 4;
        grid.nU                 = 5;
        grid.nU0                = 6;
        grid.nV                 = 7;
        grid.nT                 = 8;
        grid.nP                 = grid.nNumVars+0;
        grid.nQ0                = grid.nNumVars+1;
        grid.nDenAve            = grid.nNumVars+2;
        grid.nDCosThetaIJK      = grid.nNumVars+3;
        grid.nE                 = grid.nNumVars+4;
        grid.nKappa             = grid.nNumVars+5;
        grid.nGamma             = grid.nNumVars+6;
        grid.nQ1                = grid.nNumVars+7;
        grid.nDTheta            = grid.nNumVars+8;
        grid.nSinThetaIJK       = grid.nNumVars+9;
        grid.nSinThetaIJp1halfK = grid.nNumVars+10;
        grid.nDonorCellFrac     = grid.nNumVars+11;
      }
      else if(grid.nNumDims==3){
        grid.nNumIntVars=16;
        
        grid.nM                 = 0;
        grid.nTheta             = 1;
        grid.nPhi               = 2;
        grid.nDM                = 3;
        grid.nR                 = 4;
        grid.nD                 = 5;
        grid.nU                 = 6;
        grid.nU0                = 7;
        grid.nV                 = 8;
        grid.nW                 = 9;
        grid.nT                 = 10;
        grid.nP                 = grid.nNumVars+0;
        grid.nQ0                = grid.nNumVars+1;
        grid.nDenAve            = grid.nNumVars+2;
        grid.nDPhi              = grid.nNumVars+3;
        grid.nDCosThetaIJK      = grid.nNumVars+4;
        grid.nE                 = grid.nNumVars+5;
        grid.nKappa             = grid.nNumVars+6;
        grid.nGamma             = grid.nNumVars+7;
        grid.nQ1                = grid.nNumVars+8;
        grid.nDTheta            = grid.nNumVars+9;
        grid.nSinThetaIJK       = grid.nNumVars+10;
        grid.nSinThetaIJp1halfK = grid.nNumVars+11;
        grid.nCotThetaIJK       = grid.nNumVars+12;
        grid.nCotThetaIJp1halfK = grid.nNumVars+13;
        grid.nQ2                = grid.nNumVars+14;
        grid.nDonorCellFrac     = grid.nNumVars+15;
      }
    }
  }
  
  //set variable infos, for non-internal variables
  grid.nVariables=new int*[grid.nNumVars+grid.nNumIntVars];
  for(int n=0;n<grid.nNumVars;n++){
    grid.nVariables[n]=new int[3];//+1 because of keeping track of time info
    ifIn.read((char*)(grid.nVariables[n]),(4)*sizeof(int));
    if(grid.nNum1DZones==grid.nGlobalGridDims[0]){//there is no need to define variable in any direction other than radial
      grid.nVariables[n][1]=-1;//not defined in theta
      grid.nVariables[n][2]=-1;//not defined in phi
    }
    
    //checking that all variables are not defined in direciton 2 if only 1D or 2D
    if(grid.nNumDims<3&&grid.nVariables[n][2]!=-1){
      std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
        <<": WARNING the input file, \""<<sFileName<<"\", had variable "<<n
        <<" set as defined in direction 2 but model has no direciton 2, undefining variable in"
        <<" direciton 2\n";
      grid.nVariables[n][2]=-1;
    }
    
    //checking that all variables are not defined in direction 1 if only 1D
    if(grid.nNumDims<2&&grid.nVariables[n][1]!=-1){
      std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
        <<": WARNING the input file, \""<<sFileName<<"\", had variable "<<n
        <<" set as defined in direction 1 but model has no direciton 2, undefining variable in"
        <<" direciton 1\n";
      grid.nVariables[n][1]=-1;
    }
  }
  
  //set internal variable infos
  setInternalVarInf(grid,parameters);
  
  if(grid.nNum1DZones==grid.nGlobalGridDims[0]){/*check that there is something for all processors
    to do*/
    if(procTop.nNumProcs>1){
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": file \""<<sFileName<<"\" has "<<grid.nNum1DZones
        <<" 1D zones and an equal number of global grid zones,"<<grid.nGlobalGridDims[0]
        <<", in direction 0, and the number of processors, "<<procTop.nNumProcs
        <<", is greater than 1, nothing for other processors to do.\n";
      throw exception2(ssTemp.str(),INPUT);
    }
  }
  
  //set up data storage and processor topography
  setupLocalGrid(procTop,grid);
  
  if(procTop.nRank==0){//read in grid
    for(int n=0;n<grid.nNumVars;n++){
      
      int nGhostCellsX=1;
      
      //set global size of grid for var n, minus the 1D region, and ghostcells
      int nGlobalSize[3]={1,1,1};
      if(grid.nVariables[n][0]!=-1){
        nGlobalSize[0]=grid.nGlobalGridDims[0]-grid.nNum1DZones;/*don't need to add extra for
          interface variables since 1d region allready read extra interface*/
        if(procTop.nProcDims[0]==1){
          nGlobalSize[0]=0;//only 1 processor in this direction, so no x-zones left
        }
      }
      if(grid.nVariables[n][1]!=-1){
        nGlobalSize[1]=grid.nGlobalGridDims[1];
        if(procTop.nPeriodic[1]==0){//if not periodic add
          nGlobalSize[1]+=grid.nVariables[n][1];
        }
      }
      if(grid.nVariables[n][2]!=-1){
        nGlobalSize[2]=grid.nGlobalGridDims[2];
        if(procTop.nPeriodic[2]==0){
          nGlobalSize[2]+=grid.nVariables[n][2];
        }
      }
      if(grid.nVariables[n][0]!=-1){//only read it if variable is defined in x-direction
      
        //read in inner 1D region
        for(int i=0;i<grid.nLocalGridDims[procTop.nRank][n][0]+grid.nNumGhostCells;i++){
          for(int j=0;j<grid.nLocalGridDims[procTop.nRank][n][1];j++){
            ifIn.read((char*)(grid.dLocalGridOld[n][i][j])
              ,(grid.nLocalGridDims[procTop.nRank][n][2])*sizeof(double));
          }
        }
        
        int nSkipSize[3]={0,0,0};
        if(grid.nVariables[n][1]!=-1){//y-skip size
          if(grid.nNumDims>1){
            nSkipSize[1]=grid.nNumGhostCells;
            if(grid.nNumDims>2){
              if(grid.nVariables[n][2]!=-1){
                nSkipSize[1]=nSkipSize[1]*(2*grid.nNumGhostCells+nGlobalSize[2]);
              }
            }
          }
        }
        if(grid.nVariables[n][2]!=-1){//z-skip size
          if(grid.nNumDims>2){
            nSkipSize[2]=grid.nNumGhostCells;
          }
        }
        
        //read in outter 1D region
        for(int i=grid.nLocalGridDims[procTop.nRank][n][0]+nGhostCellsX*grid.nNumGhostCells;
          i<grid.nLocalGridDims[procTop.nRank][n][0]+2*nGhostCellsX*grid.nNumGhostCells;i++){
          
          //skip inner y-ghostcells
          ifIn.seekg(nSkipSize[1]*sizeof(double),std::ios_base::cur);
          
          for(int j=0;j<nGlobalSize[1];j++){
            
            //skip inner z-ghost cells
            ifIn.seekg(nSkipSize[2]*sizeof(double),std::ios_base::cur);
            
            ifIn.read((char*)(grid.dLocalGridOld[n][i][j]),(nGlobalSize[2])*sizeof(double));
            //may need to copy these around if the variable is not defined in y and or z directions
            //grid will be the size of the y and z processor dimensions
            
            //skip outer z-ghost cells
            ifIn.seekg(nSkipSize[2]*sizeof(double),std::ios_base::cur);
          }
          //skip outer y-ghost cells
          ifIn.seekg(nSkipSize[1]*sizeof(double),std::ios_base::cur);
        }
        
        //skip rest of grid for var n
        int nNumSkip=1;
        if(grid.nVariables[n][0]!=-1){
          nNumSkip=nGlobalSize[0];
        }
        if(grid.nVariables[n][1]!=-1){
          if(grid.nNumDims>1){
            nNumSkip=nNumSkip*(2*grid.nNumGhostCells+nGlobalSize[1]);
          }
        }
        if(grid.nVariables[n][2]!=-1){
          if(grid.nNumDims>2){
            nNumSkip=nNumSkip*(2*grid.nNumGhostCells+nGlobalSize[2]);
          }
        }
        ifIn.seekg(nNumSkip*sizeof(double),std::ios_base::cur);
      }
      else{
        //skip variable, since it isn't defined in x-direction
        int nNumSkip=1;
        if(grid.nVariables[n][1]!=-1){
          if(grid.nNumDims>1){
            nNumSkip=nNumSkip*(2*grid.nNumGhostCells+nGlobalSize[1]);
          }
        }
        if(grid.nVariables[n][2]!=-1){
          if(grid.nNumDims>2){
            nNumSkip=nNumSkip*(2*grid.nNumGhostCells+nGlobalSize[2]);
          }
        }
        ifIn.seekg(nNumSkip*sizeof(double),std::ios_base::cur);
      }
    }
  }
  else{//read in grid
    for(int n=0;n<grid.nNumVars;n++){
      
      int nPosGrid[3]={0,0,0};//holds start position of processor procTop.nRank in global grid
      
      //add any offset due to position in dimension 2
      for(int p=1;p<procTop.nRank;p++){
        if(procTop.nCoords[p][2]<procTop.nCoords[procTop.nRank][2]
          &&procTop.nCoords[p][1]==procTop.nCoords[procTop.nRank][1]
          &&procTop.nCoords[p][0]==procTop.nCoords[procTop.nRank][0]){
          if(grid.nVariables[n][2]!=-1){
            nPosGrid[2]+=grid.nLocalGridDims[p][n][2];
          }
        }
      }
      
      //Add any offset due to position in dimension 1
      for(int p=1;p<procTop.nRank;p++){
        if(procTop.nCoords[p][2]==procTop.nCoords[procTop.nRank][2]
          &&procTop.nCoords[p][1]<procTop.nCoords[procTop.nRank][1]
          &&procTop.nCoords[p][0]==procTop.nCoords[procTop.nRank][0]){
          if(grid.nVariables[n][1]!=-1){
            nPosGrid[1]+=grid.nLocalGridDims[p][n][1];
          }
        }
      }
      
      //Add any offset due to position in dimension 0
      for(int p=1;p<procTop.nRank;p++){
        if(procTop.nCoords[p][2]==procTop.nCoords[procTop.nRank][2]
          &&procTop.nCoords[p][1]==procTop.nCoords[procTop.nRank][1]
          &&procTop.nCoords[p][0]<procTop.nCoords[procTop.nRank][0]){
          if(grid.nVariables[n][0]!=-1){
            nPosGrid[0]+=grid.nLocalGridDims[p][n][0];
          }
        }
      }
      
      //calculate total global sizes
      int nSize[3]={1,1,1};//1 if variable not defined in that direction
      for(int l=0;l<3;l++){
        if(grid.nVariables[n][l]!=-1){
          if(procTop.nPeriodic[l]==0){
            nSize[l]=grid.nGlobalGridDims[l]+grid.nVariables[n][l]+2*grid.nNumGhostCells;
          }
          else{
            nSize[l]=grid.nGlobalGridDims[l]+2*grid.nNumGhostCells;
          }
        }
      }
      
      //set starting position
      int nStart=0;
      if(procTop.nCoords[procTop.nRank][0]==1){
        if(grid.nVariables[n][0]==1&&procTop.nPeriodic[0]==0){
          //adds an additional number for the inner interface
          nStart=grid.nNum1DZones+1+nPosGrid[2]+nPosGrid[1]*nSize[2]+nPosGrid[0]*nSize[2]*nSize[1];
        }
        else{
          nStart=grid.nNum1DZones+nPosGrid[2]+nPosGrid[1]*nSize[2]+nPosGrid[0]*nSize[2]*nSize[1];
        }
      }
      else{
        if(grid.nVariables[n][0]==1){//adds an additional number for the inner interface
          nStart=grid.nNum1DZones+1+grid.nNumGhostCells+nPosGrid[2]+nPosGrid[1]*nSize[2]
          +(nPosGrid[0]-grid.nNumGhostCells)*nSize[2]*nSize[1];
          //accounts for the fact that the first 2 ghost cells are in 1D region
        }
        else{
          nStart=grid.nNum1DZones+grid.nNumGhostCells+nPosGrid[2]+nPosGrid[1]*nSize[2]
          +(nPosGrid[0]-grid.nNumGhostCells)*nSize[2]*nSize[1];
          //accounts for the fact that the first 2 ghost cells are in 1D region
        }
      }
        
      //find out if we need ghost cells
      int nGhostCellsX=0;
      if(grid.nVariables[n][0]!=-1){
        nGhostCellsX=1;
      }
      int nGhostCellsY=0;
      if(grid.nVariables[n][1]!=-1&&grid.nNumDims>1){
        nGhostCellsY=1;
      }
      int nGhostCellsZ=0;
      if(grid.nVariables[n][2]!=-1&&grid.nNumDims>2){
        nGhostCellsZ=1;
      }
      
      //calculate spacings
      int nZSpacing=0;
      if(grid.nVariables[n][2]!=-1){
        nZSpacing=(nSize[2]-grid.nLocalGridDims[procTop.nRank][n][2]
          -2*grid.nNumGhostCells)*sizeof(double);
      }
      int nYSpacing=0;
      if(grid.nVariables[n][1]!=-1){
        nYSpacing=(nSize[1]-grid.nLocalGridDims[procTop.nRank][n][1]-2*grid.nNumGhostCells)
          *(nSize[2])*sizeof(double);
      }
      int nROffset=0;
      if(grid.nVariables[n][0]==1&&procTop.nPeriodic[0]==0){//if interface variable and not periodic
        nROffset=1;
      }
      if(procTop.nCoords[procTop.nRank][0]==1){//boardering 1D region
        
        //read in inner ghost cells
        ifIn.seekg(nGhostCellsX*(grid.nNum1DZones+nROffset)*sizeof(double),std::ios_base::cur);
        for(int i=0;i<nGhostCellsX*grid.nNumGhostCells;i++){
          double dTemp;
          ifIn.read((char*)(&dTemp),sizeof(double));
          
          //copy to all y and z at that x in old grid
          for(int j=0;j<grid.nLocalGridDims[procTop.nRank][n][1]
            +nGhostCellsY*2*grid.nNumGhostCells;j++){
            for(int k=0;k<grid.nLocalGridDims[procTop.nRank][n][2]
              +nGhostCellsZ*2*grid.nNumGhostCells;k++){
              grid.dLocalGridOld[n][i][j][k]=dTemp;
            }
          }
        }
        
        //move to the start for current processor, after inner ghost cells read in
        ifIn.seekg(nGhostCellsX*(nStart-(grid.nNum1DZones+nROffset))*sizeof(double)
          ,std::ios_base::cur);
        
        //read in rest of grid
        for(int i=nGhostCellsX*grid.nNumGhostCells;i<grid.nLocalGridDims[procTop.nRank][n][0]
          +nGhostCellsX*2*grid.nNumGhostCells;i++){
          for(int j=0;j<grid.nLocalGridDims[procTop.nRank][n][1]
            +nGhostCellsY*2*grid.nNumGhostCells;j++){
            ifIn.read((char*)(grid.dLocalGridOld[n][i][j])
              ,(grid.nLocalGridDims[procTop.nRank][n][2]+nGhostCellsZ*2*grid.nNumGhostCells)
              *sizeof(double));
            
            //skip some z for other processors
            ifIn.seekg(nZSpacing,std::ios_base::cur);
          }
          //skip some y for other processors
          ifIn.seekg(nYSpacing,std::ios_base::cur);
        }
      }
      else{//not boardering 1D region
        
        //move to the start for current processor
        ifIn.seekg(nGhostCellsX*nStart*sizeof(double),std::ios_base::cur);
        
        //read in rest of grid
        for(int i=0;i<grid.nLocalGridDims[procTop.nRank][n][0]
          +nGhostCellsX*2*grid.nNumGhostCells;i++){
          for(int j=0;j<grid.nLocalGridDims[procTop.nRank][n][1]
            +nGhostCellsY*2*grid.nNumGhostCells;j++){
            ifIn.read((char*)(grid.dLocalGridOld[n][i][j])
              ,(grid.nLocalGridDims[procTop.nRank][n][2]+nGhostCellsZ*2*grid.nNumGhostCells)
              *sizeof(double));
              
            //skip some z for other processors
            ifIn.seekg(nZSpacing,std::ios_base::cur);
          }
          //skip some y for other processors
          ifIn.seekg(nYSpacing,std::ios_base::cur);
        }
      }
      
      //skip rest of grid for var n
      int nTotal=nSize[2]*nSize[1]*(nSize[0]-nGhostCellsX*(grid.nNumGhostCells
        +grid.nNum1DZones+nROffset))+nGhostCellsX*(grid.nNumGhostCells+grid.nNum1DZones+nROffset);
      int nReadIn=nSize[2]*nSize[1]*(grid.nLocalGridDims[procTop.nRank][n][0]
        +2*nGhostCellsX*grid.nNumGhostCells);
      if(procTop.nCoords[procTop.nRank][0]==1){
        nReadIn=nGhostCellsX*grid.nNumGhostCells
          +nSize[2]*nSize[1]*(grid.nLocalGridDims[procTop.nRank][n][0]
          +nGhostCellsX*grid.nNumGhostCells);
      }
      int nSkip=nTotal-nGhostCellsX*nStart-nReadIn;
      ifIn.seekg(nSkip*sizeof(double),std::ios_base::cur);
    }
  }
  ifIn.close();
}
void initUpdateLocalBoundaries(ProcTop &procTop, Grid &grid, MessPass &messPass,Implicit &implicit){
  
  //create send and recieve types
  if(procTop.nRank==0){

    //count number of neighbors
    //no periodic neighbors, go straight to the total number of neighbors
    for(int l=0;l<3;l++){
      for(int p=1;p<procTop.nNumProcs;p++){
        
        //matches x
        if((procTop.nCoords[procTop.nRank][0]+l-1)==procTop.nCoords[p][0]){
          procTop.nNumNeighbors++;
        }
      }
    }
    
    //set id's of neighbors
    procTop.nNeighborRanks=new int[procTop.nNumNeighbors];
    int nIndex=0;
    for(int l=0;l<3;l++){
      for(int p=0;p<procTop.nNumProcs;p++){
        if( ((procTop.nCoords[procTop.nRank][0]+l-1)==procTop.nCoords[p][0]) &&//matches x
          p!=procTop.nRank ){//not the current processor
          procTop.nNeighborRanks[nIndex]=p;
          nIndex++;
        }
      }
    }
    
    //allocate memory for send and recive types
    messPass.typeSendNewGrid=new MPI::Datatype[procTop.nNumNeighbors];
    messPass.typeRecvOldGrid=new MPI::Datatype[procTop.nNumNeighbors];
    messPass.typeSendNewVar=new MPI::Datatype*[procTop.nNumNeighbors];
    messPass.typeRecvNewVar=new MPI::Datatype*[procTop.nNumNeighbors];
    for(int i=0;i<procTop.nNumNeighbors;i++){//allocate space for MPI::Datatypes for each variable
      messPass.typeRecvNewVar[i]=new MPI::Datatype[grid.nNumVars+grid.nNumIntVars];
      messPass.typeSendNewVar[i]=new MPI::Datatype[grid.nNumVars+grid.nNumIntVars];
    }

    //allocate memory to hold send and recieve block info
    int **nSendBlockStart=new int*[grid.nNumVars+grid.nNumIntVars];
    int **nRecvBlockStart=new int*[grid.nNumVars+grid.nNumIntVars];
    int **nSendBlockDims=new int*[grid.nNumVars+grid.nNumIntVars];
    int **nRecvBlockDims=new int*[grid.nNumVars+grid.nNumIntVars];
    for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
      nSendBlockDims[n]=new int[3];
      nRecvBlockDims[n]=new int[3];
      nSendBlockStart[n]=new int[3];
      nRecvBlockStart[n]=new int[3];
    }
    
    //set send and recieve types
    for(int p=0;p<procTop.nNumNeighbors;p++){
      
      //set start of send and recive blocks and set block dimensions
      for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
        
        //set start of send and recive blocks
        nSendBlockStart[n][0]=grid.nNum1DZones;
        nSendBlockStart[n][1]=0;
        nSendBlockStart[n][2]=0;
        nRecvBlockStart[n][0]=grid.nNum1DZones+grid.nNumGhostCells;
        nRecvBlockStart[n][1]=0;
        nRecvBlockStart[n][2]=0;
        if(grid.nVariables[n][0]==1){
          nSendBlockStart[n][0]+=grid.nCenIntOffset[0];
          nRecvBlockStart[n][0]+=grid.nCenIntOffset[0];
        }
        
        for(int j=0;j<procTop.nCoords[procTop.nNeighborRanks[p]][1];j++){
          nRecvBlockStart[n][1]+=grid.nLocalGridDims[procTop.nNeighborRanks[
            j*procTop.nProcDims[2]]][n][1];
        }
        for(int k=0;k<procTop.nCoords[procTop.nNeighborRanks[p]][2];k++){
          nRecvBlockStart[n][2]+=grid.nLocalGridDims[procTop.nNeighborRanks[k]][n][2];
        }
        
        //set block dimensions
        nRecvBlockDims[n][0]=grid.nNumGhostCells;
        nRecvBlockDims[n][1]=grid.nLocalGridDims[procTop.nNeighborRanks[p]][n][1];
        nRecvBlockDims[n][2]=grid.nLocalGridDims[procTop.nNeighborRanks[p]][n][2];
        
        nSendBlockDims[n][0]=grid.nNumGhostCells;
        nSendBlockDims[n][1]=grid.nLocalGridDims[procTop.nNeighborRanks[p]][n][1]
          +2*grid.nNumGhostCells;
        nSendBlockDims[n][2]=grid.nLocalGridDims[procTop.nNeighborRanks[p]][n][2]
          +2*grid.nNumGhostCells;
        
        //adjust for non-3D variables
        for(int l=0;l<3;l++){
          
          //if variable not defined in direction l
          if(grid.nVariables[n][l]==-1){
            
            //set to 0 since there is no other cells in direction l for variable i
            nSendBlockStart[n][l]=0;
            nRecvBlockStart[n][l]=0;
            if(l!=0){//if non-radial
              nRecvBlockStart[n][l]=procTop.nCoords[procTop.nNeighborRanks[p]][l];
              //set a unique spot for each neighbor based in their coordinates
            }
            
            //set to 1 in l direction
            nSendBlockDims[n][l]=1;
            nRecvBlockDims[n][l]=1;
            
          }
        }
        
        //adjust for non-time dependent variables
        if(grid.nVariables[n][3]==0){//don't need to update in time
          for(int l=0;l<3;l++){
            nSendBlockDims[n][l]=0;//set all block dims to 0 (nothing to be sent)
            nRecvBlockDims[n][l]=0;
          }
        }
        
        //check to see if we need to send current neighbor a message for current variable
        for(int l=0;l<1;l++){
          
          //if variable not defined in direction l don't need to send variable info to it
          if(grid.nVariables[n][l]==-1){
            //if current neighbor is not at the same coordinate in direction l
            if(procTop.nCoords[procTop.nNeighborRanks[p]][l]!=procTop.nCoords[procTop.nRank][l]
              &&procTop.nCoords[procTop.nNeighborRanks[p]][l]!=-1){
              for(int q=0;q<3;q++){
                //all block dimensions set to 0, since nothing to send
                nSendBlockDims[n][q]=0;
                nRecvBlockDims[n][q]=0;
              }
              //break;//exit testing for current neighbor
            }
          }
        }
      }
      
      //set sub-block lengths and types for send
      int nNumSubBlocksSend=0;
      for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){//convert to 1D array of blocks
        nNumSubBlocksSend+=nSendBlockDims[n][0]*nSendBlockDims[n][1]*nSendBlockDims[n][2];
      }
      int* nBlockLenSend=new int[nNumSubBlocksSend];
      MPI::Datatype *typeBaseSend=new MPI::Datatype[nNumSubBlocksSend];
      int nIter=0;
      std::vector<std::vector<int> > vecnBlockLenSendNewVar;
      std::vector<std::vector<MPI::Datatype> > vectypeBaseSendNewVar;
      for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
        vecnBlockLenSendNewVar.push_back(std::vector<int>());
        vectypeBaseSendNewVar.push_back(std::vector<MPI::Datatype>());
        for(int i=0;i<nSendBlockDims[n][0];i++){
          for(int j=0;j<nSendBlockDims[n][1];j++){
            for(int k=0;k<nSendBlockDims[n][2];k++){
              nBlockLenSend[nIter]=1;
              vecnBlockLenSendNewVar[n].push_back(1);
              typeBaseSend[nIter]=MPI::DOUBLE;
              vectypeBaseSendNewVar[n].push_back(MPI::DOUBLE);
              nIter++;
            }
          }
        }
      }
      
      //set starting send address
      MPI::Aint nStartAddressSend;
      nStartAddressSend=MPI::Get_address(grid.dLocalGridNew);
      
      //set addresses for send
      MPI::Aint *nSendAddresses=new MPI::Aint[nNumSubBlocksSend];
      int nCount=0;
      std::vector<std::vector<MPI::Aint> > vecSendNewVarAddresses;
      for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
        vecSendNewVarAddresses.push_back( std::vector<MPI::Aint>() );
        for(int i=nSendBlockStart[n][0];i<nSendBlockStart[n][0]+nSendBlockDims[n][0];i++){
          for(int j=nSendBlockStart[n][1];j<nSendBlockStart[n][1]+nSendBlockDims[n][1];j++){
            for(int k=nSendBlockStart[n][2];k<nSendBlockStart[n][2]+nSendBlockDims[n][2];k++){
              MPI::Aint nCurAddress=MPI::Get_address(
                &grid.dLocalGridNew[n][i][nSendBlockStart[n][1]][nSendBlockStart[n][2]]);
              nSendAddresses[nCount]=nCurAddress-nStartAddressSend;
              nCount++;
              vecSendNewVarAddresses[n].push_back(nCurAddress-nStartAddressSend);
            }
          }
        }
      }
      
      //set sub-block lengths and types for recv
      int nNumSubBlocksRecv=0;
      for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
        nNumSubBlocksRecv+=nRecvBlockDims[n][0]*nRecvBlockDims[n][1]*nRecvBlockDims[n][2];
      }
      int* nBlockLenRecv=new int[nNumSubBlocksRecv];
      MPI::Datatype *typeBaseRecv=new MPI::Datatype[nNumSubBlocksRecv];
      nIter=0;
      std::vector<std::vector<int> > vecnBlockLenRecvNewVar;
      std::vector<std::vector<MPI::Datatype> > vectypeBaseRecvNewVar;
      for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
        vecnBlockLenRecvNewVar.push_back(std::vector<int>());
        vectypeBaseRecvNewVar.push_back(std::vector<MPI::Datatype>());
        for(int i=0;i<nRecvBlockDims[n][0];i++){
          for(int j=0;j<nRecvBlockDims[n][1];j++){
            for(int k=0;k<nRecvBlockDims[n][2];k++){
              nBlockLenRecv[nIter]=1;
              vecnBlockLenRecvNewVar[n].push_back(1);
              typeBaseRecv[nIter]=MPI::DOUBLE;
              vectypeBaseRecvNewVar[n].push_back(MPI::DOUBLE);
              nIter++;
            }
          }
        }
      }
      
      //set starting recv address
      MPI::Aint nStartAddressRecv=MPI::Get_address(grid.dLocalGridOld);
      MPI::Aint nStartAddressRecv2=MPI::Get_address(grid.dLocalGridNew);
      
      //set addresses for recv
      MPI::Aint *nRecvAddresses=new MPI::Aint[nNumSubBlocksRecv];
      nCount=0;
      std::vector<std::vector<MPI::Aint> > vecRecvNewVarAddresses;
      for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
        vecRecvNewVarAddresses.push_back( std::vector<MPI::Aint>() );
        for(int i=nRecvBlockStart[n][0];i<nRecvBlockStart[n][0]+nRecvBlockDims[n][0];i++){
          for(int j=nRecvBlockStart[n][1];j<nRecvBlockStart[n][1]+nRecvBlockDims[n][1];j++){
            for(int k=nRecvBlockStart[n][2];k<nRecvBlockStart[n][2]+nRecvBlockDims[n][2];k++){
              MPI::Aint nCurAddress=MPI::Get_address(&grid.dLocalGridOld[n][i][j][k]);
              MPI::Aint nCurAddress2=MPI::Get_address(&grid.dLocalGridNew[n][i][j][k]);
              nRecvAddresses[nCount]=nCurAddress-nStartAddressRecv;
              nCount++;
              vecRecvNewVarAddresses[n].push_back(nCurAddress2-nStartAddressRecv2);
            }
          }
        }
      }
      
      //set send type new grid
      messPass.typeSendNewGrid[p]=MPI::Datatype::Create_struct(nNumSubBlocksSend,nBlockLenSend
        ,nSendAddresses,typeBaseSend);
      messPass.typeSendNewGrid[p].Commit();//note: must delete datatypes
      
      //set send types for new vars
      for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
        
        //copy vecnBlockLenSendNewVar[n] to 1D array
        int *nBlockLen=new int[vecnBlockLenSendNewVar[n].size()];
        std::copy(vecnBlockLenSendNewVar[n].begin(),vecnBlockLenSendNewVar[n].end(),nBlockLen);
        
        //copy vectypeBaseSendNewVar[n] to 1D array
        MPI::Datatype *nBaseType=new MPI::Datatype[vectypeBaseSendNewVar[n].size()];
        std::copy(vectypeBaseSendNewVar[n].begin(),vectypeBaseSendNewVar[n].end(),nBaseType);
        
        //copy vecSendNewVarAddresses[n] to 1D array
        MPI::Aint *nAddresses=new MPI::Aint[vecSendNewVarAddresses[n].size()];
        std::copy(vecSendNewVarAddresses[n].begin(),vecSendNewVarAddresses[n].end(),nAddresses);
        
        //set send type
        messPass.typeSendNewVar[p][n]=MPI::Datatype::Create_struct(vecSendNewVarAddresses[n].size()
          ,nBlockLen,nAddresses,nBaseType);
        messPass.typeSendNewVar[p][n].Commit();
        
        delete [] nBlockLen;
        delete [] nBaseType;
        delete [] nAddresses;
      }
      
      //set recv type
      messPass.typeRecvOldGrid[p]=MPI::Datatype::Create_struct(nNumSubBlocksRecv,nBlockLenRecv
        ,nRecvAddresses,typeBaseRecv);
      messPass.typeRecvOldGrid[p].Commit();//note: must delete datatypes
      
      //set recv types for new vars
      for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
        
        //copy vecnBlockLenSendNewVar[n] to 1D array
        int *nBlockLen=new int[vecnBlockLenRecvNewVar[n].size()];
        std::copy(vecnBlockLenRecvNewVar[n].begin(),vecnBlockLenRecvNewVar[n].end(),nBlockLen);
        
        //copy vectypeBaseSendNewVar[n] to 1D array
        MPI::Datatype *nBaseType=new MPI::Datatype[vectypeBaseRecvNewVar[n].size()];
        std::copy(vectypeBaseRecvNewVar[n].begin(),vectypeBaseRecvNewVar[n].end(),nBaseType);
        
        //copy vecSendNewVarAddresses[n] to 1D array
        MPI::Aint *nAddresses=new MPI::Aint[vecRecvNewVarAddresses[n].size()];
        std::copy(vecRecvNewVarAddresses[n].begin(),vecRecvNewVarAddresses[n].end(),nAddresses);
        
        //set send type
        messPass.typeRecvNewVar[p][n]=MPI::Datatype::Create_struct(vecRecvNewVarAddresses[n].size()
          ,nBlockLen,nAddresses,nBaseType);
        messPass.typeRecvNewVar[p][n].Commit();
        
        delete [] nBlockLen;
        delete [] nBaseType;
        delete [] nAddresses;
      }
      
      //delete dynamic memory
      delete [] nBlockLenRecv;
      delete [] nBlockLenSend;
      delete [] nSendAddresses;
      delete [] nRecvAddresses;
      delete [] typeBaseSend;
      delete [] typeBaseRecv;
    }
    
    //delete dynmic memory
    for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
      delete [] nSendBlockDims[n];
      delete [] nRecvBlockDims[n];
      delete [] nRecvBlockStart[n];
      delete [] nSendBlockStart[n];
    }
    delete [] nSendBlockDims;
    delete [] nRecvBlockDims;
    delete [] nRecvBlockStart;
    delete [] nSendBlockStart;
  }
  else{//if procTop.nRank!=0
    
    //count number of non-periodic neighbors
    int nNumNonPeriodicNeighbors=0;
    for(int i=0;i<3;i++){
      
      //look for neighbors in 3D region
      for(int j=0;j<3;j++){
        for(int k=0;k<3;k++){
          for(int p=1;p<procTop.nNumProcs;p++){
            if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                ( (procTop.nCoords[procTop.nRank][1]+j-1)==procTop.nCoords[p][1])&&//matches y
                ( (procTop.nCoords[procTop.nRank][2]+k-1)==procTop.nCoords[p][2])&&//mathces z
              p!=procTop.nRank){//not the current processor
              nNumNonPeriodicNeighbors++;
              break;//if found one, break, since only one processor will have coordinates (i,j,k)
            }
          }
        }
      }
      
      //check if processor 0 matches in x
      if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[0][0])){//not the current processor
        nNumNonPeriodicNeighbors++;
      }
    }
    
    //count number of periodic neighbors, including the current processor
    int nNumPeriodicNeighbors=0;
    if(procTop.nPeriodic[0]){//if periodic in x-direction
      //not needed, if at a later time it is needed, code will
      // have to be added here
    }
    if(procTop.nPeriodic[1]){//if periodic in y-direction
      if(procTop.nCoords[procTop.nRank][1]==(procTop.nProcDims[1]-1)){//on outer edge of processor grid
        
        //search for neighbors at inner edge
        for(int i=0;i<3;i++){//search x-direction
          for(int k=0;k<3;k++){//search z-direction
            for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
              if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                  ( (0                    )==procTop.nCoords[p][1])&&//on inner edge
                  ( (procTop.nCoords[procTop.nRank][2]+k-1)==procTop.nCoords[p][2])  //mathces z
                ){//current processor allowed as periodic neighbor
                nNumPeriodicNeighbors++;
              }
            }
          }
        }
      }
      if(procTop.nCoords[procTop.nRank][1]==0){//on inner edge of processor grid
        
        //search for neighbors at outter edge
        for(int i=0;i<3;i++){//search x-direction
          for(int k=0;k<3;k++){//search z-direction
            for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
              if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                  ( (procTop.nProcDims[1]-1       )==procTop.nCoords[p][1])&&//on inner edge
                  ( (procTop.nCoords[procTop.nRank][2]+k-1)==procTop.nCoords[p][2])  //mathces z
                ){//current processor allowed as periodic neighbor
                nNumPeriodicNeighbors++;
              }
            }
          }
        }
      }
    }
    if(procTop.nPeriodic[2]){//if periodic in z-direction
      if(procTop.nCoords[procTop.nRank][2]==(procTop.nProcDims[2]-1)){//on outer edge of processor grid
        
        //search for neighbors at inner edge
        for(int i=0;i<3;i++){//search x-direction
          for(int j=0;j<3;j++){//search y-direction
            for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
              if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                  ( (procTop.nCoords[procTop.nRank][1]+j-1)==procTop.nCoords[p][1])&&//matches y
                  ( (0                    )==procTop.nCoords[p][2])  //on inner edge
                ){//current processor allowed as periodic neighbor
                nNumPeriodicNeighbors++;
              }
            }
          }
        }
      }
      if(procTop.nCoords[procTop.nRank][2]==0){//on inner edge of processor grid
      
        //search for neighbors at outter edge
        for(int i=0;i<3;i++){//search x-direction
          for(int j=0;j<3;j++){//search y-direction
            for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
              if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                  ( (procTop.nCoords[procTop.nRank][1]+j-1)==procTop.nCoords[p][1])&&//matches y
                  ( (procTop.nProcDims[2]-1       )==procTop.nCoords[p][2])  //on outer edge
                ){//current processor allowed as periodic neighbor
                nNumPeriodicNeighbors++;
              }
            }
          }
        }
      }
    }
    if(procTop.nPeriodic[0]&&procTop.nPeriodic[1]){//if periodic in x & y -directions
      //periodicity not needed in x-direction
    }
    if(procTop.nPeriodic[0]&&procTop.nPeriodic[2]){//if periodic in x & z -directions
      //periodicity not needed in x-direction
    }
    if(procTop.nPeriodic[1]&&procTop.nPeriodic[2]){//if periodic in y & z -directions
      if( (procTop.nCoords[procTop.nRank][1]==(procTop.nProcDims[1]-1))
        &&(procTop.nCoords[procTop.nRank][2]==(procTop.nProcDims[2]-1)) ){//on outer corner in y and z direction
        
        //search for neighbors at corner 
        for(int i=0;i<3;i++){//search x-direction
          for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
            if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                ( (0                    )==procTop.nCoords[p][1])&&//on inner edge in y
                ( (0                    )==procTop.nCoords[p][2])  //on inner edge in z
              ){//current processor allowed as periodic neighbor
              nNumPeriodicNeighbors++;
            }
          }
        }
      }
      if( (procTop.nCoords[procTop.nRank][1]==(0             ))
        &&(procTop.nCoords[procTop.nRank][2]==(procTop.nProcDims[2]-1)) ){//on inner corner in y, and outer corner in z
                                                  //direction
        //search for neighbors at corner 
        for(int i=0;i<3;i++){//search x-direction
          for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
            if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                ( (procTop.nProcDims[1]-1       )==procTop.nCoords[p][1])&&//on inner edge in y
                ( (0                    )==procTop.nCoords[p][2])  //on inner edge in z
              ){//current processor allowed as periodic neighbor
              nNumPeriodicNeighbors++;
            }
          }
        }
      }
      if( (procTop.nCoords[procTop.nRank][1]==(procTop.nProcDims[1]-1))
        &&(procTop.nCoords[procTop.nRank][2]==(0             )) ){//on outer corner in y, and inner corner in z 
                                                  //direction
        //search for neighbors at corner 
        for(int i=0;i<3;i++){//search x-direction
          for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
            if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                ( (0                    )==procTop.nCoords[p][1])&&//on inner edge in y
                ( (procTop.nProcDims[2]-1       )==procTop.nCoords[p][2])  //on inner edge in z
              ){//current processor allowed as periodic neighbor
              nNumPeriodicNeighbors++;
            }
          }
        }
      }
      if( (procTop.nCoords[procTop.nRank][1]==(0             ))
        &&(procTop.nCoords[procTop.nRank][2]==(0             )) ){//on outer corner in y, and inner corner in z 
                                                  //direction
        //search for neighbors at corner 
        for(int i=0;i<3;i++){//search x-direction
          for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
            if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                ( (procTop.nProcDims[1]-1       )==procTop.nCoords[p][1])&&//on inner edge in y
                ( (procTop.nProcDims[2]-1       )==procTop.nCoords[p][2])  //on inner edge in z
              ){//current processor allowed as periodic neighbor
              nNumPeriodicNeighbors++;
            }
          }
        }
      }
    }
    if(procTop.nPeriodic[0]&&procTop.nPeriodic[1]&&procTop.nPeriodic[2]){//if periodic in x, y & z -directions
      //periodicity not needed in x-direction
    }
    
    //set id's of non-periodic neighbors
    int *nNonPeriodicNeighborRanks=new int[nNumNonPeriodicNeighbors];
    int l=0;
    for(int i=0;i<3;i++){
    
      //in 3D region
      for(int j=0;j<3;j++){
        for(int k=0;k<3;k++){
          for(int p=1;p<procTop.nNumProcs;p++){
            if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0]) &&//matches x
                ( (procTop.nCoords[procTop.nRank][1]+j-1)==procTop.nCoords[p][1]) &&//matches y
                ( (procTop.nCoords[procTop.nRank][2]+k-1)==procTop.nCoords[p][2]) &&//mathces z
              p!=procTop.nRank ){//not the current processor
              nNonPeriodicNeighborRanks[l]=p;
              l++;
            }
          }
        }
      }
      
      //in 1D region
      for(int p=0;p<1;p++){
        if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
            (  procTop.nCoords[p][1]==-1)&&
            (  procTop.nCoords[p][2]==-1)&&
          p!=procTop.nRank ){//not the current processor
          nNonPeriodicNeighborRanks[l]=p;
          l++;
        }
      }
    }
    
    //set id's and edges of periodic neighbors
    int *nPeriodicNeighborRanks=new int[nNumPeriodicNeighbors];
    int **nEdge=new int*[nNumPeriodicNeighbors];//inner=0 or outer=1 edge
    l=0;
    if(procTop.nPeriodic[0]){//if periodic in x-direction
      //not needed, if at a later time it is needed, code will
      // have to be added here
    }
    if(procTop.nPeriodic[1]){//if periodic in y-direction
      if(procTop.nCoords[procTop.nRank][1]==(procTop.nProcDims[1]-1)){//on outer edge of processor grid
        
        //search for neighbors at inner edge
        for(int i=0;i<3;i++){//search x-direction
          for(int k=0;k<3;k++){//search z-direction
            for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
              if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                  ( (0                    )==procTop.nCoords[p][1])&&//on inner edge
                  ( (procTop.nCoords[procTop.nRank][2]+k-1)==procTop.nCoords[p][2])  //mathces z
                ){//current processor allowed as periodic neighbor
                  nPeriodicNeighborRanks[l]=p;
                  nEdge[l]=new int[3];
                  nEdge[l][0]=-1;
                  nEdge[l][1]=0;
                  nEdge[l][2]=-1;
                  l++;
              }
            }
          }
        }
      }
      if(procTop.nCoords[procTop.nRank][1]==0){//on inner edge of processor grid
        
        //search for neighbors at outter edge
        for(int i=0;i<3;i++){//search x-direction
          for(int k=0;k<3;k++){//search z-direction
            for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
              if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                  ( (procTop.nProcDims[1]-1       )==procTop.nCoords[p][1])&&//on inner edge
                  ( (procTop.nCoords[procTop.nRank][2]+k-1)==procTop.nCoords[p][2])  //mathces z
                ){//current processor allowed as periodic neighbor
                  nPeriodicNeighborRanks[l]=p;
                  nEdge[l]=new int[3];
                  nEdge[l][0]=-1;
                  nEdge[l][1]=1;
                  nEdge[l][2]=-1;
                  l++;
              }
            }
          }
        }
      }
    }
    if(procTop.nPeriodic[2]){//if periodic in z-direction
      if(procTop.nCoords[procTop.nRank][2]==(procTop.nProcDims[2]-1)){//on outer edge of processor grid
        
        //search for neighbors at inner edge
        for(int i=0;i<3;i++){//search x-direction
          for(int j=0;j<3;j++){//search y-direction
            for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
              if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                  ( (procTop.nCoords[procTop.nRank][1]+j-1)==procTop.nCoords[p][1])&&//matches y
                  ( (0                    )==procTop.nCoords[p][2])  //on inner edge
                ){//current processor allowed as periodic neighbor
                  nPeriodicNeighborRanks[l]=p;
                  nEdge[l]=new int[3];
                  nEdge[l][0]=-1;
                  nEdge[l][1]=-1;
                  nEdge[l][2]=0;
                  l++;
              }
            }
          }
        }
      }
      if(procTop.nCoords[procTop.nRank][2]==0){//on inner edge of processor grid
      
        //search for neighbors at outter edge
        for(int i=0;i<3;i++){//search x-direction
          for(int j=0;j<3;j++){//search y-direction
            for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
              if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                  ( (procTop.nCoords[procTop.nRank][1]+j-1)==procTop.nCoords[p][1])&&//matches y
                  ( (procTop.nProcDims[2]-1       )==procTop.nCoords[p][2])  //on outer edge
                ){//current processor allowed as periodic neighbor
                  nPeriodicNeighborRanks[l]=p;
                  nEdge[l]=new int[3];
                  nEdge[l][0]=-1;
                  nEdge[l][1]=-1;
                  nEdge[l][2]=1;
                  l++;
              }
            }
          }
        }
      }
    }
    if(procTop.nPeriodic[0]&&procTop.nPeriodic[1]){//if periodic in x & y -directions
      //periodicity not needed in x-direction
    }
    if(procTop.nPeriodic[0]&&procTop.nPeriodic[2]){//if periodic in x & z -directions
      //periodicity not needed in x-direction
    }
    if(procTop.nPeriodic[1]&&procTop.nPeriodic[2]){//if periodic in y & z -directions
      if( (procTop.nCoords[procTop.nRank][1]==(procTop.nProcDims[1]-1))
        &&(procTop.nCoords[procTop.nRank][2]==(procTop.nProcDims[2]-1)) ){//on outer corner in y and z direction
        
        //search for neighbors at corner 
        for(int i=0;i<3;i++){//search x-direction
          for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
            if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                ( (0                    )==procTop.nCoords[p][1])&&//on inner edge in y
                ( (0                    )==procTop.nCoords[p][2])  //on inner edge in z
              ){//current processor allowed as periodic neighbor
              nPeriodicNeighborRanks[l]=p;
              nEdge[l]=new int[3];
              nEdge[l][0]=-1;
              nEdge[l][1]=0;
              nEdge[l][2]=0;
              l++;
            }
          }
        }
      }
      if( (procTop.nCoords[procTop.nRank][1]==(0             ))
        &&(procTop.nCoords[procTop.nRank][2]==(procTop.nProcDims[2]-1)) ){//on inner corner in y, and outer corner in z
                                                  //direction
        //search for neighbors at corner 
        for(int i=0;i<3;i++){//search x-direction
          for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
            if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                ( (procTop.nProcDims[1]-1       )==procTop.nCoords[p][1])&&//on outer edge in y
                ( (0                    )==procTop.nCoords[p][2])  //on inner edge in z
              ){//current processor allowed as periodic neighbor
              nPeriodicNeighborRanks[l]=p;
              nEdge[l]=new int[3];
              nEdge[l][0]=-1;
              nEdge[l][1]=1;
              nEdge[l][2]=0;
              l++;
            }
          }
        }
      }
      if( (procTop.nCoords[procTop.nRank][1]==(procTop.nProcDims[1]-1))
        &&(procTop.nCoords[procTop.nRank][2]==(0             )) ){//on outer corner in y, and inner corner in z 
                                                  //direction
        //search for neighbors at corner 
        for(int i=0;i<3;i++){//search x-direction
          for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
            if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                ( (0                    )==procTop.nCoords[p][1])&&//on inner edge in y
                ( (procTop.nProcDims[2]-1       )==procTop.nCoords[p][2])  //on outer edge in z
              ){//current processor allowed as periodic neighbor
              nPeriodicNeighborRanks[l]=p;
              nEdge[l]=new int[3];
              nEdge[l][0]=-1;
              nEdge[l][1]=0;
              nEdge[l][2]=1;
              l++;
            }
          }
        }
      }
      if( (procTop.nCoords[procTop.nRank][1]==(0             ))
        &&(procTop.nCoords[procTop.nRank][2]==(0             )) ){//on inner corner in y, and inner corner in z 
                                                  //direction
        //search for neighbors at corner 
        for(int i=0;i<3;i++){//search x-direction
          for(int p=1;p<procTop.nNumProcs;p++){//processor 0 already handeled
            if( ( (procTop.nCoords[procTop.nRank][0]+i-1)==procTop.nCoords[p][0])&&//matches x
                ( (procTop.nProcDims[1]-1       )==procTop.nCoords[p][1])&&//on outer edge in y
                ( (procTop.nProcDims[2]-1       )==procTop.nCoords[p][2])  //on outer edge in z
              ){//current processor allowed as periodic neighbor
              nPeriodicNeighborRanks[l]=p;
              nEdge[l]=new int[3];
              nEdge[l][0]=-1;
              nEdge[l][1]=1;
              nEdge[l][2]=1;
              l++;
            }
          }
        }
      }
    }
    if(procTop.nPeriodic[0]&&procTop.nPeriodic[1]&&procTop.nPeriodic[2]){//if periodic in x, y & z -directions
      //periodicity not needed in x-direction
    }
    
    //find total number of unique neighbors
    procTop.nNumNeighbors=nNumNonPeriodicNeighbors;
    for(int i=0;i<nNumPeriodicNeighbors;i++){
      bool bIsUnique=true;
      
      //already accounted for in non-periodic neighbors
      for(int j=0;j<nNumNonPeriodicNeighbors;j++){
        if(nPeriodicNeighborRanks[i]==nNonPeriodicNeighborRanks[j]){
          bIsUnique=false;//found one already in non-periodic nieghbors
          break;
        }
      }
      
      //already counted in periodic neighbors
      for(int j=0;j<i;j++){
        if(nPeriodicNeighborRanks[i]==nPeriodicNeighborRanks[j]){
          bIsUnique=false;//already account for once in periodic neighbors
          break;
        }
      }
      if(bIsUnique){
        procTop.nNumNeighbors++;
      }
    }
    
    //set ranks of unique neighbors
    procTop.nNeighborRanks=new int[procTop.nNumNeighbors];
    int nIndex=0;
    for(int i=0;i<nNumNonPeriodicNeighbors;i++){
      procTop.nNeighborRanks[nIndex]=nNonPeriodicNeighborRanks[i];
      nIndex++;
    }
    for(int i=0;i<nNumPeriodicNeighbors;i++){
      bool bIsUnique=true;
      
      //already accounted for in non-periodic neighbors
      for(int j=0;j<nNumNonPeriodicNeighbors;j++){
        if(nPeriodicNeighborRanks[i]==nNonPeriodicNeighborRanks[j]){
          bIsUnique=false;//found one already in non-periodic nieghbors
          break;
        }
      }
      
      //already counted in periodic neighbors
      for(int j=0;j<i;j++){
        if(nPeriodicNeighborRanks[i]==nPeriodicNeighborRanks[j]){
          bIsUnique=false;//already account for once in periodic neighbors
          break;
        }
      }
      if(bIsUnique){
        procTop.nNeighborRanks[nIndex]=nPeriodicNeighborRanks[i];
        nIndex++;
      }
    }
    
    //allocate memory for send and recive types
    messPass.typeSendNewGrid=new MPI::Datatype[procTop.nNumNeighbors];
    messPass.typeRecvOldGrid=new MPI::Datatype[procTop.nNumNeighbors];
    messPass.typeSendNewVar=new MPI::Datatype*[procTop.nNumNeighbors];
    messPass.typeRecvNewVar=new MPI::Datatype*[procTop.nNumNeighbors];
    for(int i=0;i<procTop.nNumNeighbors;i++){//allocate space for MPI::Datatypes for each variable
      messPass.typeRecvNewVar[i]=new MPI::Datatype[grid.nNumVars+grid.nNumIntVars];
      messPass.typeSendNewVar[i]=new MPI::Datatype[grid.nNumVars+grid.nNumIntVars];
    }

    //set send and recieve types
    for(int n=0;n<procTop.nNumNeighbors;n++){
      
      //find number of blocks to send and recieve, they are the same
      std::vector<int> nBlockTypes;//0=non-Periodic, 1=periodic
      for(int j=0;j<nNumNonPeriodicNeighbors;j++){//check non-periodic neighbors
        if(procTop.nNeighborRanks[n]==nNonPeriodicNeighborRanks[j]){
          nBlockTypes.push_back(0);
        }
      }
      for(int j=0;j<nNumPeriodicNeighbors;j++){//check periodic neighbors
        if(procTop.nNeighborRanks[n]==nPeriodicNeighborRanks[j]){
          nBlockTypes.push_back(1);
        }
      }
      
      //allocate memory to hold send and recieve block info
      int ***nSendBlockStart=new int**[nBlockTypes.size()];
      int ***nRecvBlockStart=new int**[nBlockTypes.size()];
      int ***nSendBlockDims=new int**[nBlockTypes.size()];
      int ***nRecvBlockDims=new int**[nBlockTypes.size()];
      for(unsigned int j=0;j<nBlockTypes.size();j++){
        nSendBlockStart[j]=new int*[grid.nNumVars+grid.nNumIntVars];
        nRecvBlockStart[j]=new int*[grid.nNumVars+grid.nNumIntVars];
        nSendBlockDims[j]=new int*[grid.nNumVars+grid.nNumIntVars];
        nRecvBlockDims[j]=new int*[grid.nNumVars+grid.nNumIntVars];
        for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
          nSendBlockDims[j][i]=new int[3];
          nRecvBlockDims[j][i]=new int[3];
          nSendBlockStart[j][i]=new int[3];
          nRecvBlockStart[j][i]=new int[3];
        }
      }
      if(procTop.nNeighborRanks[n]==0){//if 1D neighbor
        
        for(unsigned int j=0;j<nBlockTypes.size();j++){
          
          //set start of blocks
          for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
            
            for(int l=0;l<3;l++){
            
              //recv blocks
              nRecvBlockStart[j][i][l]=0;
              
              //send blocks
              if(procTop.nCoords[procTop.nNeighborRanks[n]][l]<=procTop.nCoords[procTop.nRank][l]
                ||procTop.nCoords[procTop.nNeighborRanks[n]][l]==-1){
                nSendBlockStart[j][i][l]=grid.nNumGhostCells;
              }
              else{
                nSendBlockStart[j][i][l]=grid.nLocalGridDims[procTop.nRank][i][l];
              }
            }
            
            for(int l=0;l<3;l++){
              if(procTop.nCoords[procTop.nNeighborRanks[n]][l]==procTop.nCoords[procTop.nRank][l]
                ||procTop.nCoords[procTop.nNeighborRanks[n]][l]==-1){
                nSendBlockDims[j][i][l]=grid.nLocalGridDims[procTop.nRank][i][l];
                nRecvBlockDims[j][i][l]=grid.nLocalGridDims[procTop.nRank][i][l];
              }
              else{
                nSendBlockDims[j][i][l]=grid.nNumGhostCells;//same for recv and send
                nRecvBlockDims[j][i][l]=grid.nNumGhostCells;//same for recv and send
              }
              if(l!=0){//if not radial direction, then add ghost cells to the recieve block
                nRecvBlockDims[j][i][l]+=2*grid.nNumGhostCells;
              }
            }
            
            //adjust for non-3D variables
            for(int l=0;l<3;l++){
              
              //if variable not defined in direction l
              if(grid.nVariables[i][l]==-1){
                
                //set to 0 since there is no other cells in direction l for variable i
                nSendBlockStart[j][i][l]=0;
                nRecvBlockStart[j][i][l]=0;
                
                //set to 1 in this direction
                nSendBlockDims[j][i][l]=1;
                nRecvBlockDims[j][i][l]=1;
                
              }
            }
            
            //adjust for non-time dependent variables
            if(grid.nVariables[i][3]==0){//don't need to update in time
              for(int l=0;l<3;l++){
                nSendBlockDims[j][i][l]=0;//set all block dims to 0 (nothing to be sent)
                nRecvBlockDims[j][i][l]=0;
              }
            }
            
            //check to see if we need to send current neighbor a message for current variable
            for(int l=0;l<1;l++){//only need to test x-direction
              
              //if variable not defined in direction l
              if(grid.nVariables[i][l]==-1){
                
                //if current neighbor is not at the same coordinate in direction l
                if(procTop.nCoords[procTop.nNeighborRanks[n]][l]!=procTop.nCoords[procTop.nRank][l]
                  &&procTop.nCoords[procTop.nRank][l]!=-1){
                  for(int l=0;l<3;l++){//all block dimensions set to 0, since nothing to send
                    nSendBlockDims[j][i][l]=0;
                    nRecvBlockDims[j][i][l]=0;
                  }
                  break;//exit testing for current neighbor
                }
              }
            }
          }
        }
        
        //set sub-block lengths and types for send
        int nNumSubBlocksSend=0;
        for(unsigned int j=0;j<nBlockTypes.size();j++){
          for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
            nNumSubBlocksSend+=nSendBlockDims[j][i][0]*nSendBlockDims[j][i][1]
              *nSendBlockDims[j][i][2];
          }
        }
        
        int* nSendBlockLen=new int[nNumSubBlocksSend];
        MPI::Datatype *typeBaseSend=new MPI::Datatype[nNumSubBlocksSend];
        int nIter=0;
        std::vector<std::vector<int> > vecnBlockLenSendNewVar;
        std::vector<std::vector<MPI::Datatype> > vectypeBaseSendNewVar;
        for(unsigned int p=0;p<nBlockTypes.size();p++){
          for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
            vecnBlockLenSendNewVar.push_back(std::vector<int>());
            vectypeBaseSendNewVar.push_back(std::vector<MPI::Datatype>());
            for(int j=0;j<nSendBlockDims[p][i][0];j++){
              for(int k=0;k<nSendBlockDims[p][i][1];k++){
                for(int l=0;l<nSendBlockDims[p][i][2];l++){
                  nSendBlockLen[nIter]=1;
                  vecnBlockLenSendNewVar[i].push_back(1);
                  typeBaseSend[nIter]=MPI::DOUBLE;
                  vectypeBaseSendNewVar[i].push_back(MPI::DOUBLE);
                  nIter++;
                }
              }
            }
          }
        }
        
        //set starting send address
        MPI::Aint nStartAddressSend=MPI::Get_address(grid.dLocalGridNew);
        
        //set addresses for send
        MPI::Aint *nSendAddresses=new MPI::Aint[nNumSubBlocksSend];
        int nCount=0;
        std::vector<std::vector<MPI::Aint> > vecSendNewVarAddresses;
        for(unsigned int p=0;p<nBlockTypes.size();p++){
          for(int l=0;l<grid.nNumVars+grid.nNumIntVars;l++){
            vecSendNewVarAddresses.push_back(std::vector<MPI::Aint>() );
            for(int i=nSendBlockStart[p][l][0];i<nSendBlockStart[p][l][0]+nSendBlockDims[p][l][0]
              ;i++){
              for(int j=nSendBlockStart[p][l][1];j<nSendBlockStart[p][l][1]+nSendBlockDims[p][l][1]
                ;j++){
                for(int k=nSendBlockStart[p][l][2];k<nSendBlockStart[p][l][2]
                  +nSendBlockDims[p][l][2];k++){
                  MPI::Aint nCurAddress=MPI::Get_address(&grid.dLocalGridNew[l][i][j][k]);
                  nSendAddresses[nCount]=nCurAddress-nStartAddressSend;
                  nCount++;
                  vecSendNewVarAddresses[l].push_back(nCurAddress-nStartAddressSend);
                }
              }
            }
          }
        }
        //set sub-block lengths and types for recv
        int nNumSubBlocksRecv=0;
        for(unsigned int j=0;j<nBlockTypes.size();j++){
          for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
            nNumSubBlocksRecv+=nRecvBlockDims[j][i][0]*nRecvBlockDims[j][i][1]
              *nRecvBlockDims[j][i][2];
          }
        }
        int* nRecvBlockLen=new int[nNumSubBlocksRecv];
        MPI::Datatype *typeBaseRecv=new MPI::Datatype[nNumSubBlocksRecv];
        nIter=0;
        std::vector<std::vector<int> > vecnBlockLenRecvNewVar;
        std::vector<std::vector<MPI::Datatype> > vectypeBaseRecvNewVar;
        for(unsigned int p=0;p<nBlockTypes.size();p++){
          for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
            vecnBlockLenRecvNewVar.push_back(std::vector<int>());
            vectypeBaseRecvNewVar.push_back(std::vector<MPI::Datatype>());
            for(int j=0;j<nRecvBlockDims[p][i][0];j++){
              for(int k=0;k<nRecvBlockDims[p][i][1];k++){
                for(int l=0;l<nRecvBlockDims[p][i][2];l++){
                  nRecvBlockLen[nIter]=1;
                  vecnBlockLenRecvNewVar[i].push_back(1);
                  typeBaseRecv[nIter]=MPI::DOUBLE;
                  vectypeBaseRecvNewVar[i].push_back(MPI::DOUBLE);
                  nIter++;
                }
              }
            }
          }
        }
        
        //set starting recv address
        MPI::Aint nStartAddressRecv=MPI::Get_address(grid.dLocalGridOld);
        MPI::Aint nStartAddressRecv2=MPI::Get_address(grid.dLocalGridNew);
        
        //set addresses for recv
        MPI::Aint *nRecvAddresses=new MPI::Aint[nNumSubBlocksRecv];
        nCount=0;
        std::vector<std::vector<MPI::Aint> > vecRecvNewVarAddresses;
        for(unsigned int p=0;p<nBlockTypes.size();p++){
          for(int l=0;l<grid.nNumVars+grid.nNumIntVars;l++){
            vecRecvNewVarAddresses.push_back(std::vector<MPI::Aint>() );
            for(int i=nRecvBlockStart[p][l][0];i<nRecvBlockStart[p][l][0]+nRecvBlockDims[p][l][0];
              i++){
              for(int j=nRecvBlockStart[p][l][1];j<nRecvBlockStart[p][l][1]+nRecvBlockDims[p][l][1];
                j++){
                for(int k=nRecvBlockStart[p][l][2];k<nRecvBlockStart[p][l][2]
                  +nRecvBlockDims[p][l][2];k++){
                  MPI::Aint nCurAddress=MPI::Get_address(&grid.dLocalGridOld[l][i][j][k]);
                  MPI::Aint nCurAddress2=MPI::Get_address(&grid.dLocalGridNew[l][i][j][k]);
                  nRecvAddresses[nCount]=nCurAddress-nStartAddressRecv;
                  nCount++;
                  vecRecvNewVarAddresses[l].push_back(nCurAddress2-nStartAddressRecv2);
                }
              }
            }
          }
        }
        
        //set send type
        messPass.typeSendNewGrid[n]=MPI::Datatype::Create_struct(nNumSubBlocksSend,nSendBlockLen
          ,nSendAddresses,typeBaseSend);
        messPass.typeSendNewGrid[n].Commit();//note: must delete datatypes
        
        //set send types for new vars
        for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
          
          //copy vecnBlockLenSendNewVar[n] to 1D array
          int *nBlockLen=new int[vecnBlockLenSendNewVar[i].size()];
          std::copy(vecnBlockLenSendNewVar[i].begin(),vecnBlockLenSendNewVar[i].end(),nBlockLen);
          
          //copy vectypeBaseSendNewVar[n] to 1D array
          MPI::Datatype *nBaseType=new MPI::Datatype[vectypeBaseSendNewVar[i].size()];
          std::copy(vectypeBaseSendNewVar[i].begin(),vectypeBaseSendNewVar[i].end(),nBaseType);
          
          //copy vecSendNewVarAddresses[i] to 1D array
          MPI::Aint *nAddresses=new MPI::Aint[vecSendNewVarAddresses[i].size()];
          std::copy(vecSendNewVarAddresses[i].begin(),vecSendNewVarAddresses[i].end(),nAddresses);
          
          //set send type
          messPass.typeSendNewVar[n][i]=MPI::Datatype::Create_struct(
            vecSendNewVarAddresses[i].size(),nBlockLen,nAddresses,nBaseType);
          messPass.typeSendNewVar[n][i].Commit();
          
          delete [] nBlockLen;
          delete [] nBaseType;
          delete [] nAddresses;
        }
        
        //set recv type
        messPass.typeRecvOldGrid[n]=MPI::Datatype::Create_struct(nNumSubBlocksRecv,nRecvBlockLen
          ,nRecvAddresses,typeBaseRecv);
        messPass.typeRecvOldGrid[n].Commit();//note: must delete datatypes
        
        //set recv types for new vars
        for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
          
          //copy vecnBlockLenSendNewVar[n] to 1D array
          int *nBlockLen=new int[vecnBlockLenRecvNewVar[i].size()];
          std::copy(vecnBlockLenRecvNewVar[i].begin(),vecnBlockLenRecvNewVar[i].end(),nBlockLen);
          
          //copy vectypeBaseSendNewVar[n] to 1D array
          MPI::Datatype *nBaseType=new MPI::Datatype[vectypeBaseRecvNewVar[i].size()];
          std::copy(vectypeBaseRecvNewVar[i].begin(),vectypeBaseRecvNewVar[i].end(),nBaseType);
          
          //copy vecSendNewVarAddresses[i] to 1D array
          MPI::Aint *nAddresses=new MPI::Aint[vecRecvNewVarAddresses[i].size()];
          std::copy(vecRecvNewVarAddresses[i].begin(),vecRecvNewVarAddresses[i].end(),nAddresses);
          
          //set send type
          messPass.typeRecvNewVar[n][i]=MPI::Datatype::Create_struct(
            vecRecvNewVarAddresses[i].size(),nBlockLen,nAddresses,nBaseType);
          messPass.typeRecvNewVar[n][i].Commit();
          
          delete [] nBlockLen;
          delete [] nBaseType;
          delete [] nAddresses;
        }
        
        delete [] nRecvBlockLen;
        delete [] nSendBlockLen;
        delete [] nSendAddresses;
        delete [] nRecvAddresses;
        delete [] typeBaseRecv;
        delete [] typeBaseSend;
      }
      else{//if normal neighbor
        
        //set send and recv block start and dimensions
        int nPeriodicBlock=0;
        for(unsigned int p=0;p<nBlockTypes.size();p++){
          if(nBlockTypes[p]==0){//not-periodic
            for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
              
              //set send and recv block stars
              for(int l=0;l<3;l++){
              
                //recv blocks start
                if(procTop.nCoords[procTop.nNeighborRanks[n]][l]<procTop.nCoords[procTop.nRank][l]){
                  nRecvBlockStart[p][i][l]=0;
                }
                else if(procTop.nCoords[procTop.nNeighborRanks[n]][l]==
                  procTop.nCoords[procTop.nRank][l]){
                  nRecvBlockStart[p][i][l]=grid.nNumGhostCells;
                }
                else{//greater than procTop.nRank
                  nRecvBlockStart[p][i][l]=grid.nLocalGridDims[procTop.nRank][i][l]
                    +grid.nNumGhostCells;
                }
                
                //send blocks start
                if(procTop.nCoords[procTop.nNeighborRanks[n]][l]
                  <=procTop.nCoords[procTop.nRank][l]){
                  nSendBlockStart[p][i][l]=grid.nNumGhostCells;
                }
                else{
                  nSendBlockStart[p][i][l]=grid.nLocalGridDims[procTop.nRank][i][l];
                }
              }
              
              //set block dimensions
              for(int l=0;l<3;l++){
                if(
                  procTop.nCoords[procTop.nNeighborRanks[n]][l]==procTop.nCoords[procTop.nRank][l]){
                  nSendBlockDims[p][i][l]=grid.nLocalGridDims[procTop.nRank][i][l];
                  nRecvBlockDims[p][i][l]=grid.nLocalGridDims[procTop.nRank][i][l];
                  if( (procTop.nCoords[procTop.nRank][l]==procTop.nProcDims[l]-1)
                    && (procTop.nPeriodic[l]==0) ){//if last processor in direction l,
                                                   // and direction l is not periodic
                    nSendBlockDims[p][i][l]+=grid.nNumGhostCells;
                    nRecvBlockDims[p][i][l]+=grid.nNumGhostCells;
                  }
                }
                else{
                  nSendBlockDims[p][i][l]=grid.nNumGhostCells;//same for recv and send
                  nRecvBlockDims[p][i][l]=grid.nNumGhostCells;//same for recv and send
                }
              }
              
              //adjust for non-3D variables
              for(int l=0;l<3;l++){
                
                //if variable not defined in direction l
                if(grid.nVariables[i][l]==-1){
                  
                  //set to 0 since there is no other cells in direction l for variable i
                  nSendBlockStart[p][i][l]=0;
                  nRecvBlockStart[p][i][l]=0;
                  
                  //set to 1 in this direction
                  nSendBlockDims[p][i][l]=1;
                  nRecvBlockDims[p][i][l]=1;
                  
                }
              }
              
              //adjust for non-time dependent variables
              if(grid.nVariables[i][3]==0){//don't need to update in time
                for(int l=0;l<3;l++){
                  nSendBlockDims[p][i][l]=0;//set all block dims to 0 (nothing to be sent)
                  nRecvBlockDims[p][i][l]=0;
                }
              }
              
              //check to see if we need to send current neighbor a message for current variable
              for(int l=0;l<3;l++){
                //if variable not defined in direction l
                if(grid.nVariables[i][l]==-1){
                
                  //if current neighbor is not at the same coordinate in direction l
                  if(
                    procTop.nCoords[procTop.nNeighborRanks[n]][l]!=procTop.nCoords[procTop.nRank][l]
                    &&procTop.nCoords[procTop.nNeighborRanks[n]][l]!=-1){
                    for(int l=0;l<3;l++){//all block dimensions set to 0, since nothing to send
                      nSendBlockDims[p][i][l]=0;
                      nRecvBlockDims[p][i][l]=0;
                    }
                    break;//exit testing for current neighbor
                  }
                }
              }
            }
          }
          else{//periodic neighbor
            
            //get edge index of block
            int nIndex;
            int nCount=0;
            for(int i=0;i<nNumPeriodicNeighbors;i++){
              if(procTop.nNeighborRanks[n]==nPeriodicNeighborRanks[i]){
                nIndex=i;
                if(nCount==nPeriodicBlock){
                  break;//if we are at the correct periodic block
                }
                nCount++;
              }
            }
            for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
              
              //set initial send and recv start
              for(int l=0;l<3;l++){
                
                //recv blocks start
                if(procTop.nCoords[procTop.nNeighborRanks[n]][l]<procTop.nCoords[procTop.nRank][l]){
                  nRecvBlockStart[p][i][l]=0;
                }
                else if(procTop.nCoords[procTop.nNeighborRanks[n]][l]
                  ==procTop.nCoords[procTop.nRank][l]){
                  nRecvBlockStart[p][i][l]=grid.nNumGhostCells;
                }
                else{//greater than procTop.nRank
                  nRecvBlockStart[p][i][l]=grid.nLocalGridDims[procTop.nRank][i][l]
                    +grid.nNumGhostCells;
                }
                
                //send blocks start
                if(procTop.nCoords[procTop.nNeighborRanks[n]][l]
                  <=procTop.nCoords[procTop.nRank][l]){
                  nSendBlockStart[p][i][l]=grid.nNumGhostCells;
                }
                else{
                  nSendBlockStart[p][i][l]=grid.nLocalGridDims[procTop.nRank][i][l];
                }
              }
              
              //set initial block dimensions
              for(int l=0;l<3;l++){//need to make adjustments for periodic neighbors
                if(procTop.nCoords[procTop.nNeighborRanks[n]][l]
                  ==procTop.nCoords[procTop.nRank][l]){
                  nSendBlockDims[p][i][l]=grid.nLocalGridDims[procTop.nRank][i][l];
                  nRecvBlockDims[p][i][l]=grid.nLocalGridDims[procTop.nRank][i][l];
                  if( (procTop.nCoords[procTop.nRank][l]==procTop.nProcDims[l]-1)
                    && (procTop.nPeriodic[l]==0) ){//if last processor in direction l,
                                                   // and direction l is not periodic
                    nSendBlockDims[p][i][l]+=grid.nNumGhostCells;
                    nRecvBlockDims[p][i][l]+=grid.nNumGhostCells;
                  }
                }
                else{
                  nSendBlockDims[p][i][l]=grid.nNumGhostCells;//same for recv and send
                  nRecvBlockDims[p][i][l]=grid.nNumGhostCells;//same for recv and send
                }
              }
              
              if(nEdge[nIndex][0]==-1&&nEdge[nIndex][1]==0&&nEdge[nIndex][2]==-1){//neighbor at y inner edge
                
                //recv blocks start
                nRecvBlockStart[p][i][1]=grid.nLocalGridDims[procTop.nRank][i][1]
                  +grid.nNumGhostCells;
                
                //send blocks start
                nSendBlockStart[p][i][1]=grid.nLocalGridDims[procTop.nRank][i][1];
                
                //set send/recv block dimensions
                nSendBlockDims[p][i][1]=grid.nNumGhostCells;
                nRecvBlockDims[p][i][1]=grid.nNumGhostCells;
                
                //required when only 1 processor wide in y-direction
                if(procTop.nCoords[procTop.nNeighborRanks[n]][1]
                  ==procTop.nCoords[procTop.nRank][1]){
                  nSendBlockStart[p][i][1]=grid.nNumGhostCells;
                }
              }
              
              if(nEdge[nIndex][0]==-1&&nEdge[nIndex][1]==1&&nEdge[nIndex][2]==-1){//neighbor at y outter edge
                //recv blocks start
                nRecvBlockStart[p][i][1]=0;//j=block #, i=variable #
                
                //send blocks start
                nSendBlockStart[p][i][1]=grid.nNumGhostCells;
                
                //set send/recv block dimensions
                nSendBlockDims[p][i][1]=grid.nNumGhostCells;
                nRecvBlockDims[p][i][1]=grid.nNumGhostCells;
                
                //required when only 1 processor wide in y-direction
                if(procTop.nCoords[procTop.nNeighborRanks[n]][1]
                  ==procTop.nCoords[procTop.nRank][1]){
                  nSendBlockStart[p][i][1]=grid.nLocalGridDims[procTop.nRank][i][1];
                }
              }
              if(nEdge[nIndex][0]==-1&&nEdge[nIndex][1]==-1&&nEdge[nIndex][2]==0){//neighbor at z inner edge
                
                //recv blocks start
                nRecvBlockStart[p][i][2]=grid.nLocalGridDims[procTop.nRank][i][2]+grid.nNumGhostCells;
                
                //send blocks start
                nSendBlockStart[p][i][2]=grid.nLocalGridDims[procTop.nRank][i][2];
                
                //set send/recv block dimensions
                nSendBlockDims[p][i][2]=grid.nNumGhostCells;
                nRecvBlockDims[p][i][2]=grid.nNumGhostCells;
                
                //required when only 1 processor wide in z-direction
                if(
                  procTop.nCoords[procTop.nNeighborRanks[n]][2]==procTop.nCoords[procTop.nRank][2]){
                  nSendBlockStart[p][i][2]=grid.nNumGhostCells;
                }
              }
              if(nEdge[nIndex][0]==-1&&nEdge[nIndex][1]==-1&&nEdge[nIndex][2]==1){//neighbor at z outter edge
                
                //recv blocks start
                nRecvBlockStart[p][i][2]=0;
                
                //send blocks start
                nSendBlockStart[p][i][2]=grid.nNumGhostCells;
                
                //set send/recv block dimensions
                nSendBlockDims[p][i][2]=grid.nNumGhostCells;
                nRecvBlockDims[p][i][2]=grid.nNumGhostCells;
                
                //required when only 1 processor wide in z-direction
                if(
                  procTop.nCoords[procTop.nNeighborRanks[n]][2]==procTop.nCoords[procTop.nRank][2]){
                  nSendBlockStart[p][i][2]=grid.nLocalGridDims[procTop.nRank][i][2];
                }
              }
              if(nEdge[nIndex][0]==-1&&nEdge[nIndex][1]==0&&nEdge[nIndex][2]==0){//neighbor at corner, inner y, inner z
                
                //recv blocks start
                nRecvBlockStart[p][i][1]=grid.nLocalGridDims[procTop.nRank][i][1]
                  +grid.nNumGhostCells;
                nRecvBlockStart[p][i][2]=grid.nLocalGridDims[procTop.nRank][i][2]
                  +grid.nNumGhostCells;
                
                //send blocks start
                nSendBlockStart[p][i][1]=grid.nLocalGridDims[procTop.nRank][i][1];
                nSendBlockStart[p][i][2]=grid.nLocalGridDims[procTop.nRank][i][2];
                
                //set send/recv block dimensions
                nSendBlockDims[p][i][1]=grid.nNumGhostCells;
                nRecvBlockDims[p][i][1]=grid.nNumGhostCells;
                nSendBlockDims[p][i][2]=grid.nNumGhostCells;
                nRecvBlockDims[p][i][2]=grid.nNumGhostCells;
                
                //required when only 1 processor wide in y-direction
                if(
                  procTop.nCoords[procTop.nNeighborRanks[n]][1]==procTop.nCoords[procTop.nRank][1]){
                  nSendBlockStart[p][i][1]=grid.nNumGhostCells;
                }
                //required when only 1 processor wide in z-direction
                if(
                  procTop.nCoords[procTop.nNeighborRanks[n]][2]==procTop.nCoords[procTop.nRank][2]){
                  nSendBlockStart[p][i][2]=grid.nNumGhostCells;
                }
              }
              if(nEdge[nIndex][0]==-1&&nEdge[nIndex][1]==1&&nEdge[nIndex][2]==0){//neighbor at corner, outer y, inner z
                
                //recv blocks start
                nRecvBlockStart[p][i][1]=0;
                nRecvBlockStart[p][i][2]=grid.nLocalGridDims[procTop.nRank][i][2]
                  +grid.nNumGhostCells;
                
                //send blocks start
                nSendBlockStart[p][i][1]=grid.nNumGhostCells;
                nSendBlockStart[p][i][2]=grid.nLocalGridDims[procTop.nRank][i][2];
                
                //set send/recv block dimensions
                nSendBlockDims[p][i][1]=grid.nNumGhostCells;
                nRecvBlockDims[p][i][1]=grid.nNumGhostCells;
                nSendBlockDims[p][i][2]=grid.nNumGhostCells;
                nRecvBlockDims[p][i][2]=grid.nNumGhostCells;
                
                //required when only 1 processor wide in y-direction
                if(
                  procTop.nCoords[procTop.nNeighborRanks[n]][1]==procTop.nCoords[procTop.nRank][1]){
                  nSendBlockStart[p][i][1]=grid.nLocalGridDims[procTop.nRank][i][1];
                }
                //required when only 1 processor wide in z-direction
                if(
                  procTop.nCoords[procTop.nNeighborRanks[n]][2]==procTop.nCoords[procTop.nRank][2]){
                  nSendBlockStart[p][i][2]=grid.nNumGhostCells;
                }
              }
              if(nEdge[nIndex][0]==-1&&nEdge[nIndex][1]==1&&nEdge[nIndex][2]==1){//neighbor at corner, outer y, outer z
                
                //recv blocks start
                nRecvBlockStart[p][i][1]=0;
                nRecvBlockStart[p][i][2]=0;
                
                //send blocks start
                nSendBlockStart[p][i][1]=grid.nNumGhostCells;
                nSendBlockStart[p][i][2]=grid.nNumGhostCells;
                
                //set send/recv block dimensions
                nSendBlockDims[p][i][1]=grid.nNumGhostCells;
                nRecvBlockDims[p][i][1]=grid.nNumGhostCells;
                nSendBlockDims[p][i][2]=grid.nNumGhostCells;
                nRecvBlockDims[p][i][2]=grid.nNumGhostCells;
                
                //required when only 1 processor wide in y-direction
                if(
                  procTop.nCoords[procTop.nNeighborRanks[n]][1]==procTop.nCoords[procTop.nRank][1]){
                  nSendBlockStart[p][i][1]=grid.nLocalGridDims[procTop.nRank][i][1];
                }
                //required when only 1 processor wide in z-direction
                if(
                  procTop.nCoords[procTop.nNeighborRanks[n]][2]==procTop.nCoords[procTop.nRank][2]){
                  nSendBlockStart[p][i][2]=grid.nLocalGridDims[procTop.nRank][i][2];
                }
              }
              if(nEdge[nIndex][0]==-1&&nEdge[nIndex][1]==0&&nEdge[nIndex][2]==1){//neighbor at corner, inner y, outer z
                
                //recv blocks start
                nRecvBlockStart[p][i][1]=grid.nLocalGridDims[procTop.nRank][i][1]
                  +grid.nNumGhostCells;
                nRecvBlockStart[p][i][2]=0;
                
                //send blocks start
                nSendBlockStart[p][i][1]=grid.nLocalGridDims[procTop.nRank][i][1];
                nSendBlockStart[p][i][2]=grid.nNumGhostCells;
                
                //set send/recv block dimensions
                nSendBlockDims[p][i][1]=grid.nNumGhostCells;
                nRecvBlockDims[p][i][1]=grid.nNumGhostCells;
                nSendBlockDims[p][i][2]=grid.nNumGhostCells;
                nRecvBlockDims[p][i][2]=grid.nNumGhostCells;
                
                //required when only 1 processor wide in y-direction
                if(
                  procTop.nCoords[procTop.nNeighborRanks[n]][1]==procTop.nCoords[procTop.nRank][1]){
                  nSendBlockStart[p][i][1]=grid.nNumGhostCells;
                }
                //required when only 1 processor wide in z-direction
                if(
                  procTop.nCoords[procTop.nNeighborRanks[n]][2]==procTop.nCoords[procTop.nRank][2]){
                  nSendBlockStart[p][i][2]=grid.nLocalGridDims[procTop.nRank][i][2];
                }
              }
              
              //adjust for non-3D variables
              for(int l=0;l<3;l++){
                
                //if variable not defined in direction l
                if(grid.nVariables[i][l]==-1){
                  
                  //set to 0 since there is no other cells in direction l for variable i
                  nSendBlockStart[p][i][l]=0;
                  nRecvBlockStart[p][i][l]=0;
                  
                  //set to 1 in this direction
                  nSendBlockDims[p][i][l]=1;
                  nRecvBlockDims[p][i][l]=1;
                  
                }
              }
              
              //adjust for non-time dependent variables
              if(grid.nVariables[i][3]==0){//don't need to update in time
                for(int l=0;l<3;l++){
                  nSendBlockDims[p][i][l]=0;//set all block dims to 0 (nothing to be sent)
                  nRecvBlockDims[p][i][l]=0;
                }
              }
              
              //check to see if we need to send current neighbor a message for current variable
              for(int l=0;l<3;l++){
                
                //if variable not defined in direction l
                if(grid.nVariables[i][l]==-1){
                
                  //if current neighbor is not at the same coordinate in direction l
                  if(
                    procTop.nCoords[procTop.nNeighborRanks[n]][l]!=procTop.nCoords[procTop.nRank][l]
                    &&procTop.nCoords[procTop.nNeighborRanks[n]][l]!=-1){
                    for(int l=0;l<3;l++){//all block dimensions set to 0, since nothing to send
                      nSendBlockDims[p][i][l]=0;
                      nRecvBlockDims[p][i][l]=0;
                    }
                    break;//exit testing for current neighbor
                  }
                }
              }
            }
            nPeriodicBlock++;//do next periodic block, if there is another one
          }
        }
        
        //set sub-block lengths and types, same for send and recieve
        int nNumSubBlocks=0;
        for(unsigned int p=0;p<nBlockTypes.size();p++){
          for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
            nNumSubBlocks+=nSendBlockDims[p][i][0]*nSendBlockDims[p][i][1]*nSendBlockDims[p][i][2];
          }
        }
        int* nBlockLen=new int[nNumSubBlocks];
        MPI::Datatype *typeBase=new MPI::Datatype[nNumSubBlocks];
        int nIter=0;
        std::vector<std::vector<int> > vecnBlockLenNewVar;
        std::vector<std::vector<MPI::Datatype> > vectypeBaseNewVar;
        for(unsigned int p=0;p<nBlockTypes.size();p++){
          for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
            vecnBlockLenNewVar.push_back(std::vector<int>());
            vectypeBaseNewVar.push_back(std::vector<MPI::Datatype>());
            for(int j=0;j<nSendBlockDims[p][i][0];j++){
              for(int k=0;k<nSendBlockDims[p][i][1];k++){
                for(int l=0;l<nSendBlockDims[p][i][2];l++){
                  vecnBlockLenNewVar[i].push_back(1);
                  nBlockLen[nIter]=1;
                  typeBase[nIter]=MPI::DOUBLE;
                  vectypeBaseNewVar[i].push_back(MPI::DOUBLE);
                  nIter++;
                }
              }
            }
          }
        }
        
        //set starting send address
        MPI::Aint nStartAddressSend;
        nStartAddressSend=MPI::Get_address(grid.dLocalGridNew);
        
        //set addresses for send
        MPI::Aint *nSendAddresses=new MPI::Aint[nNumSubBlocks];
        int nCount=0;
        std::vector<std::vector<MPI::Aint> > vecSendNewVarAddresses;
        for(unsigned int p=0;p<nBlockTypes.size();p++){
          for(int l=0;l<grid.nNumVars+grid.nNumIntVars;l++){
            vecSendNewVarAddresses.push_back(std::vector<MPI::Aint>() );
            for(int i=nSendBlockStart[p][l][0];i<nSendBlockStart[p][l][0]+nSendBlockDims[p][l][0]
              ;i++){
              for(int j=nSendBlockStart[p][l][1];j<nSendBlockStart[p][l][1]+nSendBlockDims[p][l][1]
                ;j++){
                for(int k=nSendBlockStart[p][l][2];k<nSendBlockStart[p][l][2]
                  +nSendBlockDims[p][l][2];k++){
                  MPI::Aint nCurAddress=MPI::Get_address(&grid.dLocalGridNew[l][i][j][k]);
                  nSendAddresses[nCount]=nCurAddress-nStartAddressSend;
                  vecSendNewVarAddresses[l].push_back(nCurAddress-nStartAddressSend);
                  nCount++;
                }
              }
            }
          }
        }
        
        //set starting recv address
        MPI::Aint nStartAddressRecv=MPI::Get_address(grid.dLocalGridOld);
        MPI::Aint nStartAddressRecv2=MPI::Get_address(grid.dLocalGridNew);
        
        //set addresses for recv
        MPI::Aint *nRecvAddresses=new MPI::Aint[nNumSubBlocks];
        nCount=0;
        std::vector<std::vector<MPI::Aint> > vecRecvNewVarAddresses;
        for(unsigned int p=0;p<nBlockTypes.size();p++){
          for(int l=0;l<grid.nNumVars+grid.nNumIntVars;l++){
            vecRecvNewVarAddresses.push_back(std::vector<MPI::Aint> ());
            for(int i=nRecvBlockStart[p][l][0];i<nRecvBlockStart[p][l][0]+nSendBlockDims[p][l][0]
              ;i++){
              for(int j=nRecvBlockStart[p][l][1];j<nRecvBlockStart[p][l][1]+nSendBlockDims[p][l][1]
                ;j++){
                for(int k=nRecvBlockStart[p][l][2];k<nRecvBlockStart[p][l][2]
                  +nSendBlockDims[p][l][2];k++){
                  MPI::Aint nCurAddress=MPI::Get_address(&grid.dLocalGridOld[l][i][j][k]);
                  MPI::Aint nCurAddress2=MPI::Get_address(&grid.dLocalGridNew[l][i][j][k]);
                  nRecvAddresses[nCount]=nCurAddress-nStartAddressRecv;
                  vecRecvNewVarAddresses[l].push_back(nCurAddress2-nStartAddressRecv2);
                  nCount++;
                }
              }
            }
          }
        }
        
        //set send type
        messPass.typeSendNewGrid[n]=MPI::Datatype::Create_struct(nNumSubBlocks,nBlockLen
          ,nSendAddresses,typeBase);
        messPass.typeSendNewGrid[n].Commit();//note: must delete datatypes
        
        //set send types for new vars
        for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
          
          //copy vecnBlockLenSendNewVar[n] to 1D array
          int *nBlockLen=new int[vecnBlockLenNewVar[i].size()];
          std::copy(vecnBlockLenNewVar[i].begin(),vecnBlockLenNewVar[i].end(),nBlockLen);
          
          //copy vectypeBaseSendNewVar[n] to 1D array
          MPI::Datatype *nBaseType=new MPI::Datatype[vectypeBaseNewVar[i].size()];
          std::copy(vectypeBaseNewVar[i].begin(),vectypeBaseNewVar[i].end(),nBaseType);
          
          //copy vecSendNewVarAddresses[i] to 1D array
          MPI::Aint *nAddresses=new MPI::Aint[vecSendNewVarAddresses[i].size()];
          std::copy(vecSendNewVarAddresses[i].begin(),vecSendNewVarAddresses[i].end(),nAddresses);
          
          //set send type
          messPass.typeSendNewVar[n][i]=MPI::Datatype::Create_struct(
            vecSendNewVarAddresses[i].size(),nBlockLen,nAddresses,nBaseType);
          messPass.typeSendNewVar[n][i].Commit();
          
          delete [] nBlockLen;
          delete [] nBaseType;
          delete [] nAddresses;
        }
        //set recv type
        messPass.typeRecvOldGrid[n]=MPI::Datatype::Create_struct(nNumSubBlocks,nBlockLen
          ,nRecvAddresses,typeBase);
        messPass.typeRecvOldGrid[n].Commit();//note: must delete datatypes
        
        //set recv types for new vars
        for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
          
          //copy vecnBlockLenSendNewVar[n] to 1D array
          int *nBlockLen=new int[vecnBlockLenNewVar[i].size()];
          std::copy(vecnBlockLenNewVar[i].begin(),vecnBlockLenNewVar[i].end(),nBlockLen);
          
          //copy vectypeBaseSendNewVar[n] to 1D array
          MPI::Datatype *nBaseType=new MPI::Datatype[vectypeBaseNewVar[i].size()];
          std::copy(vectypeBaseNewVar[i].begin(),vectypeBaseNewVar[i].end(),nBaseType);
          
          //copy vecSendNewVarAddresses[i] to 1D array
          MPI::Aint *nAddresses=new MPI::Aint[vecRecvNewVarAddresses[i].size()];
          std::copy(vecRecvNewVarAddresses[i].begin(),vecRecvNewVarAddresses[i].end(),nAddresses);
          
          //set send type
          messPass.typeRecvNewVar[n][i]=MPI::Datatype::Create_struct(
            vecRecvNewVarAddresses[i].size(),nBlockLen,nAddresses,nBaseType);
          messPass.typeRecvNewVar[n][i].Commit();
          
          delete [] nBlockLen;
          delete [] nBaseType;
          delete [] nAddresses;
        }
        //delete dynamic memory used for current block
        delete [] nBlockLen;
        delete [] nSendAddresses;
        delete [] nRecvAddresses;
        delete [] typeBase;
      }
      
      //delete dynamic memory used for current neighbor
      for(unsigned int j=0;j<nBlockTypes.size();j++){
        for(int i=0;i<grid.nNumVars+grid.nNumIntVars;i++){
          delete [] nSendBlockDims[j][i];
          delete [] nRecvBlockDims[j][i];
          delete [] nSendBlockStart[j][i];
          delete [] nRecvBlockStart[j][i];
        }
        delete [] nSendBlockStart[j];
        delete [] nRecvBlockStart[j];
        delete [] nSendBlockDims[j];
        delete [] nRecvBlockDims[j];
      }
      delete [] nSendBlockStart;
      delete [] nRecvBlockStart;
      delete [] nSendBlockDims;
      delete [] nRecvBlockDims;
    }
    
    //delete dynamic memory used to determine neighbors
    delete [] nNonPeriodicNeighborRanks;
    delete [] nPeriodicNeighborRanks;
    for(int i=0;i<nNumPeriodicNeighbors;i++){
      delete [] nEdge[i];
    }
    delete [] nEdge;
  }
  
  //allocate memory for send and recieve request handles
  messPass.requestSend=new MPI::Request[procTop.nNumNeighbors];
  messPass.requestRecv=new MPI::Request[procTop.nNumNeighbors];
  
  //allocate memory for send and recieve status
  messPass.statusSend=new MPI::Status[procTop.nNumNeighbors];
  messPass.statusRecv=new MPI::Status[procTop.nNumNeighbors];
  
  //determine starting points for updating old grid, and calculating ghost cell regions
  grid.nStartUpdateExplicit=new int*[grid.nNumVars+grid.nNumIntVars];
  grid.nEndUpdateExplicit=new int*[grid.nNumVars+grid.nNumIntVars];
  grid.nStartUpdateImplicit=new int*[grid.nNumVars+grid.nNumIntVars];
  grid.nEndUpdateImplicit=new int*[grid.nNumVars+grid.nNumIntVars];
  grid.nStartGhostUpdateExplicit=new int**[grid.nNumVars+grid.nNumIntVars];//for each variable
  grid.nEndGhostUpdateExplicit=new int**[grid.nNumVars+grid.nNumIntVars];
  grid.nStartGhostUpdateImplicit=new int**[grid.nNumVars+grid.nNumIntVars];//for each variable
  grid.nEndGhostUpdateImplicit=new int**[grid.nNumVars+grid.nNumIntVars];
  for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){/*works the same for both procTop.nRank==0 and 
    procTop.nRank==* */
    
    //allocate memory
    grid.nStartGhostUpdateExplicit[n]=new int*[2*3];//can have 6 different ghost cell regions
    grid.nEndGhostUpdateExplicit[n]=new int*[2*3];
    grid.nStartGhostUpdateImplicit[n]=new int*[2*3];//can have 6 different ghost cell regions
    grid.nEndGhostUpdateImplicit[n]=new int*[2*3];
    grid.nStartUpdateExplicit[n]=new int[3];
    grid.nEndUpdateExplicit[n]=new int[3];
    grid.nStartUpdateImplicit[n]=new int[3];
    grid.nEndUpdateImplicit[n]=new int[3];
    for(int i=0;i<2*3;i++){
      grid.nStartGhostUpdateExplicit[n][i]=new int[3];
      grid.nEndGhostUpdateExplicit[n][i]=new int[3];
      grid.nStartGhostUpdateImplicit[n][i]=new int[3];
      grid.nEndGhostUpdateImplicit[n][i]=new int[3];
    }
    
    /*determine global start and stop positions of local grid w.r.t. the surface in the radial 
      direction, used to determine extent of implicit region on local grid*/
    int nGlobalStart=0;
    int nGlobalEnd=0;
    for(int p=(procTop.nNumProcs-1);p>=0;p--){
      if( (procTop.nCoords[p][0]>procTop.nCoords[procTop.nRank][0])
        &&((procTop.nCoords[p][1]==procTop.nCoords[procTop.nRank][1])
          ||(procTop.nCoords[p][1]==-1)||(procTop.nCoords[procTop.nRank][1]==-1))
        &&((procTop.nCoords[p][2]==procTop.nCoords[procTop.nRank][2])
          ||(procTop.nCoords[p][2]==-1)||(procTop.nCoords[procTop.nRank][2]==-1))
        ){//add only grids of those processors inside the current processor
        nGlobalStart+=grid.nLocalGridDims[p][n][0];
      }
    }
    nGlobalStart+=grid.nNumGhostCells;
    nGlobalEnd=nGlobalStart+grid.nLocalGridDims[procTop.nRank][n][0];
    if(procTop.nCoords[procTop.nRank][0]==procTop.nProcDims[0]-1){/*adjust for surface boundary 
      conditions being only applied at 1 zone at the surface and not two like the rest of the grid*/
      nGlobalStart=nGlobalStart-grid.nNumGhostCells+1;
    }
    
    //for each direction
    for(int l=0;l<3;l++){
      if(grid.nVariables[n][l]==-1||grid.nLocalGridDims[procTop.nRank][n][l]<grid.nNumGhostCells){/* 
        if variable not defined in direction l, or if the local grid dimensions are less than the 
         number of ghost cells*/
        grid.nStartUpdateExplicit[n][l]=0;
      }
      else{
        grid.nStartUpdateExplicit[n][l]=grid.nNumGhostCells;
      }
      for(int i=0;i<2*3;i++){//most ghost cells regions don't need calculating, initilize to zero
        grid.nStartGhostUpdateExplicit[n][i][l]=0;
        grid.nEndGhostUpdateExplicit[n][i][l]=0;
      }
      
      /*If interface quantity, if not periodic in direction l, and first processor in direction l
        there is an extra inner interface which is bracketing a ghost zone so it shouldn't be 
        updated*/
      if(grid.nVariables[n][l]==1&&procTop.nPeriodic[l]==0&&procTop.nCoords[procTop.nRank][l]==0){
        grid.nStartUpdateExplicit[n][l]++;
      }
      grid.nEndUpdateExplicit[n][l]=grid.nStartUpdateExplicit[n][l]
        +grid.nLocalGridDims[procTop.nRank][n][l];
      
      /*Set outer radial ghost zones for calculating if last in radial direction and there is more 
        than 1 processor in radial direction*/
      if(procTop.nCoords[procTop.nRank][0]==procTop.nProcDims[0]-1&&procTop.nProcDims[0]>1){
        grid.nStartGhostUpdateExplicit[n][0][l]=grid.nStartUpdateExplicit[n][l];//0==outter x1 boundary
        grid.nEndGhostUpdateExplicit[n][0][l]=grid.nEndUpdateExplicit[n][l];
        if(l==0){//only one ghost cell region which is at the outer radial boundary
          grid.nEndUpdateExplicit[n][0]+=grid.nNumGhostCells-1;
          grid.nStartGhostUpdateExplicit[n][0][0]=grid.nEndUpdateExplicit[n][0];
          grid.nEndGhostUpdateExplicit[n][0][0]=grid.nEndUpdateExplicit[n][0]+1;/* one cell in outer radial 
            direction is a ghost cell calculation for processors at outter edge of grid*/
        }
      }
      
      #if SEDOV==1
        /*Set inner radial ghost zones for calculating if first in radial direction*/
        if(procTop.nCoords[procTop.nRank][0]==0){
          grid.nStartGhostUpdateExplicit[n][1][l]=grid.nStartUpdateExplicit[n][l];//0==outter x1 boundary
          grid.nEndGhostUpdateExplicit[n][1][l]=grid.nEndUpdateExplicit[n][l];
          if(l==0){//only one ghost cell region which is at the outer radial boundary
            grid.nStartUpdateExplicit[n][0]-=(grid.nNumGhostCells-1);
            grid.nEndGhostUpdateExplicit[n][1][0]=grid.nStartUpdateExplicit[n][0];
            grid.nStartGhostUpdateExplicit[n][1][0]=grid.nStartUpdateExplicit[n][0]-1;/* one cell in inner radial 
              direction is a ghost cell calculation for processors at inner edge of grid*/
          }
        }
      #endif
      
      /*Set ghost region for 1D processor, used to average the 3D boundary to the 1D boundary if 
        more than one processor in the radial direction, and 1D processor*/
      if(procTop.nProcDims[0]>1&&procTop.nRank==0){
        int nInterface=0;
        if(procTop.nPeriodic[l]!=1){//if not periodic
          if(grid.nVariables[n][l]==1){//if interface variable
            nInterface=1;
          }
        }
        if(l==0){//if in the radial direction
          if(grid.nVariables[n][l]!=-1){/* if not defined in radial direction, nothing updated no 
            matter the size in the other directions*/
            grid.nStartGhostUpdateExplicit[n][0][0]=0;
            grid.nEndGhostUpdateExplicit[n][0][0]=0;
          }
          else{//if defined in the radial direction
            grid.nStartGhostUpdateExplicit[n][0][0]=grid.nEndUpdateExplicit[n][0];
            grid.nEndGhostUpdateExplicit[n][0][0]=grid.nEndUpdateExplicit[n][0]+grid.nNumGhostCells;
          }
        }
        else{//if not the radial direction
          if(grid.nVariables[n][l]==-1){//if not defined in direction l
            grid.nStartGhostUpdateExplicit[n][0][l]=0;
            grid.nEndGhostUpdateExplicit[n][0][l]=grid.nStartGhostUpdateExplicit[n][0][l]+1;
          }
          else{
            grid.nStartGhostUpdateExplicit[n][0][l]=grid.nGlobalGridDims[l]+nInterface;
          }
        }
      }
      
      //adjust for non-time dependent variables
      if(grid.nVariables[n][3]==0){//if not dependent on time
        grid.nEndUpdateExplicit[n][l]=grid.nStartUpdateExplicit[n][l];/*update end is same as start
          i.e. nothing needs to be updated*/
        for(int nGR=0;nGR<2*3;nGR++){//for all regions
          grid.nEndGhostUpdateExplicit[n][nGR][l]=grid.nStartGhostUpdateExplicit[n][nGR][l];
        }
      }
      
      //initilize implicit regions to not update
      for(int i=0;i<2*3;i++){
        grid.nStartGhostUpdateImplicit[n][i][l]=0;
        grid.nEndGhostUpdateImplicit[n][i][l]=0;
      }
      grid.nStartUpdateImplicit[n][l]=0;
      grid.nEndUpdateImplicit[n][l]=0;
      
      //adjust for implicit calculation of variables affected by implicit calculation
      if(n==grid.nE||n==grid.nT||n==grid.nKappa||n==grid.nP){
        
        //set main grid updates
        bool bHasImplicitRegion=false;
        if(nGlobalEnd<=implicit.nNumImplicitZones){//completely implicit
          grid.nStartUpdateImplicit[n][l]=grid.nStartUpdateExplicit[n][l];
          grid.nEndUpdateImplicit[n][l]=grid.nEndUpdateExplicit[n][l];
          grid.nStartUpdateExplicit[n][l]=0;
          grid.nEndUpdateExplicit[n][l]=0;
          bHasImplicitRegion=true;
        }
        else if(nGlobalStart<implicit.nNumImplicitZones&&nGlobalEnd>implicit.nNumImplicitZones){/*
          it is partially implicit*/
          
          //theta and phi directions are the same
          grid.nEndUpdateImplicit[n][l]=grid.nEndUpdateExplicit[n][l];
          grid.nStartUpdateImplicit[n][l]=grid.nStartUpdateExplicit[n][l];
          
          //adjust radial boundaries of the explicit and implicit regions
          if(l==0){
            grid.nEndUpdateExplicit[n][l]=grid.nEndUpdateExplicit[n][l]-(implicit.nNumImplicitZones
              -nGlobalStart);
            grid.nStartUpdateImplicit[n][l]=grid.nEndUpdateExplicit[n][l];
          }
          bHasImplicitRegion=true;
        }
        
        //set ghost updates
        if( (procTop.nCoords[procTop.nRank][0]==procTop.nProcDims[0]-1)
          &&(bHasImplicitRegion||implicit.nNumImplicitZones==1)){/*if a processor at the surface and 
          has an implicit region in main grid, or only the ghost region is implicit.*/
          grid.nStartGhostUpdateImplicit[n][0][l]=grid.nStartGhostUpdateExplicit[n][0][l];
          grid.nEndGhostUpdateImplicit[n][0][l]=grid.nEndGhostUpdateExplicit[n][0][l];
          grid.nStartGhostUpdateExplicit[n][0][l]=0;
          grid.nEndGhostUpdateExplicit[n][0][l]=0;
        }
      }
    }
    
    //start in a few zones from 1D boundary when updating v and w
    if((n==grid.nV||n==grid.nW)&&procTop.nCoords[procTop.nRank][0]!=0){
      
      //get start and end of velocity zering region
      int nStartV0=grid.nNum1DZones+grid.nNumGhostCells;
      int nEndV0=grid.nNum1DZones+grid.nNumGhostCells
        +grid.nNumZones1DBoundaryZeroHorizontalVelocity;
      
      int nStartGlobal=grid.nStartUpdateExplicit[n][0]+grid.nGlobalGridPositionLocalGrid[0]
        -grid.nNumGhostCells;
      int nEndGlobal=grid.nEndUpdateExplicit[n][0]+grid.nGlobalGridPositionLocalGrid[0]
        -grid.nNumGhostCells;
      if(nEndGlobal<=nEndV0&&nStartGlobal<nEndV0){
        grid.nStartUpdateExplicit[n][0]=nEndGlobal-grid.nGlobalGridPositionLocalGrid[0]
          +grid.nNumGhostCells;
      }
      if(nEndGlobal>nEndV0&&nStartGlobal<nEndV0){
        grid.nStartUpdateExplicit[n][0]=nEndV0-grid.nGlobalGridPositionLocalGrid[0]
          +grid.nNumGhostCells;
      }
    }
  }
  
  //set radial neighbors
  procTop.nNumRadialNeighbors=0;
  for(int i=0;i<procTop.nNumNeighbors;i++){
    if( (procTop.nCoords[procTop.nNeighborRanks[i]][0]==procTop.nCoords[procTop.nRank][0]+1
        ||procTop.nCoords[procTop.nNeighborRanks[i]][0]==procTop.nCoords[procTop.nRank][0]-1)
      &&((procTop.nCoords[procTop.nNeighborRanks[i]][1]==procTop.nCoords[procTop.nRank][1])
        ||(procTop.nCoords[procTop.nNeighborRanks[i]][1]==-1)
        ||(procTop.nCoords[procTop.nRank][2]==-1))
      &&((procTop.nCoords[procTop.nNeighborRanks[i]][2]==procTop.nCoords[procTop.nRank][2])
        ||(procTop.nCoords[procTop.nNeighborRanks[i]][2]==-1)
        ||(procTop.nCoords[procTop.nRank][2]==-1))){
      procTop.nNumRadialNeighbors++;
    }
  }
  procTop.nRadialNeighborRanks=new int[procTop.nNumRadialNeighbors];
  procTop.nRadialNeighborNeighborIDs=new int[procTop.nNumRadialNeighbors];
  int nCount=0;
  for(int i=0;i<procTop.nNumNeighbors;i++){
    if( (procTop.nCoords[procTop.nNeighborRanks[i]][0]==procTop.nCoords[procTop.nRank][0]+1
        ||procTop.nCoords[procTop.nNeighborRanks[i]][0]==procTop.nCoords[procTop.nRank][0]-1)
      &&((procTop.nCoords[procTop.nNeighborRanks[i]][1]==procTop.nCoords[procTop.nRank][1])
        ||(procTop.nCoords[procTop.nNeighborRanks[i]][1]==-1)
        ||(procTop.nCoords[procTop.nRank][2]==-1))
      &&((procTop.nCoords[procTop.nNeighborRanks[i]][2]==procTop.nCoords[procTop.nRank][2])
        ||(procTop.nCoords[procTop.nNeighborRanks[i]][2]==-1)
        ||(procTop.nCoords[procTop.nRank][2]==-1))){
      procTop.nRadialNeighborRanks[nCount]=procTop.nNeighborRanks[i];
      procTop.nRadialNeighborNeighborIDs[nCount]=i;
      nCount++;
    }
  }
}
void updateLocalBoundaries(ProcTop &procTop, MessPass &messPass, Grid &grid){
  //reciev from neighbors, into old grid
  for(int i=0;i<procTop.nNumNeighbors;i++){
    messPass.requestRecv[i]=MPI::COMM_WORLD.Irecv(grid.dLocalGridOld,1,messPass.typeRecvOldGrid[i]
      ,procTop.nNeighborRanks[i],0);
  }
  
  //send to neighbors, from new grid
  for(int i=0;i<procTop.nNumNeighbors;i++){
    messPass.requestSend[i]=MPI::COMM_WORLD.Isend(grid.dLocalGridNew,1,messPass.typeSendNewGrid[i]
      ,procTop.nNeighborRanks[i],0);
  }
  
  //update old grid with new grid
  updateOldGrid(procTop,grid);
  
  //wait till all recieves complet on current processor
  MPI::Request::Waitall(procTop.nNumNeighbors,messPass.requestRecv,messPass.statusRecv);
  
  if(procTop.nRank==0){
    //average recieved values
    average3DTo1DBoundariesOld(grid);
  }
  
  //wait till all sends completed on current processor, since the send buffer can't be modified 
  //until after all sends complete.
  MPI::Request::Waitall(procTop.nNumNeighbors,messPass.requestSend,messPass.statusSend);
  
  //wait for all processors to finish
  /**\todo Shouldn't need MPI::COMM_WORLD.Barrier() may want to test out removing this at some
  point as it might produce a bit of a speed up.*/
  MPI::COMM_WORLD.Barrier();
}
void updateLocalBoundariesNewGrid(int nVar, ProcTop &procTop, MessPass &messPass,Grid &grid){
  
  //reciev from neighbors, into new grid
  for(int i=0;i<procTop.nNumNeighbors;i++){
    //irecv might not work, since sending and recieving from same grid
    messPass.requestRecv[i]=MPI::COMM_WORLD.Irecv(grid.dLocalGridNew,1
      ,messPass.typeRecvNewVar[i][nVar],procTop.nNeighborRanks[i],1);
  }
  
  //send to neighbors, from new grid
  for(int i=0;i<procTop.nNumNeighbors;i++){
    //isend might not work since sending and recieving from same grid
    MPI::COMM_WORLD.Send(grid.dLocalGridNew,1,messPass.typeSendNewVar[i][nVar]
      ,procTop.nNeighborRanks[i],1);
  }
  
  //wait till all recieves complet on current processor
  MPI::Request::Waitall(procTop.nNumNeighbors,messPass.requestRecv,messPass.statusRecv);
  
  if(procTop.nRank==0){
    //average recieved values
    average3DTo1DBoundariesNew(grid, nVar);
  }
  
  //wait till all sends completed on current processor, can't modify the send buffer until
  //after all sends complete. This shouldn't be a problem as the new buffer won't be modified until
  //the next time step, 
  /**\todo May want to do some waiting on this message at some point before the end of the timestep,
  but it doesn't need to be done in this function. It might also be that this is built into the code
  by waiting at some other point. This is something that should be checked out at somepoint, perhaps
  once the preformance starts to be analyzed. I would think that if the send buffer was being 
  modified before the send was completed, that there would be some errors poping up that would 
  likely kill the program.*/
  //MPI::Request::Waitall(procTop.nNumNeighbors,messPass.requestRecv,messPass.statusSend);
  
  //wait for all processors to finish, this prevents modification
  //MPI::COMM_WORLD.Barrier();
}
void updateOldGrid(ProcTop &procTop, Grid &grid){
  
  //update the old grid
  for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
    
    //explicit parts of the grid
    for(int i=grid.nStartUpdateExplicit[n][0];i<grid.nEndUpdateExplicit[n][0];i++){
      for(int j=grid.nStartUpdateExplicit[n][1];j<grid.nEndUpdateExplicit[n][1];j++){
        for(int k=grid.nStartUpdateExplicit[n][2];k<grid.nEndUpdateExplicit[n][2];k++){
          grid.dLocalGridOld[n][i][j][k]=grid.dLocalGridNew[n][i][j][k];
        }
      }
    }
    for(int l=0;l<2*3;l++){//for each possible ghost region, usually all zero
      for(int i=grid.nStartGhostUpdateExplicit[n][l][0];i<grid.nEndGhostUpdateExplicit[n][l][0];
        i++){
        for(int j=grid.nStartGhostUpdateExplicit[n][l][1];j<grid.nEndGhostUpdateExplicit[n][l][1];
          j++){
          for(int k=grid.nStartGhostUpdateExplicit[n][l][2];k<grid.nEndGhostUpdateExplicit[n][l][2];
            k++){
            grid.dLocalGridOld[n][i][j][k]=grid.dLocalGridNew[n][i][j][k];
          }
        }
      }
    }
    
    //implicit parts of the grid
    for(int i=grid.nStartUpdateImplicit[n][0];i<grid.nEndUpdateImplicit[n][0];i++){
      for(int j=grid.nStartUpdateImplicit[n][1];j<grid.nEndUpdateImplicit[n][1];j++){
        for(int k=grid.nStartUpdateImplicit[n][2];k<grid.nEndUpdateImplicit[n][2];k++){
          grid.dLocalGridOld[n][i][j][k]=grid.dLocalGridNew[n][i][j][k];
        }
      }
    }
    for(int l=0;l<2*3;l++){//for each possible ghost region, usually all zero
      for(int i=grid.nStartGhostUpdateImplicit[n][l][0];i<grid.nEndGhostUpdateImplicit[n][l][0];
        i++){
        for(int j=grid.nStartGhostUpdateImplicit[n][l][1];j<grid.nEndGhostUpdateImplicit[n][l][1];
          j++){
          for(int k=grid.nStartGhostUpdateImplicit[n][l][2];k<grid.nEndGhostUpdateImplicit[n][l][2];
            k++){
            grid.dLocalGridOld[n][i][j][k]=grid.dLocalGridNew[n][i][j][k];
          }
        }
      }
    }
  }
}
void updateNewGridWithOld(Grid &grid, ProcTop &procTop){
  //update the old grid
  for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
    int nStartX=0;
    int nStartY=0;
    int nStartZ=0;
    int nEndX=grid.nLocalGridDims[procTop.nRank][n][0]+2*grid.nNumGhostCells;
    int nEndY=1;
    if(grid.nNumDims>1){
      nEndY=grid.nLocalGridDims[procTop.nRank][n][1]+2*grid.nNumGhostCells;
    }
    int nEndZ=1;
    if(grid.nNumDims>2){
      nEndZ=grid.nLocalGridDims[procTop.nRank][n][2]+2*grid.nNumGhostCells;
    }
    if(procTop.nRank==0){//only has 1D
      nEndY=1;
      nEndZ=1;
    }
    if(grid.nVariables[n][0]==-1){
      nEndX=1;
    }
    if(grid.nVariables[n][0]==-1&&procTop.nRank==0){
      nEndX=nStartX;
    }
    if(grid.nVariables[n][1]==-1){
      nEndY=1;
    }
    if(grid.nVariables[n][2]==-1){
      nEndZ=1;
    }
    for(int i=nStartX;i<nEndX;i++){
      for(int j=nStartY;j<nEndY;j++){
        for(int k=nStartZ;k<nEndZ;k++){
          grid.dLocalGridNew[n][i][j][k]=grid.dLocalGridOld[n][i][j][k];
        }
      }
    }
  }
}
void average3DTo1DBoundariesOld(Grid &grid){
  //only to be called by procTop.nRank==0
  //not needed if procTop.nNumProcs=1
  //old grid has been completely updated with new grid at this point
  
  for(int n=0;n<grid.nNumVars+grid.nNumIntVars;n++){
    for(int i=grid.nStartGhostUpdateExplicit[n][0][0];i<grid.nEndGhostUpdateExplicit[n][0][0];i++){
        
      //calculate i for interface centered quantities and zone centered quantities
      //depends on weather the variable is zone or interface centered
      int nIInt=i;
      int nICen=i;
      if(grid.nVariables[n][0]==1){//interface variable
        nICen-=grid.nCenIntOffset[0];
      }
      if(grid.nVariables[n][0]==0){//centered variable
        nIInt+=grid.nCenIntOffset[0];
      }
      
      double dSum=0.0;
      double dVolume=0.0;//total volume of shell
      double dRFactor=0.33333333333333333*(pow(grid.dLocalGridOld[grid.nR][nIInt][0][0],3.0)
        -pow(grid.dLocalGridOld[grid.nR][nIInt-1][0][0],3.0));
      
      for(int j=grid.nStartGhostUpdateExplicit[n][0][1];j<grid.nEndGhostUpdateExplicit[n][0][1];j++){
        
        //need nJCen, may or may not be j, depending on if n is an interface or zone centered quantity
        int nJCen=j;
        int nJInt=j;
        if(grid.nVariables[n][1]==1){//interface variable
          nJCen-=grid.nCenIntOffset[1];
        }
        if(grid.nVariables[n][1]==0){//centered variable
          nJInt+=grid.nCenIntOffset[1];
        }
        
        for(int k=grid.nStartGhostUpdateExplicit[n][0][2];k<grid.nEndGhostUpdateExplicit[n][0][2];k++){
          
          //need nKCen, may or may not be k, depending on if n is an interface or zone centered quantity
          int nKCen=k;
          int nKInt=k;
          if(grid.nVariables[n][2]==1){//interface variable
            nKCen-=grid.nCenIntOffset[2];
          }
          if(grid.nVariables[n][2]==0){//centered variable
            nKInt+=grid.nCenIntOffset[2];
          }
          
          double dVolumeTemp=dRFactor*grid.dLocalGridOld[grid.nDCosThetaIJK][0][nJCen][0]
            *grid.dLocalGridOld[grid.nDPhi][0][0][nKCen];;
            
          dSum+=dVolumeTemp*grid.dLocalGridOld[n][i][j][k];
          dVolume+=dVolumeTemp;
        }
      }
      grid.dLocalGridOld[n][i][0][0]=dSum/dVolume;
    }
  }
}
void average3DTo1DBoundariesNew(Grid &grid, int nVar){
  //only to be called by procTop.nRank==0
  //not needed if procTop.nNumProcs=1
  for(int i=grid.nStartGhostUpdateExplicit[nVar][0][0];i<grid.nEndGhostUpdateExplicit[nVar][0][0];
    i++){
      
    //calculate i for interface centered quantities and zone centered quantities
    //depends on weather the variable is zone or interface centered
    int nIInt=i;
    int nICen=i;
    if(grid.nVariables[nVar][0]==1){//interface variable
      nICen-=grid.nCenIntOffset[0];
    }
    if(grid.nVariables[nVar][0]==0){//centered variable
      nIInt+=grid.nCenIntOffset[0];
    }
    
    double dSum=0.0;
    double dVolume=0.0;//total volume of shell
    double dRFactor=0.33333333333333333*(pow(grid.dLocalGridOld[grid.nR][nIInt][0][0],3.0)
      -pow(grid.dLocalGridOld[grid.nR][nIInt-1][0][0],3.0));
    
    for(int j=grid.nStartGhostUpdateExplicit[nVar][0][1];j<grid.nEndGhostUpdateExplicit[nVar][0][1];
      j++){
      
      //need nJCen, may or may not be j, depending on if n is an interface or zone centered quantity
      int nJCen=j;
      int nJInt=j;
      if(grid.nVariables[nVar][1]==1){//interface variable
        nJCen-=grid.nCenIntOffset[1];
      }
      if(grid.nVariables[nVar][1]==0){//centered variable
        nJInt+=grid.nCenIntOffset[1];
      }
      
      for(int k=grid.nStartGhostUpdateExplicit[nVar][0][2];
        k<grid.nEndGhostUpdateExplicit[nVar][0][2];k++){
        
        //need nKCen, may or may not be k, depending on if n is an interface or zone centered quantity
        int nKCen=k;
        int nKInt=k;
        if(grid.nVariables[nVar][2]==1){//interface variable
          nKCen-=grid.nCenIntOffset[2];
        }
        if(grid.nVariables[nVar][2]==0){//centered variable
          nKInt+=grid.nCenIntOffset[2];
        }
        
        double dVolumeTemp=dRFactor*grid.dLocalGridOld[grid.nDCosThetaIJK][0][nJCen][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][nKCen];;
          
        dSum+=dVolumeTemp*grid.dLocalGridNew[nVar][i][j][k];
        dVolume+=dVolumeTemp;
      }
    }
    grid.dLocalGridNew[nVar][i][0][0]=dSum/dVolume;
  }
}
void updateLocalBoundaryVelocitiesNewGrid_R(ProcTop &procTop,MessPass &messPass,Grid &grid){
  updateLocalBoundariesNewGrid(grid.nU,procTop,messPass,grid);
}
void updateLocalBoundaryVelocitiesNewGrid_RT(ProcTop &procTop,MessPass &messPass,Grid &grid){
  updateLocalBoundariesNewGrid(grid.nU,procTop,messPass,grid);
  updateLocalBoundariesNewGrid(grid.nV,procTop,messPass,grid);
}
void updateLocalBoundaryVelocitiesNewGrid_RTP(ProcTop &procTop,MessPass &messPass,Grid &grid){
  updateLocalBoundariesNewGrid(grid.nU,procTop,messPass,grid);
  updateLocalBoundariesNewGrid(grid.nV,procTop,messPass,grid);
  updateLocalBoundariesNewGrid(grid.nW,procTop,messPass,grid);
}
void initImplicitCalculation(Implicit &implicit, Grid &grid, ProcTop &procTop, int nNumArgs
  , char* cArgs[]){
  
  //initilize PETSc
  PetscInitialize(&nNumArgs,&cArgs,PETSC_NULL,PETSC_NULL);
  
  
  //INITIALIZE COEFFECIENT MATRIX DATA STRUCTURE
  
  //figure out distribution of rows on each processor
  int nNumGlobalRows=implicit.nNumImplicitZones*grid.nGlobalGridDims[1]*grid.nGlobalGridDims[2];
  int *nNumLocalRows=new int[procTop.nNumProcs];
  int nTemp=int(nNumGlobalRows/procTop.nNumProcs);
  int nRemainder=nNumGlobalRows%procTop.nNumProcs;
  for(int p=0;p<procTop.nNumProcs;p++){
    nNumLocalRows[p]=nTemp;
    if(p<nRemainder){
      nNumLocalRows[p]++;
    }
  }
  
  //count number of rows on other processors before first row on current processor
  int nNumRowsBefore=0;
  for(int p=0;p<procTop.nNumProcs;p++){
    if(p<procTop.nRank){
      nNumRowsBefore+=nNumLocalRows[p];
    }
  }
  
  //count number of non-zero elements per row in the diagonal and off-diagonal submatrices
  int *nNumNonzeroElementsPerRowOD=new int[nNumLocalRows[procTop.nRank]];
  int *nNumNonzeroElementsPerRowD=new int[nNumLocalRows[procTop.nRank]];
  int nNumHorizontalZones=grid.nGlobalGridDims[1]*grid.nGlobalGridDims[2];
  int nIndex=0;
  for(int i=0;i<implicit.nNumImplicitZones;i++){//i,j,k are row iterators
    for(int j=0;j<grid.nGlobalGridDims[1];j++){
      for(int k=0;k<grid.nGlobalGridDims[2];k++){
        int p=k+j*grid.nGlobalGridDims[2]+i*nNumHorizontalZones;
        if(p>=nNumRowsBefore&&p<(nNumRowsBefore+nNumLocalRows[procTop.nRank])){/* if row on local
          processor*/
          nNumNonzeroElementsPerRowD[nIndex]=0;/* initialize number of non-zero diagonal submatrix
            columns*/
          nNumNonzeroElementsPerRowOD[nIndex]=0;
          for(int l=0;l<implicit.nNumImplicitZones;l++){//l,m,n are column iterators
            for(int m=0;m<grid.nGlobalGridDims[1];m++){
              for(int n=0;n<grid.nGlobalGridDims[2];n++){
                int q=n+m*grid.nGlobalGridDims[2]+l*nNumHorizontalZones;
                
                //determine if the element is non-zero
                bool bNonZero=false;
                if(l==i&&m==j&&n==k){//i,j,k
                  bNonZero=true;
                }
                else if(l==i+1&&m==j&&n==k){//i+1,j,k
                  bNonZero=true;
                }
                else if(l==i-1&&m==j&&n==k){//i-1,j,k
                  bNonZero=true;
                }
                else if(l==i&&m==j+1&&n==k&&grid.nNumDims>1){//i,j+1,k
                  bNonZero=true;
                }
                else if(l==i&&m==j-1&&n==k&&grid.nNumDims>1){//i,j-1,k
                  bNonZero=true;
                }
                else if(l==i&&m==j&&n==k+1&&grid.nNumDims>2){//i,j,k+1
                  bNonZero=true;
                }
                else if(l==i&&m==j&&n==k-1&&grid.nNumDims>2){//i,k,k-1
                  bNonZero=true;
                }
                else if(l==i&&m==grid.nGlobalGridDims[1]-1&&j==0&&n==k&&grid.nNumDims>1){//theta inner boundary when j==0
                  bNonZero=true;
                }
                else if(l==i&&m==0&&j==grid.nGlobalGridDims[1]-1&&n==k&&grid.nNumDims>1){//theta outter boundary when j==nNumTheta-1
                  bNonZero=true;
                }
                else if(l==i&&m==j&&n==grid.nGlobalGridDims[2]-1&&k==0&&grid.nNumDims>2){//phi inner boundary when k==0
                  bNonZero=true;
                }
                else if(l==i&&m==j&&n==0&&k==grid.nGlobalGridDims[2]-1&&grid.nNumDims>2){//phi outter boundary when k==nNumPhi-1
                  bNonZero=true;
                }
                else{//element is zero
                  bNonZero=false;
                }
                
                if(bNonZero){/*if non-zero, count it in either Diagonal, or off diaginal non-zero 
                  element count*/
                  if(q>=nNumRowsBefore&&q<nNumRowsBefore+nNumLocalRows[procTop.nRank]){/*column is
                    in diagonal submatrix*/
                    nNumNonzeroElementsPerRowD[nIndex]++;
                  }
                  else{//column is not in diagonal submatrix (off-diagonal submatrix)
                    nNumNonzeroElementsPerRowOD[nIndex]++;
                  }
                }
                
              }
            }
          }
          nIndex++;//next row
        }
      }
    }
  }
  
  //initilize coeffecient matrix
  MatCreateMPIAIJ(PETSC_COMM_WORLD
    ,nNumLocalRows[procTop.nRank]//local number of rows in the rhs vector
    ,nNumLocalRows[procTop.nRank]//local number of rows in the solution vector
    ,nNumGlobalRows//global number of rows of the coeffecient matrix
    ,nNumGlobalRows//global number of columns of the coeffecient matrix
    ,0//set size of diaginal submatrix to zero
    ,nNumNonzeroElementsPerRowD//set array of diaginal submatrix rows sizes to null
    ,0//set size of off-diagonal submatrix to zero
    ,nNumNonzeroElementsPerRowOD//set array of off-diaginal submatrix rows sizes to null
    ,&implicit.matCoeff);
  
  //initialize rhs vector
  VecCreateMPI(PETSC_COMM_WORLD,nNumLocalRows[procTop.nRank],nNumGlobalRows,&implicit.vecRHS);
  
  //initialize solution vector
  VecCreateMPI(PETSC_COMM_WORLD,nNumLocalRows[procTop.nRank],nNumGlobalRows
    ,&implicit.vecTCorrections);
  
  //SET DERVIATIVES INFOS
  
  //find start and end position of local grid in the global grid
  int nLocalGridStart[3]={0,0,0};//holds start position of processor procTop.nRank in global grid
  int nLocalGridEnd[3]={0,0,0};//holds end position of processor procTop.nRank in global grid
  if(grid.nNumDims>0){
    nLocalGridStart[0]=grid.nNumGhostCells;
  }
  for(int p=0;p<procTop.nRank;p++){
    if(procTop.nCoords[p][2]<procTop.nCoords[procTop.nRank][2]
      &&procTop.nCoords[p][1]==procTop.nCoords[procTop.nRank][1]
      &&procTop.nCoords[p][0]==procTop.nCoords[procTop.nRank][0]){
      
      //add any offset due to position in dimension 2
      nLocalGridStart[2]+=grid.nLocalGridDims[p][grid.nT][2];
    }
  }
  nLocalGridEnd[2]=nLocalGridStart[2]+grid.nLocalGridDims[procTop.nRank][grid.nT][2];
  for(int p=0;p<procTop.nRank;p++){
    if(procTop.nCoords[p][2]==procTop.nCoords[procTop.nRank][2]
      &&procTop.nCoords[p][1]<procTop.nCoords[procTop.nRank][1]
      &&procTop.nCoords[p][0]==procTop.nCoords[procTop.nRank][0]){
      
      //Add any offset due to position in dimension 1
      nLocalGridStart[1]+=grid.nLocalGridDims[p][grid.nT][1];
    }
  }
  nLocalGridEnd[1]=nLocalGridStart[1]+grid.nLocalGridDims[procTop.nRank][grid.nT][1];
  for(int p=0;p<procTop.nRank;p++){
    if( (procTop.nCoords[p][2]==procTop.nCoords[procTop.nRank][2]||procTop.nCoords[p][2]==-1)
      &&(procTop.nCoords[p][1]==procTop.nCoords[procTop.nRank][1]||procTop.nCoords[p][2]==-1)
      &&procTop.nCoords[p][0]<procTop.nCoords[procTop.nRank][0]){
      
      //Add any offset due to position in dimension 0
      nLocalGridStart[0]+=grid.nLocalGridDims[p][grid.nT][0];
    }
  }
  nLocalGridEnd[0]=nLocalGridStart[0]+grid.nLocalGridDims[procTop.nRank][grid.nT][0];
  int nNumRadialZonesInSB=0;
  if(procTop.nCoords[procTop.nRank][0]==(procTop.nProcDims[0]-1)){/*if at the surface, add ghost 
    cells*/
    nLocalGridEnd[0]+=grid.nNumGhostCells-1;
    nNumRadialZonesInSB=1;
  }
  
  //count number of rows on local processor
  implicit.nNumRowsALocal=0;
  implicit.nNumRowsALocalSB=0;
  int nStartGlobal0=grid.nGlobalGridDims[0]+2*grid.nNumGhostCells-implicit.nNumImplicitZones;
  for(int i=0;i<implicit.nNumImplicitZones;i++){//i,j,k are row iterators
    for(int j=0;j<grid.nGlobalGridDims[1];j++){
      for(int k=0;k<grid.nGlobalGridDims[2];k++){
        if( (nStartGlobal0+i>=nLocalGridStart[0]&&nStartGlobal0+i<nLocalGridEnd[0])
          &&(j>=nLocalGridStart[1]&&j<nLocalGridEnd[1])
          &&(k>=nLocalGridStart[2]&&k<nLocalGridEnd[2])){//if on local grid
          implicit.nNumRowsALocal++;
        }
        if( 
          (nStartGlobal0+i>=nLocalGridEnd[0]&&nStartGlobal0+i<nLocalGridEnd[0]+nNumRadialZonesInSB)
          &&(j>=nLocalGridStart[1]&&j<nLocalGridEnd[1])
          &&(k>=nLocalGridStart[2]&&k<nLocalGridEnd[2])){//if on local grid
          implicit.nNumRowsALocalSB++;
        }
      }
    }
  }
  
  //count number of deriviatives per row on local processor
  implicit.nNumDerPerRow=new int[implicit.nNumRowsALocal+implicit.nNumRowsALocalSB];
  nIndex=0;
  for(int i=0;i<implicit.nNumImplicitZones;i++){//i,j,k are row iterators
    for(int j=0;j<grid.nGlobalGridDims[1];j++){
      for(int k=0;k<grid.nGlobalGridDims[2];k++){
        if((nStartGlobal0+i>=nLocalGridStart[0]&&nStartGlobal0+i<nLocalGridEnd[0]
          +nNumRadialZonesInSB)
          &&(j>=nLocalGridStart[1]&&j<nLocalGridEnd[1])
          &&(k>=nLocalGridStart[2]&&k<nLocalGridEnd[2])){//if on local grid
          implicit.nNumDerPerRow[nIndex]=0;//initialize to zero
          for(int l=0;l<implicit.nNumImplicitZones;l++){//l,m,n are column iterators
            for(int m=0;m<grid.nGlobalGridDims[1];m++){
              for(int n=0;n<grid.nGlobalGridDims[2];n++){
                
                //determine if the element is non-zero
                if(l==i&&m==j&&n==k){//i,j,k
                  implicit.nNumDerPerRow[nIndex]++;
                }
                else if(l==i+1&&m==j&&n==k){//i+1,j,k
                  implicit.nNumDerPerRow[nIndex]++;
                }
                else if(l==i-1&&m==j&&n==k){//i-1,j,k
                  implicit.nNumDerPerRow[nIndex]++;
                }
                else if(l==i&&m==j+1&&n==k&&grid.nNumDims>1){//i,j+1,k
                  implicit.nNumDerPerRow[nIndex]++;
                }
                else if(l==i&&m==j-1&&n==k&&grid.nNumDims>1){//i,j-1,k
                  implicit.nNumDerPerRow[nIndex]++;
                }
                else if(l==i&&m==j&&n==k+1&&grid.nNumDims>2){//i,j,k+1
                  implicit.nNumDerPerRow[nIndex]++;
                }
                else if(l==i&&m==j&&n==k-1&&grid.nNumDims>2){//i,k,k-1
                  implicit.nNumDerPerRow[nIndex]++;
                }
                else if(l==i&&m==grid.nGlobalGridDims[1]-1&&j==0&&n==k&&grid.nNumDims>1){//theta inner boundary when j==0
                  implicit.nNumDerPerRow[nIndex]++;
                }
                else if(l==i&&m==0&&j==grid.nGlobalGridDims[1]-1&&n==k&&grid.nNumDims>1){//theta outter boundary when j==nNumTheta-1
                  implicit.nNumDerPerRow[nIndex]++;
                }
                else if(l==i&&m==j&&n==grid.nGlobalGridDims[2]-1&&k==0&&grid.nNumDims>2){//phi inner boundary when k==0
                  implicit.nNumDerPerRow[nIndex]++;
                }
                else if(l==i&&m==j&&n==0&&k==grid.nGlobalGridDims[2]-1&&grid.nNumDims>2){//phi outter boundary when k==nNumPhi-1
                  implicit.nNumDerPerRow[nIndex]++;
                }
                else{//element is zero
                }
              }
            }
          }
          nIndex++;//next row
        }
      }
    }
  }
    
  //set global coeffecient matrix row and column indices, set local grid indicies and derivative type
  implicit.nTypeDer=new int*[implicit.nNumRowsALocal+implicit.nNumRowsALocalSB];
  implicit.nLocFun=new int*[implicit.nNumRowsALocal+implicit.nNumRowsALocalSB];
  implicit.nLocDer=new int**[implicit.nNumRowsALocal+implicit.nNumRowsALocalSB];
  nIndex=0;//start back at frist row
  int *nFromIndex=new int[implicit.nNumRowsALocal+implicit.nNumRowsALocalSB];/*global index in the
    distributed vecTCorrections vector from which the local rows are needed.*/
  int *nToIndex=new int[implicit.nNumRowsALocal+implicit.nNumRowsALocalSB];/*local index in the
    local vecTCorrectionsLocal vector into which to put the temperature corrections.*/
  for(int i=0;i<implicit.nNumImplicitZones;i++){//i,j,k are row iterators
    for(int j=0;j<grid.nGlobalGridDims[1];j++){
      for(int k=0;k<grid.nGlobalGridDims[2];k++){
        int p=k+j*grid.nGlobalGridDims[2]+i*nNumHorizontalZones;
          if((nStartGlobal0+i>=nLocalGridStart[0]&&nStartGlobal0+i<nLocalGridEnd[0]
            +nNumRadialZonesInSB)
            && (j>=nLocalGridStart[1]&&j<nLocalGridEnd[1])
            && (k>=nLocalGridStart[2]&&k<nLocalGridEnd[2])
            ){//if on local grid
            nFromIndex[nIndex]=p;
            nToIndex[nIndex]=nIndex;
            implicit.nTypeDer[nIndex]=new int[implicit.nNumDerPerRow[nIndex]];
            implicit.nLocDer[nIndex]=new int*[2];
            implicit.nLocDer[nIndex][0]=new int[implicit.nNumDerPerRow[nIndex]];
            implicit.nLocDer[nIndex][1]=new int[implicit.nNumDerPerRow[nIndex]];
            implicit.nLocFun[nIndex]=new int[3];
            implicit.nLocFun[nIndex][0]=nStartGlobal0+i-nLocalGridStart[0]+grid.nNumGhostCells;
            implicit.nLocFun[nIndex][1]=j-nLocalGridStart[1];
            if(grid.nNumDims>1){//add in ghost cells
              implicit.nLocFun[nIndex][1]+=grid.nNumGhostCells;
            }
            implicit.nLocFun[nIndex][2]=k-nLocalGridStart[2];
            if(grid.nNumDims>2){//add in ghost cells
              implicit.nLocFun[nIndex][2]+=grid.nNumGhostCells;
            }
            int nIndex1=0;
            for(int l=0;l<implicit.nNumImplicitZones;l++){//l,m,n are column iterators
              for(int m=0;m<grid.nGlobalGridDims[1];m++){
                for(int n=0;n<grid.nGlobalGridDims[2];n++){
                  int q=n+m*grid.nGlobalGridDims[2]+l*nNumHorizontalZones;
                  
                  //determine if the element is non-zero
                  if(l==i&&m==j&&n==k){//i,j,k
                    implicit.nLocDer[nIndex][0][nIndex1]=p;
                    implicit.nLocDer[nIndex][1][nIndex1]=q;
                    implicit.nTypeDer[nIndex][nIndex1]=0;
                    nIndex1++;
                  }
                  else if(l==i+1&&m==j&&n==k){//i+1,j,k
                    implicit.nLocDer[nIndex][0][nIndex1]=p;
                    implicit.nLocDer[nIndex][1][nIndex1]=q;
                    implicit.nTypeDer[nIndex][nIndex1]=1;
                    nIndex1++;
                  }
                  else if(l==i-1&&m==j&&n==k){//i-1,j,k
                    implicit.nLocDer[nIndex][0][nIndex1]=p;
                    implicit.nLocDer[nIndex][1][nIndex1]=q;
                    implicit.nTypeDer[nIndex][nIndex1]=2;
                    nIndex1++;
                  }
                  else if(l==i&&m==j+1&&n==k&&grid.nNumDims>1){//i,j+1,k
                    implicit.nLocDer[nIndex][0][nIndex1]=p;
                    implicit.nLocDer[nIndex][1][nIndex1]=q;
                    implicit.nTypeDer[nIndex][nIndex1]=3;
                    if(m==grid.nGlobalGridDims[1]-1&&j==0){//theta inner boundary when j==0
                      implicit.nLocDer[nIndex][0][nIndex1]=p;
                      implicit.nLocDer[nIndex][1][nIndex1]=q;
                      implicit.nTypeDer[nIndex][nIndex1]=34;//really just j+1 and j-1
                    }
                    nIndex1++;
                  }
                  else if(l==i&&m==j-1&&n==k&&grid.nNumDims>1){//i,j-1,k
                    implicit.nLocDer[nIndex][0][nIndex1]=p;
                    implicit.nLocDer[nIndex][1][nIndex1]=q;
                    implicit.nTypeDer[nIndex][nIndex1]=4;
                    if(j==grid.nGlobalGridDims[1]-1&&m==0){//theta outter boundary when j==nNumTheta-1. m==0
                      implicit.nLocDer[nIndex][0][nIndex1]=p;
                      implicit.nLocDer[nIndex][1][nIndex1]=q;
                      implicit.nTypeDer[nIndex][nIndex1]=34;//really just j-1 and j+1
                    }
                    nIndex1++;
                  }
                  else if(l==i&&m==j&&n==k+1&&grid.nNumDims>2){//i,j,k+1
                    implicit.nLocDer[nIndex][0][nIndex1]=p;
                    implicit.nLocDer[nIndex][1][nIndex1]=q;
                    implicit.nTypeDer[nIndex][nIndex1]=5;
                    if(n==grid.nGlobalGridDims[2]-1&&k==0){//phi inner boundary when k==0
                      implicit.nLocDer[nIndex][0][nIndex1]=p;
                      implicit.nLocDer[nIndex][1][nIndex1]=q;
                      implicit.nTypeDer[nIndex][nIndex1]=56;//really just k-1 and k+1
                    }
                    nIndex1++;
                  }
                  else if(l==i&&m==j&&n==k-1&&grid.nNumDims>2){//i,k,k-1
                    implicit.nLocDer[nIndex][0][nIndex1]=p;
                    implicit.nLocDer[nIndex][1][nIndex1]=q;
                    implicit.nTypeDer[nIndex][nIndex1]=6;
                    if(n==0&&k==grid.nGlobalGridDims[2]-1){//phi outter boundary when k==nNumPhi-1
                      implicit.nLocDer[nIndex][0][nIndex1]=p;
                      implicit.nLocDer[nIndex][1][nIndex1]=q;
                      implicit.nTypeDer[nIndex][nIndex1]=56;//really just k+1 and k-1
                    }
                    nIndex1++;
                  }
                  else if(l==i&&m==grid.nGlobalGridDims[1]-1&&j==0&&n==k&&grid.nNumDims>1){//theta inner boundary when j==0
                    implicit.nLocDer[nIndex][0][nIndex1]=p;
                    implicit.nLocDer[nIndex][1][nIndex1]=q;
                    implicit.nTypeDer[nIndex][nIndex1]=4;//really just j-1
                    nIndex1++;
                  }
                  else if(l==i&&m==0&&j==grid.nGlobalGridDims[1]-1&&n==k&&grid.nNumDims>1){//theta outter boundary when j==nNumTheta-1
                    implicit.nLocDer[nIndex][0][nIndex1]=p;
                    implicit.nLocDer[nIndex][1][nIndex1]=q;
                    implicit.nTypeDer[nIndex][nIndex1]=3;//really just j+1
                    nIndex1++;
                  }
                  else if(l==i&&m==j&&n==grid.nGlobalGridDims[2]-1&&k==0&&grid.nNumDims>2){//phi inner boundary when k==0
                    implicit.nLocDer[nIndex][0][nIndex1]=p;
                    implicit.nLocDer[nIndex][1][nIndex1]=q;
                    implicit.nTypeDer[nIndex][nIndex1]=6;//really just k-1
                    nIndex1++;
                  }
                  else if(l==i&&m==j&&n==0&&k==grid.nGlobalGridDims[2]-1&&grid.nNumDims>2){//phi outter boundary when k==nNumPhi-1
                    implicit.nLocDer[nIndex][0][nIndex1]=p;
                    implicit.nLocDer[nIndex][1][nIndex1]=q;
                    implicit.nTypeDer[nIndex][nIndex1]=5;//really just k+1
                    nIndex1++;
                  }
                  else{//element is zero
                  }
                }
              }
            }
            nIndex++;
          }
        }
     }
  }
  
  //create solver context
  KSPCreate(PETSC_COMM_WORLD,&implicit.kspContext);
  int ierr;
  //PC pcPreconditioner;
  //ierr = KSPGetPC(implicit.kspContext,&pcPreconditioner);
  //ierr = PCSetType(pcPreconditioner,PCLU);
  //ierr=KSPSetType(implicit.kspContext,KSPCG);
  KSPSetFromOptions(implicit.kspContext);//set from command line options
  ierr = KSPSetTolerances(implicit.kspContext,implicit.dTolerance,PETSC_DEFAULT,PETSC_DEFAULT
    ,implicit.nMaxNumIterations);
  
  //initialize the local vector to hold the corrections
  VecCreateSeq(PETSC_COMM_SELF,implicit.nNumRowsALocal+implicit.nNumRowsALocalSB
    ,&implicit.vecTCorrectionsLocal);//maybe of zero size
  
  //use VecScatterCreate to create the pathway to send information back to the local processors
  IS isFrom;
  IS isTo;/*not sure if these need to be declared in a larger scope or not, it might cause problems
    if they go out of scope. If it doesn't cause problems, maybe I can just destroy them here in 
    this scope when I am done with them.*/

  ISCreateGeneral(PETSC_COMM_SELF,implicit.nNumRowsALocal+implicit.nNumRowsALocalSB,nFromIndex
    ,&isFrom);
  ISCreateGeneral(PETSC_COMM_SELF,implicit.nNumRowsALocal+implicit.nNumRowsALocalSB,nToIndex,&isTo);
  VecScatterCreate(implicit.vecTCorrections,isFrom,implicit.vecTCorrectionsLocal,isTo
    ,&implicit.vecscatTCorrections);
  
  /**\todo isFrom, isTo, matCoeff,vecTCorrections, vecTCorrections,vecRHS,vecTCorrectionsLocal
  ,kspContext,vecscatTCorrections all need to be destroyed before program finishes.*/
}
