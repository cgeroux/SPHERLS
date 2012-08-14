void setMainFunctions(Functions &functions,ProcTop &procTop,Parameters &parameters,Grid &grid
  ,Time &time,Implicit &implicit){
  
  //set some defaults
  functions.fpCalculateNewEddyVisc=&calNewEddyVisc_None;
  functions.fpImplicitSolve=&implicitSolve_None;
  functions.fpImplicitEnergyFunction=&dImplicitEnergyFunction_None;
  functions.fpImplicitEnergyFunction_SB=&dImplicitEnergyFunction_None;
  
  //rank 0 will be 1D, so always want to use 1D version of these equations
  if(procTop.nRank==0){// proc 1 always uses 1D
    
    functions.fpCalculateNewGridVelocities=&calNewU0_R;
    functions.fpCalculateNewRadii=&calNewR;
    functions.fpCalculateNewDensities=&calNewD_R;
    functions.fpCalculateAveDensities=&calNewDenave_R;
    functions.fpUpdateLocalBoundaryVelocitiesNewGrid=&updateLocalBoundaryVelocitiesNewGrid_R;
    if(parameters.bEOSGammaLaw){//use gamma law gas
      functions.fpCalculateDeltat=&calDelt_R_GL;
      functions.fpCalculateNewEOSVars=&calNewP_GL;
      functions.fpCalculateNewAV=&calNewQ0_R_GL;
      functions.fpModelWrite=&modelWrite_GL;
      functions.fpWriteWatchZones=&writeWatchZones_R_GL;
    }
    else{//use tabulated equation of state
      functions.fpCalculateDeltat=&calDelt_R_TEOS;
      functions.fpCalculateNewEOSVars=&calNewTPKappaGamma_TEOS;
      functions.fpCalculateNewAV=&calNewQ0_R_TEOS;
      functions.fpModelWrite=&modelWrite_TEOS;
      functions.fpWriteWatchZones=&writeWatchZones_R_TEOS;
    }
    
    //velocity equation
    functions.fpCalculateNewVelocities=&calNewVelocities_R;
    
    //energy equation, eddy viscosity
    if(parameters.bAdiabatic){//adiabatic
      functions.fpCalculateNewEnergies=&calNewE_R_AD;
    }
    else{//non-adaibatic
      if(!parameters.bEOSGammaLaw){//needs a tabulated equation of state
        functions.fpCalculateNewEnergies=&calNewE_R_NA;
        functions.fpCalculateNewEddyVisc=&calNewEddyVisc_None;
        if(implicit.nNumImplicitZones>0){//implicit, requires non-adiabatic
          functions.fpImplicitSolve=&implicitSolve_R;
          functions.fpImplicitEnergyFunction=&dImplicitEnergyFunction_R;
          functions.fpImplicitEnergyFunction_SB=&dImplicitEnergyFunction_R_SB;
        }
      }
      else{//can't do a non-adiabatic calculation, with a gamma-law gas
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
          <<": User selected to do a non-adiabatic calculation but starting model uses a gamma-law"
          <<" gas. Starting model must use a tabulated equation of state in order to perform a"
          <<" non-adiabatic calculation.\n";
        throw exception2(ssTemp.str(),CALCULATION);
      }
    }
    
    /*Processor 0 must update all the velocities that the other processors do, so that they don't
    get stuck on mpi::waitall calls in the update routines*/
    
    //set density average function, and grid update function
    if(grid.nNumDims==1){
      functions.fpCalculateAveDensities=&calNewDenave_None;
      functions.fpUpdateLocalBoundaryVelocitiesNewGrid=&updateLocalBoundaryVelocitiesNewGrid_R;
    }
    else if(grid.nNumDims==2){
      functions.fpUpdateLocalBoundaryVelocitiesNewGrid=&updateLocalBoundaryVelocitiesNewGrid_RT;
    }
    else if(grid.nNumDims==3){
      functions.fpUpdateLocalBoundaryVelocitiesNewGrid=&updateLocalBoundaryVelocitiesNewGrid_RTP;
    }
  }
  else{
    if(grid.nNumDims==3){//use 3D
      functions.fpCalculateNewGridVelocities=&calNewU0_RTP;
      functions.fpCalculateNewRadii=&calNewR;
      functions.fpCalculateNewDensities=&calNewD_RTP;
      functions.fpCalculateAveDensities=&calNewDenave_RTP;
      functions.fpUpdateLocalBoundaryVelocitiesNewGrid=&updateLocalBoundaryVelocitiesNewGrid_RTP;
      if(parameters.bEOSGammaLaw){//use gamma law gas
        functions.fpCalculateDeltat=&calDelt_RTP_GL;
        functions.fpCalculateNewEOSVars=&calNewP_GL;
        functions.fpCalculateNewAV=&calNewQ0Q1Q2_RTP_GL;
        functions.fpModelWrite=&modelWrite_GL;
        functions.fpWriteWatchZones=&writeWatchZones_RTP_GL;
      }
      else{//use tabulated equation of state
        functions.fpCalculateDeltat=&calDelt_RTP_TEOS;
        functions.fpCalculateNewEOSVars=&calNewTPKappaGamma_TEOS;
        functions.fpCalculateNewAV=&calNewQ0Q1Q2_RTP_TEOS;
        functions.fpModelWrite=&modelWrite_TEOS;
        functions.fpWriteWatchZones=&writeWatchZones_RTP_TEOS;
      }
      
      //velocity equation
      if(parameters.nTypeTurbulanceMod>0){
        functions.fpCalculateNewVelocities=&calNewVelocities_RTP_LES;
      }
      else{
        functions.fpCalculateNewVelocities=&calNewVelocities_RTP;
      }
      
      if(parameters.nTypeTurbulanceMod==1){
        functions.fpCalculateNewEddyVisc=&calNewEddyVisc_RTP_CN;
      }
      else if(parameters.nTypeTurbulanceMod==2){
        functions.fpCalculateNewEddyVisc=&calNewEddyVisc_RTP_SM;
      }
      
      //energy equation, eddy viscosity
      if(parameters.bAdiabatic){//adiabatic
        functions.fpCalculateNewEnergies=&calNewE_RTP_AD;
      }
      else{//non-adaibatic
        if(!parameters.bEOSGammaLaw){//needs a tabulated equation of state
          if(parameters.nTypeTurbulanceMod>0){//with turbulance model
            functions.fpCalculateNewEnergies=&calNewE_RTP_NA_LES;
          }
          else{
            functions.fpCalculateNewEnergies=&calNewE_RTP_NA;
          }
          if(implicit.nNumImplicitZones>0){//implicit, requires non-adiabatic
            functions.fpImplicitSolve=&implicitSolve_RTP;
            if(parameters.nTypeTurbulanceMod>0){
              functions.fpImplicitEnergyFunction=&dImplicitEnergyFunction_RTP_LES;
              functions.fpImplicitEnergyFunction_SB=&dImplicitEnergyFunction_RTP_LES_SB;
            }
            else{
              functions.fpImplicitEnergyFunction=&dImplicitEnergyFunction_RTP;
              functions.fpImplicitEnergyFunction_SB=&dImplicitEnergyFunction_RTP_SB;
            }
          }
        }
        else{//can't do a non-adiabatic calculation, with a gamma-law gas
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
            <<": User selected to do a non-adiabatic calculation but starting model uses a "
            <<"gamma-law gas. Starting model must use a tabulated equation of state in order to "
            <<"perform a non-adiabatic calculation.\n";
          throw exception2(ssTemp.str(),CALCULATION);
        }
      }
    }
    else if(grid.nNumDims==2){//use 2D
      functions.fpCalculateNewGridVelocities=&calNewU0_RT;
      functions.fpCalculateNewRadii=&calNewR;
      functions.fpCalculateNewDensities=&calNewD_RT;
      functions.fpCalculateAveDensities=&calNewDenave_RT;
      functions.fpUpdateLocalBoundaryVelocitiesNewGrid=&updateLocalBoundaryVelocitiesNewGrid_RT;
      if(parameters.bEOSGammaLaw){//use gamma law gas
        functions.fpCalculateDeltat=&calDelt_RT_GL;
        functions.fpCalculateNewEOSVars=&calNewP_GL;
        functions.fpCalculateNewAV=&calNewQ0Q1_RT_GL;
        functions.fpModelWrite=&modelWrite_GL;
        functions.fpWriteWatchZones=&writeWatchZones_RT_GL;
      }
      else{//use tabulated equation of state
        functions.fpCalculateDeltat=&calDelt_RT_TEOS;
        functions.fpCalculateNewEOSVars=&calNewTPKappaGamma_TEOS;
        functions.fpCalculateNewAV=&calNewQ0Q1_RT_TEOS;
        functions.fpModelWrite=&modelWrite_TEOS;
        functions.fpWriteWatchZones=&writeWatchZones_RT_TEOS;
      }
      
      //velocity equation
      if(parameters.nTypeTurbulanceMod>0){
        functions.fpCalculateNewVelocities=&calNewVelocities_RT_LES;
      }
      else{
        functions.fpCalculateNewVelocities=&calNewVelocities_RT;
      }
      
      if(parameters.nTypeTurbulanceMod==1){
        functions.fpCalculateNewEddyVisc=&calNewEddyVisc_RT_CN;
      }
      else if(parameters.nTypeTurbulanceMod==2){
        functions.fpCalculateNewEddyVisc=&calNewEddyVisc_RT_SM;
      }
      
      //energy equation, eddy viscosity
      if(parameters.bAdiabatic){//adiabatic
        functions.fpCalculateNewEnergies=&calNewE_RT_AD;
      }
      else{//non-adaibatic
        if(!parameters.bEOSGammaLaw){//needs a tabulated equation of state
          if(parameters.nTypeTurbulanceMod>0){//with turbulance model
            functions.fpCalculateNewEnergies=&calNewE_RT_NA_LES;
          }
          else{
            functions.fpCalculateNewEnergies=&calNewE_RT_NA;
          }
          if(implicit.nNumImplicitZones>0){//implicit, requires non-adiabatic
            functions.fpImplicitSolve=&implicitSolve_RT;
            if(parameters.nTypeTurbulanceMod>0){
              functions.fpImplicitEnergyFunction=&dImplicitEnergyFunction_RT_LES;
              functions.fpImplicitEnergyFunction_SB=&dImplicitEnergyFunction_RT_LES_SB;
            }
            else{
              functions.fpImplicitEnergyFunction=&dImplicitEnergyFunction_RT;
              functions.fpImplicitEnergyFunction_SB=&dImplicitEnergyFunction_RT_SB;
            }
          }
        }
        else{//can't do a non-adiabatic calculation, with a gamma-law gas
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
            <<": User selected to do a non-adiabatic calculation but starting model uses a "
            <<"gamma-law gas. Starting model must use a tabulated equation of state in order to "
            <<"perform a non-adiabatic calculation.\n";
          throw exception2(ssTemp.str(),CALCULATION);
        }
      }
    }
    else if(grid.nNumDims==1){//use 1D for all
      functions.fpCalculateNewGridVelocities=&calNewU0_R;
      functions.fpCalculateNewRadii=&calNewR;
      functions.fpCalculateNewDensities=&calNewD_R;
      functions.fpCalculateAveDensities=&calNewDenave_None;
      functions.fpUpdateLocalBoundaryVelocitiesNewGrid=&updateLocalBoundaryVelocitiesNewGrid_R;
      if(parameters.bEOSGammaLaw){//use gamma law gas
        functions.fpCalculateDeltat=&calDelt_R_GL;
        functions.fpCalculateNewEOSVars=&calNewP_GL;
        functions.fpCalculateNewAV=&calNewQ0_R_GL;
        functions.fpModelWrite=&modelWrite_GL;
        functions.fpWriteWatchZones=&writeWatchZones_R_GL;
      }
      else{//use tabulated equation of state
        functions.fpCalculateDeltat=&calDelt_R_TEOS;
        functions.fpCalculateNewEOSVars=&calNewTPKappaGamma_TEOS;
        functions.fpCalculateNewAV=&calNewQ0_R_TEOS;
        functions.fpModelWrite=&modelWrite_TEOS;
        functions.fpWriteWatchZones=&writeWatchZones_R_TEOS;
      }
      
      //velocity equation
      functions.fpCalculateNewVelocities=&calNewVelocities_R;
      
      //energy equation, eddy viscosity
      if(parameters.bAdiabatic){//adiabatic
        functions.fpCalculateNewEnergies=&calNewE_R_AD;
      }
      else{//non-adaibatic
        if(!parameters.bEOSGammaLaw){//needs a tabulated equation of state
          functions.fpCalculateNewEnergies=&calNewE_R_NA;
          functions.fpCalculateNewEddyVisc=&calNewEddyVisc_None;
          if(implicit.nNumImplicitZones>0){//implicit, requires non-adiabatic
            functions.fpImplicitSolve=&implicitSolve_R;
            functions.fpImplicitEnergyFunction=&dImplicitEnergyFunction_R;
            functions.fpImplicitEnergyFunction_SB=&dImplicitEnergyFunction_R_SB;
          }
        }
        else{//can't do a non-adiabatic calculation, with a gamma-law gas
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
            <<": User selected to do a non-adiabatic calculation but starting model uses a "
            <<"gamma-law gas. Starting model must use a tabulated equation of state in order to "
            <<"perform a non-adiabatic calculation.\n";
          throw exception2(ssTemp.str(),CALCULATION);
        }
      }
    }
  }
  if(!time.bVariableTimeStep){//if not a variable time step use a constant one
    functions.fpCalculateDeltat=&calDelt_CONST;
  }
}
