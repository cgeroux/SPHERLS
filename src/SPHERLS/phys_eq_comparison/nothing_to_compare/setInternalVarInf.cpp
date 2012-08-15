void setInternalVarInf(Grid &grid, Parameters &parameters){
  
  //allocate space for internal variable infos
  for(int n=grid.nNumVars;n<grid.nNumVars+grid.nNumIntVars;n++){
    grid.nVariables[n]=new int[4];//+1 because of keeping track of time info
  }
  
  //VARIABLE INFOS.
  
  //P
  grid.nVariables[grid.nP][0]=0;//r centered
  grid.nVariables[grid.nP][1]=0;//centered in theta
  grid.nVariables[grid.nP][2]=0;//centered in phi
  grid.nVariables[grid.nP][3]=1;//updated with time
  
  //Q0
  grid.nVariables[grid.nQ0][0]=0;//r centered
  grid.nVariables[grid.nQ0][1]=0;//theta centered
  grid.nVariables[grid.nQ0][2]=0;//phi centered
  grid.nVariables[grid.nQ0][3]=1;//updated with time
  
  //nDonorCellFrac
  grid.nVariables[grid.nDonorCellFrac][0]=0;//r centered
  grid.nVariables[grid.nDonorCellFrac][1]=-1;//theta centered
  grid.nVariables[grid.nDonorCellFrac][2]=-1;//phi centered
  grid.nVariables[grid.nDonorCellFrac][3]=1;//updated with time
  
  if(!parameters.bEOSGammaLaw){  //If using TEOS these are extra internal variables
    
    //E
    grid.nVariables[grid.nE][0]=0;//r centered
    grid.nVariables[grid.nE][1]=0;//theta centered
    grid.nVariables[grid.nE][2]=0;//phi centered
    grid.nVariables[grid.nE][3]=1;//updated with time
    
    //KAPPA
    grid.nVariables[grid.nKappa][0]=0;//r centered
    grid.nVariables[grid.nKappa][1]=0;//centered in theta
    grid.nVariables[grid.nKappa][2]=0;//centered in phi
    grid.nVariables[grid.nKappa][3]=1;//updated with time
  
    //GAMMA
    grid.nVariables[grid.nGamma][0]=0;//r centered
    grid.nVariables[grid.nGamma][1]=0;//centered in theta
    grid.nVariables[grid.nGamma][2]=0;//centered in phi
    grid.nVariables[grid.nGamma][3]=1;//updated with time
    
    //Eddy viscosity
    if(parameters.nTypeTurbulanceMod>0){
      grid.nVariables[grid.nEddyVisc][0]=0;//r centered
      grid.nVariables[grid.nEddyVisc][1]=0;//centered in theta
      grid.nVariables[grid.nEddyVisc][2]=0;//centered in phi
      grid.nVariables[grid.nEddyVisc][3]=1;//updated with time
    }
  }
  if(grid.nNumDims>1){//not defined for 1D
    
    //DENAVE
    grid.nVariables[grid.nDenAve][0]=0;//r centered
    grid.nVariables[grid.nDenAve][1]=-1;//not defined in theta
    grid.nVariables[grid.nDenAve][2]=-1;//not defined in phi
    grid.nVariables[grid.nDenAve][3]=1;//updated with time
    
    //DCOSTHETAIJK
    grid.nVariables[grid.nDCosThetaIJK][0]=-1;//not defined in r
    grid.nVariables[grid.nDCosThetaIJK][1]=0;//theta centered
    grid.nVariables[grid.nDCosThetaIJK][2]=-1;//not defined in phi
    grid.nVariables[grid.nDCosThetaIJK][3]=0;//not updated with time
    
    //Q1
    grid.nVariables[grid.nQ1][0]=0;//r centered
    grid.nVariables[grid.nQ1][1]=0;//theta centered
    grid.nVariables[grid.nQ1][2]=0;//phi centered
    grid.nVariables[grid.nQ1][3]=1;//updated with time
    
    //DTHETA
    grid.nVariables[grid.nDTheta][0]=-1;//not defined in r
    grid.nVariables[grid.nDTheta][1]=0;//theta centered
    grid.nVariables[grid.nDTheta][2]=-1;//not defined in phi
    grid.nVariables[grid.nDTheta][3]=0;//not updated with time
    
    //SINTHETAIJK
    grid.nVariables[grid.nSinThetaIJK][0]=-1;//not defined in r
    grid.nVariables[grid.nSinThetaIJK][1]=0;//theta centered
    grid.nVariables[grid.nSinThetaIJK][2]=-1;//not defined in phi
    grid.nVariables[grid.nSinThetaIJK][3]=0;//not updated with time
    
    //SINTHETAIJP1HALFK
    grid.nVariables[grid.nSinThetaIJp1halfK][0]=-1;//not defined in r
    grid.nVariables[grid.nSinThetaIJp1halfK][1]=1;//theta interface
    grid.nVariables[grid.nSinThetaIJp1halfK][2]=-1;//not defined in phi
    grid.nVariables[grid.nSinThetaIJp1halfK][3]=0;//not updated with time
  }
  if(grid.nNumDims>2){//not defined for 1D or 2D
    
    //DPHI
    grid.nVariables[grid.nDPhi][0]=-1;//not defined in r
    grid.nVariables[grid.nDPhi][1]=-1;//not defined in theta
    grid.nVariables[grid.nDPhi][2]=0;//phi centered
    grid.nVariables[grid.nDPhi][3]=0;//not updated with time
    
    //Q2
    grid.nVariables[grid.nQ2][0]=0;//r centered
    grid.nVariables[grid.nQ2][1]=0;//theta centered
    grid.nVariables[grid.nQ2][2]=0;//phi centered
    grid.nVariables[grid.nQ2][3]=1;//updated with time
  }
  if(grid.nNumDims>2||(grid.nNumDims>1&&parameters.nTypeTurbulanceMod>0)){/* Need these for 3D 
    calculations and 2D calculations that use a turbulance model*/
    
    //COTTHETAIJP1HALFK
    grid.nVariables[grid.nCotThetaIJp1halfK][0]=-1;//not defined in r
    grid.nVariables[grid.nCotThetaIJp1halfK][1]=1;//theta interface
    grid.nVariables[grid.nCotThetaIJp1halfK][2]=-1;//not defined in phi
    grid.nVariables[grid.nCotThetaIJp1halfK][3]=0;//not updated with time
    
    //COTTHETAIJK
    grid.nVariables[grid.nCotThetaIJK][0]=-1;//not defined in r
    grid.nVariables[grid.nCotThetaIJK][1]=0;//theta center
    grid.nVariables[grid.nCotThetaIJK][2]=-1;//not defined in phi
    grid.nVariables[grid.nCotThetaIJK][3]=0;//not updated with time
  }
  
  //adjust based on number of dimensions
  if(grid.nNumDims<3){
    for(int n=grid.nNumVars;n<grid.nNumVars+grid.nNumIntVars;n++){
      grid.nVariables[n][2]=-1;//not defined in phi direction if not 3D
    }
  }
  if(grid.nNumDims<2){
    for(int n=grid.nNumVars;n<grid.nNumVars+grid.nNumIntVars;n++){
      grid.nVariables[n][1]=-1;//not defined in theta direction if not 2D
    }
  }
}
