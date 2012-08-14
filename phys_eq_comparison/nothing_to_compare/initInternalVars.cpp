void initInternalVars(Grid &grid, ProcTop &procTop, Parameters &parameters){
  /** \warning \f$\Delta \theta \f$, \f$\Delta \phi \f$, \f$\sin \theta_{i,j,k} \f$,
      \f$\Delta \cos\theta_{i,j,k} \f$, all don't have the first zone calculated. At the moment
      this is a ghost cell that doesn't matter, but it may become a problem if calculations require
      this quantity. This is an issue for quantities that aren't updated in time, as those that are
      will have boundary cells updated with periodic boundary conditions.
  */
  
  //set values from equation of state, they don't care if 1D, 2D or 3D
  if(parameters.bEOSGammaLaw){
    calOldP_GL(grid,parameters);//set pressure
  }
  else{
    //initialize P, E, Kappa and Gamma, suitable for both 1D and 3D regions
    calOldPEKappaGamma_TEOS(grid, parameters);//set pressure, energy, kappa, and gamma
  }
  
  if(procTop.nRank!=0){
    if(grid.nNumDims>1){//both 2D, and 3D
      
      //initialize DCOSTHETAIJK
      for(int j=1;j<grid.nLocalGridDims[procTop.nRank][grid.nDCosThetaIJK][grid.nTheta]
        +2*grid.nNumGhostCells;j++){
        
        //calculate j for interface centered quantities
        int nJInt=j+grid.nCenIntOffset[1];
        
        grid.dLocalGridOld[grid.nDCosThetaIJK][0][j][0]
          =cos(grid.dLocalGridOld[grid.nTheta][0][nJInt-1][0])
          -cos(grid.dLocalGridOld[grid.nTheta][0][nJInt][0]);
      }
      
      //initialize DTHETA
      for(int j=1;j<grid.nLocalGridDims[procTop.nRank][grid.nDTheta][grid.nTheta]
        +2*grid.nNumGhostCells;j++){
        
        //calculate j for interface centered quantities
        int nJInt=j+grid.nCenIntOffset[1];
        
        grid.dLocalGridOld[grid.nDTheta][0][j][0]=grid.dLocalGridOld[grid.nTheta][0][nJInt][0]
          -grid.dLocalGridOld[grid.nTheta][0][nJInt-1][0];
      }
      
      //initialize SINTHETAIJK
      for(int j=1;j<grid.nLocalGridDims[procTop.nRank][grid.nSinThetaIJK][grid.nTheta]
        +2*grid.nNumGhostCells;j++){
        
        //calculate j for interface centered quantities
        int nJInt=j+grid.nCenIntOffset[1];
        
        grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          =sin((grid.dLocalGridOld[grid.nTheta][0][nJInt][0]
          +grid.dLocalGridOld[grid.nTheta][0][nJInt-1][0])*0.5);
      }
    
      //initialize SINTHETAIJP1HALFK
      for(int j=0;j<grid.nLocalGridDims[procTop.nRank][grid.nSinThetaIJp1halfK][grid.nTheta]
        +2*grid.nNumGhostCells;j++){
        
        grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]
          =sin(grid.dLocalGridOld[grid.nTheta][0][j][0]);
      }
    }
    if(grid.nNumDims==2){//only 2D
      
      //initialize DENAVE
      calOldDenave_RT(grid);
    }
    if(grid.nNumDims==3){//only 3D
      
      //initialize DPHI
      for(int k=1;k<grid.nLocalGridDims[procTop.nRank][grid.nDPhi][grid.nPhi]
        +2*grid.nNumGhostCells;k++){
        
        //calculate k for interface centered quantities
        int nKInt=k+grid.nCenIntOffset[2];
        
        grid.dLocalGridOld[grid.nDPhi][0][0][k]=grid.dLocalGridOld[grid.nPhi][0][0][nKInt]
          -grid.dLocalGridOld[grid.nPhi][0][0][nKInt-1];
      }
      
      //initialize DENAVE
      calOldDenave_RTP(grid);
    }
    if(grid.nNumDims>2||(grid.nNumDims>1&&parameters.nTypeTurbulanceMod>0)){/* Need these for 3D and
      2D calculations that use a turbulance model*/
      
      //initialize COTTHETAIJP1HALFK
      for(int j=0;j<grid.nLocalGridDims[procTop.nRank][grid.nCotThetaIJp1halfK][grid.nTheta]
        +2*grid.nNumGhostCells;j++){
        
        grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]
          =1.0/tan(grid.dLocalGridOld[grid.nTheta][0][j][0]);
      }
      
      //initialize COTTHETAIJK
      for(int j=1;j<grid.nLocalGridDims[procTop.nRank][grid.nCotThetaIJK][grid.nTheta]
        +2*grid.nNumGhostCells;j++){
        
        //calculate j for interface centered quantities
        int nJInt=j+grid.nCenIntOffset[1];
        double dTheta_ijk=(grid.dLocalGridOld[grid.nTheta][0][nJInt][0]
          +grid.dLocalGridOld[grid.nTheta][0][nJInt-1][0])*0.5;
        grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]=1.0/tan(dTheta_ijk);
      }
    }
      
    //initialize Q (Artificial Viscosity), donor fraction, and maximum convective velocity
    if(parameters.bEOSGammaLaw){//Gamma law
      if(grid.nNumDims==1){
        calOldQ0_R_GL(grid,parameters);
        initDonorFracAndMaxConVel_R_GL(grid,parameters);
      }
      if(grid.nNumDims==2){
        calOldQ0Q1_RT_GL(grid,parameters);
        initDonorFracAndMaxConVel_RT_GL(grid,parameters);
      }
      if(grid.nNumDims==3){
        calOldQ0Q1Q2_RTP_GL(grid,parameters);
        initDonorFracAndMaxConVel_RTP_GL(grid,parameters);
      }
    }
    else{//tabulated equation of state
      if(grid.nNumDims==1){
        calOldQ0_R_TEOS(grid,parameters);
        initDonorFracAndMaxConVel_R_TEOS(grid,parameters);
      }
      if(grid.nNumDims==2){
        calOldQ0Q1_RT_TEOS(grid,parameters);
        initDonorFracAndMaxConVel_RT_TEOS(grid,parameters);
      }
      if(grid.nNumDims==3){
        calOldQ0Q1Q2_RTP_TEOS(grid,parameters);
        initDonorFracAndMaxConVel_RTP_TEOS(grid,parameters);
      }
    }
    
    //if using a turblance model, initilize the eddy viscosity
    if(parameters.nTypeTurbulanceMod==1){
      if(grid.nNumDims==1){
        calOldEddyVisc_R_CN(grid,parameters);
      }
      if(grid.nNumDims==2){
        calOldEddyVisc_RT_CN(grid,parameters);
      }
      if(grid.nNumDims==3){
        calOldEddyVisc_RTP_CN(grid,parameters);
      }
    }
    if(parameters.nTypeTurbulanceMod==2){
      if(grid.nNumDims==1){
        calOldEddyVisc_R_SM(grid,parameters);
      }
      if(grid.nNumDims==2){
        calOldEddyVisc_RT_SM(grid,parameters);
      }
      if(grid.nNumDims==3){
        calOldEddyVisc_RTP_SM(grid,parameters);
      }
    }
  }
  else{//processor 0, always 1D
    
    //initialize DENAVE only if number of dimensions greater than one, this will allow the grid in
    //the 1D region to be compatible with the grid in the 3D region for message passing perposes
    if(grid.nNumDims>1){
      calOldDenave_R(grid);
    }
    
    //initialize Q (Artificial Viscosity), donor fraction, and maximum convective velocity
    if(parameters.bEOSGammaLaw){
      calOldQ0_R_GL(grid,parameters);
      initDonorFracAndMaxConVel_R_GL(grid,parameters);
    }
    else{
      calOldQ0_R_TEOS(grid,parameters);
      initDonorFracAndMaxConVel_R_TEOS(grid,parameters);
    }
  }
}
