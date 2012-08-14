void calOldEddyVisc_RT_CN(Grid &grid, Parameters &parameters){
  
  int nIInt;//i+1/2 index
  double dR_ip1half_np1half;
  double dR_im1half_np1half;
  double dR_i_np1half;
  double dDelR_i_np1half;
  double dLengthScaleSq;
  double dConstant=parameters.dEddyViscosity;
  
  //main grid explicit
  for(int i=grid.nStartUpdateExplicit[grid.nEddyVisc][0];
    i<grid.nEndUpdateExplicit[grid.nEddyVisc][0];i++){
      
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_np1half=grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dR_im1half_np1half=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_i_np1half=(dR_ip1half_np1half+dR_im1half_np1half)*0.5;
    dDelR_i_np1half=dR_ip1half_np1half-dR_im1half_np1half;
      
    for(int j=grid.nStartUpdateExplicit[grid.nEddyVisc][1];
      j<grid.nEndUpdateExplicit[grid.nEddyVisc][1];j++){
      for(int k=grid.nStartUpdateExplicit[grid.nEddyVisc][2];
        k<grid.nEndUpdateExplicit[grid.nEddyVisc][2];k++){
        dLengthScaleSq=dDelR_i_np1half*dR_i_np1half*grid.dLocalGridOld[grid.nDTheta][0][j][0];
        grid.dLocalGridNew[grid.nEddyVisc][i][j][k]=dLengthScaleSq*dConstant
          *parameters.dMaxConvectiveVelocity/1.0e6;
      }
    }
  }
  
  //inner radial ghost cells, set to zero
  for(int i=0;i<grid.nNumGhostCells;i++){
    for(int j=grid.nStartUpdateExplicit[grid.nEddyVisc][1];
      j<grid.nEndUpdateExplicit[grid.nEddyVisc][1];j++){
      for(int k=grid.nStartUpdateExplicit[grid.nEddyVisc][2];
        k<grid.nEndUpdateExplicit[grid.nEddyVisc][2];k++){
        grid.dLocalGridOld[grid.nEddyVisc][i][j][k]=0.0;
      }
    }
  }
  
  //outter radial ghost cells,explicit
  for(int i=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][0];i++){
      
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_np1half=grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dR_im1half_np1half=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_i_np1half=(dR_ip1half_np1half+dR_im1half_np1half)*0.5;
    dDelR_i_np1half=dR_ip1half_np1half-dR_im1half_np1half;
    for(int j=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][1];j++){
      for(int k=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][2];k++){
        dLengthScaleSq=dDelR_i_np1half*dR_i_np1half*grid.dLocalGridOld[grid.nDTheta][0][j][0];
        grid.dLocalGridNew[grid.nEddyVisc][i][j][k]=dLengthScaleSq*dConstant
          *parameters.dMaxConvectiveVelocity/1.0e6;
      }
    }
  }
}
