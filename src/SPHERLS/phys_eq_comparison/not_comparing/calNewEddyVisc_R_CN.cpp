void calNewEddyVisc_R_CN(Grid &grid, Parameters &parameters){
  
  int i;
  int j;
  int k;
  int nIInt;//i+1/2 index
  double dR_ip1half_np1half;
  double dR_im1half_np1half;
  double dDelR_i_np1half;
  double dConstant=parameters.dEddyViscosity;
  double dLengthScaleSq;
  
  //main grid explicit
  for(i=grid.nStartUpdateExplicit[grid.nEddyVisc][0];
    i<grid.nEndUpdateExplicit[grid.nEddyVisc][0];i++){
      
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_np1half=grid.dLocalGridNew[grid.nR][nIInt][0][0];
    dR_im1half_np1half=grid.dLocalGridNew[grid.nR][nIInt-1][0][0];
    dDelR_i_np1half=dR_ip1half_np1half-dR_im1half_np1half;
    
    for(j=grid.nStartUpdateExplicit[grid.nEddyVisc][1];
      j<grid.nEndUpdateExplicit[grid.nEddyVisc][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nEddyVisc][2];
        k<grid.nEndUpdateExplicit[grid.nEddyVisc][2];k++){
        dLengthScaleSq=dDelR_i_np1half*dDelR_i_np1half;
        grid.dLocalGridNew[grid.nEddyVisc][i][j][k]=dLengthScaleSq*dConstant
          *parameters.dMaxConvectiveVelocity/1.0e6;
      }
    }
  }
  
  //outter radial ghost cells,explicit
  for(i=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][0];i++){
      
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_np1half=grid.dLocalGridNew[grid.nR][nIInt][0][0];
    dR_im1half_np1half=grid.dLocalGridNew[grid.nR][nIInt-1][0][0];
    dDelR_i_np1half=dR_ip1half_np1half-dR_im1half_np1half;
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][1];j++){
      for(k=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][2];k++){
        dLengthScaleSq=dDelR_i_np1half*dDelR_i_np1half;
        grid.dLocalGridNew[grid.nEddyVisc][i][j][k]=dLengthScaleSq*dConstant
          *parameters.dMaxConvectiveVelocity/1.0e6;
      }
    }
  }
}
