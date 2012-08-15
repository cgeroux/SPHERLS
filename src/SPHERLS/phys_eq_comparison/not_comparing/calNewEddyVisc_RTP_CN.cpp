void calNewEddyVisc_RTP_CN(Grid &grid, Parameters &parameters){
  int i;
  int j;
  int k;
  int nIInt;//i+1/2 index
  double dR_ip1half_np1half;
  double dR_im1half_np1half;
  double dR_i_np1half;
  double dR_i_np1half_Sq;
  double dDelR_i_np1half;
  double dConstant=parameters.dEddyViscosity;
  double dLengthScaleSq;
  
  //main grid
  for(i=grid.nStartUpdateExplicit[grid.nEddyVisc][0];
    i<grid.nEndUpdateExplicit[grid.nEddyVisc][0];i++){
      
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_np1half=grid.dLocalGridNew[grid.nR][nIInt][0][0];
    dR_im1half_np1half=grid.dLocalGridNew[grid.nR][nIInt-1][0][0];
    dR_i_np1half=(dR_ip1half_np1half+dR_im1half_np1half)*0.5;
    dR_i_np1half_Sq=dR_i_np1half*dR_i_np1half;
    dDelR_i_np1half=dR_ip1half_np1half-dR_im1half_np1half;
    
    for(j=grid.nStartUpdateExplicit[grid.nEddyVisc][1];
      j<grid.nEndUpdateExplicit[grid.nEddyVisc][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nEddyVisc][2];
        k<grid.nEndUpdateExplicit[grid.nEddyVisc][2];k++){
        dLengthScaleSq=dR_i_np1half_Sq*dDelR_i_np1half*grid.dLocalGridOld[grid.nDTheta][0][j][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dLengthScaleSq=pow(dLengthScaleSq,0.666666666666666);
        grid.dLocalGridNew[grid.nEddyVisc][i][j][k]=dConstant*dLengthScaleSq
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
    dR_i_np1half=(dR_ip1half_np1half+dR_im1half_np1half)*0.5;
    dDelR_i_np1half=dR_ip1half_np1half-dR_im1half_np1half;
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][1];j++){
      for(k=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][2];k++){
        dLengthScaleSq=dR_i_np1half_Sq*dDelR_i_np1half*grid.dLocalGridOld[grid.nDTheta][0][j][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dLengthScaleSq=pow(dLengthScaleSq,0.666666666666666);
        grid.dLocalGridNew[grid.nEddyVisc][i][j][k]=dConstant*dLengthScaleSq
          *parameters.dMaxConvectiveVelocity/1.0e6;
      }
    }
  }
}
