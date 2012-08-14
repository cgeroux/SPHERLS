void calNewEddyVisc_R_SM(Grid &grid, Parameters &parameters){
  
  int i;
  int j;
  int k;
  int nIInt;//i+1/2 index
  double d1;
  double dA;
  double dTerms;
  double dR_ip1half_np1half;
  double dR_im1half_np1half;
  double dR_i_np1half;
  double dU_ip1halfjk_np1half;
  double dU_im1halfjk_np1half;
  double dDelR_i_np1half;
  double dConstantSq=parameters.dEddyViscosity*parameters.dEddyViscosity/pow(2.0,0.5);
  double dLengthScaleSq;
  
  //main grid explicit
  for(i=grid.nStartUpdateExplicit[grid.nEddyVisc][0];
    i<grid.nEndUpdateExplicit[grid.nEddyVisc][0];i++){
      
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_np1half=(grid.dLocalGridNew[grid.nR][nIInt][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt][0][0])*0.5;
    dR_im1half_np1half=(grid.dLocalGridNew[grid.nR][nIInt-1][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt-1][0][0])*0.5;
    dR_i_np1half=(dR_ip1half_np1half+dR_im1half_np1half)*0.5;
    dDelR_i_np1half=dR_ip1half_np1half-dR_im1half_np1half;
    
    for(j=grid.nStartUpdateExplicit[grid.nEddyVisc][1];
      j<grid.nEndUpdateExplicit[grid.nEddyVisc][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nEddyVisc][2];
        k<grid.nEndUpdateExplicit[grid.nEddyVisc][2];k++){
        
        dLengthScaleSq=dDelR_i_np1half*dDelR_i_np1half;
        
        //interpolate
        dU_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k];
        dU_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k];
        
        //term 1
        d1=((dU_ip1halfjk_np1half-grid.dLocalGridNew[grid.nU0][nIInt][0][0])
          -(dU_im1halfjk_np1half-grid.dLocalGridNew[grid.nU0][nIInt-1][0][0]))
          /(dR_ip1half_np1half-dR_im1half_np1half);
        dA=2.0*d1*d1;
        
        dTerms=dA;
        grid.dLocalGridNew[grid.nEddyVisc][i][j][k]=dConstantSq*dLengthScaleSq
          *grid.dLocalGridNew[grid.nD][i][j][k]*pow(dTerms,0.5);
      }
    }
  }
  
  //outter radial ghost cells,explicit
  for(i=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_np1half=(grid.dLocalGridNew[grid.nR][nIInt][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt][0][0])*0.5;
    dR_im1half_np1half=(grid.dLocalGridNew[grid.nR][nIInt-1][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt-1][0][0])*0.5;
    dR_i_np1half=(dR_ip1half_np1half+dR_im1half_np1half)*0.5;
    dDelR_i_np1half=dR_ip1half_np1half-dR_im1half_np1half;
    for(j=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][1];j++){
      for(k=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][2];k++){
        
        dLengthScaleSq=dDelR_i_np1half*dDelR_i_np1half;
        
        //interpolate
        dU_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k];
        dU_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k];
        
        //term 1
        d1=((dU_ip1halfjk_np1half-grid.dLocalGridNew[grid.nU0][nIInt][0][0])
          -(dU_im1halfjk_np1half-grid.dLocalGridNew[grid.nU0][nIInt-1][0][0]))
          /(dR_ip1half_np1half-dR_im1half_np1half);
        dA=2.0*d1*d1;
        
        dTerms=dA;
        grid.dLocalGridNew[grid.nEddyVisc][i][j][k]=dConstantSq*dLengthScaleSq
          *grid.dLocalGridNew[grid.nD][i][j][k]*pow(dTerms,0.5);
      }
    }
  }
}
