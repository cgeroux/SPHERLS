void calOldEddyVisc_R_SM(Grid &grid, Parameters &parameters){
  
  int nIInt;//i+1/2 index
  double d1;//term 1
  double dA;
  double dR_ip1half_np1half;
  double dR_im1half_np1half;
  double dR_i_np1half;
  double dU_ip1halfjk_np1half;
  double dU_im1halfjk_np1half;
  double dDelR_i_np1half;
  double dConstantSq=parameters.dEddyViscosity*parameters.dEddyViscosity/pow(2.0,0.5);
  double dLengthScaleSq;
  double dTerms;
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
        
        dLengthScaleSq=dDelR_i_np1half*dDelR_i_np1half;
        
        //interpolate
        dU_ip1halfjk_np1half=grid.dLocalGridOld[grid.nU][nIInt][j][k];
        dU_im1halfjk_np1half=grid.dLocalGridOld[grid.nU][nIInt-1][j][k];
        
        //calculate dA term
        d1=(dU_ip1halfjk_np1half-dU_im1halfjk_np1half)/(dR_ip1half_np1half-dR_im1half_np1half);
        dA=2.0*d1*d1;
        
        dTerms=dA;
        grid.dLocalGridOld[grid.nEddyVisc][i][j][k]=dLengthScaleSq*dConstantSq
          *grid.dLocalGridOld[grid.nD][i][j][k]*pow(dTerms,0.5);
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
        
        dLengthScaleSq=dDelR_i_np1half*dDelR_i_np1half;
        
        //interpolate
        dU_ip1halfjk_np1half=grid.dLocalGridOld[grid.nU][nIInt][j][k];
        dU_im1halfjk_np1half=grid.dLocalGridOld[grid.nU][nIInt-1][j][k];
        
        //calculate dA term
        d1=(dU_ip1halfjk_np1half-dU_im1halfjk_np1half)/(dR_ip1half_np1half-dR_im1half_np1half);
        dA=2.0*d1*d1;
        
        dTerms=dA;
        grid.dLocalGridOld[grid.nEddyVisc][i][j][k]=dLengthScaleSq*dConstantSq
          *grid.dLocalGridOld[grid.nD][i][j][k]*pow(dTerms,0.5);
      }
    }
  }
}
