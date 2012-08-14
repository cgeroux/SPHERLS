void calNewQ0_R_TEOS(Grid& grid,Parameters &parameters){
  
  double dA_sq=parameters.dA*parameters.dA;
  double dDVDt;
  double dA_ip1half;
  double dA_im1half;
  double dR_i;
  int nIInt;
  double dC;
  double dDVDtThreshold;
  double dR_i_sq;
  double dDVDt_mthreshold;
  int i;
  
  //in the main grid
  for(i=grid.nStartUpdateExplicit[grid.nQ0][0];i<grid.nEndUpdateExplicit[grid.nQ0][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    
    //calculate volume change
    dA_ip1half=grid.dLocalGridNew[grid.nR][nIInt][0][0]*grid.dLocalGridNew[grid.nR][nIInt][0][0];
    dA_im1half=grid.dLocalGridNew[grid.nR][nIInt-1][0][0]
      *grid.dLocalGridNew[grid.nR][nIInt-1][0][0];
    dR_i=(grid.dLocalGridNew[grid.nR][nIInt][0][0]+grid.dLocalGridNew[grid.nR][nIInt-1][0][0])*0.5;
    dR_i_sq=dR_i*dR_i;
    dDVDt=(dA_ip1half*grid.dLocalGridNew[grid.nU][nIInt][0][0]
      -dA_im1half*grid.dLocalGridNew[grid.nU][nIInt-1][0][0])/dR_i_sq;
    
    //calculate threshold compression to turn viscosity on at
    dC=sqrt(grid.dLocalGridNew[grid.nGamma][i][0][0]*(grid.dLocalGridNew[grid.nP][i][0][0])
      /grid.dLocalGridNew[grid.nD][i][0][0]);
    dDVDtThreshold=parameters.dAVThreshold*dC;
    
    if(dDVDt<-1.0*dDVDtThreshold){//being compressed
      dDVDt_mthreshold=dDVDt+dDVDtThreshold;
      grid.dLocalGridNew[grid.nQ0][i][0][0]=dA_sq*grid.dLocalGridNew[grid.nD][i][0][0]
        *dDVDt_mthreshold*dDVDt_mthreshold;
    }
    else{
      grid.dLocalGridNew[grid.nQ0][i][0][0]=0.0;
    }
  }
  
  //outter ghost region
  for(i=grid.nStartGhostUpdateExplicit[grid.nQ0][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nQ0][0][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    
    //calculate volume change
    dA_ip1half=grid.dLocalGridNew[grid.nR][nIInt][0][0]*grid.dLocalGridNew[grid.nR][nIInt][0][0];
    dA_im1half=grid.dLocalGridNew[grid.nR][nIInt-1][0][0]
      *grid.dLocalGridNew[grid.nR][nIInt-1][0][0];
    dR_i=(grid.dLocalGridNew[grid.nR][nIInt][0][0]+grid.dLocalGridNew[grid.nR][nIInt-1][0][0])*0.5;
    dR_i_sq=dR_i*dR_i;
    dDVDt=(dA_ip1half*grid.dLocalGridNew[grid.nU][nIInt][0][0]
      -dA_im1half*grid.dLocalGridNew[grid.nU][nIInt-1][0][0])/dR_i_sq;
    
    //calculate threshold compression to turn viscosity on at
    dC=sqrt(grid.dLocalGridNew[grid.nGamma][i][0][0]*(grid.dLocalGridNew[grid.nP][i][0][0])
      /grid.dLocalGridNew[grid.nD][i][0][0]);
    dDVDtThreshold=parameters.dAVThreshold*dC;
    
    if(dDVDt<-1.0*dDVDtThreshold){//being compressed
      dDVDt_mthreshold=dDVDt+dDVDtThreshold;
      grid.dLocalGridNew[grid.nQ0][i][0][0]=dA_sq*grid.dLocalGridNew[grid.nD][i][0][0]
        *dDVDt_mthreshold*dDVDt_mthreshold;
    }
    else{
      grid.dLocalGridNew[grid.nQ0][i][0][0]=0.0;
    }
  }
}
