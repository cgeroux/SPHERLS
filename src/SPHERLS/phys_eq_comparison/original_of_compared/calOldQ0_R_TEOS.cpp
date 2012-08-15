void calOldQ0_R_TEOS(Grid& grid, Parameters &parameters){
  
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
    dA_ip1half=grid.dLocalGridOld[grid.nR][nIInt][0][0]*grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dA_im1half=grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      *grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_i=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])*0.5;
    dR_i_sq=dR_i*dR_i;
    dDVDt=(dA_ip1half*grid.dLocalGridOld[grid.nU][nIInt][0][0]
      -dA_im1half*grid.dLocalGridOld[grid.nU][nIInt-1][0][0])/dR_i_sq;
    
    //calculate threshold compression to turn viscosity on at
    dC=sqrt(grid.dLocalGridOld[grid.nGamma][i][0][0]*(grid.dLocalGridOld[grid.nP][i][0][0])
      /grid.dLocalGridOld[grid.nD][i][0][0]);
    dDVDtThreshold=parameters.dAVThreshold*dC;
    
    if(dDVDt<-1.0*dDVDtThreshold){//being compressed
      dDVDt_mthreshold=dDVDt+dDVDtThreshold;
      grid.dLocalGridOld[grid.nQ0][i][0][0]=dA_sq*grid.dLocalGridOld[grid.nD][i][0][0]
        *dDVDt_mthreshold*dDVDt_mthreshold;
    }
    else{
      grid.dLocalGridOld[grid.nQ0][i][0][0]=0.0;
    }
  }
  
  //outter ghost region
  for(i=grid.nStartGhostUpdateExplicit[grid.nQ0][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nQ0][0][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    
    //calculate volume change
    dA_ip1half=grid.dLocalGridOld[grid.nR][nIInt][0][0]*grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dA_im1half=grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      *grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_i=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])*0.5;
    dR_i_sq=dR_i*dR_i;
    dDVDt=(dA_ip1half*grid.dLocalGridOld[grid.nU][nIInt][0][0]
      -dA_im1half*grid.dLocalGridOld[grid.nU][nIInt-1][0][0])/dR_i_sq;
    
    //calculate threshold compression to turn viscosity on at
    dC=sqrt(grid.dLocalGridOld[grid.nGamma][i][0][0]*(grid.dLocalGridOld[grid.nP][i][0][0])
      /grid.dLocalGridOld[grid.nD][i][0][0]);
    dDVDtThreshold=parameters.dAVThreshold*dC;
    
    if(dDVDt<-1.0*dDVDtThreshold){//being compressed
      dDVDt_mthreshold=dDVDt+dDVDtThreshold;
      grid.dLocalGridOld[grid.nQ0][i][0][0]=dA_sq*grid.dLocalGridOld[grid.nD][i][0][0]
        *dDVDt_mthreshold*dDVDt_mthreshold;
    }
    else{
      grid.dLocalGridOld[grid.nQ0][i][0][0]=0.0;
    }
  }
}
