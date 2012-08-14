void calNewQ0_R_TEOS(Grid& grid,Parameters &parameters){
  
  double dA_sq=parameters.dA*parameters.dA;
  double dDVDt;
  double dA_ip1half;
  double dA_im1half;
  double dR_i_n;
  int nIInt;
  double dC;
  double dDVDtThreshold;
  double dRSq_i_n;
  double dDVDt_mthreshold;
  int i;
  
  //in the main grid
  for(i=grid.nStartUpdateExplicit[grid.nQ0][0];i<grid.nEndUpdateExplicit[grid.nQ0][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_i_n=(grid.dLocalGridNew[grid.nR][nIInt][0][0]+grid.dLocalGridNew[grid.nR][nIInt-1][0][0])*0.5;
    dRSq_i_n=dR_i_n*dR_i_n;
    dA_ip1half=grid.dLocalGridNew[grid.nR][nIInt][0][0]*grid.dLocalGridNew[grid.nR][nIInt][0][0];
    dA_im1half=grid.dLocalGridNew[grid.nR][nIInt-1][0][0]
      *grid.dLocalGridNew[grid.nR][nIInt-1][0][0];
    
    for(j=grid.nStartUpdateExplicit[grid.nQ0][1];j<grid.nEndUpdateExplicit[grid.nQ0][1];j++){
      
      for(k=grid.nStartUpdateExplicit[grid.nQ0][2];k<grid.nEndUpdateExplicit[grid.nQ0][2];k++){
        
        //calculate threshold compression to turn viscosity on at
        dC=sqrt(grid.dLocalGridNew[grid.nGamma][i][j][k]*(grid.dLocalGridNew[grid.nP][i][j][k])
          /grid.dLocalGridNew[grid.nD][i][j][k]);
        dDVDtThreshold=parameters.dAVThreshold*dC;
        
        //calculate Q0
        dDVDt=(dA_ip1half*grid.dLocalGridNew[grid.nU][nIInt][j][k]
          -dA_im1half*grid.dLocalGridNew[grid.nU][nIInt-1][j][k])/dRSq_i_n;
        if(dDVDt<-1.0*dDVDtThreshold){//being compressed
          dDVDt_mthreshold=dDVDt+dDVDtThreshold;
          grid.dLocalGridNew[grid.nQ0][i][j][k]=dA_sq*grid.dLocalGridNew[grid.nD][i][j][k]
            *dDVDt_mthreshold*dDVDt_mthreshold;
        }
        else{
          grid.dLocalGridNew[grid.nQ0][i][j][k]=0.0;
        }
      }
    }
  }
  
  //outter ghost region
  for(i=grid.nStartGhostUpdateExplicit[grid.nQ0][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nQ0][0][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_i_n=(grid.dLocalGridNew[grid.nR][nIInt][0][0]+grid.dLocalGridNew[grid.nR][nIInt-1][0][0])*0.5;
    dRSq_i_n=dR_i_n*dR_i_n;
    dA_ip1half=grid.dLocalGridNew[grid.nR][nIInt][0][0]*grid.dLocalGridNew[grid.nR][nIInt][0][0];
    dA_im1half=grid.dLocalGridNew[grid.nR][nIInt-1][0][0]
      *grid.dLocalGridNew[grid.nR][nIInt-1][0][0];
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nQ0][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nQ0][0][1];j++){
      for(k=grid.nStartGhostUpdateExplicit[grid.nQ0][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nQ0][0][2];k++){
        
        //calculate threshold compression to turn viscosity on at
        dC=sqrt(grid.dLocalGridNew[grid.nGamma][i][j][k]*(grid.dLocalGridNew[grid.nP][i][j][k])
          /grid.dLocalGridNew[grid.nD][i][j][k]);
        dDVDtThreshold=parameters.dAVThreshold*dC;
        
        //calculate Q0
        dDVDt=(dA_ip1half*grid.dLocalGridNew[grid.nU][nIInt][j][k]
          -dA_im1half*grid.dLocalGridNew[grid.nU][nIInt-1][j][k])/dRSq_i_n;
        if(dDVDt<-1.0*dDVDtThreshold){//being compressed
          dDVDt_mthreshold=dDVDt+dDVDtThreshold;
          grid.dLocalGridNew[grid.nQ0][i][j][k]=dA_sq*grid.dLocalGridNew[grid.nD][i][j][k]
            *dDVDt_mthreshold*dDVDt_mthreshold;
        }
        else{
          grid.dLocalGridNew[grid.nQ0][i][j][k]=0.0;
        }
      }
    }
  }
}
