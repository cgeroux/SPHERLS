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
    dR_i=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])*0.5;
    dR_i_sq=dR_i*dR_i;
    dA_ip1half=grid.dLocalGridOld[grid.nR][nIInt][0][0]*grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dA_im1half=grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      *grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    
    for(j=grid.nStartUpdateExplicit[grid.nQ0][1];j<grid.nEndUpdateExplicit[grid.nQ0][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nQ0][2];k<grid.nEndUpdateExplicit[grid.nQ0][2];k++){
        
        //calculate threshold compression to turn viscosity on at
        dC=sqrt(grid.dLocalGridOld[grid.nGamma][i][j][k]*(grid.dLocalGridOld[grid.nP][i][j][k])
          /grid.dLocalGridOld[grid.nD][i][j][k]);
        dDVDtThreshold=parameters.dAVThreshold*dC;
        
        //calculate Q0
        dDVDt=(dA_ip1half*grid.dLocalGridOld[grid.nU][nIInt][j][k]
          -dA_im1half*grid.dLocalGridOld[grid.nU][nIInt-1][j][k])/dR_i_sq;
        if(dDVDt<-1.0*dDVDtThreshold){//being compressed
          dDVDt_mthreshold=dDVDt+dDVDtThreshold;
          grid.dLocalGridOld[grid.nQ0][i][j][k]=dA_sq*grid.dLocalGridOld[grid.nD][i][j][k]
            *dDVDt_mthreshold*dDVDt_mthreshold;
        }
        else{
          grid.dLocalGridOld[grid.nQ0][i][j][k]=0.0;
        }
      }
    }
  }
  
  //outter ghost region
  for(i=grid.nStartGhostUpdateExplicit[grid.nQ0][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nQ0][0][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_i=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])*0.5;
    dR_i_sq=dR_i*dR_i;
    dA_ip1half=grid.dLocalGridOld[grid.nR][nIInt][0][0]*grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dA_im1half=grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      *grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nQ0][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nQ0][0][1];j++){
      
      for(k=grid.nStartGhostUpdateExplicit[grid.nQ0][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nQ0][0][2];k++){
        
        //calculate threshold compression to turn viscosity on at
        dC=sqrt(grid.dLocalGridOld[grid.nGamma][i][j][k]*(grid.dLocalGridOld[grid.nP][i][j][k])
          /grid.dLocalGridOld[grid.nD][i][j][k]);
        dDVDtThreshold=parameters.dAVThreshold*dC;
        
        //calculate Q0
        dDVDt=(dA_ip1half*grid.dLocalGridOld[grid.nU][nIInt][j][k]
          -dA_im1half*grid.dLocalGridOld[grid.nU][nIInt-1][j][k])/dR_i_sq;
        if(dDVDt<-1.0*dDVDtThreshold){//being compressed
          dDVDt_mthreshold=dDVDt+dDVDtThreshold;
          grid.dLocalGridOld[grid.nQ0][i][j][k]=dA_sq*grid.dLocalGridOld[grid.nD][i][j][k]
            *dDVDt_mthreshold*dDVDt_mthreshold;
        }
        else{
          grid.dLocalGridOld[grid.nQ0][i][j][k]=0.0;
        }
      }
    }
  }
  
  //inner ghost region
  #if SEDOV==1
    for(i=grid.nStartGhostUpdateExplicit[grid.nQ0][1][0];
      i<grid.nEndGhostUpdateExplicit[grid.nQ0][1][0];i++){
      
      //calculate i for interface centered quantities
      nIInt=i+grid.nCenIntOffset[0];
      dR_i=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])*0.5;
      dR_i_sq=dR_i*dR_i;
      dA_ip1half=grid.dLocalGridOld[grid.nR][nIInt][0][0]*grid.dLocalGridOld[grid.nR][nIInt][0][0];
      dA_im1half=grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
        *grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
      
      for(j=grid.nStartGhostUpdateExplicit[grid.nQ0][1][1];
        j<grid.nEndGhostUpdateExplicit[grid.nQ0][1][1];j++){
        
        for(int k=grid.nStartGhostUpdateExplicit[grid.nQ0][1][2];
          k<grid.nEndGhostUpdateExplicit[grid.nQ0][1][2];k++){
          
          //calculate threshold compression to turn viscosity on at
          dC=sqrt(grid.dLocalGridOld[grid.nGamma][i][j][k]*(grid.dLocalGridOld[grid.nP][i][j][k])
            /grid.dLocalGridOld[grid.nD][i][j][k]);
          dDVDtThreshold=parameters.dAVThreshold*dC;
          
          //calculate Q0
          dDVDt=(dA_ip1half*grid.dLocalGridOld[grid.nU][nIInt][j][k]
            -dA_im1half*grid.dLocalGridOld[grid.nU][nIInt-1][j][k])/dR_i_sq;
          if(dDVDt<-1.0*dDVDtThreshold){//being compressed
            dDVDt_mthreshold=dDVDt+dDVDtThreshold;
            grid.dLocalGridOld[grid.nQ0][i][j][k]=dA_sq*grid.dLocalGridOld[grid.nD][i][j][k]
              *dDVDt_mthreshold*dDVDt_mthreshold;
          }
          else{
            grid.dLocalGridOld[grid.nQ0][i][j][k]=0.0;
          }
        }
      }
    }
  #endif
}
