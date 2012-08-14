void calOldDenave_RTP(Grid &grid){
  
  //EXPLICIT REGION
  //most ghost cell, since we don't have R at outer interface. This should be ok in most cases
  for(int i=grid.nStartUpdateExplicit[grid.nDenAve][0];i<grid.nEndUpdateExplicit[grid.nDenAve][0];
    i++){
    
    //calculate i for interface centered quantities
    int nIInt=i+grid.nCenIntOffset[0];
    
    double dSum=0.0;
    double dVolume=0.0;
    double dRFactor=0.33333333333333333*(pow(grid.dLocalGridOld[grid.nR][nIInt][0][0],3.0)
      -pow(grid.dLocalGridOld[grid.nR][nIInt-1][0][0],3.0));
    for(int j=grid.nStartUpdateExplicit[grid.nD][1];j<grid.nEndUpdateExplicit[grid.nD][1];j++){
      for(int k=grid.nStartUpdateExplicit[grid.nD][2];k<grid.nEndUpdateExplicit[grid.nD][2];k++){
        double dVolumeTemp=dRFactor*grid.dLocalGridOld[grid.nDCosThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dSum+=dVolumeTemp*grid.dLocalGridOld[grid.nD][i][j][k];
        dVolume+=dVolumeTemp;
      }
    }
    grid.dLocalGridOld[grid.nDenAve][i][0][0]=dSum/dVolume;
  }
  //ghost region 0, outter most ghost region in x1 direction
  for(int i=grid.nStartGhostUpdateExplicit[grid.nDenAve][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nDenAve][0][0];i++){
    double dSum=0.0;
    double dVolume=0.0;
    double dRFactor=0.33333333333333333*(pow(grid.dLocalGridOld[grid.nR][i][0][0],3.0)
      -pow(grid.dLocalGridOld[grid.nR-1][i][0][0],3.0));
    for(int j=grid.nStartUpdateExplicit[grid.nD][1];j<grid.nEndUpdateExplicit[grid.nD][1];j++){
      for(int k=grid.nStartUpdateExplicit[grid.nD][2];k<grid.nEndUpdateExplicit[grid.nD][2];k++){
        double dVolumeTemp=dRFactor*grid.dLocalGridOld[grid.nDCosThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dSum+=dVolumeTemp*grid.dLocalGridOld[grid.nD][i][j][k];
        dVolume+=dVolumeTemp;
      }
    }
    grid.dLocalGridOld[grid.nDenAve][i][0][0]=dSum/dVolume;
  }
  
  //IMPLICT REGION
  //most ghost cell, since we don't have R at outer interface. This should be ok in most cases
  for(int i=grid.nStartUpdateImplicit[grid.nDenAve][0];i<grid.nEndUpdateImplicit[grid.nDenAve][0];
    i++){
    
    //calculate i for interface centered quantities
    int nIInt=i+grid.nCenIntOffset[0];
    
    double dSum=0.0;
    double dVolume=0.0;
    double dRFactor=0.33333333333333333*(pow(grid.dLocalGridOld[grid.nR][nIInt][0][0],3.0)
      -pow(grid.dLocalGridOld[grid.nR][nIInt-1][0][0],3.0));
    for(int j=grid.nStartUpdateImplicit[grid.nD][1];j<grid.nEndUpdateImplicit[grid.nD][1];j++){
      for(int k=grid.nStartUpdateImplicit[grid.nD][2];k<grid.nEndUpdateImplicit[grid.nD][2];k++){
        double dVolumeTemp=dRFactor*grid.dLocalGridOld[grid.nDCosThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dSum+=dVolumeTemp*grid.dLocalGridOld[grid.nD][i][j][k];
        dVolume+=dVolumeTemp;
      }
    }
    grid.dLocalGridOld[grid.nDenAve][i][0][0]=dSum/dVolume;
  }
  //ghost region 0, outter most ghost region in x1 direction
  for(int i=grid.nStartGhostUpdateImplicit[grid.nDenAve][0][0];
    i<grid.nEndGhostUpdateImplicit[grid.nDenAve][0][0];i++){
    double dSum=0.0;
    double dVolume=0.0;
    double dRFactor=0.33333333333333333*(pow(grid.dLocalGridOld[grid.nR][i][0][0],3.0)
      -pow(grid.dLocalGridOld[grid.nR-1][i][0][0],3.0));
    for(int j=grid.nStartUpdateExplicit[grid.nD][1];j<grid.nEndUpdateExplicit[grid.nD][1];j++){
      for(int k=grid.nStartUpdateExplicit[grid.nD][2];k<grid.nEndUpdateExplicit[grid.nD][2];k++){
        double dVolumeTemp=dRFactor*grid.dLocalGridOld[grid.nDCosThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dSum+=dVolumeTemp*grid.dLocalGridOld[grid.nD][i][j][k];
        dVolume+=dVolumeTemp;
      }
    }
    grid.dLocalGridOld[grid.nDenAve][i][0][0]=dSum/dVolume;
  }

}
