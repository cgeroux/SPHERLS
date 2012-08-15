void calNewDenave_RTP(Grid &grid){
  int i;
  int j;
  int k;
  int nIInt;
  double dSum;
  double dVolume;
  double dRFactor;
  double dVolumeTemp;
  for(i=grid.nStartUpdateExplicit[grid.nDenAve][0];i<grid.nEndUpdateExplicit[grid.nDenAve][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    
    dSum=0.0;
    dVolume=0.0;
    dRFactor=0.33333333333333333*(pow(grid.dLocalGridNew[grid.nR][nIInt][0][0],3.0)
      -pow(grid.dLocalGridNew[grid.nR][nIInt-1][0][0],3.0));
    for(j=grid.nStartUpdateExplicit[grid.nD][1];j<grid.nEndUpdateExplicit[grid.nD][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nD][2];k<grid.nEndUpdateExplicit[grid.nD][2];k++){
        dVolumeTemp=dRFactor*grid.dLocalGridOld[grid.nDCosThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dSum+=dVolumeTemp*grid.dLocalGridNew[grid.nD][i][j][k];
        dVolume+=dVolumeTemp;
      }
    }
    grid.dLocalGridNew[grid.nDenAve][i][0][0]=dSum/dVolume;
  }
  
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nDenAve][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nDenAve][0][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    
    dSum=0.0;
    dVolume=0.0;
    dRFactor=0.33333333333333333*(pow(grid.dLocalGridNew[grid.nR][nIInt][0][0],3.0)
      -pow(grid.dLocalGridNew[grid.nR][nIInt-1][0][0],3.0));
    for(j=grid.nStartUpdateExplicit[grid.nD][1];j<grid.nEndUpdateExplicit[grid.nD][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nD][2];k<grid.nEndUpdateExplicit[grid.nD][2];k++){
        dVolumeTemp=dRFactor*grid.dLocalGridOld[grid.nDCosThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dSum+=dVolumeTemp*grid.dLocalGridNew[grid.nD][i][j][k];
        dVolume+=dVolumeTemp;
      }
    }
    grid.dLocalGridNew[grid.nDenAve][i][0][0]=dSum/dVolume;
  }
}
