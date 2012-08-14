void calOldDenave_R(Grid &grid){//resume, moving decalerations
  
  //explicit region
  for(int i=grid.nStartUpdateExplicit[grid.nDenAve][0];
    i<grid.nEndUpdateExplicit[grid.nDenAve][0];i++){
    grid.dLocalGridOld[grid.nDenAve][i][0][0]=grid.dLocalGridOld[grid.nD][i][0][0];
  }
  for(int i=grid.nStartGhostUpdateExplicit[grid.nDenAve][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nDenAve][0][0];i++){
    grid.dLocalGridOld[grid.nDenAve][i][0][0]=grid.dLocalGridOld[grid.nD][i][0][0];
  }
  
  //implicit region
  for(int i=grid.nStartUpdateImplicit[grid.nDenAve][0];i<grid.nEndUpdateImplicit[grid.nDenAve][0];
    i++){
    grid.dLocalGridOld[grid.nDenAve][i][0][0]=grid.dLocalGridOld[grid.nD][i][0][0];
  }
  for(int i=grid.nStartGhostUpdateImplicit[grid.nDenAve][0][0];
    i<grid.nEndGhostUpdateImplicit[grid.nDenAve][0][0];i++){
    grid.dLocalGridOld[grid.nDenAve][i][0][0]=grid.dLocalGridOld[grid.nD][i][0][0];
  }
}
