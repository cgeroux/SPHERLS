void calNewDenave_R(Grid &grid){
  int i;
  //most ghost cell, since we don't have R at outer interface. This should be ok in most cases
  for(i=grid.nStartUpdateExplicit[grid.nDenAve][0];i<grid.nEndUpdateExplicit[grid.nDenAve][0];i++){
    grid.dLocalGridNew[grid.nDenAve][i][0][0]=grid.dLocalGridNew[grid.nD][i][0][0];
  }
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nDenAve][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nDenAve][0][0];i++){
    grid.dLocalGridNew[grid.nDenAve][i][0][0]=grid.dLocalGridNew[grid.nD][i][0][0];
  }
}
