void calNewR(Grid &grid, Time &time){
  int i;
  int l;
  for(i=grid.nStartUpdateExplicit[grid.nR][0];i<grid.nEndUpdateExplicit[grid.nR][0];i++){//nR needs to be 1D
    grid.dLocalGridNew[grid.nR][i][0][0]=grid.dLocalGridOld[grid.nR][i][0][0]
      +time.dDeltat_np1half*grid.dLocalGridNew[grid.nU0][i][0][0];
  }
  for(l=0;l<6;l++){
    for(i=grid.nStartGhostUpdateExplicit[grid.nR][l][0];i<grid.nEndGhostUpdateExplicit[grid.nR][l][0];i++){
      grid.dLocalGridNew[grid.nR][i][0][0]=grid.dLocalGridOld[grid.nR][i][0][0]
        +time.dDeltat_np1half*grid.dLocalGridNew[grid.nU0][i][0][0];
    }
  }
}
