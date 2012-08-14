void calNewP_GL(Grid& grid,Parameters &parameters){
  int i;
  int j;
  int k;
  for(i=grid.nStartUpdateExplicit[grid.nP][0];i<grid.nEndUpdateExplicit[grid.nP][0];i++){
    for(j=grid.nStartUpdateExplicit[grid.nP][1];j<grid.nEndUpdateExplicit[grid.nP][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nP][2];k<grid.nEndUpdateExplicit[grid.nP][2];k++){
        grid.dLocalGridNew[grid.nP][i][j][k]=dEOS_GL(grid.dLocalGridNew[grid.nD][i][j][k]
          ,grid.dLocalGridNew[grid.nE][i][j][k],parameters);
      }
    }
  }
  for(i=grid.nStartGhostUpdateExplicit[grid.nP][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nP][0][0];i++){
    for(j=grid.nStartGhostUpdateExplicit[grid.nP][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nP][0][1];j++){
      for(k=grid.nStartGhostUpdateExplicit[grid.nP][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nP][0][2];k++){
        grid.dLocalGridNew[grid.nP][i][j][k]=dEOS_GL(grid.dLocalGridNew[grid.nD][i][j][k]
          ,grid.dLocalGridNew[grid.nE][i][j][k],parameters);
      }
    }
  }
  #if SEDOV==1 //use zero P, E, and rho gradients
    for(i=grid.nStartGhostUpdateExplicit[grid.nP][1][0];
      i<grid.nEndGhostUpdateExplicit[grid.nP][1][0];i++){
      for(j=grid.nStartGhostUpdateExplicit[grid.nP][1][1];
        j<grid.nEndGhostUpdateExplicit[grid.nP][1][1];j++){
        for(k=grid.nStartGhostUpdateExplicit[grid.nP][1][2];
          k<grid.nEndGhostUpdateExplicit[grid.nP][1][2];k++){
          grid.dLocalGridNew[grid.nP][i][j][k]=dEOS_GL(grid.dLocalGridNew[grid.nD][i][j][k]
            ,grid.dLocalGridNew[grid.nE][i][j][k],parameters);
        }
      }
    }
  #endif
}
