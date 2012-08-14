void calOldP_GL(Grid& grid,Parameters &parameters){
  for(int i=grid.nStartUpdateExplicit[grid.nP][0];i<grid.nEndUpdateExplicit[grid.nP][0];i++){
    for(int j=grid.nStartUpdateExplicit[grid.nP][1];j<grid.nEndUpdateExplicit[grid.nP][1];j++){
      for(int k=grid.nStartUpdateExplicit[grid.nP][2];k<grid.nEndUpdateExplicit[grid.nP][2];k++){
        grid.dLocalGridOld[grid.nP][i][j][k]=dEOS_GL(grid.dLocalGridOld[grid.nD][i][j][k]
          ,grid.dLocalGridOld[grid.nE][i][j][k],parameters);
      }
    }
  }
  for(int i=grid.nStartGhostUpdateExplicit[grid.nP][0][0];i<grid.nEndGhostUpdateExplicit[grid.nP][0][0];i++){
    for(int j=grid.nStartGhostUpdateExplicit[grid.nP][0][1];j<grid.nEndGhostUpdateExplicit[grid.nP][0][1];j++){
      for(int k=grid.nStartGhostUpdateExplicit[grid.nP][0][2];k<grid.nEndGhostUpdateExplicit[grid.nP][0][2];k++){
        grid.dLocalGridOld[grid.nP][i][j][k]=dEOS_GL(grid.dLocalGridOld[grid.nD][i][j][k]
          ,grid.dLocalGridOld[grid.nE][i][j][k],parameters);
      }
    }
  }
  for(int i=grid.nStartGhostUpdateExplicit[grid.nP][1][0];i<grid.nEndGhostUpdateExplicit[grid.nP][1][0];i++){
    for(int j=grid.nStartGhostUpdateExplicit[grid.nP][1][1];j<grid.nEndGhostUpdateExplicit[grid.nP][1][1];j++){
      for(int k=grid.nStartGhostUpdateExplicit[grid.nP][1][2];k<grid.nEndGhostUpdateExplicit[grid.nP][1][2];k++){
        grid.dLocalGridOld[grid.nP][i][j][k]=dEOS_GL(grid.dLocalGridOld[grid.nD][i][j][k]
          ,grid.dLocalGridOld[grid.nE][i][j][k],parameters);
      }
    }
  }
  
  //NO IMPLICIT REGION FOR GAMMA LAWA GAS
}
