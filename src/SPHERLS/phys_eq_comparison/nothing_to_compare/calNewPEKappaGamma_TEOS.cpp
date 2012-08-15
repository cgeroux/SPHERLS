void calNewPEKappaGamma_TEOS(Grid& grid,Parameters &parameters){
  int i;
  int j;
  int k;
  
  //P, T, Kappa, and Gamma are all cenetered quantities, so bounds of any will be the same
  for(i=grid.nStartUpdateImplicit[grid.nP][0];i<grid.nEndUpdateImplicit[grid.nP][0];i++){
    for(j=grid.nStartUpdateImplicit[grid.nP][1];j<grid.nEndUpdateImplicit[grid.nP][1];j++){
      for(k=grid.nStartUpdateImplicit[grid.nP][2];k<grid.nEndUpdateImplicit[grid.nP][2];k++){
        
        parameters.eosTable.getPEKappaGamma(grid.dLocalGridNew[grid.nT][i][j][k]
          ,grid.dLocalGridNew[grid.nD][i][j][k],grid.dLocalGridNew[grid.nP][i][j][k]
          ,grid.dLocalGridNew[grid.nE][i][j][k],grid.dLocalGridNew[grid.nKappa][i][j][k]
          ,grid.dLocalGridNew[grid.nGamma][i][j][k]);
      }
    }
  }
  for(i=grid.nStartGhostUpdateImplicit[grid.nP][0][0];
    i<grid.nEndGhostUpdateImplicit[grid.nP][0][0];i++){
    for(j=grid.nStartGhostUpdateImplicit[grid.nP][0][1];
      j<grid.nEndGhostUpdateImplicit[grid.nP][0][1];j++){
      for(k=grid.nStartGhostUpdateImplicit[grid.nP][0][2];
        k<grid.nEndGhostUpdateImplicit[grid.nP][0][2];k++){
        
        parameters.eosTable.getPEKappaGamma(grid.dLocalGridNew[grid.nT][i][j][k]
          ,grid.dLocalGridNew[grid.nD][i][j][k],grid.dLocalGridNew[grid.nP][i][j][k]
          ,grid.dLocalGridNew[grid.nE][i][j][k],grid.dLocalGridNew[grid.nKappa][i][j][k]
          ,grid.dLocalGridNew[grid.nGamma][i][j][k]);
      }
    }
  }
}
