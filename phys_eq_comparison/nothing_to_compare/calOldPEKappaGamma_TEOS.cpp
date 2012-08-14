void calOldPEKappaGamma_TEOS(Grid& grid,Parameters &parameters){
  
  //main grid
  for(int i=grid.nStartUpdateExplicit[grid.nP][0];i<grid.nEndUpdateExplicit[grid.nP][0];i++){
    for(int j=grid.nStartUpdateExplicit[grid.nP][1];j<grid.nEndUpdateExplicit[grid.nP][1];j++){
      for(int k=grid.nStartUpdateExplicit[grid.nP][2];k<grid.nEndUpdateExplicit[grid.nP][2];k++){
        parameters.eosTable.getPEKappaGamma(grid.dLocalGridOld[grid.nT][i][j][k]
          ,grid.dLocalGridOld[grid.nD][i][j][k],grid.dLocalGridOld[grid.nP][i][j][k]
          ,grid.dLocalGridOld[grid.nE][i][j][k],grid.dLocalGridOld[grid.nKappa][i][j][k]
          ,grid.dLocalGridOld[grid.nGamma][i][j][k]);
      }
    }
  }
  
  //inner radial ghost cells
  for(int i=0;i<grid.nNumGhostCells;i++){
    for(int j=grid.nStartUpdateExplicit[grid.nP][1];j<grid.nEndUpdateExplicit[grid.nP][1];j++){
      for(int k=grid.nStartUpdateExplicit[grid.nP][2];k<grid.nEndUpdateExplicit[grid.nP][2];k++){
        parameters.eosTable.getPEKappaGamma(grid.dLocalGridOld[grid.nT][i][j][k]
          ,grid.dLocalGridOld[grid.nD][i][j][k],grid.dLocalGridOld[grid.nP][i][j][k]
          ,grid.dLocalGridOld[grid.nE][i][j][k],grid.dLocalGridOld[grid.nKappa][i][j][k]
          ,grid.dLocalGridOld[grid.nGamma][i][j][k]);
      }
    }
  }
  
  //outter radial ghost cells
  for(int i=grid.nStartGhostUpdateExplicit[grid.nP][0][0];i<grid.nEndGhostUpdateExplicit[grid.nP][0][0];i++){
    for(int j=grid.nStartGhostUpdateExplicit[grid.nP][0][1];j<grid.nEndGhostUpdateExplicit[grid.nP][0][1];j++){
      for(int k=grid.nStartGhostUpdateExplicit[grid.nP][0][2];k<grid.nEndGhostUpdateExplicit[grid.nP][0][2];k++){
        parameters.eosTable.getPEKappaGamma(grid.dLocalGridOld[grid.nT][i][j][k]
          ,grid.dLocalGridOld[grid.nD][i][j][k],grid.dLocalGridOld[grid.nP][i][j][k]
          ,grid.dLocalGridOld[grid.nE][i][j][k],grid.dLocalGridOld[grid.nKappa][i][j][k]
          ,grid.dLocalGridOld[grid.nGamma][i][j][k]);
      }
    }
  }
  
  //main grid
  for(int i=grid.nStartUpdateImplicit[grid.nP][0];i<grid.nEndUpdateImplicit[grid.nP][0];i++){
    for(int j=grid.nStartUpdateImplicit[grid.nP][1];j<grid.nEndUpdateImplicit[grid.nP][1];j++){
      for(int k=grid.nStartUpdateImplicit[grid.nP][2];k<grid.nEndUpdateImplicit[grid.nP][2];k++){
        parameters.eosTable.getPEKappaGamma(grid.dLocalGridOld[grid.nT][i][j][k]
          ,grid.dLocalGridOld[grid.nD][i][j][k],grid.dLocalGridOld[grid.nP][i][j][k]
          ,grid.dLocalGridOld[grid.nE][i][j][k],grid.dLocalGridOld[grid.nKappa][i][j][k]
          ,grid.dLocalGridOld[grid.nGamma][i][j][k]);
      }
    }
  }
  
  //inner radial ghost cells
  for(int i=0;i<grid.nNumGhostCells;i++){
    for(int j=grid.nStartUpdateImplicit[grid.nP][1];j<grid.nEndUpdateImplicit[grid.nP][1];j++){
      for(int k=grid.nStartUpdateImplicit[grid.nP][2];k<grid.nEndUpdateImplicit[grid.nP][2];k++){
        parameters.eosTable.getPEKappaGamma(grid.dLocalGridOld[grid.nT][i][j][k]
          ,grid.dLocalGridOld[grid.nD][i][j][k],grid.dLocalGridOld[grid.nP][i][j][k]
          ,grid.dLocalGridOld[grid.nE][i][j][k],grid.dLocalGridOld[grid.nKappa][i][j][k]
          ,grid.dLocalGridOld[grid.nGamma][i][j][k]);
      }
    }
  }
  
  //outter radial ghost cells
  for(int i=grid.nStartGhostUpdateImplicit[grid.nP][0][0];i<grid.nEndGhostUpdateImplicit[grid.nP][0][0];i++){
    for(int j=grid.nStartGhostUpdateImplicit[grid.nP][0][1];j<grid.nEndGhostUpdateImplicit[grid.nP][0][1];j++){
      for(int k=grid.nStartGhostUpdateImplicit[grid.nP][0][2];k<grid.nEndGhostUpdateImplicit[grid.nP][0][2];k++){
        parameters.eosTable.getPEKappaGamma(grid.dLocalGridOld[grid.nT][i][j][k]
          ,grid.dLocalGridOld[grid.nD][i][j][k],grid.dLocalGridOld[grid.nP][i][j][k]
          ,grid.dLocalGridOld[grid.nE][i][j][k],grid.dLocalGridOld[grid.nKappa][i][j][k]
          ,grid.dLocalGridOld[grid.nGamma][i][j][k]);
      }
    }
  }
}
