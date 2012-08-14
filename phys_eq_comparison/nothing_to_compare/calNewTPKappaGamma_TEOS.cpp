void calNewTPKappaGamma_TEOS(Grid& grid,Parameters &parameters){
  int i;
  int j;
  int k;
  int nCount;
  double dError;
  double dDTDE;
  double dT;
  double dE;
  double dDelE;
  
  //P, T, Kappa, and Gamma are all cenetered quantities, so bounds of any will be the same
  for(i=grid.nStartUpdateExplicit[grid.nP][0];i<grid.nEndUpdateExplicit[grid.nP][0];i++){
    for(j=grid.nStartUpdateExplicit[grid.nP][1];j<grid.nEndUpdateExplicit[grid.nP][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nP][2];k<grid.nEndUpdateExplicit[grid.nP][2];k++){
        
        //calculate new temperature
        dError=std::numeric_limits<double>::max();
        dT=grid.dLocalGridOld[grid.nT][i][j][k];
        nCount=0;
        while(dError>parameters.dTolerance&&nCount<parameters.nMaxIterations){
          parameters.eosTable.getEAndDTDE(dT,grid.dLocalGridNew[grid.nD][i][j][k],dE,dDTDE);
          
          //correct temperature
          dDelE=grid.dLocalGridNew[grid.nE][i][j][k]-dE;
          dT=dDelE*dDTDE+dT;
          
          //how far off was the energy
          dError=fabs(dDelE)/grid.dLocalGridNew[grid.nE][i][j][k];
          nCount++;
        }
        if(nCount>=parameters.nMaxIterations){
          std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": The maximum number of iteration"
          <<" for converging temperature in explicit region from equation of state ("
          <<parameters.nMaxIterations<<") has been exceeded with a maximum relative error in "
          <<"matching the energy of "<<dError<<std::endl;
        }
        grid.dLocalGridNew[grid.nT][i][j][k]=dT;
        
        //get P, Kappa, Gamma
        parameters.eosTable.getPKappaGamma(grid.dLocalGridNew[grid.nT][i][j][k]
          ,grid.dLocalGridNew[grid.nD][i][j][k],grid.dLocalGridNew[grid.nP][i][j][k]
          ,grid.dLocalGridNew[grid.nKappa][i][j][k],grid.dLocalGridNew[grid.nGamma][i][j][k]);
      }
    }
  }
  for(i=grid.nStartGhostUpdateExplicit[grid.nP][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nP][0][0];i++){
    for(j=grid.nStartGhostUpdateExplicit[grid.nP][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nP][0][1];j++){
      for(k=grid.nStartGhostUpdateExplicit[grid.nP][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nP][0][2];k++){
        
        //calculate new temperature
        dError=std::numeric_limits<double>::max();
        dT=grid.dLocalGridOld[grid.nT][i][j][k];
        nCount=0;
        while(dError>parameters.dTolerance&&nCount<parameters.nMaxIterations){
          parameters.eosTable.getEAndDTDE(dT,grid.dLocalGridNew[grid.nD][i][j][k],dE,dDTDE);
          
          //correct temperature
          dDelE=grid.dLocalGridNew[grid.nE][i][j][k]-dE;
          dT=dDelE*dDTDE+dT;
          
          //how far off was the energy
          dError=fabs(dDelE)/grid.dLocalGridNew[grid.nE][i][j][k];
          nCount++;
        }
        if(nCount>=parameters.nMaxIterations){
          std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": The maximum number of iteration"
          <<" for converging temperature in explicit region from equation of state ("
          <<parameters.nMaxIterations<<") has been exceeded with a maximum relative error in "
          <<"matching the energy of "<<dError<<std::endl;
        }
        grid.dLocalGridNew[grid.nT][i][j][k]=dT;
        
        //get P, Kappa, Gamma
        parameters.eosTable.getPKappaGamma(grid.dLocalGridNew[grid.nT][i][j][k]
          ,grid.dLocalGridNew[grid.nD][i][j][k],grid.dLocalGridNew[grid.nP][i][j][k]
          ,grid.dLocalGridNew[grid.nKappa][i][j][k],grid.dLocalGridNew[grid.nGamma][i][j][k]);
      }
    }
  }
}
