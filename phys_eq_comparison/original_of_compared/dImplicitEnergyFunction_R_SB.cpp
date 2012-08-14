double dImplicitEnergyFunction_R_SB(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k){
  
  double dT_ijk_np1=dTemps[0];
  double dT_im1jk_np1=dTemps[1];
  int nIInt=i+grid.nCenIntOffset[0];
  
  //Calculate interpolated quantities
  double dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
    +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
  double dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
    +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
  double dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]
    +grid.dLocalGridOld[grid.nR][nIInt-1][0][0])*0.5;
  double dR_im1half_n=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
  double dR_ip1half_n=grid.dLocalGridOld[grid.nR][nIInt][0][0];
  double dRSq_i_n=dR_i_n*dR_i_n;
  double dRSq_ip1half=dR_ip1half_n*dR_ip1half_n;
  double dR_im1half_sq=dR_im1half_n*dR_im1half_n;
  double dR_im1half_4=dR_im1half_sq*dR_im1half_sq;
  double dRhoAve_im1half=(grid.dLocalGridOld[grid.nD][i][0][0]
    +grid.dLocalGridOld[grid.nD][i-1][0][0])*0.5;
  double dRho_im1halfjk=(grid.dLocalGridOld[grid.nD][i][j][k]
    +grid.dLocalGridOld[grid.nD][i-1][j][k])*0.5;
  
  double dT_ijk_np1half=(dT_ijk_np1+grid.dLocalGridOld[grid.nT][i][j][k])*0.5;
  double dT_ijk_np1half_sq=dT_ijk_np1half*dT_ijk_np1half;
  double dT_ijk_np1half_4=dT_ijk_np1half_sq*dT_ijk_np1half_sq;
  
  double dT_im1jk_np1half=(dT_im1jk_np1+grid.dLocalGridOld[grid.nT][i-1][j][k])*0.5;
  double dT_im1jk_np1half_sq=dT_im1jk_np1half*dT_im1jk_np1half;
  double dT_im1jk_np1half_4=dT_im1jk_np1half_sq*dT_im1jk_np1half_sq;
  
  double dE_ijk_np1=parameters.eosTable.dGetEnergy(dT_ijk_np1
    ,grid.dLocalGridNew[grid.nD][i][j][k]);
  double dE_ijk_np1half=parameters.eosTable.dGetEnergy(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  double dE_im1jk_np1half=parameters.eosTable.dGetEnergy(dT_im1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  
  double dE_ip1halfjk_np1half=dE_ijk_np1half;/**\BC Missing
    grid.dLocalGridOld[grid.nE][i+1][j][k] in calculation of \f$E_{i+1/2,j,k}\f$ 
    setting it equal to value at i.*/
  double dE_im1halfjk_np1half=(dE_ijk_np1half+dE_im1jk_np1half)*0.5;
  
  double dKappa_ijk_np1half=parameters.eosTable.dGetOpacity(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  double dKappa_im1jk_np1half=parameters.eosTable.dGetOpacity(dT_im1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i-1][j][k]);
  
  double dKappa_im1halfjk_n_np1half=(dT_im1jk_np1half_4+dT_ijk_np1half_4)
    /(dT_ijk_np1half_4/dKappa_ijk_np1half+dT_im1jk_np1half_4/dKappa_im1jk_np1half);
    
  double dP_ijk_np1half=parameters.eosTable.dGetPressure(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  #if VISCOUS_ENERGY_EQ==1
    dP_ijk_np1half=dP_ijk_np1half+grid.dLocalGridOld[grid.nQ0][i][j][k];
  #endif
  
  //Calcuate dA1
  double dA1CenGrad=(dE_ip1halfjk_np1half-dE_im1halfjk_np1half)
    /grid.dLocalGridOld[grid.nDM][i][0][0];
  double dA1UpWindGrad=0.0;
  double dU_U0_Diff=(dU_ijk_np1half-dU0_i_np1half);
  if(dU_U0_Diff<0.0){//moving in the negative direction
    dA1UpWindGrad=0.0;/**\BC A1 upwind set to zero as no material is flowing into the star*/
  }
  else{//moving in the postive direction*/
    dA1UpWindGrad=(dE_ijk_np1half-dE_im1jk_np1half)/(grid.dLocalGridOld[grid.nDM][i][0][0]
      +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
  }
  double dA1=dU_U0_Diff*dRSq_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA1CenGrad
    +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
  
  //Calculate dS4
  double dTGrad_im1half_np1half=(dT_ijk_np1half_4-dT_im1jk_np1half_4)
    /(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
  double dGrad_ip1half_np1half=-3.0*dRSq_ip1half*dT_ijk_np1half_4/(8.0*parameters.dPi);/**\BC
    Missing grid.dLocalGridOld[grid.nT][i+1][0][0]*/
  double dGrad_im1half_np1half=dRhoAve_im1half*dR_im1half_4/(dKappa_im1halfjk_n_np1half
    *dRho_im1halfjk)*dTGrad_im1half_np1half;
  double dS4=16.0*parameters.dPi*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]
    *(dGrad_ip1half_np1half-dGrad_im1half_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
  
  //calculate dS1
  double dUR2_im1half_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dR_im1half_n
    *dR_im1half_n;
  double dUR2_ip1half_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dR_ip1half_n
    *dR_ip1half_n;
  double dS1=dP_ijk_np1half/grid.dLocalGridOld[grid.nD][i][j][k]
    *(dUR2_ip1half_np1half-dUR2_im1half_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
  
  #if DEBUG_EQUATIONS==1
  if(parameters.bSetThisCall){
    
    //add E
    parameters.profileDataDebug.setMaxAbs("E"
      ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
      ,grid.dLocalGridOld[grid.nE][i][j][k]);
    
    //add A1
    parameters.profileDataDebug.setMaxAbs("E_A1"
      ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
      ,-4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1));
    
    //add S1
    parameters.profileDataDebug.setMaxAbs("E_S1"
      ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
      ,-4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dS1));
    
    //add S4
    parameters.profileDataDebug.setMaxAbs("E_S4"
      ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
      ,4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS4));
  }
  #endif

  //calculate new energy
  return (dE_ijk_np1-grid.dLocalGridOld[grid.nE][i][j][k])/time.dDeltat_np1half
    +4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]*(dA1+dS1)
    -4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*dS4;
  
}
