double dImplicitEnergyFunction_R_LES(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k){
  
  /**\todo this funciton should probably be turffed, the LES terms aren't needed in 1D. keeping it 
  for now though.*/
  
  double dT_ijk_np1=dTemps[0];
  double dT_ip1jk_np1=dTemps[1];
  double dT_im1jk_np1=dTemps[2];
  
  //Calculate interpolated quantities
  int nIInt=i+grid.nCenIntOffset[0];
  double dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
    +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
  double dU_im1jk_np1half=(grid.dLocalGridNew[grid.nU][nIInt-1][j][k]
    +grid.dLocalGridNew[grid.nU][nIInt-2][j][k])*0.5;
  double dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
    +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
  double dR_i_np1half=(grid.dLocalGridOld[grid.nR][nIInt][0][0]
    +grid.dLocalGridOld[grid.nR][nIInt-1][0][0]+grid.dLocalGridNew[grid.nR][nIInt][0][0]
    +grid.dLocalGridNew[grid.nR][nIInt-1][0][0])*0.25;
  double dR_im1_np1half=(grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
    +grid.dLocalGridOld[grid.nR][nIInt-2][0][0]+grid.dLocalGridNew[grid.nR][nIInt-1][0][0]
    +grid.dLocalGridNew[grid.nR][nIInt-2][0][0])*0.25;
  double dRSq_im1_np1half=dR_im1_np1half*dR_im1_np1half;
  double dR_im1half_np1half=(grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
    +grid.dLocalGridNew[grid.nR][nIInt-1][0][0])*0.5;
  double dR_ip1half_np1half=(grid.dLocalGridOld[grid.nR][nIInt][0][0]
    +grid.dLocalGridNew[grid.nR][nIInt][0][0])*0.5;
  double dRSq_i_np1half=dR_i_np1half*dR_i_np1half;
  double dRSq_ip1half_np1half=dR_ip1half_np1half*dR_ip1half_np1half;
  double dR4_ip1half_np1half=dRSq_ip1half_np1half*dRSq_ip1half_np1half;
  double dRSq_im1half_np1half=dR_im1half_np1half*dR_im1half_np1half;
  double dR4_im1half_np1half=dRSq_im1half_np1half*dRSq_im1half_np1half;
  double dRhoAve_ip1half_np1half=(grid.dLocalGridOld[grid.nD][i+1][0][0]
    +grid.dLocalGridOld[grid.nD][i][0][0])*0.5;
  double dRhoAve_im1half_np1half=(grid.dLocalGridOld[grid.nD][i][0][0]
    +grid.dLocalGridOld[grid.nD][i-1][0][0])*0.5;
  double dRho_ip1halfjk_np1half=(grid.dLocalGridOld[grid.nD][i+1][j][k]
    +grid.dLocalGridOld[grid.nD][i][j][k])*0.5;
  double dRho_im1halfjk_np1half=(grid.dLocalGridOld[grid.nD][i][j][k]
    +grid.dLocalGridOld[grid.nD][i-1][j][k])*0.5;
  double dDM_ip1half=(grid.dLocalGridOld[grid.nDM][i+1][0][0]+grid.dLocalGridOld[grid.nDM][i][0][0])
    *0.5;
  double dDM_im1half=(grid.dLocalGridOld[grid.nDM][i-1][0][0]+grid.dLocalGridOld[grid.nDM][i][0][0])
    *0.5;
        
  double dEddyVisc_ip1half_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j][k]
    +grid.dLocalGridNew[grid.nEddyVisc][i+1][j][k])*0.5;
  double dEddyVisc_im1half_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j][k]
    +grid.dLocalGridNew[grid.nEddyVisc][i-1][j][k])*0.5;
  
  double dT_ip1jk_np1half=(dT_ip1jk_np1+grid.dLocalGridOld[grid.nT][i+1][j][k])*0.5;
  double dT_ip1jk_np1half_sq=dT_ip1jk_np1half*dT_ip1jk_np1half;
  double dT_ip1jk_np1half_4=dT_ip1jk_np1half_sq*dT_ip1jk_np1half_sq;
  
  double dT_ijk_np1half=(dT_ijk_np1+grid.dLocalGridOld[grid.nT][i][j][k])*0.5;
  double dT_ijk_np1half_sq=dT_ijk_np1half*dT_ijk_np1half;
  double dT_ijk_np1half_4=dT_ijk_np1half_sq*dT_ijk_np1half_sq;
  
  double dT_im1jk_np1half=(dT_im1jk_np1+grid.dLocalGridOld[grid.nT][i-1][j][k])*0.5;
  double dT_im1jk_np1half_sq=dT_im1jk_np1half*dT_im1jk_np1half;
  double dT_im1jk_np1half_4=dT_im1jk_np1half_sq*dT_im1jk_np1half_sq;
  
  double dE_ijk_np1=parameters.eosTable.dGetEnergy(dT_ijk_np1
    ,grid.dLocalGridNew[grid.nD][i][j][k]);
  double dE_ip1jk_np1half=parameters.eosTable.dGetEnergy(dT_ip1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i+1][j][k]);
  double dE_ijk_np1half=parameters.eosTable.dGetEnergy(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  double dE_im1jk_np1half=parameters.eosTable.dGetEnergy(dT_im1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i-1][j][k]);
  
  double dE_ip1halfjk_np1half=(dE_ip1jk_np1half+dE_ijk_np1half)*0.5;
  double dE_im1halfjk_np1half=(dE_ijk_np1half+dE_im1jk_np1half)*0.5;
  
  double dP_ijk_np1half=parameters.eosTable.dGetPressure(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  #if VISCOUS_ENERGY_EQ==1
    dP_ijk_np1half=dP_ijk_np1half+grid.dLocalGridOld[grid.nQ0][i][j][k];
  #endif
  
  double dKappa_ip1jk_np1half=parameters.eosTable.dGetOpacity(dT_ip1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i+1][j][k]);
  double dKappa_ijk_np1half=parameters.eosTable.dGetOpacity(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  double dKappa_im1jk_np1half=parameters.eosTable.dGetOpacity(dT_im1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i-1][j][k]);
  
  double dKappa_ip1halfjk_n_np1half=(dT_ip1jk_np1half_4+dT_ijk_np1half_4)
    /(dT_ijk_np1half_4/dKappa_ijk_np1half+dT_ip1jk_np1half_4/dKappa_ip1jk_np1half);
  double dKappa_im1halfjk_n_np1half=(dT_im1jk_np1half_4+dT_ijk_np1half_4)
    /(dT_ijk_np1half_4/dKappa_ijk_np1half+dT_im1jk_np1half_4/dKappa_im1jk_np1half);
    
  //Calcuate dA1
  double dA1CenGrad=(dE_ip1halfjk_np1half-dE_im1halfjk_np1half)
    /grid.dLocalGridOld[grid.nDM][i][0][0];
  double dA1UpWindGrad=0.0;
  double dU_U0_Diff=(dU_ijk_np1half-dU0_i_np1half);
  if(dU_U0_Diff<0.0){//moving in the negative direction
    dA1UpWindGrad=(dE_ip1jk_np1half-dE_ijk_np1half)/(grid.dLocalGridOld[grid.nDM][i+1][0][0]
      +grid.dLocalGridOld[grid.nDM][i][0][0])*2.0;
  }
  else{//moving in the postive direction
    dA1UpWindGrad=(dE_ijk_np1half-dE_im1jk_np1half)/(grid.dLocalGridOld[grid.nDM][i][0][0]
      +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
  }
  double dA1=dU_U0_Diff*dRSq_i_np1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
    *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
  
  //calculate dUR2_im1half_np1half
  double dUR2_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dRSq_im1half_np1half;
  
  //calculate dR2_ip1half_np1half
  double dUR2_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dRSq_ip1half_np1half;
  
  //calculate dURsq_ijk_np1half
  double dUR2_ijk_np1half=dU_ijk_np1half*dRSq_i_np1half;
  
  //calculate dURsq_im1jk_np1half
  double dUR2_im1jk_np1half=dU_im1jk_np1half*dRSq_im1_np1half;
  
  //calculate dS1
  double dS1=dP_ijk_np1half/grid.dLocalGridOld[grid.nD][i][j][k]
    *(dUR2_ip1halfjk_np1half-dUR2_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
  
  //Calculate dS4
  double dTGrad_ip1half_np1half=(dT_ip1jk_np1half_4-dT_ijk_np1half_4)
    /(grid.dLocalGridOld[grid.nDM][i+1][0][0]+grid.dLocalGridOld[grid.nDM][i][0][0])*2.0;
  double dTGrad_im1half_np1half=(dT_ijk_np1half_4-dT_im1jk_np1half_4)
    /(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
  double dGrad_ip1half_np1half=dRhoAve_ip1half_np1half*dR4_ip1half_np1half
    /(dKappa_ip1halfjk_n_np1half*dRho_ip1halfjk_np1half)*dTGrad_ip1half_np1half;
  double dGrad_im1half_np1half=dRhoAve_im1half_np1half*dR4_im1half_np1half
    /(dKappa_im1halfjk_n_np1half*dRho_im1halfjk_np1half)*dTGrad_im1half_np1half;
  double dS4=16.0*parameters.dPi*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]
    *(dGrad_ip1half_np1half-dGrad_im1half_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
       
  //calculate dDivU_ip1halfjk_np1half
  double dDivU_ip1halfjk_np1half=4.0*parameters.dPi*dRhoAve_ip1half_np1half*(dUR2_ip1halfjk_np1half
    -dUR2_ijk_np1half)/dDM_ip1half;
  
  //calculate dDivU_im1halfjk_np1half
  double dDivU_im1halfjk_np1half=4.0*parameters.dPi*dRhoAve_im1half_np1half*(dUR2_ijk_np1half
    -dUR2_im1jk_np1half)/dDM_im1half;
  
  //calculate dTau_rr_ip1halfjk_np1half
  double dTau_rr_ip1halfjk_np1half=2.0*dEddyVisc_ip1half_np1half*(4.0*parameters.dPi
    *dRSq_ip1half_np1half*dRhoAve_ip1half_np1half*(grid.dLocalGridNew[grid.nU][nIInt][j][k]
    -dU_ijk_np1half)/dDM_ip1half-0.333333333333333*dDivU_ip1halfjk_np1half);
  
  //calculate dTau_rr_im1halfjk_np1half
  double dTau_rr_im1halfjk_np1half=2.0*dEddyVisc_im1half_np1half*(4.0*parameters.dPi
    *dRSq_im1half_np1half*dRhoAve_im1half_np1half*(dU_ijk_np1half-dU_im1jk_np1half)
    /dDM_im1half-0.333333333333333*dDivU_im1halfjk_np1half);
  
  //calculate dTauVR2_ip1halfjk_np1half
  double dTauVR2_ip1halfjk_np1half=dRSq_ip1half_np1half*dTau_rr_ip1halfjk_np1half
    *grid.dLocalGridNew[grid.nU][nIInt][j][k];
  
  //calculate dTauVR2_im1halfjk_np1half
  double dTauVR2_im1halfjk_np1half=dRSq_im1half_np1half*dTau_rr_im1halfjk_np1half
    *grid.dLocalGridNew[grid.nU][nIInt-1][j][k];
  
  //calculate dT1
  double dT1=(dTauVR2_ip1halfjk_np1half-dTauVR2_im1halfjk_np1half)
    /(grid.dLocalGridOld[grid.nDM][i][0][0]*grid.dLocalGridOld[grid.nD][i][j][k]);
  
  //calculate PDV
  double dPDV=-1.0*dP_ijk_np1half*(1.0/grid.dLocalGridNew[grid.nD][i][j][k]
    -1.0/grid.dLocalGridOld[grid.nD][i][j][k]);
  
  #if DEBUG_EQUATIONS==1
  if(parameters.bSetThisCall){
    
    //add E
    parameters.profileDataDebug.setMaxAbs("E"
      ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
      ,grid.dLocalGridOld[grid.nE][i][j][k]);
    
    //add A1
    parameters.profileDataDebug.setMaxAbs("E_RTP_A1"
      ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
      ,-4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1));
    
    //add PDV
    parameters.profileDataDebug.setMaxAbs("E_RTP_PDV"
      ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
      ,dPDV);
    
    //add S4
    parameters.profileDataDebug.setMaxAbs("E_RTP_S4"
      ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
      ,4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS4));
  }
  #endif
  
  //calculate new energy
  return (dE_ijk_np1-grid.dLocalGridOld[grid.nE][i][j][k]-dPDV)/time.dDeltat_np1half
    +4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]*(dA1-dT1)
    -4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS4);
}
