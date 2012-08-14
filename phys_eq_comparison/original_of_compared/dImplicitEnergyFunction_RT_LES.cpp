double dImplicitEnergyFunction_RT_LES(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k){
  
  double dT_ijk_np1=dTemps[0];
  double dT_ip1jk_np1=dTemps[1];
  double dT_im1jk_np1=dTemps[2];
  double dT_ijp1k_np1=dTemps[3];
  double dT_ijm1k_np1=dTemps[4];
  
  double dPiSq=parameters.dPi*parameters.dPi;
  
  //calculate i,j,k for interface centered quantities
  int nIInt=i+grid.nCenIntOffset[0];
  int nJInt=j+grid.nCenIntOffset[1];
  
  //Calculate interpolated quantities
  double dDM_ip1half=(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i+1][0][0])
    *0.5;
  double dDM_im1half=(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i-1][0][0])
    *0.5;
  double dDelTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][j][0]
    +grid.dLocalGridOld[grid.nDTheta][0][j+1][0])*0.5;
  double dDelTheta_jm1half=(grid.dLocalGridOld[grid.nDTheta][0][j][0]
    +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*0.5;
  double dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
    +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
  double dR_ip1half_np1half=grid.dLocalGridOld[grid.nR][nIInt][0][0];
  double dR_im1half_np1half=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
  double dR_i_np1half=(dR_ip1half_np1half+dR_im1half_np1half)*0.5;
  double dRSq_i_np1half=dR_i_np1half*dR_i_np1half;
  double dRSq_ip1half_np1half=dR_ip1half_np1half*dR_ip1half_np1half;
  double dR4_ip1half_np1half=dRSq_ip1half_np1half*dRSq_ip1half_np1half;
  double dRSq_im1half_np1half=dR_im1half_np1half*dR_im1half_np1half;
  double dR4_im1half_np1half=dRSq_im1half_np1half*dRSq_im1half_np1half;
  double dRhoAve_ip1half_n=(grid.dLocalGridOld[grid.nDenAve][i+1][0][0]
    +grid.dLocalGridOld[grid.nDenAve][i][0][0])*0.5;
  double dRhoAve_im1half_n=(grid.dLocalGridOld[grid.nDenAve][i][0][0]
    +grid.dLocalGridOld[grid.nDenAve][i-1][0][0])*0.5;
  double dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][i+1][j][k]
    +grid.dLocalGridOld[grid.nD][i][j][k])*0.5;
  double dRho_im1halfjk_n=(grid.dLocalGridOld[grid.nD][i][j][k]
    +grid.dLocalGridOld[grid.nD][i-1][j][k])*0.5;
  double dRho_ijp1halfk_n=(grid.dLocalGridOld[grid.nD][i][j+1][k]
    +grid.dLocalGridOld[grid.nD][i][j][k])*0.5;
  double dRho_ijm1halfk_n=(grid.dLocalGridOld[grid.nD][i][j][k]
    +grid.dLocalGridOld[grid.nD][i][j-1][k])*0.5;
  double dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
    +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
  double dU_ijp1halfk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j+1][k]
    +grid.dLocalGridNew[grid.nU][nIInt-1][j+1][k]+grid.dLocalGridNew[grid.nU][nIInt][j][k]
    +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.25;
  double dU_ijm1halfk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j-1][k]
    +grid.dLocalGridNew[grid.nU][nIInt-1][j-1][k]+grid.dLocalGridNew[grid.nU][nIInt][j][k]
    +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.25;
  double dV_ijk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
    +grid.dLocalGridNew[grid.nV][i][nJInt-1][k])*0.5;
  double dV_ip1halfjk_np1half=(grid.dLocalGridNew[grid.nV][i+1][nJInt][k]
    +grid.dLocalGridNew[grid.nV][i+1][nJInt-1][k]+grid.dLocalGridNew[grid.nV][i][nJInt][k]
    +grid.dLocalGridNew[grid.nV][i][nJInt-1][k])*0.25;
  double dV_im1halfjk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
    +grid.dLocalGridNew[grid.nV][i][nJInt-1][k]+grid.dLocalGridNew[grid.nV][i-1][nJInt][k]
    +grid.dLocalGridNew[grid.nV][i-1][nJInt-1][k])*0.25;
  double dVSinTheta_ijp1halfk_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
    *grid.dLocalGridNew[grid.nV][i][nJInt][k];
  double dVSinTheta_ijm1halfk_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
    *grid.dLocalGridNew[grid.nV][i][nJInt-1][k];
  double dEddyVisc_ip1halfjk_n=(grid.dLocalGridNew[grid.nEddyVisc][i+1][j][k]
    +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
  double dEddyVisc_im1halfjk_n=(grid.dLocalGridNew[grid.nEddyVisc][i-1][j][k]
    +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
  double dEddyVisc_ijp1halfk_n=(grid.dLocalGridNew[grid.nEddyVisc][i][j+1][k]
  +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
  double dEddyVisc_ijm1halfk_n=(grid.dLocalGridNew[grid.nEddyVisc][i][j-1][k]
    +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
  
  double dT_ip1jk_np1half=(dT_ip1jk_np1+grid.dLocalGridOld[grid.nT][i+1][j][k])*0.5;
  double dTSq_ip1jk_np1half=dT_ip1jk_np1half*dT_ip1jk_np1half;
  double dT4_ip1jk_np1half=dTSq_ip1jk_np1half*dTSq_ip1jk_np1half;
  
  double dT_ijk_np1half=(dT_ijk_np1+grid.dLocalGridOld[grid.nT][i][j][k])*0.5;
  double dTSq_ijk_np1half=dT_ijk_np1half*dT_ijk_np1half;
  double dT4_ijk_np1half=dTSq_ijk_np1half*dTSq_ijk_np1half;
  
  double dT_im1jk_np1half=(dT_im1jk_np1+grid.dLocalGridOld[grid.nT][i-1][j][k])*0.5;
  double dTSq_im1jk_np1half=dT_im1jk_np1half*dT_im1jk_np1half;
  double dT4_im1jk_np1half = dTSq_im1jk_np1half*dTSq_im1jk_np1half;
  
  double dT_ijp1k_np1half=(dT_ijp1k_np1+grid.dLocalGridOld[grid.nT][i][j+1][k])*0.5;
  double dTSq_ijp1k_np1half=dT_ijp1k_np1half*dT_ijp1k_np1half;
  double dT4_ijp1k_np1half=dTSq_ijp1k_np1half*dTSq_ijp1k_np1half;
  
  double dT_ijm1k_np1half=(dT_ijm1k_np1+grid.dLocalGridOld[grid.nT][i][j-1][k])*0.5;
  double dTSq_ijm1k_np1half=dT_ijm1k_np1half*dT_ijm1k_np1half;
  double dT4_ijm1k_np1half=dTSq_ijm1k_np1half*dTSq_ijm1k_np1half;
  
  double dE_ijk_np1=parameters.eosTable.dGetEnergy(dT_ijk_np1
    ,grid.dLocalGridNew[grid.nD][i][j][k]);
  double dE_ip1jk_np1half=parameters.eosTable.dGetEnergy(dT_ip1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i+1][j][k]);
  double dE_ijk_np1half=parameters.eosTable.dGetEnergy(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  double dE_im1jk_np1half=parameters.eosTable.dGetEnergy(dT_im1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i-1][j][k]);
  double dE_ijp1k_np1half=parameters.eosTable.dGetEnergy(dT_ijp1k_np1half
    ,grid.dLocalGridOld[grid.nD][i][j+1][k]);
  double dE_ijm1k_np1half=parameters.eosTable.dGetEnergy(dT_ijm1k_np1half
    ,grid.dLocalGridOld[grid.nD][i][j-1][k]);
  
  double dE_ip1halfjk_np1half=(dE_ip1jk_np1half+dE_ijk_np1half)*0.5;
  double dE_im1halfjk_np1half=(dE_im1jk_np1half+dE_ijk_np1half)*0.5;
  double dE_ijp1halfk_np1half=(dE_ijp1k_np1half+dE_ijk_np1half)*0.5;
  double dE_ijm1halfk_np1half=(dE_ijm1k_np1half+dE_ijk_np1half)*0.5;
  
  double dP_ijk_np1half=parameters.eosTable.dGetPressure(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  #if VISCOUS_ENERGY_EQ==1
    dP_ijk_np1half=dP_ijk_np1half+grid.dLocalGridOld[grid.nQ0][i][j][k]
      +grid.dLocalGridOld[grid.nQ1][i][j][k];
  #endif
  
  double dKappa_ip1jk_np1half=parameters.eosTable.dGetOpacity(dT_ip1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i+1][j][k]);
  double dKappa_ijk_np1half=parameters.eosTable.dGetOpacity(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  double dKappa_im1jk_np1half=parameters.eosTable.dGetOpacity(dT_im1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i-1][j][k]);
  double dKappa_ijp1k_np1half=parameters.eosTable.dGetOpacity(dT_ijp1k_np1half
    ,grid.dLocalGridOld[grid.nD][i][j+1][k]);
  double dKappa_ijm1k_np1half=parameters.eosTable.dGetOpacity(dT_ijm1k_np1half
    ,grid.dLocalGridOld[grid.nD][i][j-1][k]);
  
  double dKappa_ip1halfjk_np1half=(dT4_ip1jk_np1half+dT4_ijk_np1half)/(dT4_ijk_np1half
    /dKappa_ijk_np1half+dT4_ip1jk_np1half/dKappa_ip1jk_np1half);
  double dKappa_im1halfjk_np1half=(dT4_im1jk_np1half+dT4_ijk_np1half)/(dT4_ijk_np1half
    /dKappa_ijk_np1half+dT4_im1jk_np1half/dKappa_im1jk_np1half);
  double dKappa_ijp1halfk_np1half=(dT4_ijp1k_np1half+dT4_ijk_np1half)/(dT4_ijk_np1half
    /dKappa_ijk_np1half+dT4_ijp1k_np1half/dKappa_ijp1k_np1half);
  double dKappa_ijm1halfk_np1half=(dT4_ijm1k_np1half+dT4_ijk_np1half)/(dT4_ijk_np1half
    /dKappa_ijk_np1half+dT4_ijm1k_np1half/dKappa_ijm1k_np1half);
  
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
  
  double dDEDM=((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
    *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
  
  //apply DEDM clamp if set, and above the required mass
  if(parameters.bDEDMClamp){
    if(parameters.dDEDMClampMr!=-1.0){//clamp has been set
      if(grid.dLocalGridOld[grid.nM][nIInt][0][0]>=parameters.dDEDMClampMr){
        dDEDM=parameters.dDEDMClampValue;
      }
    }
    else{//if the clamp hasn't been set lets see if we can set it
      
      //this should only be set on the first, static and spherically symetric model
      if(grid.dLocalGridOld[grid.nT][i][0][0]<=parameters.dEDMClampTemperature){
        parameters.dDEDMClampMr=grid.dLocalGridOld[grid.nM][nIInt][0][0];
        parameters.dDEDMClampValue=dDEDM;
      }
    }
  }
  
  double dA1=dU_U0_Diff*dRSq_i_np1half*dDEDM;
  
  //calculate dS1
  double dUR2_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dRSq_im1half_np1half;
  double dUR2_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dRSq_ip1half_np1half;
  double dS1=dP_ijk_np1half/grid.dLocalGridOld[grid.nD][i][j][k]
    *(dUR2_ip1halfjk_np1half-dUR2_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
  
  //Calcualte dA2
  double dA2CenGrad=(dE_ijp1halfk_np1half-dE_ijm1halfk_np1half)
    /grid.dLocalGridOld[grid.nDTheta][0][j][0];
  double dA2UpWindGrad=0.0;
  if(dV_ijk_np1half<0.0){//moving in the negative direction
    dA2UpWindGrad=(dE_ijp1k_np1half-dE_ijk_np1half)/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
      +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
  }
  else{//moving in the positive direction
    dA2UpWindGrad=(dE_ijk_np1half-dE_ijm1k_np1half)/(grid.dLocalGridOld[grid.nDTheta][0][j][0]
      +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;
  }
  double dA2=dV_ijk_np1half/dR_i_np1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
    *dA2CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
    
  //Calcualte dS2
  double dS2=dP_ijk_np1half/(grid.dLocalGridOld[grid.nD][i][j][k]*dR_i_np1half
    *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDTheta][0][j][0])
    *(dVSinTheta_ijp1halfk_np1half-dVSinTheta_ijm1halfk_np1half);
  
  //Calculate dS4
  double dTGrad_ip1half_np1half=(dT4_ip1jk_np1half-dT4_ijk_np1half)
    /(grid.dLocalGridOld[grid.nDM][i+1][0][0]+grid.dLocalGridOld[grid.nDM][i][0][0])*2.0;
  double dTGrad_im1half_np1half=(dT4_ijk_np1half-dT4_im1jk_np1half)
    /(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
  double dGrad_ip1half_np1half=dRhoAve_ip1half_n*dR4_ip1half_np1half/(dKappa_ip1halfjk_np1half
    *dRho_ip1halfjk_n)*dTGrad_ip1half_np1half;
  double dGrad_im1half_np1half=dRhoAve_im1half_n*dR4_im1half_np1half/(dKappa_im1halfjk_np1half
    *dRho_im1halfjk_n)*dTGrad_im1half_np1half;
  double dS4=16.0*parameters.dPi*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
    *(dGrad_ip1half_np1half-dGrad_im1half_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
  
  //Calculate dS5
  double dTGrad_jp1half_np1half=(dT4_ijp1k_np1half-dT4_ijk_np1half)
    /(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]+grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
  double dTGrad_jm1half_np1half=(dT4_ijk_np1half-dT4_ijm1k_np1half)
    /(grid.dLocalGridOld[grid.nDTheta][0][j][0]+grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;
  double dGrad_jp1half_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
    /(dKappa_ijp1halfk_np1half*dRho_ijp1halfk_n)*dTGrad_jp1half_np1half;
  double dGrad_jm1half_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
    /(dKappa_ijm1halfk_np1half*dRho_ijm1halfk_n)*dTGrad_jm1half_np1half;
  double dS5=(dGrad_jp1half_np1half-dGrad_jm1half_np1half)
    /(grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
    *dRSq_i_np1half*grid.dLocalGridOld[grid.nDTheta][0][j][0]);
  
  //calculate dT1
  double dEGrad_ip1halfjk_np1half=dR4_ip1half_np1half*dEddyVisc_ip1halfjk_n*dRhoAve_ip1half_n
    *(dE_ip1jk_np1half-dE_ijk_np1half)/(dRho_ip1halfjk_n*dDM_ip1half);
  double dEGrad_im1halfjk_np1half=dR4_im1half_np1half*dEddyVisc_im1halfjk_n*dRhoAve_im1half_n
    *(dE_ijk_np1half-dE_im1jk_np1half)/(dRho_im1halfjk_n*dDM_im1half);
  double dT1=16.0*dPiSq*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dEGrad_ip1halfjk_np1half
    -dEGrad_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
  
  //calculate dT2
  double dEGrad_ijp1halfk_np1half=dEddyVisc_ijp1halfk_n
    *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
    *(dE_ijp1k_np1half-dE_ijk_np1half)/(dRho_ijp1halfk_n*dR_i_np1half*dDelTheta_jp1half);
  double dEGrad_ijm1halfk_np1half=dEddyVisc_ijm1halfk_n
    *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
    *(dE_ijk_np1half-dE_ijm1k_np1half)/(dRho_ijm1halfk_n*dR_i_np1half*dDelTheta_jm1half);
  double dT2=(dEGrad_ijp1halfk_np1half-dEGrad_ijm1halfk_np1half)/(dR_i_np1half
    *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
    *grid.dLocalGridOld[grid.nDTheta][0][j][0]);
  
  //eddy viscosity terms
  double dEddyViscosityTerms=(dT1+dT2)/parameters.dPrt;
  
  //calculate energy equation discrepancy
  return (dE_ijk_np1-grid.dLocalGridOld[grid.nE][i][j][k])/time.dDeltat_np1half
    +4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1+dS1)+dA2+dS2
    -4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS4+dS5)-dEddyViscosityTerms;
}
