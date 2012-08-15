double dImplicitEnergyFunction_RTP_SB(Grid &grid,Parameters &parameters,Time &time,double dTemps[]
  ,int i,int j,int k){
  
  double dT_ijk_np1=dTemps[0];
  double dT_im1jk_np1=dTemps[1];
  double dT_ijp1k_np1=dTemps[2];
  double dT_ijm1k_np1=dTemps[3];
  double dT_ijkp1_np1=dTemps[4];
  double dT_ijkm1_np1=dTemps[5];
  
  //calculate i,j,k for interface centered quantities
  int nIInt=i+grid.nCenIntOffset[0];
  int nJInt=j+grid.nCenIntOffset[1];
  int nKInt=k+grid.nCenIntOffset[2];
    
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
  double dV_ijk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
    +grid.dLocalGridNew[grid.nV][i][nJInt-1][k])*0.5;
  double dW_ijk_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt]
    +grid.dLocalGridNew[grid.nW][i][j][nKInt-1])*0.5;
  double dW_ijkp1half_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt]);
  double dW_ijkm1half_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt-1]);
  double dVSinTheta_ijp1halfk_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
    *grid.dLocalGridNew[grid.nV][i][nJInt][k];
  double dVSinTheta_ijm1halfk_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
    *grid.dLocalGridNew[grid.nV][i][nJInt-1][k];
  double dRhoAve_im1half=(grid.dLocalGridOld[grid.nDenAve][i][0][0]
    +grid.dLocalGridOld[grid.nDenAve][i-1][0][0])*0.5;
  double dRho_im1halfjk=(grid.dLocalGridOld[grid.nD][i][j][k]
    +grid.dLocalGridOld[grid.nD][i-1][j][k])*0.5;
  double dRho_ijp1halfk=(grid.dLocalGridOld[grid.nD][i][j+1][k]
    +grid.dLocalGridOld[grid.nD][i][j][k])*0.5;
  double dRho_ijm1halfk=(grid.dLocalGridOld[grid.nD][i][j][k]
    +grid.dLocalGridOld[grid.nD][i][j-1][k])*0.5;
  double dRho_ijkp1half=(grid.dLocalGridOld[grid.nD][i][j][k+1]
    +grid.dLocalGridOld[grid.nD][i][j][k])*0.5;
  double dRho_ijkm1half=(grid.dLocalGridOld[grid.nD][i][j][k]
    +grid.dLocalGridOld[grid.nD][i][j][k-1])*0.5;
  
  double dT_ijk_np1half   =(dT_ijk_np1+grid.dLocalGridOld[grid.nT][i][j][k])*0.5;
  double dT_ijk_np1half_sq= dT_ijk_np1half   *dT_ijk_np1half;
  double dT_ijk_np1half_4 = dT_ijk_np1half_sq*dT_ijk_np1half_sq;
  
  double dT_im1jk_np1half   =(dT_im1jk_np1+grid.dLocalGridOld[grid.nT][i-1][j][k])*0.5;
  double dT_im1jk_np1half_sq= dT_im1jk_np1half   *dT_im1jk_np1half;
  double dT_im1jk_np1half_4 = dT_im1jk_np1half_sq*dT_im1jk_np1half_sq;
  
  double dT_ijp1k_np1half   =(dT_ijp1k_np1+grid.dLocalGridOld[grid.nT][i][j+1][k])*0.5;
  double dT_ijp1k_np1half_sq= dT_ijp1k_np1half   *dT_ijp1k_np1half;
  double dT_ijp1k_np1half_4 = dT_ijp1k_np1half_sq*dT_ijp1k_np1half_sq;
  
  double dT_ijm1k_np1half   =(dT_ijm1k_np1+grid.dLocalGridOld[grid.nT][i][j-1][k])*0.5;
  double dT_ijm1k_np1half_sq= dT_ijm1k_np1half   *dT_ijm1k_np1half;
  double dT_ijm1k_np1half_4 = dT_ijm1k_np1half_sq*dT_ijm1k_np1half_sq;
  
  double dT_ijkp1_np1half   =(dT_ijkp1_np1+grid.dLocalGridOld[grid.nT][i][j][k+1])*0.5;
  double dT_ijkp1_np1half_sq= dT_ijkp1_np1half   *dT_ijkp1_np1half;
  double dT_ijkp1_np1half_4 = dT_ijkp1_np1half_sq*dT_ijkp1_np1half_sq;
  
  double dT_ijkm1_np1half   =(dT_ijkm1_np1+grid.dLocalGridOld[grid.nT][i][j][k-1])*0.5;
  double dT_ijkm1_np1half_sq= dT_ijkm1_np1half   *dT_ijkm1_np1half;
  double dT_ijkm1_np1half_4 = dT_ijkm1_np1half_sq*dT_ijkm1_np1half_sq;
  
  double dE_ijk_np1=parameters.eosTable.dGetEnergy(dT_ijk_np1,grid.dLocalGridNew[grid.nD][i][j][k]);
  double dE_ijk_np1half=parameters.eosTable.dGetEnergy(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  double dE_im1jk_np1half=parameters.eosTable.dGetEnergy(dT_im1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i-1][j][k]);
  double dE_ijp1k_np1half=parameters.eosTable.dGetEnergy(dT_ijp1k_np1half
    ,grid.dLocalGridOld[grid.nD][i][j+1][k]);
  double dE_ijm1k_np1half=parameters.eosTable.dGetEnergy(dT_ijm1k_np1half
    ,grid.dLocalGridOld[grid.nD][i][j-1][k]);
  double dE_ijkp1_np1half=parameters.eosTable.dGetEnergy(dT_ijkp1_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k+1]);
  double dE_ijkm1_np1half=parameters.eosTable.dGetEnergy(dT_ijkm1_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k-1]);
  
  double dE_ip1halfjk_np1half=dE_ijk_np1half;/**\BC Using $E_{i,j,k}^{n+1/2}$ for 
    $E_{i+1/2,j,k}^{n+1/2}$*/
  double dE_im1halfjk_np1half=(dE_im1jk_np1half+dE_ijk_np1half)*0.5;
  double dE_ijp1halfk_np1half=(dE_ijp1k_np1half+dE_ijk_np1half)*0.5;
  double dE_ijm1halfk_np1half=(dE_ijm1k_np1half+dE_ijk_np1half)*0.5;
  double dE_ijkp1half_np1half=(dE_ijkp1_np1half+dE_ijk_np1half)*0.5;
  double dE_ijkm1half_np1half=(dE_ijkm1_np1half+dE_ijk_np1half)*0.5;
  
  double dP_ijk_np1half=parameters.eosTable.dGetPressure(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  #if VISCOUS_ENERGY_EQ==1
    dP_ijk_np1half=dP_ijk_np1half+grid.dLocalGridOld[grid.nQ0][i][j][k]
      +grid.dLocalGridOld[grid.nQ1][i][j][k]+grid.dLocalGridOld[grid.nQ2][i][j][k];
  #endif
  
  double dKappa_ijk_np1half=parameters.eosTable.dGetOpacity(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  double dKappa_im1jk_np1half=parameters.eosTable.dGetOpacity(dT_im1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i-1][j][k]);
  double dKappa_ijp1k_np1half=parameters.eosTable.dGetOpacity(dT_ijp1k_np1half
    ,grid.dLocalGridOld[grid.nD][i][j+1][k]);
  double dKappa_ijm1k_np1half=parameters.eosTable.dGetOpacity(dT_ijm1k_np1half
    ,grid.dLocalGridOld[grid.nD][i][j-1][k]);
  double dKappa_ijkp1_np1half=parameters.eosTable.dGetOpacity(dT_ijkp1_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k+1]);
  double dKappa_ijkm1_np1half=parameters.eosTable.dGetOpacity(dT_ijkm1_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k-1]);
  
  double dKappa_im1halfjk_n_np1half=(dT_im1jk_np1half_4+dT_ijk_np1half_4)/(dT_ijk_np1half_4
    /dKappa_ijk_np1half+dT_im1jk_np1half_4/dKappa_im1jk_np1half);
  double dKappa_ijp1halfk_n_np1half=(dT_ijp1k_np1half_4+dT_ijk_np1half_4)/(dT_ijk_np1half_4
    /dKappa_ijk_np1half+dT_ijp1k_np1half_4/dKappa_ijp1k_np1half);
  double dKappa_ijm1halfk_n_np1half=(dT_ijm1k_np1half_4+dT_ijk_np1half_4)/(dT_ijk_np1half_4
    /dKappa_ijk_np1half+dT_ijm1k_np1half_4/dKappa_ijm1k_np1half);
  double dKappa_ijkp1half_np1half=(dT_ijkp1_np1half_4+dT_ijk_np1half_4)/(dT_ijk_np1half_4
    /dKappa_ijk_np1half+dT_ijkp1_np1half_4/dKappa_ijkp1_np1half);
  double dKappa_ijkm1half_np1half=(dT_ijkm1_np1half_4+dT_ijk_np1half_4)/(dT_ijk_np1half_4
    /dKappa_ijk_np1half+dT_ijkm1_np1half_4/dKappa_ijkm1_np1half);
  
  //Calcuate dA1
  double dA1CenGrad=(dE_ip1halfjk_np1half-dE_im1halfjk_np1half)
    /grid.dLocalGridOld[grid.nDM][i][0][0];
  double dA1UpWindGrad=0.0;
  double dU_U0_Diff=(dU_ijk_np1half-dU0_i_np1half);
  if(dU_U0_Diff<0.0){//moving in the negative direction
    dA1UpWindGrad=dA1UpWindGrad;/**\BC Using centered gradient for upwind gradient when motion is 
      into the star at the surface*/
  }
  else{//moving in the postive direction
    dA1UpWindGrad=(dE_ijk_np1half-dE_im1jk_np1half)/(grid.dLocalGridOld[grid.nDM][i][0][0]
      +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
  }
  double dA1=dU_U0_Diff*dRSq_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA1CenGrad
    +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
  
  //calculate dS1
  double dUR2_im1half_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]
    *dR_im1half_n*dR_im1half_n;
  double dUR2_ip1half_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dR_ip1half_n*dR_ip1half_n;
  double dS1=dP_ijk_np1half/grid.dLocalGridOld[grid.nD][i][j][k]
    *(dUR2_ip1half_np1half-dUR2_im1half_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
  
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
  double dA2=dV_ijk_np1half/dR_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
    *dA2CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
    
  //Calcualte dS2
  double dS2=dP_ijk_np1half/(grid.dLocalGridOld[grid.nD][i][j][k]*dR_i_n
    *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDTheta][0][j][0])
    *(dVSinTheta_ijp1halfk_np1half-dVSinTheta_ijm1halfk_np1half);
  
  //Calcualte dA3
  double dA3CenGrad=(dE_ijkp1half_np1half-dE_ijkm1half_np1half)
    /grid.dLocalGridOld[grid.nDPhi][0][0][k];
  double dA3UpWindGrad=0.0;
  if(dW_ijk_np1half<0.0){//moving in the negative direction
    dA3UpWindGrad=(dE_ijkp1_np1half-dE_ijk_np1half)/(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
      +grid.dLocalGridOld[grid.nDPhi][0][0][k])*2.0;
  }
  else{//moving in the positive direction
    dA3UpWindGrad=(dE_ijk_np1half-dE_ijkm1_np1half)/(grid.dLocalGridOld[grid.nDPhi][0][0][k]
      +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*2.0;
  }
  double dA3=dW_ijk_np1half/(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0])*
    ((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA3CenGrad
    +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA3UpWindGrad);
  
  //Calcualte dS3
  double dS3=dP_ijk_np1half/(grid.dLocalGridOld[grid.nD][i][j][k]*dR_i_n
    *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k])
    *(dW_ijkp1half_np1half-dW_ijkm1half_np1half);
  
  //Calculate dS4
  double dTGrad_im1half_np1half=(dT_ijk_np1half_4-dT_im1jk_np1half_4)
    /(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
  double dGrad_ip1half_np1half=-3.0*dRSq_ip1half*dT_ijk_np1half_4/(8.0*parameters.dPi);/**\BC 
    Missing grid.dLocalGridOld[grid.nT][i+1][0][0] using flux equals \f$2\sigma T^4\f$ at surface.*/
  double dGrad_im1half_np1half=dRhoAve_im1half*dR_im1half_4/(dKappa_im1halfjk_n_np1half
    *dRho_im1halfjk)*dTGrad_im1half_np1half;
  double dS4=16.0*parameters.dPi*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
    *(dGrad_ip1half_np1half-dGrad_im1half_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
  
  //Calculate dS5
  double dTGrad_jp1half_np1half=(dT_ijp1k_np1half_4-dT_ijk_np1half_4)
    /(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]+grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
  double dTGrad_jm1half_np1half=(dT_ijk_np1half_4-dT_ijm1k_np1half_4)
    /(grid.dLocalGridOld[grid.nDTheta][0][j][0]+grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;
  double dGrad_jp1half_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
    /(dKappa_ijp1halfk_n_np1half*dRho_ijp1halfk*dR_i_n)*dTGrad_jp1half_np1half;
  double dGrad_jm1half_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
    /(dKappa_ijm1halfk_n_np1half*dRho_ijm1halfk*dR_i_n)*dTGrad_jm1half_np1half;
  double dS5=(dGrad_jp1half_np1half-dGrad_jm1half_np1half)
    /(grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*dR_i_n
    *grid.dLocalGridOld[grid.nDTheta][0][j][0]);
  
  //Calculate dS6
  double dTGrad_kp1half_np1half=(dT_ijkp1_np1half_4-dT_ijk_np1half_4)
    /(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]+grid.dLocalGridOld[grid.nDPhi][0][0][k])*2.0;
  double dTGrad_km1half_np1half=(dT_ijk_np1half_4-dT_ijkm1_np1half_4)
    /(grid.dLocalGridOld[grid.nDPhi][0][0][k]+grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*2.0;
  double dGrad_kp1half_np1half=dTGrad_kp1half_np1half/(dKappa_ijkp1half_np1half*dRho_ijkp1half
    *dR_i_n);
  double dGrad_km1half_np1half=dTGrad_km1half_np1half/(dKappa_ijkm1half_np1half*dRho_ijkm1half
    *dR_i_n);
  double dS6=(dGrad_kp1half_np1half-dGrad_km1half_np1half)/(dR_i_n
    *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k]);
  
  //calculate new energy
  return (dE_ijk_np1-grid.dLocalGridOld[grid.nE][i][j][k])/time.dDeltat_np1half
    +4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1+dS1)+dA2+dS2+dA3+dS3
    -4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS4+dS5+dS6);
}
