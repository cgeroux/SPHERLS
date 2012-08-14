void calNewE_RTP_NA_LES(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop){
  int i;
  int j;
  int k;
  int nIInt;
  int nJInt;
  int nKInt;
  double dDM_ip1half;
  double dDM_im1half;
  double dDelTheta_jp1half;
  double dDelTheta_jm1half;
  double dDelPhi_kp1half;
  double dDelPhi_km1half;
  double dU_ijk_np1half;
  double dU_ijp1halfk_np1half;
  double dU_ijm1halfk_np1half;
  double dU_ijkp1half_np1half;
  double dU_ijkm1half_np1half;
  double dV_ijk_np1half;
  double dV_ip1halfjk_np1half;
  double dV_im1halfjk_np1half;
  double dV_ijkp1half_np1half;
  double dV_ijkm1half_np1half;
  double dW_ijk_np1half;
  double dW_ijkp1half_np1half;
  double dW_ijkm1half_np1half;
  double dW_ip1halfjk_np1half;
  double dW_im1halfjk_np1half;
  double dW_ijp1halfk_np1half;
  double dW_ijm1halfk_np1half;
  double dU0_i_np1half;
  double dE_ip1halfjk_n;
  double dE_im1halfjk_n;
  double dR_ip1_n;
  double dRSq_ip1_n;
  double dR_i_n;
  double dR_im1half_n;
  double dR_ip1half_n;
  double dRSq_i_n;
  double dR_im1_n;
  double dRSq_im1_n;
  double dRSq_ip1half_n;
  double dR4_ip1half_n;
  double dRSq_im1half_n;
  double dR4_im1half_n;
  double dE_ijp1halfk_n;
  double dE_ijm1halfk_n;
  double dUmU0_ijk_np1half;
  double dVSinTheta_ijp1halfk_np1half;
  double dVSinTheta_ijm1halfk_np1half;
  double dUR2_im1halfjk_np1half;
  double dUR2_ip1halfjk_np1half;
  double dE_ijkp1half_n;
  double dE_ijkm1half_n;
  double dRhoAve_ip1half_n;
  double dRhoAve_im1half_n;
  double dRho_ip1halfjk_n;
  double dRho_im1halfjk_n;
  double dRho_ijp1halfk_n;
  double dRho_ijm1halfk_n;
  double dRho_ijkp1half_n;
  double dRho_ijkm1half_n;
  double dTSq_ip1jk_n;
  double dT4_ip1jk_n;
  double dTSq_ijk_n;
  double dT4_ijk_n;
  double dTSq_im1jk_n;
  double dT4_im1jk_n;
  double dTSq_ijp1k_n;
  double dT4_ijp1k_n;
  double dTSq_ijm1k_n;
  double dT4_ijm1k_n;
  double dTSq_ijkp1_n;
  double dT4_ijkp1_n;
  double dTSq_ijkm1_n;
  double dT4_ijkm1_n;
  double dKappa_ip1halfjk_n;
  double dKappa_im1halfjk_n;
  double dKappa_ijp1halfk_n;
  double dKappa_ijm1halfk_n;
  double dKappa_ijkp1half_n;
  double dKappa_ijkm1half_n;
  double dA1CenGrad;
  double dA1UpWindGrad;
  double dU_U0_Diff;
  double dA1;
  double dP_ijk_n;
  double dA2CenGrad;
  double dA2UpWindGrad;
  double dA2;
  double dA3CenGrad;
  double dA3UpWindGrad;
  double dA3;
  double dTGrad_ip1half;
  double dTGrad_im1half;
  double dGrad_ip1half;
  double dGrad_im1half;
  double dS4;
  double dTGrad_jp1half;
  double dTGrad_jm1half;
  double dGrad_jp1half;
  double dGrad_jm1half;
  double dS5;
  double dTGrad_kp1half;
  double dTGrad_km1half;
  double dGrad_kp1half;
  double dGrad_km1half;
  double dS6;
  double dT1;
  double dT2;
  double dT3;
  double dEGrad_ip1halfjk_np1half;
  double dEGrad_im1halfjk_np1half;
  double dEGrad_ijp1halfk_np1half;
  double dEGrad_ijm1halfk_np1half;
  double dEGrad_ijkp1half_np1half;
  double dEGrad_ijkm1half_np1half;
  double dEddyVisc_ip1halfjk_np1half;
  double dEddyVisc_im1halfjk_np1half;
  double dEddyVisc_ijp1halfk_np1half;
  double dEddyVisc_ijm1halfk_np1half;
  double dEddyVisc_ijkm1half_np1half;
  double dEddyVisc_ijkp1half_np1half;
  double VSinTheta_ijp1halfk_np1half;
  double dEddyViscosityTerms;
  double dPiSq=parameters.dPi*parameters.dPi;
  double dS1;
  double dS2;
  double dS3;
  
  for(i=grid.nStartUpdateExplicit[grid.nE][0];i<grid.nEndUpdateExplicit[grid.nE][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_n=grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dR_im1half_n=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_ip1_n=(grid.dLocalGridOld[grid.nR][nIInt+1][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt][0][0])*0.5;
    dRSq_ip1_n=dR_ip1_n*dR_ip1_n;
    dR_im1_n=(grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt-2][0][0])*0.5;
    dRSq_im1_n=dR_im1_n*dR_im1_n;
    dR_i_n=(dR_ip1half_n+dR_im1half_n)*0.5;
    dRSq_i_n=dR_i_n*dR_i_n;
    dRSq_ip1half_n=dR_ip1half_n*dR_ip1half_n;
    dR4_ip1half_n=dRSq_ip1half_n*dRSq_ip1half_n;
    dRSq_im1half_n=dR_im1half_n*dR_im1half_n;
    dR4_im1half_n=dRSq_im1half_n*dRSq_im1half_n;
    dRhoAve_ip1half_n=(grid.dLocalGridOld[grid.nDenAve][i][0][0]
      +grid.dLocalGridOld[grid.nDenAve][i+1][0][0])*0.5;
    dRhoAve_im1half_n=(grid.dLocalGridOld[grid.nDenAve][i][0][0]
      +grid.dLocalGridOld[grid.nDenAve][i-1][0][0])*0.5;
    dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
      +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
    dDM_ip1half=(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i+1][0][0])*0.5;
    dDM_im1half=(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i-1][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nE][1];j<grid.nEndUpdateExplicit[grid.nE][1];j++){
      
      //calculate j for interface centered quantities
      nJInt=j+grid.nCenIntOffset[1];
      dDelTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][j][0]
        +grid.dLocalGridOld[grid.nDTheta][0][j+1][0])*0.5;
      dDelTheta_jm1half=(grid.dLocalGridOld[grid.nDTheta][0][j][0]
        +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*0.5;
      
      for(k=grid.nStartUpdateExplicit[grid.nE][2];k<grid.nEndUpdateExplicit[grid.nE][2];k++){
        
        //calculate k for interface centered quantities
        nKInt=k+grid.nCenIntOffset[2];
        dDelPhi_kp1half=(grid.dLocalGridOld[grid.nDPhi][0][0][k]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k+1])*0.5;
        dDelPhi_km1half=(grid.dLocalGridOld[grid.nDPhi][0][0][k]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*0.5;
        
        //Calculate interpolated quantities
        dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dU_ijp1halfk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j+1][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j+1][k]+grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.25;
        dU_ijm1halfk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j-1][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j-1][k]+grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.25;
        dU_ijkp1half_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt][j][k+1]+grid.dLocalGridNew[grid.nU][nIInt-1][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k+1])*0.25;
        dU_ijkm1half_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt][j][k-1]+grid.dLocalGridNew[grid.nU][nIInt-1][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k-1])*0.25;
        dV_ijk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k])*0.5;
        dV_ip1halfjk_np1half=(grid.dLocalGridNew[grid.nV][i+1][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i+1][nJInt-1][k]+grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k])*0.25;
        dV_im1halfjk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k]+grid.dLocalGridNew[grid.nV][i-1][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i-1][nJInt-1][k])*0.25;
        dV_ijkp1half_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k+1]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k+1]+grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k])*0.25;
        dV_ijkm1half_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k]+grid.dLocalGridNew[grid.nV][i][nJInt][k-1]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k-1])*0.25;
        dW_ijk_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j][nKInt-1])*0.5;
        dW_ijkp1half_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt]);
        dW_ijkm1half_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt-1]);
        dW_ip1halfjk_np1half=(grid.dLocalGridNew[grid.nW][i+1][j][nKInt]
          +grid.dLocalGridNew[grid.nW][i+1][j][nKInt-1]+grid.dLocalGridNew[grid.nW][i][j][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j][nKInt-1])*0.25;
        dW_im1halfjk_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j][nKInt-1]+grid.dLocalGridNew[grid.nW][i-1][j][nKInt]
          +grid.dLocalGridNew[grid.nW][i-1][j][nKInt-1])*0.25;
        dW_ijp1halfk_np1half=(grid.dLocalGridNew[grid.nW][i][j+1][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j+1][nKInt-1]+grid.dLocalGridNew[grid.nW][i][j][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j][nKInt-1])*0.25;
        dW_ijm1halfk_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j][nKInt-1]+grid.dLocalGridNew[grid.nW][i][j-1][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j-1][nKInt-1])*0.25;
        dE_ip1halfjk_n=(grid.dLocalGridOld[grid.nE][i+1][j][k]
          +grid.dLocalGridOld[grid.nE][i][j][k])*0.5;
        dE_im1halfjk_n=(grid.dLocalGridOld[grid.nE][i][j][k]
          +grid.dLocalGridOld[grid.nE][i-1][j][k])*0.5;
        dE_ijp1halfk_n=(grid.dLocalGridOld[grid.nE][i][j+1][k]
          +grid.dLocalGridOld[grid.nE][i][j][k])*0.5;
        dE_ijm1halfk_n=(grid.dLocalGridOld[grid.nE][i][j][k]
          +grid.dLocalGridOld[grid.nE][i][j-1][k])*0.5;
        dE_ijkp1half_n=(grid.dLocalGridOld[grid.nE][i][j][k+1]+grid.dLocalGridOld[grid.nE][i][j][k])
          *0.5;
        dE_ijkm1half_n=(grid.dLocalGridOld[grid.nE][i][j][k-1]+grid.dLocalGridOld[grid.nE][i][j][k])
          *0.5;
        dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][i+1][j][k]
          +grid.dLocalGridOld[grid.nD][i][j][k])*0.5;
        dRho_im1halfjk_n=(grid.dLocalGridOld[grid.nD][i][j][k]
          +grid.dLocalGridOld[grid.nD][i-1][j][k])*0.5;
        dRho_ijp1halfk_n=(grid.dLocalGridOld[grid.nD][i][j+1][k]
          +grid.dLocalGridOld[grid.nD][i][j][k])*0.5;
        dRho_ijm1halfk_n=(grid.dLocalGridOld[grid.nD][i][j][k]
          +grid.dLocalGridOld[grid.nD][i][j-1][k])*0.5;
        dRho_ijkp1half_n=(grid.dLocalGridOld[grid.nD][i][j][k+1]
          +grid.dLocalGridOld[grid.nD][i][j][k])*0.5;
        dRho_ijkm1half_n=(grid.dLocalGridOld[grid.nD][i][j][k]
          +grid.dLocalGridOld[grid.nD][i][j][k-1])*0.5;
        dEddyVisc_ip1halfjk_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i+1][j][k]
          +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
        dEddyVisc_im1halfjk_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i-1][j][k]
          +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
        dEddyVisc_ijp1halfk_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j+1][k]
        +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
        dEddyVisc_ijm1halfk_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j-1][k]
          +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
        dEddyVisc_ijkp1half_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j][k+1]
          +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
        dEddyVisc_ijkm1half_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j][k-1]
          +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
        
        //calculate derived quantities
        dVSinTheta_ijp1halfk_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          *grid.dLocalGridNew[grid.nV][i][nJInt][k];
        dVSinTheta_ijm1halfk_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
          *grid.dLocalGridNew[grid.nV][i][nJInt-1][k];
        dUR2_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dRSq_im1half_n;
        dUR2_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dRSq_ip1half_n;
        dTSq_ip1jk_n=grid.dLocalGridOld[grid.nT][i+1][j][k]*grid.dLocalGridOld[grid.nT][i+1][j][k];
        dT4_ip1jk_n=dTSq_ip1jk_n*dTSq_ip1jk_n;
        dTSq_ijk_n=grid.dLocalGridOld[grid.nT][i][j][k]*grid.dLocalGridOld[grid.nT][i][j][k];
        dT4_ijk_n=dTSq_ijk_n*dTSq_ijk_n;
        dTSq_im1jk_n=grid.dLocalGridOld[grid.nT][i-1][j][k]*grid.dLocalGridOld[grid.nT][i-1][j][k];
        dT4_im1jk_n=dTSq_im1jk_n*dTSq_im1jk_n;
        dTSq_ijp1k_n=grid.dLocalGridOld[grid.nT][i][j+1][k]*grid.dLocalGridOld[grid.nT][i][j+1][k];
        dT4_ijp1k_n=dTSq_ijp1k_n*dTSq_ijp1k_n;
        dTSq_ijm1k_n=grid.dLocalGridOld[grid.nT][i][j-1][k]*grid.dLocalGridOld[grid.nT][i][j-1][k];
        dT4_ijm1k_n=dTSq_ijm1k_n*dTSq_ijm1k_n;
        dTSq_ijkp1_n=grid.dLocalGridOld[grid.nT][i][j][k+1]*grid.dLocalGridOld[grid.nT][i][j][k+1];
        dT4_ijkp1_n=dTSq_ijkp1_n*dTSq_ijkp1_n;
        dTSq_ijkm1_n=grid.dLocalGridOld[grid.nT][i][j][k-1]*grid.dLocalGridOld[grid.nT][i][j][k-1];
        dT4_ijkm1_n=dTSq_ijkm1_n*dTSq_ijkm1_n;
        dKappa_ip1halfjk_n=(dT4_ip1jk_n+dT4_ijk_n)/(dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]+dT4_ip1jk_n
          /grid.dLocalGridOld[grid.nKappa][i+1][j][k]);
        dKappa_im1halfjk_n=(dT4_im1jk_n+dT4_ijk_n)/(dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]+dT4_im1jk_n
          /grid.dLocalGridOld[grid.nKappa][i-1][j][k]);
        dKappa_ijp1halfk_n=(dT4_ijp1k_n+dT4_ijk_n)/(dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]+dT4_ijp1k_n
          /grid.dLocalGridOld[grid.nKappa][i][j+1][k]);
        dKappa_ijm1halfk_n=(dT4_ijm1k_n+dT4_ijk_n)/(dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]+dT4_ijm1k_n
          /grid.dLocalGridOld[grid.nKappa][i][j-1][k]);
        dKappa_ijkp1half_n=(dT4_ijkp1_n+dT4_ijk_n)/(dT4_ijkp1_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k+1]+dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]);
        dKappa_ijkm1half_n=(dT4_ijkm1_n+dT4_ijk_n)/(dT4_ijkm1_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k-1]+dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]);
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP_ijk_n=dP_ijk_n+grid.dLocalGridOld[grid.nQ0][i][j][k]
            +grid.dLocalGridOld[grid.nQ1][i][j][k]+grid.dLocalGridOld[grid.nQ2][i][j][k];
        #endif
        
        //Calcuate dA1
        dA1CenGrad=(dE_ip1halfjk_n-dE_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0];
        dA1UpWindGrad=0.0;
        dUmU0_ijk_np1half=(dU_ijk_np1half-dU0_i_np1half);
        if(dUmU0_ijk_np1half<0.0){//moving in the negative direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nE][i+1][j][k]
            -grid.dLocalGridOld[grid.nE][i][j][k])/(grid.dLocalGridOld[grid.nDM][i+1][0][0]
            +grid.dLocalGridOld[grid.nDM][i][0][0])*2.0;
        }
        else{//moving in the postive direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k]
            -grid.dLocalGridOld[grid.nE][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=dUmU0_ijk_np1half*dRSq_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate dS1
        dS1=dP_ijk_n/grid.dLocalGridOld[grid.nD][i][j][k]
          *(dUR2_ip1halfjk_np1half-dUR2_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //Calcualte dA2
        dA2CenGrad=(dE_ijp1halfk_n-dE_ijm1halfk_n)/grid.dLocalGridOld[grid.nDTheta][0][j][0];
        if(dV_ijk_np1half<0.0){//moving in the negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j+1][k]
            -grid.dLocalGridOld[grid.nE][i][j][k])/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        }
        else{//moving in the positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k]
            -grid.dLocalGridOld[grid.nE][i][j-1][k])/(grid.dLocalGridOld[grid.nDTheta][0][j][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;
        }
        dA2=dV_ijk_np1half/dR_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA2CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
          
        //Calcualte dS2
        dS2=dP_ijk_n/(grid.dLocalGridOld[grid.nD][i][j][k]*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDTheta][0][j][0])
          *(dVSinTheta_ijp1halfk_np1half-dVSinTheta_ijm1halfk_np1half);
        
        //Calcualte dA3
        dA3CenGrad=(dE_ijkp1half_n-dE_ijkm1half_n)/grid.dLocalGridOld[grid.nDPhi][0][0][k];
        if(dW_ijk_np1half<0.0){//moving in the negative direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k+1]
            -grid.dLocalGridOld[grid.nE][i][j][k])/
            (grid.dLocalGridOld[grid.nDPhi][0][0][k+1]+grid.dLocalGridOld[grid.nDPhi][0][0][k])*2.0;
        }
        else{//moving in the positive direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k]
            -grid.dLocalGridOld[grid.nE][i][j][k-1])/
            (grid.dLocalGridOld[grid.nDPhi][0][0][k]+grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*2.0;
        }
        dA3=dW_ijk_np1half/(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0])*
          ((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA3CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA3UpWindGrad);
        
        //Calcualte dS3
        dS3=dP_ijk_n/(grid.dLocalGridOld[grid.nD][i][j][k]*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k])
          *(dW_ijkp1half_np1half-dW_ijkm1half_np1half);
        
        //Calculate dS4
        dTGrad_ip1half=(dT4_ip1jk_n-dT4_ijk_n)/(grid.dLocalGridOld[grid.nDM][i+1][0][0]
          +grid.dLocalGridOld[grid.nDM][i][0][0])*2.0;
        dTGrad_im1half=(dT4_ijk_n-dT4_im1jk_n)/(grid.dLocalGridOld[grid.nDM][i][0][0]
          +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        dGrad_ip1half=dRhoAve_ip1half_n*dR4_ip1half_n/(dKappa_ip1halfjk_n
          *dRho_ip1halfjk_n)*dTGrad_ip1half;
        dGrad_im1half=dRhoAve_im1half_n*dR4_im1half_n/(dKappa_im1halfjk_n
          *dRho_im1halfjk_n)*dTGrad_im1half;
        dS4=16.0*dPiSq*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dGrad_ip1half-dGrad_im1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //Calculate dS5
        dTGrad_jp1half=(dT4_ijp1k_n-dT4_ijk_n)/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
          +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        dTGrad_jm1half=(dT4_ijk_n-dT4_ijm1k_n)/(grid.dLocalGridOld[grid.nDTheta][0][j][0]
          +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;;
        dGrad_jp1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          /(dKappa_ijp1halfk_n*dRho_ijp1halfk_n)*dTGrad_jp1half;
        dGrad_jm1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
          /(dKappa_ijm1halfk_n*dRho_ijm1halfk_n)*dTGrad_jm1half;
        dS5=(dGrad_jp1half-dGrad_jm1half)/(grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *dRSq_i_n*grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //Calculate dS6
        dTGrad_kp1half=(dT4_ijkp1_n-dT4_ijk_n)/(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k])*2.0;
        dTGrad_km1half=(dT4_ijk_n-dT4_ijkm1_n)/(grid.dLocalGridOld[grid.nDPhi][0][0][k]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*2.0;;
        dGrad_kp1half=dTGrad_kp1half/(dKappa_ijkp1half_n*dRho_ijkp1half_n);
        dGrad_km1half=dTGrad_km1half/(dKappa_ijkm1half_n*dRho_ijkm1half_n);
        dS6=(dGrad_kp1half-dGrad_km1half)/(dRSq_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //calculate dT1
        dEGrad_ip1halfjk_np1half=dR4_ip1half_n*dEddyVisc_ip1halfjk_np1half*dRhoAve_ip1half_n
          *(grid.dLocalGridOld[grid.nE][i+1][j][k]-grid.dLocalGridOld[grid.nE][i][j][k])
          /(dRho_ip1halfjk_n*dDM_ip1half);
        dEGrad_im1halfjk_np1half=dR4_im1half_n*dEddyVisc_im1halfjk_np1half*dRhoAve_im1half_n
          *(grid.dLocalGridOld[grid.nE][i][j][k]-grid.dLocalGridOld[grid.nE][i-1][j][k])
          /(dRho_im1halfjk_n*dDM_im1half);
        dT1=16.0*dPiSq*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dEGrad_ip1halfjk_np1half
          -dEGrad_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //calculate dT2
        dEGrad_ijp1halfk_np1half=dEddyVisc_ijp1halfk_np1half
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          *(grid.dLocalGridOld[grid.nE][i][j+1][k]-grid.dLocalGridOld[grid.nE][i][j][k])
          /(dRho_ijp1halfk_n*dR_i_n*dDelTheta_jp1half);
        dEGrad_ijm1halfk_np1half=dEddyVisc_ijm1halfk_np1half
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
          *(grid.dLocalGridOld[grid.nE][i][j][k]-grid.dLocalGridOld[grid.nE][i][j-1][k])
          /(dRho_ijm1halfk_n*dR_i_n*dDelTheta_jm1half);
        dT2=(dEGrad_ijp1halfk_np1half-dEGrad_ijm1halfk_np1half)/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //calculate dT3
        dEGrad_ijkp1half_np1half=dEddyVisc_ijkp1half_np1half*(grid.dLocalGridOld[grid.nE][i][j][k+1]
          -grid.dLocalGridOld[grid.nE][i][j][k])/(dRho_ijkp1half_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*dR_i_n*dDelPhi_kp1half);
        dEGrad_ijkm1half_np1half=dEddyVisc_ijkm1half_np1half*(grid.dLocalGridOld[grid.nE][i][j][k]
          -grid.dLocalGridOld[grid.nE][i][j][k-1])/(dRho_ijkm1half_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*dR_i_n*dDelPhi_km1half);
        dT3=(dEGrad_ijkp1half_np1half-dEGrad_ijkm1half_np1half)/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //eddy viscosity terms
        dEddyViscosityTerms=(dT1+dT2+dT3)/parameters.dPrt;
        
        //calculate new energy
        grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridOld[grid.nE][i][j][k]
          -time.dDeltat_np1half*(4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dA1+dS1)+dA2+dS2+dA3+dS3-4.0*parameters.dSigma/(3.0
          *grid.dLocalGridOld[grid.nD][i][j][k])*(dS4+dS5+dS6)-dEddyViscosityTerms);
        
        #if DEBUG_EQUATIONS==1
        
        
        //if we don't want zone by zone, set ssEnd.str("")
        std::stringstream ssName;
        std::stringstream ssEnd;
        if(parameters.bEveryJK){
          ssEnd<<"_"<<j<<"_"<<k;
        }
        else{
          ssEnd.str("");
        }
        
        //add E
        ssName.str("");
        ssName<<"E"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,grid.dLocalGridOld[grid.nE][i][j][k]);
        
        //add A1
        ssName.str("");
        ssName<<"E_A1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1));
        
        //add A2
        ssName.str("");
        ssName<<"E_A2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dA2);
        
        //add A3
        ssName.str("");
        ssName<<"E_A3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dA3);
        
        //add S1
        ssName.str("");
        ssName<<"E_S1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dS1));
          
        //add S2
        ssName.str("");
        ssName<<"E_S2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS2);
          
        //add S3
        ssName.str("");
        ssName<<"E_S3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS3);
        
        //add S4
        ssName.str("");
        ssName<<"E_S4"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS4));
        
        //add S5
        ssName.str("");
        ssName<<"E_S5"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS5));
        
        //add S6
        ssName.str("");
        ssName<<"E_S6"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS6));
        
        //add EV
        ssName.str("");
        ssName<<"E_EV"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,dEddyViscosityTerms);
        
        //add DEDt
        ssName.str("");
        ssName<<"E_DEDt"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,(grid.dLocalGridNew[grid.nE][i][j][k]-grid.dLocalGridOld[grid.nE][i][j][k])
          /time.dDeltat_np1half);
        #endif
        
        if(grid.dLocalGridNew[grid.nE][i][j][k]<0.0){
          
          #if SIGNEGENG==1
            raise(SIGINT);
          #endif
          
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
            <<": negative energy calculated in , ("<<i<<","<<j<<","<<k<<")\n";
          throw exception2(ssTemp.str(),CALCULATION);
          
        }
      }
    }
  }
  
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nE][0][0];i<grid.nEndGhostUpdateExplicit[grid.nE][0][0]
    ;i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_n=grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dR_im1half_n=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_im1_n=(grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt-2][0][0])*0.5;
    dRSq_im1_n=dR_im1_n*dR_im1_n;
    dR_i_n=(dR_ip1half_n+dR_im1half_n)*0.5;
    dRSq_i_n=dR_i_n*dR_i_n;
    dRSq_ip1half_n=dR_ip1half_n*dR_ip1half_n;
    dR4_ip1half_n=dRSq_ip1half_n*dRSq_ip1half_n;
    dRSq_im1half_n=dR_im1half_n*dR_im1half_n;
    dR4_im1half_n=dRSq_im1half_n*dRSq_im1half_n;
    dRhoAve_ip1half_n=(grid.dLocalGridOld[grid.nDenAve][i][0][0])*0.5;/*\BC missing average density
      outside model setting it to zero*/
    dRhoAve_im1half_n=(grid.dLocalGridOld[grid.nDenAve][i][0][0]
      +grid.dLocalGridOld[grid.nDenAve][i-1][0][0])*0.5;
    dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
      +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
    dDM_ip1half=(grid.dLocalGridOld[grid.nDM][i][0][0])*(0.5+parameters.dAlpha
      +parameters.dAlphaExtra);/**\BC Missing \f$\Delta M_r\f$ outside model using 
      \ref Parameters.dAlpha times \f$\Delta M_r\f$ in the last zone instead.*/
    dDM_im1half=(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i-1][0][0])*0.5;
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nE][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nE][0][1];j++){
      
      //calculate i for interface centered quantities
      nJInt=j+grid.nCenIntOffset[1];
      dDelTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][j][0]
        +grid.dLocalGridOld[grid.nDTheta][0][j+1][0])*0.5;
      dDelTheta_jm1half=(grid.dLocalGridOld[grid.nDTheta][0][j][0]
        +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*0.5;
      
      for(k=grid.nStartGhostUpdateExplicit[grid.nE][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nE][0][2];k++){
        
        nKInt=k+grid.nCenIntOffset[2];
        dDelPhi_kp1half=(grid.dLocalGridOld[grid.nDPhi][0][0][k]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k+1])*0.5;
        dDelPhi_km1half=(grid.dLocalGridOld[grid.nDPhi][0][0][k]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*0.5;
        
        //Calculate interpolated quantities
        dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dU_ijp1halfk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j+1][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j+1][k]+grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.25;
        dU_ijm1halfk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j-1][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j-1][k]+grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.25;
        dU_ijkp1half_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt][j][k+1]+grid.dLocalGridNew[grid.nU][nIInt-1][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k+1])*0.25;
        dU_ijkm1half_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt][j][k-1]+grid.dLocalGridNew[grid.nU][nIInt-1][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k-1])*0.25;
        dV_ijk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k])*0.5;
        dV_ip1halfjk_np1half=dV_ijk_np1half;
        dV_im1halfjk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k]+grid.dLocalGridNew[grid.nV][i-1][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i-1][nJInt-1][k])*0.25;
        dV_ijkp1half_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k+1]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k+1]+grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k])*0.25;
        dV_ijkm1half_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k]+grid.dLocalGridNew[grid.nV][i][nJInt][k-1]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k-1])*0.25;
        dW_ijk_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j][nKInt-1])*0.5;
        dW_ijkp1half_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt]);
        dW_ijkm1half_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt-1]);
        dW_ip1halfjk_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j][nKInt-1])*0.5;/**\BC Missing W at i+1, assuming the 
          same as at i*/
        dW_im1halfjk_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j][nKInt-1]+grid.dLocalGridNew[grid.nW][i-1][j][nKInt]
          +grid.dLocalGridNew[grid.nW][i-1][j][nKInt-1])*0.25;
        dW_ijp1halfk_np1half=(grid.dLocalGridNew[grid.nW][i][j+1][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j+1][nKInt-1]+grid.dLocalGridNew[grid.nW][i][j][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j][nKInt-1])*0.25;
        dW_ijm1halfk_np1half=(grid.dLocalGridNew[grid.nW][i][j][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j][nKInt-1]+grid.dLocalGridNew[grid.nW][i][j-1][nKInt]
          +grid.dLocalGridNew[grid.nW][i][j-1][nKInt-1])*0.25;
        dE_ip1halfjk_n=(grid.dLocalGridOld[grid.nE][i][j][k]);/**\BC Missing
          grid.dLocalGridOld[grid.nE][i+1][j][k] in calculation of \f$E_{i+1/2,j,k}\f$ 
          setting it equal to the value at i.*/
        dE_im1halfjk_n=(grid.dLocalGridOld[grid.nE][i][j][k]+grid.dLocalGridOld[grid.nE][i-1][j][k])
          *0.5;
        dE_ijp1halfk_n=(grid.dLocalGridOld[grid.nE][i][j+1][k]+grid.dLocalGridOld[grid.nE][i][j][k])
          *0.5;
        dE_ijm1halfk_n=(grid.dLocalGridOld[grid.nE][i][j][k]+grid.dLocalGridOld[grid.nE][i][j-1][k])
          *0.5;
        dE_ijkp1half_n=(grid.dLocalGridOld[grid.nE][i][j][k+1]+grid.dLocalGridOld[grid.nE][i][j][k])
          *0.5;
        dE_ijkm1half_n=(grid.dLocalGridOld[grid.nE][i][j][k-1]+grid.dLocalGridOld[grid.nE][i][j][k])
          *0.5;
        dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][i+1][j][k])*0.5;/**\BC missing density outside
          model, setting it to zero*/
        dRho_im1halfjk_n=(grid.dLocalGridOld[grid.nD][i][j][k]
          +grid.dLocalGridOld[grid.nD][i-1][j][k])*0.5;
        dRho_ijp1halfk_n=(grid.dLocalGridOld[grid.nD][i][j+1][k]
          +grid.dLocalGridOld[grid.nD][i][j][k])*0.5;
        dRho_ijm1halfk_n=(grid.dLocalGridOld[grid.nD][i][j][k]
          +grid.dLocalGridOld[grid.nD][i][j-1][k])*0.5;
        dRho_ijkp1half_n=(grid.dLocalGridOld[grid.nD][i][j][k+1]
          +grid.dLocalGridOld[grid.nD][i][j][k])*0.5;
        dRho_ijkm1half_n=(grid.dLocalGridOld[grid.nD][i][j][k]
          +grid.dLocalGridOld[grid.nD][i][j][k-1])*0.5;
        dEddyVisc_ip1halfjk_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;/**\BC missing 
          eddy viscosity outside the model setting it to zero*/
        dEddyVisc_im1halfjk_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i-1][j][k]
          +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
        dEddyVisc_ijp1halfk_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j+1][k]
          +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
        dEddyVisc_ijm1halfk_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j-1][k]
          +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
        dEddyVisc_ijkp1half_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j][k+1]
          +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
        dEddyVisc_ijkm1half_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j][k-1]
          +grid.dLocalGridNew[grid.nEddyVisc][i][j][k])*0.5;
        
        //calculate derived quantities
        dVSinTheta_ijp1halfk_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          *grid.dLocalGridNew[grid.nV][i][nJInt][k];
        dVSinTheta_ijm1halfk_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
          *grid.dLocalGridNew[grid.nV][i][nJInt-1][k];
        dUR2_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dRSq_im1half_n;
        dUR2_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dRSq_ip1half_n;
        dTSq_ijk_n=grid.dLocalGridOld[grid.nT][i][j][k]*grid.dLocalGridOld[grid.nT][i][j][k];
        dT4_ijk_n=dTSq_ijk_n*dTSq_ijk_n;
        dTSq_im1jk_n=grid.dLocalGridOld[grid.nT][i-1][j][k]*grid.dLocalGridOld[grid.nT][i-1][j][k];
        dT4_im1jk_n=dTSq_im1jk_n*dTSq_im1jk_n;
        dTSq_ijp1k_n=grid.dLocalGridOld[grid.nT][i][j+1][k]*grid.dLocalGridOld[grid.nT][i][j+1][k];
        dT4_ijp1k_n=dTSq_ijp1k_n*dTSq_ijp1k_n;
        dTSq_ijm1k_n=grid.dLocalGridOld[grid.nT][i][j-1][k]*grid.dLocalGridOld[grid.nT][i][j-1][k];
        dT4_ijm1k_n=dTSq_ijm1k_n*dTSq_ijm1k_n;
        dTSq_ijkp1_n=grid.dLocalGridOld[grid.nT][i][j][k+1]*grid.dLocalGridOld[grid.nT][i][j][k+1];
        dT4_ijkp1_n=dTSq_ijkp1_n*dTSq_ijkp1_n;
        dTSq_ijkm1_n=grid.dLocalGridOld[grid.nT][i][j][k-1]*grid.dLocalGridOld[grid.nT][i][j][k-1];
        dT4_ijkm1_n=dTSq_ijkm1_n*dTSq_ijkm1_n;
        dKappa_im1halfjk_n=(dT4_im1jk_n+dT4_ijk_n)/(dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]+dT4_im1jk_n
          /grid.dLocalGridOld[grid.nKappa][i-1][j][k]);
        dKappa_ijp1halfk_n=(dT4_ijp1k_n+dT4_ijk_n)/(dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]+dT4_ijp1k_n
          /grid.dLocalGridOld[grid.nKappa][i][j+1][k]);
        dKappa_ijm1halfk_n=(dT4_ijm1k_n+dT4_ijk_n)/(dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]+dT4_ijm1k_n
          /grid.dLocalGridOld[grid.nKappa][i][j-1][k]);
        dKappa_ijkp1half_n=(dT4_ijkp1_n+dT4_ijk_n)/(dT4_ijkp1_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k+1]+dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]);
        dKappa_ijkm1half_n=(dT4_ijkm1_n+dT4_ijk_n)/(dT4_ijkm1_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k-1]+dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]);
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP_ijk_n=dP_ijk_n+grid.dLocalGridOld[grid.nQ0][i][j][k]
            +grid.dLocalGridOld[grid.nQ1][i][j][k]+grid.dLocalGridOld[grid.nQ2][i][j][k];
        #endif
        
        //Calcuate dA1
        dA1CenGrad=(dE_ip1halfjk_n-dE_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0];
        dU_U0_Diff=(dU_ijk_np1half-dU0_i_np1half);
        if(dU_U0_Diff<0.0){//moving in the negative direction
          dA1UpWindGrad=dA1CenGrad;/**\BC grid.dLocalGridOld[grid.nDM][i+1][0][0] and 
            grid.dLocalGridOld[grid.nE][i+1][j][k] missing in the calculation of upwind gradient in 
            dA1. Using the centered gradient instead.*/
        }
        else{//moving in the postive direction*/
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k]
            -grid.dLocalGridOld[grid.nE][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=dU_U0_Diff*dRSq_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate dS1
        dS1=dP_ijk_n/grid.dLocalGridOld[grid.nD][i][j][k]
          *(dUR2_ip1halfjk_np1half-dUR2_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //Calcualte dA2
        dA2CenGrad=(dE_ijp1halfk_n-dE_ijm1halfk_n)/grid.dLocalGridOld[grid.nDTheta][0][j][0];
        if(dV_ijk_np1half<0.0){//moving in the negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j+1][k]
            -grid.dLocalGridOld[grid.nE][i][j][k])/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        }
        else{//moving in the positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k]
            -grid.dLocalGridOld[grid.nE][i][j-1][k])/(grid.dLocalGridOld[grid.nDTheta][0][j][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;
        }
        dA2=dV_ijk_np1half/dR_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA2CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
        
        //Calcualte dS2
        dS2=dP_ijk_n/(grid.dLocalGridOld[grid.nD][i][j][k]*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDTheta][0][j][0])
          *(dVSinTheta_ijp1halfk_np1half-dVSinTheta_ijm1halfk_np1half);
        
        //Calcualte dA3
        dA3CenGrad=(dE_ijkp1half_n-dE_ijkm1half_n)/grid.dLocalGridOld[grid.nDPhi][0][0][k];
        if(dW_ijk_np1half<0.0){//moving in the negative direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k+1]
            -grid.dLocalGridOld[grid.nE][i][j][k])/(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k])*2.0;
        }
        else{//moving in the positive direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k]
            -grid.dLocalGridOld[grid.nE][i][j][k-1])/(grid.dLocalGridOld[grid.nDPhi][0][0][k]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*2.0;
        }
        dA3=dW_ijk_np1half/(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0])*
          ((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA3CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA3UpWindGrad);
        
        //Calcualte dS3
        dS3=dP_ijk_n/(grid.dLocalGridOld[grid.nD][i][j][k]*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k])
          *(dW_ijkp1half_np1half-dW_ijkm1half_np1half);
        
        //Calculate dS4
        dTGrad_im1half=(dT4_ijk_n-dT4_im1jk_n)/(grid.dLocalGridOld[grid.nDM][i][0][0]
          +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        dGrad_ip1half=-3.0*dRSq_ip1half_n*dT4_ijk_n/(8.0*parameters.dPi);/**\BC
          Missing grid.dLocalGridOld[grid.nT][i+1][0][0]*/
        dGrad_im1half=dRhoAve_im1half_n*dR4_im1half_n/(dKappa_im1halfjk_n*dRho_im1halfjk_n)
          *dTGrad_im1half;
        dS4=16.0*parameters.dPi*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dGrad_ip1half-dGrad_im1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //Calculate dS5
        dTGrad_jp1half=(dT4_ijp1k_n-dT4_ijk_n)/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
          +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        dTGrad_jm1half=(dT4_ijk_n-dT4_ijm1k_n)/(grid.dLocalGridOld[grid.nDTheta][0][j][0]
          +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;;
        dGrad_jp1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          /(dKappa_ijp1halfk_n*dRho_ijp1halfk_n)*dTGrad_jp1half;
        dGrad_jm1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
          /(dKappa_ijm1halfk_n*dRho_ijm1halfk_n)*dTGrad_jm1half;;
        dS5=(dGrad_jp1half-dGrad_jm1half)/(grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *dRSq_i_n*grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //Calculate dS6
        dTGrad_kp1half=(dT4_ijkp1_n-dT4_ijk_n)/(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k])*2.0;
        dTGrad_km1half=(dT4_ijk_n-dT4_ijkm1_n)/(grid.dLocalGridOld[grid.nDPhi][0][0][k]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*2.0;;
        dGrad_kp1half=dTGrad_kp1half/(dKappa_ijkp1half_n*dRho_ijkp1half_n);
        dGrad_km1half=dTGrad_km1half/(dKappa_ijkm1half_n*dRho_ijkm1half_n);
        dS6=(dGrad_kp1half-dGrad_km1half)/(dRSq_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //calculate dT1
        dEGrad_ip1halfjk_np1half=dR4_ip1half_n*dEddyVisc_ip1halfjk_np1half*dRhoAve_ip1half_n
          *(grid.dLocalGridOld[grid.nE][i+1][j][k]-grid.dLocalGridOld[grid.nE][i][j][k])
          /(dRho_ip1halfjk_n*dDM_ip1half);
        dEGrad_im1halfjk_np1half=dR4_im1half_n*dEddyVisc_im1halfjk_np1half*dRhoAve_im1half_n
          *(grid.dLocalGridOld[grid.nE][i][j][k]-grid.dLocalGridOld[grid.nE][i-1][j][k])
          /(dRho_im1halfjk_n*dDM_im1half);
        dT1=16.0*dPiSq*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dEGrad_ip1halfjk_np1half
          -dEGrad_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //calculate dT2
        dEGrad_ijp1halfk_np1half=dEddyVisc_ijp1halfk_np1half
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          *(grid.dLocalGridOld[grid.nE][i][j+1][k]-grid.dLocalGridOld[grid.nE][i][j][k])
          /(dRho_ijp1halfk_n*dR_i_n*dDelTheta_jp1half);
        dEGrad_ijm1halfk_np1half=dEddyVisc_ijm1halfk_np1half
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
          *(grid.dLocalGridOld[grid.nE][i][j][k]-grid.dLocalGridOld[grid.nE][i][j-1][k])
          /(dRho_ijm1halfk_n*dR_i_n*dDelTheta_jm1half);
        dT2=(dEGrad_ijp1halfk_np1half-dEGrad_ijm1halfk_np1half)/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //calculate dT3
        dEGrad_ijkp1half_np1half=dEddyVisc_ijkp1half_np1half*(grid.dLocalGridOld[grid.nE][i][j][k+1]
          -grid.dLocalGridOld[grid.nE][i][j][k])/(dRho_ijkp1half_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*dR_i_n*dDelPhi_kp1half);
        dEGrad_ijkm1half_np1half=dEddyVisc_ijkm1half_np1half*(grid.dLocalGridOld[grid.nE][i][j][k]
          -grid.dLocalGridOld[grid.nE][i][j][k-1])/(dRho_ijkm1half_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*dR_i_n*dDelPhi_km1half);
        dT3=(dEGrad_ijkp1half_np1half-dEGrad_ijkm1half_np1half)/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //eddy viscosity terms
        dEddyViscosityTerms=(dT1+dT2+dT3)/parameters.dPrt;
        
        //calculate new energy
        grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridOld[grid.nE][i][j][k]
          -time.dDeltat_np1half*(4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1
          +dS1)+dA2+dS2+dA3+dS3-4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])
          *(dS4+dS5+dS6)-dEddyViscosityTerms);
        
        if(grid.dLocalGridNew[grid.nE][i][j][k]<0.0){
          
          #if SIGNEGENG==1
            raise(SIGINT);
          #endif
          
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
            <<": negative energy calculated in , ("<<i<<","<<j<<","<<k<<")\n";
          throw exception2(ssTemp.str(),CALCULATION);
          
        }
      }
    }
  }
}
