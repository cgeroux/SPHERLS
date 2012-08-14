void calNewU_RTP_LES(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop){
  int i;
  int j;
  int k;
  int nICen;
  int nJInt;
  int nKInt;
  double dRho_ip1halfjk_n;
  double dP_ip1jk_n;
  double dP_ijk_n;
  double dA1CenGrad;
  double dA1UpWindGrad;
  double dA1;
  double dS1;
  double dA2CenGrad;
  double dA2UpWindGrad;
  double dA3CenGrad;
  double dA3UpWindGrad;
  double dA2;
  double dS2;
  double dA3;
  double dS3;
  double dS4;
  double dTA1;
  double dTS1;
  double dTA2;
  double dTS2;
  double dTA3;
  double dTS3;
  double dTS4;
  double dDivU_ijk_n;
  double dDivU_ip1jk_n;
  double dTau_rr_ip1jk_n;
  double dTau_rr_ijk_n;
  double dTau_rt_ip1halfjp1halfk_n;
  double dTau_rt_ip1halfjm1halfk_n;
  double dTau_rp_ip1halfjkp1half_n;
  double dTau_rp_ip1halfjkm1half_n;
  double dR_i_n;
  double dR_ip1_n;
  double dU_ip1jk_nm1half;
  double dU_ijk_nm1half;
  double dU_ip1halfjp1halfk_nm1half;
  double dU_ip1halfjm1halfk_nm1half;
  double dU_ip1halfjkp1half_nm1half;
  double dU_ip1halfjkm1half_nm1half;
  double dU0_ip1_nm1half;
  double dU0_i_nm1half;
  double dUmU0_ip1halfjk_nm1half;
  double dV_ip1halfjk_nm1half;
  double dV_ip1halfjp1halfk_nm1half;
  double dV_ip1halfjm1halfk_nm1half;
  double dV_ip1jk_nm1half;
  double dV_ijk_nm1half;
  double dW_ip1halfjk_nm1half;
  double dW_ip1halfjkp1half_nm1half;
  double dW_ip1halfjkm1half_nm1half;
  double dV_R_ip1jk_n;
  double dV_R_ip1jp1halfk_n;
  double dV_R_ip1jm1halfk_n;
  double dV_R_ijp1halfk_n;
  double dV_R_ijm1halfk_n;
  double dV_R_ijk_n;
  double dW_R_ip1jkp1half_n;
  double dW_R_ijkp1half_n;
  double dW_R_ip1jkm1half_n;
  double dW_R_ijkm1half_n;
  double dRSq_ip1half_n;
  double dRSq_im1half_n;
  double dRSq_ip1_n;
  double dRSq_i_n;
  double dRSq_ip3half_n;
  double dRCu_ip1half_n;
  double dRSqUmU0_ip3halfjk_n;
  double dRSqUmU0_ip1halfjk_n;
  double dRSqUmU0_im1halfjk_n;
  double dRSqUmU0_ijk_n;//needed at surface boundary
  double dRhoR_ip1halfjk_n;
  double dDM_ip1half;
  double dDTheta_jp1half;
  double dDTheta_jm1half;
  double dDPhi_kp1half;
  double dDPhi_km1half;
  double dEddyVisc_ip1halfjk_n;
  double dEddyVisc_ip1halfjp1halfk_n;
  double dEddyVisc_ip1halfjm1halfk_n;
  double dEddyVisc_ip1halfjkp1half_n;
  double dEddyVisc_ip1halfjkm1half_n;
  double dRhoAve_ip1half_n;
  double dEddyViscosityTerms;
  double dDonorFrac_ip1half;
  
  //calculate new u
  for(i=grid.nStartUpdateExplicit[grid.nU][0];i<grid.nEndUpdateExplicit[grid.nU][0];i++){
    
    //calculate i of centered quantities
    nICen=i-grid.nCenIntOffset[0];
    
    //calculate quantities that vary only with radius
    dR_ip1_n=(grid.dLocalGridOld[grid.nR][i+1][0][0]+grid.dLocalGridOld[grid.nR][i][0][0])*0.5;
    dR_i_n=(grid.dLocalGridOld[grid.nR][i][0][0]+grid.dLocalGridOld[grid.nR][i-1][0][0])*0.5;
    dRSq_ip1_n=dR_ip1_n*dR_ip1_n;
    dRSq_i_n=dR_i_n*dR_i_n;
    dRSq_ip1half_n=grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nR][i][0][0];
    dRSq_im1half_n=grid.dLocalGridOld[grid.nR][i-1][0][0]*grid.dLocalGridOld[grid.nR][i-1][0][0];
    dRSq_ip3half_n=grid.dLocalGridOld[grid.nR][i+1][0][0]*grid.dLocalGridOld[grid.nR][i+1][0][0];
    dRCu_ip1half_n=dRSq_ip1half_n*grid.dLocalGridOld[grid.nR][i][0][0];
    dDM_ip1half=(grid.dLocalGridOld[grid.nDM][nICen+1][0][0]
      +grid.dLocalGridOld[grid.nDM][nICen][0][0])*0.5;
    dRhoAve_ip1half_n=(grid.dLocalGridOld[grid.nDenAve][nICen+1][0][0]
      +grid.dLocalGridOld[grid.nDenAve][nICen][0][0])*0.5;
    dU0_ip1_nm1half=(grid.dLocalGridOld[grid.nU0][i+1][0][0]
      +grid.dLocalGridOld[grid.nU0][i][0][0])*0.5;
    dU0_i_nm1half=(grid.dLocalGridOld[grid.nU0][i][0][0]
      +grid.dLocalGridOld[grid.nU0][i-1][0][0])*0.5;
    dDonorFrac_ip1half=(grid.dLocalGridOld[grid.nDonorCellFrac][nICen+1][0][0]
          +grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nU][1];j<grid.nEndUpdateExplicit[grid.nU][1];j++){
      
      //calculate j of interface quantities
      nJInt=j+grid.nCenIntOffset[1];
      
      //calculating quantities that vary only with theta, and perhaps radius
      dDTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
        +grid.dLocalGridOld[grid.nDTheta][0][j][0])*0.5;
      dDTheta_jm1half=(grid.dLocalGridOld[grid.nDTheta][0][j-1][0]
        +grid.dLocalGridOld[grid.nDTheta][0][j][0])*0.5;
      
      for(k=grid.nStartUpdateExplicit[grid.nU][2];k<grid.nEndUpdateExplicit[grid.nU][2];k++){
        
        //calculate k of interface quantities
        nKInt=k+grid.nCenIntOffset[2];
        dDPhi_kp1half=(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k])*0.5;
        dDPhi_km1half=(grid.dLocalGridOld[grid.nDPhi][0][0][k]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*0.5;
        
        //CALCULATE INTERPOLATED QUANTITIES
        dU_ip1jk_nm1half=(grid.dLocalGridOld[grid.nU][i+1][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
        dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i-1][j][k])*0.5;
        dUmU0_ip1halfjk_nm1half=grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0];
        dU_ip1halfjp1halfk_nm1half=(grid.dLocalGridOld[grid.nU][i][j+1][k]
          +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
        dU_ip1halfjm1halfk_nm1half=(grid.dLocalGridOld[grid.nU][i][j-1][k]
          +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
        dU_ip1halfjkp1half_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k+1])*0.5;
        dU_ip1halfjkm1half_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k-1])*0.5;
        dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
        dV_ip1halfjk_nm1half=0.25*(grid.dLocalGridOld[grid.nV][nICen+1][nJInt][k]
          +grid.dLocalGridOld[grid.nV][nICen+1][nJInt-1][k]
          +grid.dLocalGridOld[grid.nV][nICen][nJInt][k]
          +grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k]);
        dV_ip1halfjp1halfk_nm1half=(grid.dLocalGridOld[grid.nV][nICen+1][nJInt][k]
          +grid.dLocalGridOld[grid.nV][nICen][nJInt][k])*0.5;
        dV_ip1halfjm1halfk_nm1half=(grid.dLocalGridOld[grid.nV][nICen+1][nJInt-1][k]
          +grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k])*0.5;
        dV_ip1jk_nm1half=(grid.dLocalGridOld[grid.nV][nICen+1][nJInt][k]
          +grid.dLocalGridOld[grid.nV][nICen+1][nJInt-1][k])*0.5;
        dV_ijk_nm1half=(grid.dLocalGridOld[grid.nV][nICen][nJInt][k]
          +grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k])*0.5;
        dW_ip1halfjk_nm1half=(grid.dLocalGridOld[grid.nW][nICen+1][j][nKInt]
          +grid.dLocalGridOld[grid.nW][nICen+1][j][nKInt-1]
          +grid.dLocalGridOld[grid.nW][nICen][j][nKInt]
          +grid.dLocalGridOld[grid.nW][nICen][j][nKInt-1])*0.25;
        dW_ip1halfjkp1half_nm1half=(grid.dLocalGridOld[grid.nW][nICen+1][j][nKInt]
          +grid.dLocalGridOld[grid.nW][nICen][j][nKInt])*0.5;
        dW_ip1halfjkm1half_nm1half=(grid.dLocalGridOld[grid.nW][nICen+1][j][nKInt-1]
          +grid.dLocalGridOld[grid.nW][nICen][j][nKInt-1])*0.5;
        dP_ip1jk_n=grid.dLocalGridOld[grid.nP][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nQ0][nICen+1][j][k]+grid.dLocalGridOld[grid.nQ1][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nQ2][nICen+1][j][k];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][nICen][j][k]
          +grid.dLocalGridOld[grid.nQ0][nICen][j][k]+grid.dLocalGridOld[grid.nQ1][nICen][j][k]
          +grid.dLocalGridOld[grid.nQ2][nICen][j][k];
        dEddyVisc_ip1halfjk_n=(grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen+1][j][k])*0.5;
        dEddyVisc_ip1halfjp1halfk_n=(grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen][j+1][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen+1][j+1][k])*0.25;
        dEddyVisc_ip1halfjm1halfk_n=(grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen][j-1][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen+1][j-1][k])*0.25;
        dEddyVisc_ip1halfjkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k+1]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen+1][j][k+1]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen+1][j][k])*0.25;
        dEddyVisc_ip1halfjkm1half_n=(grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k-1]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen+1][j][k-1])*0.25;
        
        //calculate derived quantities
        dRSqUmU0_ip3halfjk_n=dRSq_ip3half_n*(grid.dLocalGridOld[grid.nU][i+1][j][k]
          -grid.dLocalGridOld[grid.nU0][i+1][0][0]);
        dRSqUmU0_ip1halfjk_n=dRSq_ip1half_n*(grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0]);
        dRSqUmU0_im1halfjk_n=dRSq_im1half_n*(grid.dLocalGridOld[grid.nU][i-1][j][k]
          -grid.dLocalGridOld[grid.nU0][i-1][0][0]);
        dV_R_ip1jk_n=dV_ip1jk_nm1half/dR_ip1_n;
        dV_R_ip1jp1halfk_n=grid.dLocalGridOld[grid.nV][nICen+1][nJInt][k]/dR_ip1_n;
        dV_R_ip1jm1halfk_n=grid.dLocalGridOld[grid.nV][nICen+1][nJInt-1][k]/dR_ip1_n;
        dV_R_ijp1halfk_n=grid.dLocalGridOld[grid.nV][nICen][nJInt][k]/dR_i_n;
        dV_R_ijm1halfk_n=grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k]/dR_i_n;
        dV_R_ijk_n=dV_ijk_nm1half/dR_i_n;
        dW_R_ip1jkp1half_n=grid.dLocalGridOld[grid.nW][nICen+1][j][nKInt]/dR_ip1_n;
        dW_R_ijkp1half_n=grid.dLocalGridOld[grid.nW][nICen][j][nKInt]/dR_i_n;
        dW_R_ip1jkm1half_n=grid.dLocalGridOld[grid.nW][nICen+1][j][nKInt-1]/dR_ip1_n;
        dW_R_ijkm1half_n=grid.dLocalGridOld[grid.nW][nICen][j][nKInt-1]/dR_i_n;
        dRhoR_ip1halfjk_n=dRho_ip1halfjk_n*grid.dLocalGridOld[grid.nR][i][0][0];
        
        //Calculate dA1
        dA1CenGrad=(dU_ip1jk_nm1half-dU_ijk_nm1half)
          /(grid.dLocalGridOld[grid.nDM][nICen+1][0][0]
          +grid.dLocalGridOld[grid.nDM][nICen][0][0])*2.0;
        dA1UpWindGrad=0.0;
        if(dUmU0_ip1halfjk_nm1half<0.0){//moving from outside in
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nU][i+1][j][k]
            -grid.dLocalGridOld[grid.nU][i][j][k])
            /grid.dLocalGridOld[grid.nDM][nICen+1][0][0];
        }
        else{//moving from inside out
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU][i-1][j][k])
            /grid.dLocalGridOld[grid.nDM][nICen][0][0];
        }
        dA1=dUmU0_ip1halfjk_nm1half*((1.0-dDonorFrac_ip1half)*dA1CenGrad+dDonorFrac_ip1half
          *dA1UpWindGrad);
        
        //calculate dS1
        dS1=(dP_ip1jk_n-dP_ijk_n)/(dDM_ip1half*dRho_ip1halfjk_n);
        
        //Calculate dS4
        dS4=parameters.dG*grid.dLocalGridOld[grid.nM][i][0][0]/dRSq_ip1half_n;
        
        //Calculate dA2
        dA2CenGrad=(dU_ip1halfjp1halfk_nm1half-dU_ip1halfjm1halfk_nm1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        dA2UpWindGrad=0.0;
        if(dV_ip1halfjk_nm1half>0.0){//moving in positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU][i][j-1][k])
            /(grid.dLocalGridOld[grid.nDTheta][0][j][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;
        }
        else{//moving in negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j+1][k]
            -grid.dLocalGridOld[grid.nU][i][j][k])
            /(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        }
        dA2=dV_ip1halfjk_nm1half*((1.0-dDonorFrac_ip1half)*dA2CenGrad
          +dDonorFrac_ip1half*dA2UpWindGrad)/grid.dLocalGridOld[grid.nR][i][0][0];
        
        //Calculate dS2
        dS2=dV_ip1halfjk_nm1half*dV_ip1halfjk_nm1half
          /grid.dLocalGridOld[grid.nR][i][0][0];
        
        //Calculate dA3
        dA3CenGrad=(dU_ip1halfjkp1half_nm1half-dU_ip1halfjkm1half_nm1half)
          /grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dA3UpWindGrad=0.0;
        if(dW_ip1halfjk_nm1half>0.0){//moving in positive direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU][i][j][k-1])
            /(grid.dLocalGridOld[grid.nDPhi][0][0][k]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*2.0;
        }
        else{//moving in negative direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k+1]
            -grid.dLocalGridOld[grid.nU][i][j][k])
            /(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k])*2.0;
        }
        dA3=dW_ip1halfjk_nm1half*((1.0-dDonorFrac_ip1half)*dA3CenGrad+dDonorFrac_ip1half
          *dA3UpWindGrad)/(grid.dLocalGridOld[grid.nR][i][0][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
        
        //Calculate dS3
        dS3=dW_ip1halfjk_nm1half*dW_ip1halfjk_nm1half/grid.dLocalGridOld[grid.nR][i][0][0];
        
        //cal DivU_ip1jk_n
        dDivU_ip1jk_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][nICen+1][0][0]
          *(dRSqUmU0_ip3halfjk_n-dRSqUmU0_ip1halfjk_n)/grid.dLocalGridOld[grid.nDM][nICen+1][0][0]
          +(grid.dLocalGridOld[grid.nV][nICen+1][nJInt][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          -grid.dLocalGridOld[grid.nV][nICen+1][nJInt-1][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0])
          /(grid.dLocalGridOld[grid.nDTheta][0][j][0]*dR_ip1_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0])
          +(grid.dLocalGridOld[grid.nW][nICen+1][j][nKInt]
          -grid.dLocalGridOld[grid.nW][nICen+1][j][nKInt-1])
          /(grid.dLocalGridOld[grid.nDPhi][0][0][k]*dR_ip1_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
        
        //cal DivU_ijk_n
        dDivU_ijk_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][nICen][0][0]
          *(dRSqUmU0_ip1halfjk_n-dRSqUmU0_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][nICen][0][0]
          +(grid.dLocalGridOld[grid.nV][nICen][nJInt][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          -grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0])
          /(grid.dLocalGridOld[grid.nDTheta][0][j][0]*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0])
          +(grid.dLocalGridOld[grid.nW][nICen][j][nKInt]
          -grid.dLocalGridOld[grid.nW][nICen][j][nKInt-1])
          /(grid.dLocalGridOld[grid.nDPhi][0][0][k]*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
        
        //cal Tau_rr_ip1jk_n
        dTau_rr_ip1jk_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][nICen+1][j][k]*(4.0*parameters.dPi
          *dRSq_ip1_n*grid.dLocalGridOld[grid.nDenAve][nICen+1][0][0]
          *((grid.dLocalGridOld[grid.nU][i+1][j][k]-grid.dLocalGridOld[grid.nU0][i+1][0][0])
          -(grid.dLocalGridOld[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU0][i][0][0]))
          /grid.dLocalGridOld[grid.nDM][nICen+1][0][0]-0.3333333333333333*dDivU_ip1jk_n);
        
        //cal Tau_rr_ijk_n
        dTau_rr_ijk_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]*(4.0*parameters.dPi
          *dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][nICen][0][0]
          *((grid.dLocalGridOld[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU0][i][0][0])
          -(grid.dLocalGridOld[grid.nU][i-1][j][k]-grid.dLocalGridOld[grid.nU0][i-1][0][0]))
          /grid.dLocalGridOld[grid.nDM][nICen][0][0]-0.3333333333333333*dDivU_ijk_n);
        
        //calculate dTau_rt_ip1halfjp1halfk_n
        dTau_rt_ip1halfjp1halfk_n=dEddyVisc_ip1halfjp1halfk_n*(4.0*parameters.dPi*dRCu_ip1half_n
          *dRhoAve_ip1half_n*(dV_R_ip1jp1halfk_n-dV_R_ijp1halfk_n)/dDM_ip1half
          +((grid.dLocalGridOld[grid.nU][i][j+1][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0])-(grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0]))/(dDTheta_jp1half
          *grid.dLocalGridOld[grid.nR][i][0][0]));
        
        //calculate dTau_rt_ip1halfjm1halfk_n
        dTau_rt_ip1halfjm1halfk_n=dEddyVisc_ip1halfjm1halfk_n*(4.0*parameters.dPi*dRCu_ip1half_n
          *dRhoAve_ip1half_n*(dV_R_ip1jm1halfk_n-dV_R_ijm1halfk_n)/dDM_ip1half
          +((grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0])-(grid.dLocalGridOld[grid.nU][i][j-1][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0]))/(grid.dLocalGridOld[grid.nR][i][0][0]
          *dDTheta_jm1half));
        
        //calculate dTau_rp_ip1halfjkp1half_n
        dTau_rp_ip1halfjkp1half_n=dEddyVisc_ip1halfjkp1half_n*(4.0*parameters.dPi*dRCu_ip1half_n
          *dRhoAve_ip1half_n*(dW_R_ip1jkp1half_n-dW_R_ijkp1half_n)/dDM_ip1half
          +((grid.dLocalGridOld[grid.nU][i][j][k+1]-grid.dLocalGridOld[grid.nU0][i][0][0])
          -(grid.dLocalGridOld[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU0][i][0][0]))
          /(grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nR][i][0][0]
          *dDPhi_kp1half));
          
        //calculate dTau_rp_im1halfjkm1half_n
        dTau_rp_ip1halfjkm1half_n=dEddyVisc_ip1halfjkm1half_n*(4.0*parameters.dPi*dRCu_ip1half_n
          *dRhoAve_ip1half_n*(dW_R_ip1jkm1half_n-dW_R_ijkm1half_n)/dDM_ip1half
          +((grid.dLocalGridOld[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU0][i][0][0])
          -(grid.dLocalGridOld[grid.nU][i][j][k-1]-grid.dLocalGridOld[grid.nU0][i][0][0]))
          /(grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nR][i][0][0]
          *dDPhi_km1half));
        
        //cal dTA1
        dTA1=(dTau_rr_ip1jk_n-dTau_rr_ijk_n)/(dDM_ip1half*dRho_ip1halfjk_n);
        
        //cal dTS1
        dTS1=dEddyVisc_ip1halfjk_n/dRhoR_ip1halfjk_n*(4.0
          *((dU_ip1jk_nm1half-dU0_ip1_nm1half)-(dU_ijk_nm1half-dU0_i_nm1half))/dDM_ip1half
          +grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]
          *(dV_R_ip1jk_n-dV_R_ijk_n)/dDM_ip1half);
        
        //calculate dTA2
        dTA2=(dTau_rt_ip1halfjp1halfk_n-dTau_rt_ip1halfjm1halfk_n)
          /(grid.dLocalGridOld[grid.nDTheta][0][j][0]*dRhoR_ip1halfjk_n);
        
        //calculate dTS2
        dTS2=(2.0*(dV_ip1halfjp1halfk_nm1half-dV_ip1halfjm1halfk_nm1half)
          -grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]*((dU_ip1halfjp1halfk_nm1half
          -grid.dLocalGridOld[grid.nU0][i][0][0])-(dU_ip1halfjm1halfk_nm1half
          -grid.dLocalGridOld[grid.nU0][i][0][0])))/(grid.dLocalGridOld[grid.nR][i][0][0]
          *grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //calculate dTA3
        dTA3=(dTau_rp_ip1halfjkp1half_n-dTau_rp_ip1halfjkm1half_n)/(dRho_ip1halfjk_n
          *grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //calculate dTS3
        dTS3=2.0*(dW_ip1halfjkp1half_nm1half-dW_ip1halfjkm1half_nm1half)
          /(grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //cal dTS4
        dTS4=(4.0*(grid.dLocalGridOld[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU0][i][0][0])
          +2.0*grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]*dV_ip1halfjk_nm1half)
          /grid.dLocalGridOld[grid.nR][i][0][0];
        
        dEddyViscosityTerms=-4.0*parameters.dPi*dRhoAve_ip1half_n*dRSq_ip1half_n*(dTA1+dTS1)-dTA2
          -dTA3+dEddyVisc_ip1halfjk_n/dRhoR_ip1halfjk_n*(dTS2+dTS3+dTS4);
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nU][i][j][k]=grid.dLocalGridOld[grid.nU][i][j][k]
          -time.dDeltat_n*(4.0*parameters.dPi*dRhoAve_ip1half_n*dRSq_ip1half_n*(dA1+dS1)
          +dA2-dS2+dA3-dS3+dS4+dEddyViscosityTerms);
        
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
        
        //add M_r
        ssName.str("");
        ssName<<"M_r";
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,grid.dLocalGridOld[grid.nM][i][0][0]);
        
        //add A1
        ssName.str("");
        ssName<<"U_A1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-4.0*parameters.dPi*dRhoAve_ip1half_n*dRSq_ip1half_n*(dA1));
          
        //add S1
        ssName.str("");
        ssName<<"U_S1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-4.0*parameters.dPi*dRhoAve_ip1half_n*dRSq_ip1half_n*(dS1));
        
        //add A2
        ssName.str("");
        ssName<<"U_A2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-1.0*dA2);
          
        //add S2
        ssName.str("");
        ssName<<"U_S2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,dS2);
          
        //add A3
        ssName.str("");
        ssName<<"U_A3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-1.0*dA3);
          
        //add S3
        ssName.str("");
        ssName<<"U_S3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,dS3);
          
        //add S4
        ssName.str("");
        ssName<<"U_S4"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-1.0*dS4);
          
        //add dEddyViscosityTerms
        ssName.str("");
        ssName<<"U_EV"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
        ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-1.0*dEddyViscosityTerms);
        
        //add DuDt
        ssName.str("");
        ssName<<"U_DuDt"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,(grid.dLocalGridNew[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU][i][j][k])
          /time.dDeltat_n);
        #endif
      }
    }
  }
  
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nU][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nU][0][0];i++){
    
    //calculate i of centered quantities
    nICen=i-grid.nCenIntOffset[0];
    
    //calculate quantities that vary only with radius
    dR_i_n=(grid.dLocalGridOld[grid.nR][i][0][0]+grid.dLocalGridOld[grid.nR][i-1][0][0])*0.5;
    dRSq_i_n=dR_i_n*dR_i_n;
    dRSq_ip1half_n=grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nR][i][0][0];
    dRCu_ip1half_n=grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nR][i][0][0]
      *grid.dLocalGridOld[grid.nR][i][0][0];
    dDM_ip1half=(grid.dLocalGridOld[grid.nDM][nICen][0][0])*(0.5+parameters.dAlpha
      +parameters.dAlphaExtra);/**\BC Missing grid.dLocalGridOld[grid.nDM][i+1][0][0] in calculation
      of \f$S_1\f$ using \ref Parameters.dAlpha *grid.dLocalGridOld[grid.nDM][nICen][0][0] instead.
      */
    dRhoAve_ip1half_n=(grid.dLocalGridOld[grid.nDenAve][nICen][0][0])*0.5;/**\BC Missing density
      outside of surface, setting it to zero.*/
    dU0_i_nm1half=(grid.dLocalGridOld[grid.nU0][i][0][0]+grid.dLocalGridOld[grid.nU0][i-1][0][0])
      *0.5;
    dR_ip1_n=grid.dLocalGridOld[grid.nR][i][0][0];
    dDonorFrac_ip1half=grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0];
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nU][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nU][0][1];j++){
      
      //calculate j of interface quantities
      nJInt=j+grid.nCenIntOffset[1];
      
      //calculating quantities that vary only with theta, and perhaps radius
      dDTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
        +grid.dLocalGridOld[grid.nDTheta][0][j][0])*0.5;
      dDTheta_jm1half=(grid.dLocalGridOld[grid.nDTheta][0][j-1][0]
        +grid.dLocalGridOld[grid.nDTheta][0][j][0])*0.5;
      
      for(k=grid.nStartGhostUpdateExplicit[grid.nU][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nU][0][2];k++){
        
        //calculate k of interface quantities
        nKInt=k+grid.nCenIntOffset[2];
        dDPhi_kp1half=(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
        +grid.dLocalGridOld[grid.nDPhi][0][0][k])*0.5;
        dDPhi_km1half=(grid.dLocalGridOld[grid.nDPhi][0][0][k]
        +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*0.5;
        
        //CALCULATE INTERPOLATED QUANTITIES
        dU_ip1jk_nm1half=grid.dLocalGridOld[grid.nU][i][j][k];
        dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i-1][j][k])*0.5;
        dUmU0_ip1halfjk_nm1half=grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0];
        dU_ip1halfjp1halfk_nm1half=(grid.dLocalGridOld[grid.nU][i][j+1][k]
          +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
        dU_ip1halfjm1halfk_nm1half=(grid.dLocalGridOld[grid.nU][i][j-1][k]
          +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
        dU_ip1halfjkp1half_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k+1])*0.5;
        dU_ip1halfjkm1half_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k-1])*0.5;
        dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;/**\BC Missing density 
          outside model, setting it to zero. */
        dV_ip1halfjk_nm1half=0.5*(grid.dLocalGridOld[grid.nV][nICen][nJInt][k]
          +grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k]);/**\BC assuming theta and phi velocity 
          same outside star as inside.*/
        dV_ip1halfjp1halfk_nm1half=grid.dLocalGridOld[grid.nV][nICen][nJInt][k];/**\BC Assuming 
          theta velocities are constant across surface.*/
        dV_ip1halfjm1halfk_nm1half=grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k];
        dV_ip1jk_nm1half=(grid.dLocalGridOld[grid.nV][nICen][nJInt][k]
          +grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k])*0.5;/**\BC assuming that $V$ at
          $i+1$ is equal to $v$ at $i$.*/
        dV_ijk_nm1half=(grid.dLocalGridOld[grid.nV][nICen][nJInt][k]
          +grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k])*0.5;
        dW_ip1halfjk_nm1half=(grid.dLocalGridOld[grid.nW][nICen][j][nKInt]
          +grid.dLocalGridOld[grid.nW][nICen][j][nKInt-1])*0.5;
        dW_ip1halfjkp1half_nm1half=grid.dLocalGridOld[grid.nW][nICen][j][nKInt];
        dW_ip1halfjkm1half_nm1half=grid.dLocalGridOld[grid.nW][nICen][j][nKInt-1];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][nICen][j][k]
          +grid.dLocalGridOld[grid.nQ0][nICen][j][k]+grid.dLocalGridOld[grid.nQ1][nICen][j][k];
        dP_ip1jk_n=-1.0*dP_ijk_n;/**\BC Missing pressure outside surface setting it equal to 
          negative pressure in the center of the first cell so that it will be zero at surface.*/
        dEddyVisc_ip1halfjk_n=grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]*0.5;/**\BC assume 
          viscosity is zero outside the star.*/
        dEddyVisc_ip1halfjp1halfk_n=(grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen][j+1][k])*0.25;
        dEddyVisc_ip1halfjm1halfk_n=(grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen][j-1][k])*0.25;
        dEddyVisc_ip1halfjkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k+1]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k])*0.25;
        dEddyVisc_ip1halfjkm1half_n=(grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k-1])*0.25;
        
        //calculate derived quantities
        dRSqUmU0_ijk_n=dRSq_i_n*(dU_ijk_nm1half-dU0_i_nm1half);
        dRSqUmU0_ip1halfjk_n=dRSq_ip1half_n*(grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0]);
        dRSqUmU0_im1halfjk_n=dRSq_im1half_n*(grid.dLocalGridOld[grid.nU][i-1][j][k]
          -grid.dLocalGridOld[grid.nU0][i-1][0][0]);
        dV_R_ip1jk_n=dV_ip1jk_nm1half/dR_ip1_n;
        dV_R_ip1jp1halfk_n=grid.dLocalGridOld[grid.nV][nICen][nJInt][k]/dR_ip1_n;
        dV_R_ip1jm1halfk_n=grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k]/dR_ip1_n;
        dV_R_ijp1halfk_n=grid.dLocalGridOld[grid.nV][nICen][nJInt][k]/dR_i_n;
        dV_R_ijm1halfk_n=grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k]/dR_i_n;
        dV_R_ijk_n=dV_ijk_nm1half/dR_i_n;
        dW_R_ip1jkp1half_n=grid.dLocalGridOld[grid.nW][nICen][j][nKInt]/dR_ip1_n;
        dW_R_ijkp1half_n=grid.dLocalGridOld[grid.nW][nICen][j][nKInt]/dR_i_n;
        dW_R_ip1jkm1half_n=grid.dLocalGridOld[grid.nW][nICen][j][nKInt-1]/dR_ip1_n;
        dW_R_ijkm1half_n=grid.dLocalGridOld[grid.nW][nICen][j][nKInt-1]/dR_i_n;
        dRhoR_ip1halfjk_n=dRho_ip1halfjk_n*grid.dLocalGridOld[grid.nR][i][0][0];
        
        //Calculate dA1
        dA1CenGrad=(dU_ip1jk_nm1half-dU_ijk_nm1half)
          /(grid.dLocalGridOld[grid.nDM][nICen][0][0]*0.5);/**\BC Missing mass outside model,
          setting it to zero.*/
        dA1UpWindGrad=0.0;
        if(dUmU0_ip1halfjk_nm1half<0.0){//moving from outside in
          dA1UpWindGrad=dA1CenGrad;/**\BC Missing grid.dLocalGridOld[grid.nU][i+1][j][k] and 
            grid.dLocalGridOld[grid.nDM][nICen+1][0][0] in calculation of upwind gradient, when 
            moving inward. Using centered gradient instead.*/
        }
        else{//moving from inside out
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU][i-1][j][k])/grid.dLocalGridOld[grid.nDM][nICen][0][0];
        }
        dA1=dUmU0_ip1halfjk_nm1half*((1.0-dDonorFrac_ip1half)*dA1CenGrad+dDonorFrac_ip1half
          *dA1UpWindGrad);
        
        //calculate dS1
        dS1=(dP_ip1jk_n-dP_ijk_n)/(dDM_ip1half*dRho_ip1halfjk_n);
        
        //Calculate dS4
        dS4=parameters.dG*grid.dLocalGridOld[grid.nM][i][0][0]/dRSq_ip1half_n;
        
        //Calculate dA2
        dA2CenGrad=(dU_ip1halfjp1halfk_nm1half-dU_ip1halfjm1halfk_nm1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        dA2UpWindGrad=0.0;
        if(dV_ip1halfjk_nm1half>0.0){//moving in positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU][i][j-1][k])/(grid.dLocalGridOld[grid.nDTheta][0][j][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;
        }
        else{//moving in negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j+1][k]
            -grid.dLocalGridOld[grid.nU][i][j][k])/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        }
        dA2=dV_ip1halfjk_nm1half*((1.0-dDonorFrac_ip1half)*dA2CenGrad
          +dDonorFrac_ip1half*dA2UpWindGrad)/grid.dLocalGridOld[grid.nR][i][0][0];
        
        //Calculate dS2
        dS2=dV_ip1halfjk_nm1half*dV_ip1halfjk_nm1half/grid.dLocalGridOld[grid.nR][i][0][0];
        
        //Calculate dA3
        dA3CenGrad=(dU_ip1halfjkp1half_nm1half-dU_ip1halfjkm1half_nm1half)
          /grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dA3UpWindGrad=0.0;
        if(dW_ip1halfjk_nm1half>0.0){//moving in positive direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU][i][j][k-1])
            /(grid.dLocalGridOld[grid.nDPhi][0][0][k]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*2.0;
        }
        else{//moving in negative direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k+1]
            -grid.dLocalGridOld[grid.nU][i][j][k])
            /(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k])*2.0;
        }
        dA3=dW_ip1halfjk_nm1half*((1.0-dDonorFrac_ip1half)*dA3CenGrad+dDonorFrac_ip1half
          *dA3UpWindGrad)/(grid.dLocalGridOld[grid.nR][i][0][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
        
        //Calculate dS3
        dS3=dW_ip1halfjk_nm1half*dW_ip1halfjk_nm1half/grid.dLocalGridOld[grid.nR][i][0][0];
        
        //cal DivU_ip1jk_n
        dDivU_ip1jk_n=4.0*parameters.dPi*dRhoAve_ip1half_n
          *(dRSqUmU0_ip1halfjk_n-dRSqUmU0_ijk_n)/grid.dLocalGridOld[grid.nDM][nICen][0][0]*2.0
          +(grid.dLocalGridOld[grid.nV][nICen][nJInt][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          -grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0])
          /(grid.dLocalGridOld[grid.nDTheta][0][j][0]*dR_ip1_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0])
          +(grid.dLocalGridOld[grid.nW][nICen][j][nKInt]
          -grid.dLocalGridOld[grid.nW][nICen][j][nKInt-1])
          /(grid.dLocalGridOld[grid.nDPhi][0][0][k]*dR_ip1_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
        
        //cal DivU_ijk_n
        dDivU_ijk_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][nICen][0][0]
          *(dRSqUmU0_ip1halfjk_n-dRSqUmU0_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][nICen][0][0]
          +(grid.dLocalGridOld[grid.nV][nICen][nJInt][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          -grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0])
          /(grid.dLocalGridOld[grid.nDTheta][0][j][0]*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0])
          +(grid.dLocalGridOld[grid.nW][nICen][j][nKInt]
          -grid.dLocalGridOld[grid.nW][nICen][j][nKInt-1])
          /(grid.dLocalGridOld[grid.nDPhi][0][0][k]*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
        
        //cal Tau_rr_ip1jk_n
        dTau_rr_ip1jk_n=2.0*dEddyVisc_ip1halfjk_n*(4.0*parameters.dPi*dRSq_ip1half_n
          *dRhoAve_ip1half_n*((grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0])-(dU_ijk_nm1half-dU0_i_nm1half))
          /grid.dLocalGridOld[grid.nDM][nICen][0][0]*2.0-0.3333333333333333*dDivU_ip1jk_n);
        
        //cal Tau_rr_ijk_n
        dTau_rr_ijk_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]*(4.0*parameters.dPi
          *dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][nICen][0][0]
          *((grid.dLocalGridOld[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU0][i][0][0])
          -(grid.dLocalGridOld[grid.nU][i-1][j][k]-grid.dLocalGridOld[grid.nU0][i-1][0][0]))
          /grid.dLocalGridOld[grid.nDM][nICen][0][0]-0.3333333333333333*dDivU_ijk_n);
        
        //calculate dTau_rt_ip1halfjp1halfk_n
        dTau_rt_ip1halfjp1halfk_n=dEddyVisc_ip1halfjp1halfk_n*(4.0*parameters.dPi*dRCu_ip1half_n
          *dRhoAve_ip1half_n*(dV_R_ip1jp1halfk_n-dV_R_ijp1halfk_n)/dDM_ip1half
          +((grid.dLocalGridOld[grid.nU][i][j+1][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0])-(grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0]))/(dDTheta_jp1half
          *grid.dLocalGridOld[grid.nR][i][0][0]));
        
        //calculate dTau_rt_ip1halfjm1halfk_n
        dTau_rt_ip1halfjm1halfk_n=dEddyVisc_ip1halfjm1halfk_n*(4.0*parameters.dPi*dRCu_ip1half_n
          *dRhoAve_ip1half_n*(dV_R_ip1jm1halfk_n-dV_R_ijm1halfk_n)/dDM_ip1half
          +((grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0])-(grid.dLocalGridOld[grid.nU][i][j-1][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0]))/(grid.dLocalGridOld[grid.nR][i][0][0]
          *dDTheta_jm1half));
        
        //calculate dTau_rp_ip1halfjkp1half_n
        dTau_rp_ip1halfjkp1half_n=dEddyVisc_ip1halfjkp1half_n*(4.0*parameters.dPi*dRCu_ip1half_n
          *dRhoAve_ip1half_n*(dW_R_ip1jkp1half_n-dW_R_ijkp1half_n)/dDM_ip1half
          +((grid.dLocalGridOld[grid.nU][i][j][k+1]-grid.dLocalGridOld[grid.nU0][i][0][0])
          -(grid.dLocalGridOld[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU0][i][0][0]))
          /(grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *dDPhi_kp1half));
          
        //calculate dTau_rp_im1halfjkm1half_n
        dTau_rp_ip1halfjkm1half_n=dEddyVisc_ip1halfjkm1half_n*(4.0*parameters.dPi*dRCu_ip1half_n
          *dRhoAve_ip1half_n*(dW_R_ip1jkm1half_n-dW_R_ijkm1half_n)/dDM_ip1half
          +((grid.dLocalGridOld[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU0][i][0][0])
          -(grid.dLocalGridOld[grid.nU][i][j][k-1]-grid.dLocalGridOld[grid.nU0][i][0][0]))
          /(grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *dDPhi_km1half));
        
        //cal dTA1
        dTA1=(dTau_rr_ip1jk_n-dTau_rr_ijk_n)/(dDM_ip1half*dRho_ip1halfjk_n);
        
        //cal dTS1
        dTS1=dEddyVisc_ip1halfjk_n/dRhoR_ip1halfjk_n*(4.0
          *((dU_ip1jk_nm1half-dU0_ip1_nm1half)-(dU_ijk_nm1half-dU0_i_nm1half))/dDM_ip1half
          +grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]
          *(dV_R_ip1jk_n-dV_R_ijk_n)/dDM_ip1half);
        
        //calculate dTA2
        dTA2=(dTau_rt_ip1halfjp1halfk_n-dTau_rt_ip1halfjm1halfk_n)
          /(grid.dLocalGridOld[grid.nDTheta][0][j][0]*dRhoR_ip1halfjk_n);
        
        //calculate dTS2
        dTS2=(2.0*(dV_ip1halfjp1halfk_nm1half-dV_ip1halfjm1halfk_nm1half)
          -grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]*((dU_ip1halfjp1halfk_nm1half
          -grid.dLocalGridOld[grid.nU0][i][0][0])-(dU_ip1halfjm1halfk_nm1half
          -grid.dLocalGridOld[grid.nU0][i][0][0])))/(grid.dLocalGridOld[grid.nR][i][0][0]
          *grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //cal dTS4
        dTS4=(4.0*(grid.dLocalGridOld[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU0][i][0][0])
          +2.0*grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]*dV_ip1halfjk_nm1half)
          /grid.dLocalGridOld[grid.nR][i][0][0];
        
        dEddyViscosityTerms=-4.0*parameters.dPi*dRhoAve_ip1half_n*dRSq_ip1half_n*(dTA1+dTS1)-dTA2
          -dTA3+dEddyVisc_ip1halfjk_n/dRhoR_ip1halfjk_n*(dTS2+dTS3+dTS4);
        //dEddyViscosityTerms=0.0;
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nU][i][j][k]=grid.dLocalGridOld[grid.nU][i][j][k]
          -time.dDeltat_n*(4.0*parameters.dPi*dRhoAve_ip1half_n*dRSq_ip1half_n*(dA1+dS1)
          +dA2-dS2+dA3-dS3+dS4+dEddyViscosityTerms);
          
          
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
        
        //add M_r
        ssName.str("");
        ssName<<"M_r";
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,grid.dLocalGridOld[grid.nM][i][0][0]);
        
        //add A1
        ssName.str("");
        ssName<<"U_A1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-4.0*parameters.dPi*dRhoAve_ip1half_n*dRSq_ip1half_n*(dA1));
          
        //add S1
        ssName.str("");
        ssName<<"U_S1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-4.0*parameters.dPi*dRhoAve_ip1half_n*dRSq_ip1half_n*(dS1));
        
        //add A2
        ssName.str("");
        ssName<<"U_A2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-1.0*dA2);
          
        //add S2
        ssName.str("");
        ssName<<"U_S2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,dS2);
          
        //add A3
        ssName.str("");
        ssName<<"U_A3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-1.0*dA3);
          
        //add S3
        ssName.str("");
        ssName<<"U_S3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,dS3);
          
        //add S4
        ssName.str("");
        ssName<<"U_S4"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-1.0*dS4);
          
        //add dEddyViscosityTerms
        ssName.str("");
        ssName<<"U_EV"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
        ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-1.0*dEddyViscosityTerms);
        
        //add DuDt
        ssName.str("");
        ssName<<"U_DuDt"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,(grid.dLocalGridNew[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU][i][j][k])
          /time.dDeltat_n);
        #endif
      }
    }
  }
}
