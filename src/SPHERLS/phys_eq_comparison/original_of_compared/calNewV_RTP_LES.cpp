void calNewV_RTP_LES(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop){
  int i;
  int j;
  int k;
  int nIInt;
  int nJCen;
  int nKInt;
  double dR_i_n;
  double dR_ip1_n;
  double dR_im1_n;
  double dRCu_ip1half;
  double dRCu_im1half;
  double dU_ijp1halfk_nm1half;
  double dU_ijp1k_nm1half;
  double dU_im1halfjp1halfk_nm1half;
  double dU_im1jp1halfk_nm1half;
  double dU_ijk_nm1half;
  double dU0_i_nm1half;
  double dV_ip1halfjp1halfk_nm1half;
  double dV_ijp1k_nm1half;
  double dV_ijk_nm1half;
  double dV_im1halfjp1halfk_nm1half;
  double dV_ijp1halfkp1half_nm1half;
  double dV_ijp1halfkm1half_nm1half;
  double dW_ijp1halfk_nm1half;
  double dW_ijp1halfkp1half_nm1half;
  double dW_ijp1halfkm1half_nm1half;
  double dDTheta_jp1half;
  double dDPhi_kp1half;
  double dDPhi_km1half;
  double dRho_ijp1halfk_n;
  double dRhoAve_ip1half_n;
  double dRhoAve_im1half_n;
  double dP_ijp1k_n;
  double dP_ijk_n;
  double dDM_ip1half;
  double dDM_im1half;
  double dEddyVisc_ip1halfjp1halfk_n;
  double dEddyVisc_im1halfjp1halfk_n;
  double dEddyVisc_ijp1halfk_n;
  double dEddyVisc_ijp1halfkp1half_n;
  double dEddyVisc_ijp1halfkm1half_n;
  double dRSq_ip1half_n;
  double dRSq_im1half_n;
  double dRSq_i_n;
  double dRSqUmU0_ip1halfjp1k_n;
  double dRSqUmU0_im1halfjp1k_n;
  double dRSqUmU0_ip1halfjk_n;
  double dRSqUmU0_im1halfjk_n;
  double dV_R_ip1jp1halfk_n;
  double dV_R_ijp1halfk_n;
  double dV_R_im1jp1halfk_n;
  double dV_R_ip1halfjp1halfk_n;
  double dV_R_im1halfjp1halfk_n;
  double dW_SinTheta_ijp1kp1half_n;
  double dW_SinTheta_ijkp1half_n;
  double dW_SinTheta_ijp1km1half_n;
  double dW_SinTheta_ijkm1half_n;
  double dU_U0_Diff_ijp1halfk_nm1half;
  double dA1CenGrad;
  double dA1UpWindGrad;
  double dA1;
  double dS1;
  double dA2CenGrad;
  double dA2UpWindGrad;
  double dA2;
  double dS2;
  double dA3CenGrad;
  double dA3UpWindGrad;
  double dA3;
  double dS3;
  double dTau_rt_ip1halfjp1halfk_n;
  double dTau_rt_im1halfjp1halfk_n;
  double dDivU_ijp1k_n;
  double dDivU_ijk_n;
  double dTau_tt_ijp1k_n;
  double dTau_tt_ijk_n;
  double dTau_tp_ijp1halfkp1half_n;
  double dTau_tp_ijp1halfkm1half_n;
  double dTA1;
  double dTS1;
  double dTA2;
  double dTS2;
  double dTA3;
  double dTS3;
  double dTS4;
  double dEddyViscosityTerms;
  
  //calculate new v
  for(i=grid.nStartUpdateExplicit[grid.nV][0];i<grid.nEndUpdateExplicit[grid.nV][0];i++){
    
    //calculate i of interface quantities
    nIInt=i+grid.nCenIntOffset[0];
    
    //calculate quantities that only vary radially
    dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])
      *0.5;
    dR_ip1_n=(grid.dLocalGridOld[grid.nR][nIInt+1][0][0]+grid.dLocalGridOld[grid.nR][nIInt][0][0])
      *0.5;
    dR_im1_n=(grid.dLocalGridOld[grid.nR][nIInt-1][0][0]+grid.dLocalGridOld[grid.nR][nIInt-2][0][0])
      *0.5;
    dRSq_i_n=dR_i_n*dR_i_n;
    dRSq_ip1half_n=grid.dLocalGridOld[grid.nR][nIInt][0][0]
      *grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dRSq_im1half_n=grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      *grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dRCu_ip1half=dRSq_ip1half_n*grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dRCu_im1half=dRSq_im1half_n*grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dU0_i_nm1half=(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0])*0.5;
    dRhoAve_ip1half_n=(grid.dLocalGridOld[grid.nDenAve][i+1][0][0]
      +grid.dLocalGridOld[grid.nDenAve][i][0][0])*0.5;
    dRhoAve_im1half_n=(grid.dLocalGridOld[grid.nDenAve][i][0][0]
      +grid.dLocalGridOld[grid.nDenAve][i-1][0][0])*0.5;
    dDM_ip1half=(grid.dLocalGridOld[grid.nDM][i+1][0][0]+grid.dLocalGridOld[grid.nDM][i][0][0])*0.5;
    dDM_im1half=(grid.dLocalGridOld[grid.nDM][i-1][0][0]+grid.dLocalGridOld[grid.nDM][i][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nV][1];j<grid.nEndUpdateExplicit[grid.nV][1];j++){
      
      //calculate j of centered quantities
      nJCen=j-grid.nCenIntOffset[1];
      
      //calculate quantities that only vary with theta and or radius
      dDTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0]
        +grid.dLocalGridOld[grid.nDTheta][0][nJCen][0])*0.5;
      
      for(k=grid.nStartUpdateExplicit[grid.nV][2];k<grid.nEndUpdateExplicit[grid.nV][2];k++){
        
        //calculate k of interface quantities
        nKInt=k+grid.nCenIntOffset[2];
        
        dDPhi_kp1half=(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k])*0.5;
        dDPhi_km1half=(grid.dLocalGridOld[grid.nDPhi][0][0][k]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*0.5;
        
        //Calculate interpolated quantities
        dU_ijp1halfk_nm1half=0.25*(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]);
        dU_ijp1k_nm1half=(grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k])*0.5;
        dU_im1halfjp1halfk_nm1half=(grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k])*0.5;
        dU_im1jp1halfk_nm1half=0.25*(grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-2][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt-2][nJCen+1][k]);
        dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k])*0.5;
        dV_ip1halfjp1halfk_nm1half=0.5*(grid.dLocalGridOld[grid.nV][i+1][j][k]
          +grid.dLocalGridOld[grid.nV][i][j][k]);
        dV_im1halfjp1halfk_nm1half=0.5*(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i-1][j][k]);
        dV_ijp1k_nm1half=(grid.dLocalGridOld[grid.nV][i][j+1][k]
          +grid.dLocalGridOld[grid.nV][i][j][k])*0.5;
        dV_ijk_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i][j-1][k])*0.5;
        dV_ijp1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k+1]
          +grid.dLocalGridOld[grid.nV][i][j][k])*0.5;
        dV_ijp1halfkm1half_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i][j][k-1])*0.5;
        dW_ijp1halfk_nm1half=0.25*(grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt]
          +grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt-1]
          +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt]
          +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt-1]);
        dW_ijp1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt]
          +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt])*0.5;
        dW_ijp1halfkm1half_nm1half=(grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt-1]
          +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt-1])*0.5;
        dRho_ijp1halfk_n=(grid.dLocalGridOld[grid.nD][i][nJCen][k]
          +grid.dLocalGridOld[grid.nD][i][nJCen+1][k])*0.5;
        dP_ijp1k_n=grid.dLocalGridOld[grid.nP][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nQ0][i][nJCen+1][k]+grid.dLocalGridOld[grid.nQ1][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nQ2][i][nJCen+1][k];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][nJCen][k]+grid.dLocalGridOld[grid.nQ0][i][nJCen][k]
          +grid.dLocalGridOld[grid.nQ1][i][nJCen][k]+grid.dLocalGridOld[grid.nQ2][i][nJCen][k];
        dEddyVisc_ip1halfjp1halfk_n=(grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i+1][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i+1][nJCen+1][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k])*0.25;
        dEddyVisc_im1halfjp1halfk_n=(grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i-1][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i-1][nJCen+1][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k])*0.25;
        dEddyVisc_ijp1halfk_n=(grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k])*0.5;
        dEddyVisc_ijp1halfkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k])*0.25;
        dEddyVisc_ijp1halfkm1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k-1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k-1])*0.25;
        
        //calculate derived quantities
        dU_U0_Diff_ijp1halfk_nm1half=dU_ijp1halfk_nm1half-dU0_i_nm1half;
        dRSqUmU0_ip1halfjp1k_n=dRSq_ip1half_n*(grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0]);
        dRSqUmU0_im1halfjp1k_n=dRSq_im1half_n*(grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
        dRSqUmU0_ip1halfjk_n=dRSq_ip1half_n*(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0]);
        dRSqUmU0_im1halfjk_n=dRSq_im1half_n*(grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
        dU_U0_Diff_ijp1halfk_nm1half=dU_ijp1halfk_nm1half-dU0_i_nm1half;
        dV_R_ip1jp1halfk_n=grid.dLocalGridOld[grid.nV][i+1][j][k]/dR_ip1_n;
        dV_R_ijp1halfk_n=grid.dLocalGridOld[grid.nV][i][j][k]/dR_i_n;
        dV_R_im1jp1halfk_n=grid.dLocalGridOld[grid.nV][i-1][j][k]/dR_im1_n;
        dV_R_ip1halfjp1halfk_n=dV_ip1halfjp1halfk_nm1half/grid.dLocalGridOld[grid.nR][nIInt][0][0];
        dV_R_im1halfjp1halfk_n=dV_im1halfjp1halfk_nm1half
          /grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
        dW_SinTheta_ijp1kp1half_n=grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen+1][0];
        dW_SinTheta_ijkp1half_n=grid.dLocalGridOld[grid.nW][i][nJCen][nKInt]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen][0];
        dW_SinTheta_ijp1km1half_n=grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt-1]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen+1][0];
        dW_SinTheta_ijkm1half_n=grid.dLocalGridOld[grid.nW][i][nJCen][nKInt-1]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen][0];
        
        //calculate A1
        dA1CenGrad=(dV_ip1halfjp1halfk_nm1half-dV_im1halfjp1halfk_nm1half)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        dA1UpWindGrad=0.0;
        if(dU_U0_Diff_ijp1halfk_nm1half<0.0){//moving in a negative direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nV][i+1][j][k]
            -grid.dLocalGridOld[grid.nV][i][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i+1][0][0])*2.0;
        }
        else{//moving in a positive direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
            -grid.dLocalGridOld[grid.nV][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=dU_U0_Diff_ijp1halfk_nm1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate S1
        dS1=dU_ijp1halfk_nm1half*grid.dLocalGridOld[grid.nV][i][j][k]/dR_i_n;
        
        //calculate dA2
        dA2CenGrad=(dV_ijp1k_nm1half-dV_ijk_nm1half)/dDTheta_jp1half;
        dA2UpWindGrad=0.0;
        if(grid.dLocalGridOld[grid.nV][i][j][k]<0.0){//moning in a negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j+1][k]
            -grid.dLocalGridOld[grid.nV][i][j][k])/grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0];
        }
        else{//moving in a positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
            -grid.dLocalGridOld[grid.nV][i][j-1][k])/grid.dLocalGridOld[grid.nDTheta][0][nJCen][0];
        }
        dA2=grid.dLocalGridOld[grid.nV][i][j][k]/dR_i_n
          *((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA2CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
        
        //calculate S2
        dS2=(dP_ijp1k_n-dP_ijk_n)/(dDTheta_jp1half*dRho_ijp1halfk_n*dR_i_n);
        
        //calculate A3
        dA3CenGrad=(dV_ijp1halfkp1half_nm1half-dV_ijp1halfkm1half_nm1half)
          /grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dA3UpWindGrad=0.0;
        if(dW_ijp1halfk_nm1half<0.0){//moving in a negative direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k+1]
            -grid.dLocalGridOld[grid.nV][i][j][k])/(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k])*2.0;
        }
        else{//moving in a positive direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
            -grid.dLocalGridOld[grid.nV][i][j][k-1])/(grid.dLocalGridOld[grid.nDPhi][0][0][k]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*2.0;
        }
        dA3=dW_ijp1halfk_nm1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA3CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA3UpWindGrad)/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]);
        
        //calculate S3
        dS3=-1.0*dW_ijp1halfk_nm1half*dW_ijp1halfk_nm1half
          *grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]/dR_i_n;
        
        //calculate Tau_rt_ip1halfjp1halfk_n
        dTau_rt_ip1halfjp1halfk_n=dEddyVisc_ip1halfjp1halfk_n*(4.0*parameters.dPi*dRCu_ip1half
          *dRhoAve_ip1half_n*(dV_R_ip1jp1halfk_n-dV_R_ijp1halfk_n)/dDM_ip1half+1.0
          /grid.dLocalGridOld[grid.nR][nIInt][0][0]*((grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0])-(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0]))/dDTheta_jp1half);
        
        //calculate Tau_rt_im1halfjp1halfk_n
        dTau_rt_im1halfjp1halfk_n=dEddyVisc_im1halfjp1halfk_n*(4.0*parameters.dPi*dRCu_im1half
          *dRhoAve_im1half_n*(dV_R_ijp1halfk_n-dV_R_im1jp1halfk_n)/dDM_im1half+1.0
          /grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
          *((grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0])
          -(grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]))/dDTheta_jp1half);
        
        //calculate DivU_ijp1k_n
        dDivU_ijp1k_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dRSqUmU0_ip1halfjp1k_n-dRSqUmU0_im1halfjp1k_n)/grid.dLocalGridOld[grid.nDM][i][0][0]
          +((grid.dLocalGridOld[grid.nV][i][j+1][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j+1][0]
          -grid.dLocalGridOld[grid.nV][i][j][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0])
          /grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0]
          +(grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt]
          -grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt-1])
          /(grid.dLocalGridOld[grid.nDPhi][0][0][k]))
          /(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen+1][0]);
        
        //calculate DivU_ijk_n
        dDivU_ijk_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dRSqUmU0_ip1halfjk_n-dRSqUmU0_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0]
          +((grid.dLocalGridOld[grid.nV][i][j][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]
          -grid.dLocalGridOld[grid.nV][i][j-1][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j-1][0])
          /grid.dLocalGridOld[grid.nDTheta][0][nJCen][0]
          +(grid.dLocalGridOld[grid.nW][i][nJCen][nKInt]
          -grid.dLocalGridOld[grid.nW][i][nJCen][nKInt-1])
          /(grid.dLocalGridOld[grid.nDPhi][0][0][k]))
          /(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen][0]);
        
        //calculate Tau_tt_ijp1k_n
        dTau_tt_ijp1k_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k]
          *((grid.dLocalGridOld[grid.nV][i][j+1][k]-grid.dLocalGridOld[grid.nV][i][j][k])
          /(dR_i_n*grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0])
          +(dU_ijp1k_nm1half-dU0_i_nm1half)/dR_i_n-0.333333333333333*dDivU_ijp1k_n);
        
        //calculate Tau_tt_ijk_n
        dTau_tt_ijk_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k]
          *((grid.dLocalGridOld[grid.nV][i][j][k]-grid.dLocalGridOld[grid.nV][i][j-1][k])
          /(grid.dLocalGridOld[grid.nDTheta][0][nJCen][0]*dR_i_n)
          +(dU_ijk_nm1half-dU0_i_nm1half)/dR_i_n-0.333333333333333*dDivU_ijk_n);
        
        //calculate dTau_tp_ijp1halfkp1half_n
        dTau_tp_ijp1halfkp1half_n=dEddyVisc_ijp1halfkp1half_n
          *(grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]*(dW_SinTheta_ijp1kp1half_n
          -dW_SinTheta_ijkp1half_n)/(dR_i_n*dDTheta_jp1half)+(grid.dLocalGridOld[grid.nV][i][j][k+1]
          -grid.dLocalGridOld[grid.nV][i][j][k])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]*dDPhi_kp1half));
        
        //calculate dTau_tp_ijp1halfmp1half_n
        dTau_tp_ijp1halfkm1half_n=dEddyVisc_ijp1halfkm1half_n
          *(grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]*(dW_SinTheta_ijp1km1half_n
          -dW_SinTheta_ijkm1half_n)/(dR_i_n*dDTheta_jp1half)+(grid.dLocalGridOld[grid.nV][i][j][k]
          -grid.dLocalGridOld[grid.nV][i][j][k-1])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]*dDPhi_km1half));
        
        //calculate TA1
        dTA1=(dTau_rt_ip1halfjp1halfk_n-dTau_rt_im1halfjp1halfk_n)
          /(grid.dLocalGridOld[grid.nDM][i][0][0]*dRho_ijp1halfk_n);
        
        //calculate TS1
        dTS1=3.0*dEddyVisc_ijp1halfk_n*(dV_R_ip1halfjp1halfk_n-dV_R_im1halfjp1halfk_n)
          /(grid.dLocalGridOld[grid.nDM][i][0][0]*dRho_ijp1halfk_n);
        
        //calculate TA2
        dTA2=(dTau_tt_ijp1k_n-dTau_tt_ijk_n)/(dRho_ijp1halfk_n*dR_i_n*dDTheta_jp1half);
        
        //calculate TS2
        dTS2=(2.0*grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]*(dV_ijp1k_nm1half
          -dV_ijk_nm1half)+3.0*((dU_ijp1k_nm1half-dU0_i_nm1half)
          -(dU_ijk_nm1half-dU0_i_nm1half)))/(dR_i_n*dDTheta_jp1half);
        
        //calculate dTA3
        dTA3=(dTau_tp_ijp1halfkp1half_n-dTau_tp_ijp1halfkm1half_n)/(dRho_ijp1halfk_n*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //calculate dTS3
        dTS3=2.0*grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]*(dW_ijp1halfkp1half_nm1half
          -dW_ijp1halfkm1half_nm1half)/(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //calculate TS4
        dTS4=2.0*grid.dLocalGridOld[grid.nV][i][j][k]
          *grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]
          *grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]/dR_i_n;
        
        dEddyViscosityTerms=-4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dTA1+dTS1)-dTA2-dTA3-dEddyVisc_ijp1halfk_n/(dRho_ijp1halfk_n*dR_i_n)*(dTS2-dTS3-dTS4);
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nV][i][j][k]=grid.dLocalGridOld[grid.nV][i][j][k]
          -time.dDeltat_n*(4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dA1)+dS1+dA2+dS2+dA3+dS3+dEddyViscosityTerms);
          
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
        
        //add A1
        ssName.str("");
        ssName<<"V_A1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1));
        
        //add S1
        ssName.str("");
        ssName<<"V_S1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS1);
        
        //add A2
        ssName.str("");
        ssName<<"V_A2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dA2);
        
        //add S2
        ssName.str("");
        ssName<<"V_S2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS2);
        
        //add A3
        ssName.str("");
        ssName<<"V_A3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dA3);
        
        //add S3
        ssName.str("");
        ssName<<"V_S3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS3);
        
        //add EV
        ssName.str("");
        ssName<<"V_EV"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dEddyViscosityTerms);
        
        //add A3
        ssName.str("");
        ssName<<"V_DvDt"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,(grid.dLocalGridNew[grid.nV][i][j][k]-grid.dLocalGridOld[grid.nV][i][j][k])
          /time.dDeltat_n);
        #endif
        
      }
    }
  }
  
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nV][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nV][0][0];i++){
    
    //calculate i of interface quantities
    nIInt=i+grid.nCenIntOffset[0];
    
    //calculate quantities that only vary radially
    dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])
      *0.5;
    dR_ip1_n=grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dR_im1_n=(grid.dLocalGridOld[grid.nR][nIInt-1][0][0]+grid.dLocalGridOld[grid.nR][nIInt-2][0][0])
      *0.5;
    dRSq_i_n=dR_i_n*dR_i_n;
    dRSq_ip1half_n=grid.dLocalGridOld[grid.nR][nIInt][0][0]
      *grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dRSq_im1half_n=grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      *grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dRCu_ip1half=dRSq_ip1half_n*grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dRCu_im1half=dRSq_im1half_n*grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dU0_i_nm1half=(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0])*0.5;
    dRhoAve_ip1half_n=grid.dLocalGridOld[grid.nDenAve][i][0][0]*0.5;/**\BC Assuming density outside
      star is zero*/
    dRhoAve_im1half_n=(grid.dLocalGridOld[grid.nDenAve][i][0][0]
      +grid.dLocalGridOld[grid.nDenAve][i-1][0][0])*0.5;
    dDM_ip1half=grid.dLocalGridOld[grid.nDM][i][0][0]*0.5;
    dDM_im1half=(grid.dLocalGridOld[grid.nDM][i-1][0][0]+grid.dLocalGridOld[grid.nDM][i][0][0])*0.5;
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nV][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nV][0][1];j++){
      
      //calculate j of centered quantities
      nJCen=j-grid.nCenIntOffset[1];
      
      //calculate quantities that only vary with theta and or radius
      dDTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0]
        +grid.dLocalGridOld[grid.nDTheta][0][nJCen][0])*0.5;
      
      for(k=grid.nStartGhostUpdateExplicit[grid.nV][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nV][0][2];k++){
        
        //calculate k of interface quantities
        nKInt=k+grid.nCenIntOffset[2];
        dDPhi_kp1half=(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k])*0.5;
        dDPhi_km1half=(grid.dLocalGridOld[grid.nDPhi][0][0][k]
          +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*0.5;
        
        //Calculate interpolated quantities
        dU_ijp1halfk_nm1half=0.25*(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]);
        dU_ijp1k_nm1half=(grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k])*0.5;
        dU_im1halfjp1halfk_nm1half=(grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k])*0.5;
        dU_im1jp1halfk_nm1half=0.25*(grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-2][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt-2][nJCen+1][k]);
        dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k])*0.5;
        dV_ip1halfjp1halfk_nm1half=grid.dLocalGridOld[grid.nV][i][j][k];/**\BC Assuming theta 
          velocity is constant across surface.*/
        dV_im1halfjp1halfk_nm1half=0.5*(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i-1][j][k]);
        dV_ijp1k_nm1half=(grid.dLocalGridOld[grid.nV][i][j+1][k]
          +grid.dLocalGridOld[grid.nV][i][j][k])*0.5;
        dV_ijk_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i][j-1][k])*0.5;
        dV_ijp1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k+1]
          +grid.dLocalGridOld[grid.nV][i][j][k])*0.5;
        dV_ijp1halfkm1half_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i][j][k-1])*0.5;
        dW_ijp1halfk_nm1half=0.25*(grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt]
          +grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt-1]
          +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt]
          +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt-1]);
        dW_ijp1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt]
          +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt])*0.5;
        dW_ijp1halfkm1half_nm1half=(grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt-1]
          +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt-1])*0.5;
        dRho_ijp1halfk_n=(grid.dLocalGridOld[grid.nD][i][nJCen][k]
          +grid.dLocalGridOld[grid.nD][i][nJCen+1][k])*0.5;
        dP_ijp1k_n=grid.dLocalGridOld[grid.nP][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nQ0][i][nJCen+1][k]+grid.dLocalGridOld[grid.nQ1][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nQ2][i][nJCen+1][k];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][nJCen][k]+grid.dLocalGridOld[grid.nQ0][i][nJCen][k]
          +grid.dLocalGridOld[grid.nQ1][i][nJCen][k]+grid.dLocalGridOld[grid.nQ2][i][nJCen][k];
        dEddyVisc_ip1halfjp1halfk_n=(grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k])*0.25;/**\BC Assuming eddy viscosity is
          zero at surface.*/
        dEddyVisc_im1halfjp1halfk_n=(grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i-1][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i-1][nJCen+1][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k])*0.25;
        dEddyVisc_ijp1halfk_n=(grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k])*0.5;
        dEddyVisc_ijp1halfkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k])*0.25;
        dEddyVisc_ijp1halfkm1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k-1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k-1])*0.25;
        
        //calculate derived quantities
        dU_U0_Diff_ijp1halfk_nm1half=dU_ijp1halfk_nm1half-dU0_i_nm1half;
        dRSqUmU0_ip1halfjp1k_n=dRSq_ip1half_n*(grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0]);
        dRSqUmU0_im1halfjp1k_n=dRSq_im1half_n*(grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
        dRSqUmU0_ip1halfjk_n=dRSq_ip1half_n*(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0]);
        dRSqUmU0_im1halfjk_n=dRSq_im1half_n*(grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
        dV_R_ip1jp1halfk_n=grid.dLocalGridOld[grid.nV][i][j][k]/dR_ip1_n;
        dV_R_ijp1halfk_n=grid.dLocalGridOld[grid.nV][i][j][k]/dR_i_n;
        dV_R_im1jp1halfk_n=grid.dLocalGridOld[grid.nV][i-1][j][k]/dR_im1_n;
        dV_R_ip1halfjp1halfk_n=dV_ip1halfjp1halfk_nm1half/grid.dLocalGridOld[grid.nR][nIInt][0][0];
        dV_R_im1halfjp1halfk_n=dV_im1halfjp1halfk_nm1half
          /grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
        dW_SinTheta_ijp1kp1half_n=grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen+1][0];
        dW_SinTheta_ijkp1half_n=grid.dLocalGridOld[grid.nW][i][nJCen][nKInt]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen][0];
        dW_SinTheta_ijp1km1half_n=grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt-1]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen+1][0];
        dW_SinTheta_ijkm1half_n=grid.dLocalGridOld[grid.nW][i][nJCen][nKInt-1]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen][0];
        
        //calculate A1
        dA1CenGrad=(dV_ip1halfjp1halfk_nm1half-dV_im1halfjp1halfk_nm1half)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        dA1UpWindGrad=0.0;
        if(dU_U0_Diff_ijp1halfk_nm1half<0.0){//moving in a negative direction
          dA1UpWindGrad=dA1CenGrad;
        }
        else{//moving in a positive direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
            -grid.dLocalGridOld[grid.nV][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=dU_U0_Diff_ijp1halfk_nm1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate S1
        dS1=dU_ijp1halfk_nm1half*grid.dLocalGridOld[grid.nV][i][j][k]/dR_i_n;
        
        //calculate dA2
        dA2CenGrad=(dV_ijp1k_nm1half-dV_ijk_nm1half)/dDTheta_jp1half;
        dA2UpWindGrad=0.0;
        if(grid.dLocalGridOld[grid.nV][i][j][k]<0.0){//moving in a negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j+1][k]
            -grid.dLocalGridOld[grid.nV][i][j][k])/grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0];
        }
        else{//moving in a positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
            -grid.dLocalGridOld[grid.nV][i][j-1][k])/grid.dLocalGridOld[grid.nDTheta][0][nJCen][0];
        }
        dA2=grid.dLocalGridOld[grid.nV][i][j][k]/dR_i_n
          *((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA2CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
        
        //calculate S2
        dS2=(dP_ijp1k_n-dP_ijk_n)/(dDTheta_jp1half*dRho_ijp1halfk_n*dR_i_n);
        
        //calculate A3
        dA3CenGrad=(dV_ijp1halfkp1half_nm1half-dV_ijp1halfkm1half_nm1half)
          /grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dA3UpWindGrad=0.0;
        if(dW_ijp1halfk_nm1half<0.0){//moving in a negative direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k+1]
            -grid.dLocalGridOld[grid.nV][i][j][k])/(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k])*2.0;
        }
        else{//moving in a positive direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
            -grid.dLocalGridOld[grid.nV][i][j][k-1])/(grid.dLocalGridOld[grid.nDPhi][0][0][k]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*2.0;
        }
        dA3=dW_ijp1halfk_nm1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA3CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA3UpWindGrad)/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]);
        
        //calculate S3
        dS3=-1.0*dW_ijp1halfk_nm1half*dW_ijp1halfk_nm1half
          *grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]/dR_i_n;
        
        //calculate Tau_rt_ip1halfjp1halfk_n
        dTau_rt_ip1halfjp1halfk_n=dEddyVisc_ip1halfjp1halfk_n*(4.0*parameters.dPi*dRCu_ip1half
          *dRhoAve_ip1half_n*(dV_R_ip1jp1halfk_n-dV_R_ijp1halfk_n)/dDM_ip1half+1.0
          /grid.dLocalGridOld[grid.nR][nIInt][0][0]*((grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0])-(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0]))/dDTheta_jp1half);
        
        //calculate Tau_rt_im1halfjp1halfk_n
        dTau_rt_im1halfjp1halfk_n=dEddyVisc_im1halfjp1halfk_n*(4.0*parameters.dPi*dRCu_im1half
          *dRhoAve_im1half_n*(dV_R_ijp1halfk_n-dV_R_im1jp1halfk_n)/dDM_im1half+1.0
          /grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
          *((grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0])
          -(grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]))/dDTheta_jp1half);
        
        //calculate DivU_ijp1k_n
        dDivU_ijp1k_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dRSqUmU0_ip1halfjp1k_n-dRSqUmU0_im1halfjp1k_n)/grid.dLocalGridOld[grid.nDM][i][0][0]
          +((grid.dLocalGridOld[grid.nV][i][j+1][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j+1][0]
          -grid.dLocalGridOld[grid.nV][i][j][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0])
          /grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0]
          +(grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt]
          -grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt-1])
          /(grid.dLocalGridOld[grid.nDPhi][0][0][k]))
          /(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen+1][0]);
        
        //calculate DivU_ijk_n
        dDivU_ijk_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dRSqUmU0_ip1halfjk_n-dRSqUmU0_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0]
          +((grid.dLocalGridOld[grid.nV][i][j][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]
          -grid.dLocalGridOld[grid.nV][i][j-1][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j-1][0])
          /grid.dLocalGridOld[grid.nDTheta][0][nJCen][0]
          +(grid.dLocalGridOld[grid.nW][i][nJCen][nKInt]
          -grid.dLocalGridOld[grid.nW][i][nJCen][nKInt-1])
          /(grid.dLocalGridOld[grid.nDPhi][0][0][k]))
          /(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen][0]);
        
        //calculate Tau_tt_ijp1k_n
        dTau_tt_ijp1k_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k]
          *((grid.dLocalGridOld[grid.nV][i][j+1][k]-grid.dLocalGridOld[grid.nV][i][j][k])
          /(dR_i_n*grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0])
          +(dU_ijp1k_nm1half-dU0_i_nm1half)/dR_i_n-0.333333333333333*dDivU_ijp1k_n);
        
        //calculate Tau_tt_ijk_n
        dTau_tt_ijk_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k]
          *((grid.dLocalGridOld[grid.nV][i][j][k]-grid.dLocalGridOld[grid.nV][i][j-1][k])
          /(grid.dLocalGridOld[grid.nDTheta][0][nJCen][0]*dR_i_n)
          +(dU_ijk_nm1half-dU0_i_nm1half)/dR_i_n-0.333333333333333*dDivU_ijk_n);
        
        //calculate dTau_tp_ijp1halfkp1half_n
        dTau_tp_ijp1halfkp1half_n=dEddyVisc_ijp1halfkp1half_n
          *(grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]*(dW_SinTheta_ijp1kp1half_n
          -dW_SinTheta_ijkp1half_n)/(dR_i_n*dDTheta_jp1half)+(grid.dLocalGridOld[grid.nV][i][j][k+1]
          -grid.dLocalGridOld[grid.nV][i][j][k])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]*dDPhi_kp1half));
        
        //calculate dTau_tp_ijp1halfmp1half_n
        dTau_tp_ijp1halfkm1half_n=dEddyVisc_ijp1halfkm1half_n
          *(grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]*(dW_SinTheta_ijp1km1half_n
          -dW_SinTheta_ijkm1half_n)/(dR_i_n*dDTheta_jp1half)+(grid.dLocalGridOld[grid.nV][i][j][k]
          -grid.dLocalGridOld[grid.nV][i][j][k-1])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]*dDPhi_km1half));
        
        //calculate TA1
        dTA1=(dTau_rt_ip1halfjp1halfk_n-dTau_rt_im1halfjp1halfk_n)
          /(grid.dLocalGridOld[grid.nDM][i][0][0]*dRho_ijp1halfk_n);
        
        //calculate TS1
        dTS1=3.0*dEddyVisc_ijp1halfk_n*(dV_R_ip1halfjp1halfk_n-dV_R_im1halfjp1halfk_n)
          /(grid.dLocalGridOld[grid.nDM][i][0][0]*dRho_ijp1halfk_n);
        
        //calculate TA2
        dTA2=(dTau_tt_ijp1k_n-dTau_tt_ijk_n)/(dRho_ijp1halfk_n*dR_i_n*dDTheta_jp1half);
        
        //calculate TS2
        dTS2=(2.0*grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]*(dV_ijp1k_nm1half
          -dV_ijk_nm1half)+3.0*((dU_ijp1k_nm1half-dU0_i_nm1half)
          -(dU_ijk_nm1half-dU0_i_nm1half)))/(dR_i_n*dDTheta_jp1half);
        
        //calculate dTA3
        dTA3=(dTau_tp_ijp1halfkp1half_n-dTau_tp_ijp1halfkm1half_n)/(dRho_ijp1halfk_n*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //calculate dTS3
        dTS3=2.0*grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]*(dW_ijp1halfkp1half_nm1half
          -dW_ijp1halfkm1half_nm1half)/(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //calculate TS4
        dTS4=2.0*grid.dLocalGridOld[grid.nV][i][j][k]
          *grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]
          *grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]/dR_i_n;
        
        dEddyViscosityTerms=-4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dTA1+dTS1)-dTA2-dTA3-dEddyVisc_ijp1halfk_n/(dRho_ijp1halfk_n*dR_i_n)*(dTS2-dTS3-dTS4);
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nV][i][j][k]=grid.dLocalGridOld[grid.nV][i][j][k]
          -time.dDeltat_n*(4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dA1)+dS1+dA2+dS2+dA3+dS3+dEddyViscosityTerms);
          
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
        
        //add A1
        ssName.str("");
        ssName<<"V_A1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1));
        
        //add S1
        ssName.str("");
        ssName<<"V_S1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS1);
        
        //add A2
        ssName.str("");
        ssName<<"V_A2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dA2);
        
        //add S2
        ssName.str("");
        ssName<<"V_S2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS2);
        
        //add A3
        ssName.str("");
        ssName<<"V_A3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dA3);
        
        //add S3
        ssName.str("");
        ssName<<"V_S3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS3);
        
        //add EV
        ssName.str("");
        ssName<<"V_EV"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dEddyViscosityTerms);
        
        //add A3
        ssName.str("");
        ssName<<"V_DvDt"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,(grid.dLocalGridNew[grid.nV][i][j][k]-grid.dLocalGridOld[grid.nV][i][j][k])
          /time.dDeltat_n);
        #endif
      }
    }
  }
}
