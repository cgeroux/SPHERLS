void calNewW_RTP_LES(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop){
  int i;
  int j;
  int k;
  int nIInt;
  int nJInt;
  int nKCen;
  double dR_i_n;
  double dR_ip1_n;
  double dR_im1_n;
  double dRSq_i_n;
  double dRSq_ip1half_n;
  double dRSq_im1half_n;
  double dR3_ip1half_n;
  double dR3_im1half_n;
  double dRhoAve_ip1half_n;
  double dRhoAve_im1half_n;
  double dDM_ip1half;
  double dDM_im1half;
  double dDTheta_jp1half;
  double dDTheta_jm1half;
  double dDPhi_kp1half;
  double dDPhi_km1half;
  double dU0_i_nm1half;
  double dU_ijkp1half_nm1half;
  double dU_ijkp1_nm1half;
  double dU_ijk_nm1half;
  double dV_ijk_nm1half;
  double dV_ijkp1_nm1half;
  double dV_ijkp1half_nm1half;
  double dV_ijm1halfkp1half_nm1half;
  double dV_ijm1halfkm1half_nm1half;
  double dW_ijkp1half_nm1half;
  double dW_ijp1halfkp1half_nm1half;
  double dW_ijm1halfkp1half_nm1half;
  double dW_ip1halfjkp1half_nm1half;
  double dW_im1halfjkp1half_nm1half;
  double dW_ijkp1_nm1half;
  double dW_ijk_nm1half;
  double dRho_ijkp1half_n;
  double dP_ijkp1_n;
  double dP_ijk_n;
  double dEddyVisc_ip1halfjkp1half_n;
  double dEddyVisc_im1halfjkp1half_n;
  double dEddyVisc_ijp1halfkp1half_n;
  double dEddyVisc_ijm1halfkp1half_n;
  double dEddyVisc_ijkp1half_n;
  double dUmU0_ijkp1half_nm1half;
  double d1_rhoDM_ijkp1half_n;
  double dRSq_UmU0_ip1halfjkp1_n;
  double dRSq_UmU0_im1halfjkp1_n;
  double dRSq_UmU0_ip1halfjk_n;
  double dRSq_UmU0_im1halfjk_n;
  double dV_SinTheta_ijp1halfkp1_n;
  double dV_SinTheta_ijm1halfkp1_n;
  double dV_SinTheta_ijp1halfk_n;
  double dV_SinTheta_ijm1halfk_n;
  double dW_R_ip1jkp1half_n;
  double dW_R_im1jkp1half_n;
  double dW_R_ijkp1half_n;
  double dW_R_ip1halfjkp1half_n;
  double dW_R_im1halfjkp1half_n;
  double dW_SinTheta_ijp1kp1half_n;
  double dW_SinTheta_ijm1kp1half_n;
  double dW_SinTheta_ijkp1half_n;
  double dW_SinTheta_ijp1halfkp1half_n;
  double dW_SinTheta_ijm1halfkp1half_n;
  double dRRho_ijkp1half_n;
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
  double dDivU_ijkp1_n;
  double dDivU_ijk_n;
  double dTau_rp_ip1halfjkp1half_n;
  double dTau_rp_im1halfjkp1half_n;
  double dTau_tp_ijp1halfkp1half_n;
  double dTau_tp_ijm1halfkp1half_n;
  double dTau_pp_ijkp1_n;
  double dTau_pp_ijk_n;
  double dTA1;
  double dTS1;
  double dTA2;
  double dTS2;
  double dTA3;
  double dTS3;
  double dEddyViscosityTerms;
  
  //calculate new w
  for(i=grid.nStartUpdateExplicit[grid.nW][0];i<grid.nEndUpdateExplicit[grid.nW][0];i++){
    
    //calculate j of interface quantities
    nIInt=i+grid.nCenIntOffset[0];
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
    dR3_ip1half_n=dRSq_ip1half_n*grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dR3_im1half_n=dRSq_im1half_n*grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dU0_i_nm1half=(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0])*0.5;
    dRhoAve_ip1half_n=(grid.dLocalGridOld[grid.nDenAve][i+1][0][0]
      +grid.dLocalGridOld[grid.nDenAve][i][0][0])*0.5;
    dRhoAve_im1half_n=(grid.dLocalGridOld[grid.nDenAve][i][0][0]
      +grid.dLocalGridOld[grid.nDenAve][i-1][0][0])*0.5;
    dDM_ip1half=(grid.dLocalGridOld[grid.nDM][i+1][0][0]+grid.dLocalGridOld[grid.nDM][i][0][0])*0.5;
    dDM_im1half=(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i-1][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nW][1];j<grid.nEndUpdateExplicit[grid.nW][1];j++){
      
      //calculate j of centered quantities
      nJInt=j+grid.nCenIntOffset[1];
      dDTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
        +grid.dLocalGridOld[grid.nDTheta][0][j][0])*0.5;
      dDTheta_jm1half=(grid.dLocalGridOld[grid.nDTheta][0][j][0]
        +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*0.5;
      
      for(k=grid.nStartUpdateExplicit[grid.nW][2];k<grid.nEndUpdateExplicit[grid.nW][2];k++){
        
        //calculate k of interface quantities
        nKCen=k-grid.nCenIntOffset[2];
        
        //Calculate interpolated quantities
        dDPhi_kp1half=(grid.dLocalGridOld[grid.nDPhi][0][0][nKCen]
          +grid.dLocalGridOld[grid.nDPhi][0][0][nKCen+1])*0.5;
        dDPhi_km1half=(grid.dLocalGridOld[grid.nDPhi][0][0][nKCen]
          +grid.dLocalGridOld[grid.nDPhi][0][0][nKCen-1])*0.5;
        dU_ijkp1half_nm1half=(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen+1]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen+1]
          +grid.dLocalGridOld[grid.nU][nIInt][j][nKCen]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen])*0.25;
        dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen])*0.5;
        dU_ijkp1_nm1half=(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen+1]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen+1])*0.5;
        dV_ijk_nm1half=(grid.dLocalGridOld[grid.nV][i][nJInt][nKCen]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen])*0.5;
        dV_ijkp1_nm1half=(grid.dLocalGridOld[grid.nV][i][nJInt][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen+1])*0.5;
        dV_ijkp1half_nm1half=(grid.dLocalGridOld[grid.nV][i][nJInt][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt][nKCen]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen])*0.25;
        dV_ijm1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen])*0.5;
        dV_ijm1halfkm1half_nm1half=(grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen-1])*0.5;
        dW_ijkp1half_nm1half=grid.dLocalGridOld[grid.nW][i][j][k];
        dW_ijp1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nW][i][j+1][k]
          +grid.dLocalGridOld[grid.nW][i][j][k])*0.5;
        dW_ijm1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nW][i][j-1][k]
          +grid.dLocalGridOld[grid.nW][i][j][k])*0.5;
        dW_ip1halfjkp1half_nm1half=(grid.dLocalGridOld[grid.nW][i+1][j][k]
          +grid.dLocalGridOld[grid.nW][i][j][k])*0.5;
        dW_im1halfjkp1half_nm1half=(grid.dLocalGridOld[grid.nW][i-1][j][k]
          +grid.dLocalGridOld[grid.nW][i][j][k])*0.5;
        dW_ijkp1_nm1half=(grid.dLocalGridOld[grid.nW][i][j][k]
          +grid.dLocalGridOld[grid.nW][i][j][k+1])*0.5;
        dW_ijk_nm1half=(grid.dLocalGridOld[grid.nW][i][j][k]
          +grid.dLocalGridOld[grid.nW][i][j][k-1])*0.5;
        dRho_ijkp1half_n=(grid.dLocalGridOld[grid.nD][i][j][nKCen]
          +grid.dLocalGridOld[grid.nD][i][j][nKCen+1])*0.5;
        dP_ijkp1_n=grid.dLocalGridOld[grid.nP][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nQ0][i][j][nKCen+1]+grid.dLocalGridOld[grid.nQ1][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nQ2][i][j][nKCen+1];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][nKCen]+grid.dLocalGridOld[grid.nQ0][i][j][nKCen]
          +grid.dLocalGridOld[grid.nQ1][i][j][nKCen]
          +grid.dLocalGridOld[grid.nQ2][i][j][nKCen];
        dEddyVisc_ip1halfjkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i+1][j][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i+1][j][nKCen]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen])*0.25;
        dEddyVisc_im1halfjkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen]
          +grid.dLocalGridOld[grid.nEddyVisc][i-1][j][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i-1][j][nKCen])*0.25;
        dEddyVisc_ijp1halfkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i][j+1][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j+1][nKCen]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen])*0.25;
        dEddyVisc_ijm1halfkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j-1][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j-1][nKCen])*0.25;
        dEddyVisc_ijkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen])*0.5;
        
        //calculate derived quantities
        dUmU0_ijkp1half_nm1half=dU_ijkp1half_nm1half-dU0_i_nm1half;
        d1_rhoDM_ijkp1half_n=1.0/(dRho_ijkp1half_n*grid.dLocalGridOld[grid.nDM][i][0][0]);
        dRRho_ijkp1half_n=dR_i_n*dRho_ijkp1half_n;
        dRSq_UmU0_ip1halfjkp1_n=dRSq_ip1half_n*(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen+1]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0]);
        dRSq_UmU0_im1halfjkp1_n=dRSq_im1half_n*(grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen+1]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
        dRSq_UmU0_ip1halfjk_n=dRSq_ip1half_n*(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0]);
        dRSq_UmU0_im1halfjk_n=dRSq_im1half_n*(grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
        dV_SinTheta_ijp1halfkp1_n=grid.dLocalGridOld[grid.nV][i][nJInt][nKCen+1]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0];
        dV_SinTheta_ijm1halfkp1_n=grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen+1]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0];
        dV_SinTheta_ijp1halfk_n=grid.dLocalGridOld[grid.nV][i][nJInt][nKCen]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0];
        dV_SinTheta_ijm1halfk_n=grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0];
        dW_R_ip1jkp1half_n=grid.dLocalGridOld[grid.nW][i+1][j][k]/dR_ip1_n;
        dW_R_im1jkp1half_n=grid.dLocalGridOld[grid.nW][i-1][j][k]/dR_im1_n;
        dW_R_ijkp1half_n=grid.dLocalGridOld[grid.nW][i][j][k]/dR_i_n;
        dW_R_ip1halfjkp1half_n=dW_ip1halfjkp1half_nm1half/grid.dLocalGridOld[grid.nR][nIInt][0][0];
        dW_R_im1halfjkp1half_n=dW_im1halfjkp1half_nm1half
          /grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
        dW_SinTheta_ijp1kp1half_n=grid.dLocalGridOld[grid.nW][i][j+1][k]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][j+1][0];
        dW_SinTheta_ijm1kp1half_n=grid.dLocalGridOld[grid.nW][i][j-1][k]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][j-1][0];
        dW_SinTheta_ijkp1half_n=grid.dLocalGridOld[grid.nW][i][j][k]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0];
        dW_SinTheta_ijp1halfkp1half_n=dW_ijp1halfkp1half_nm1half
          /grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0];
        dW_SinTheta_ijm1halfkp1half_n=dW_ijm1halfkp1half_nm1half
          /grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0];
        
        //calculate A1
        dA1CenGrad=(dW_ip1halfjkp1half_nm1half-dW_im1halfjkp1half_nm1half)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        if(dUmU0_ijkp1half_nm1half<0.0){//moving in a negative direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nW][i+1][j][k]
            -grid.dLocalGridOld[grid.nW][i][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i+1][0][0])*2.0;
        }
        else{//moving in a positive direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j][k]
            -grid.dLocalGridOld[grid.nW][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=dUmU0_ijkp1half_nm1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate S1
        dS1=dU_ijkp1half_nm1half*dW_ijkp1half_nm1half/dR_i_n;
        
        //calculate dA2
        dA2CenGrad=(dW_ijp1halfkp1half_nm1half-dW_ijm1halfkp1half_nm1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        if(dV_ijkp1half_nm1half<0.0){//moving in a negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j+1][k]
            -grid.dLocalGridOld[grid.nW][i][j][k])/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        }
        else{//moving in a positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j][k]
            -grid.dLocalGridOld[grid.nW][i][j-1][k])/(grid.dLocalGridOld[grid.nDTheta][0][j-1][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        }
        dA2=dV_ijkp1half_nm1half/dR_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA2CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
        
        //calculate S2
        dS2=dV_ijkp1half_nm1half*grid.dLocalGridOld[grid.nW][i][j][k]
          *grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]/dR_i_n;
        
        //calculate A3
        dA3CenGrad=(dW_ijkp1_nm1half-dW_ijk_nm1half)/dDPhi_kp1half;
        if(dW_ijkp1half_nm1half<0.0){//moving in a negative direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j][k+1]
            -grid.dLocalGridOld[grid.nW][i][j][k])/grid.dLocalGridOld[grid.nDPhi][0][0][nKCen+1];
        }
        else{//moving in a positive direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j][k]
            -grid.dLocalGridOld[grid.nW][i][j][k-1])/grid.dLocalGridOld[grid.nDPhi][0][0][nKCen];
        }
        dA3=dW_ijkp1half_nm1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA3CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA3UpWindGrad)/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
        
        //calculate S3
        dS3=(dP_ijkp1_n-dP_ijk_n)/(dRho_ijkp1half_n*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*dDPhi_kp1half);
        
        //calculate dDivU_ijkp1_n
        dDivU_ijkp1_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dRSq_UmU0_ip1halfjkp1_n-dRSq_UmU0_im1halfjkp1_n)/grid.dLocalGridOld[grid.nDM][i][0][0]
          +(dV_SinTheta_ijp1halfkp1_n-dV_SinTheta_ijm1halfkp1_n)/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDTheta][0][j][0])
          +(grid.dLocalGridOld[grid.nW][i][j][k+1]-grid.dLocalGridOld[grid.nW][i][j][k])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][nKCen+1]);
        
        //calculate dDivU_ijk_n
        dDivU_ijk_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dRSq_UmU0_ip1halfjk_n-dRSq_UmU0_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0]
          +(dV_SinTheta_ijp1halfk_n-dV_SinTheta_ijm1halfk_n)/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDTheta][0][j][0])
          +(grid.dLocalGridOld[grid.nW][i][j][k]-grid.dLocalGridOld[grid.nW][i][j][k-1])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][nKCen]);
        
        //calculate dTau_rp_ip1halfjkp1half_n
        dTau_rp_ip1halfjkp1half_n=dEddyVisc_ip1halfjkp1half_n*(4.0*parameters.dPi*dR3_ip1half_n
          *dRhoAve_ip1half_n*(dW_R_ip1jkp1half_n-dW_R_ijkp1half_n)/dDM_ip1half
          +(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen+1]
          -grid.dLocalGridOld[grid.nU][nIInt][j][nKCen])
          /(dDPhi_kp1half*grid.dLocalGridOld[grid.nR][nIInt][0][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]));
        
        //calculate dTau_rp_im1halfjkp1half_n
        dTau_rp_im1halfjkp1half_n=dEddyVisc_im1halfjkp1half_n*(4.0*parameters.dPi*dR3_im1half_n
          *dRhoAve_im1half_n*(dW_R_ijkp1half_n-dW_R_im1jkp1half_n)/dDM_im1half
          +(grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen+1]
          -grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen])
          /(dDPhi_kp1half*grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]));
        
        //calculate dTau_tp_ijp1halfkp1half_n
        dTau_tp_ijp1halfkp1half_n=dEddyVisc_ijp1halfkp1half_n
          *(grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]*(dW_SinTheta_ijp1kp1half_n
          -dW_SinTheta_ijkp1half_n)/(dR_i_n*dDTheta_jp1half)
          +(grid.dLocalGridOld[grid.nV][i][nJInt][nKCen+1]
          -grid.dLocalGridOld[grid.nV][i][nJInt][nKCen])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]*dDPhi_kp1half));
        
        //calculate dTau_tp_ijm1halfkp1half_n
        dTau_tp_ijm1halfkp1half_n=dEddyVisc_ijm1halfkp1half_n
          *(grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]*(dW_SinTheta_ijkp1half_n
          -dW_SinTheta_ijm1kp1half_n)/(dR_i_n*dDTheta_jm1half)
          +(grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen+1]
          -grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]*dDPhi_kp1half));
        
        //calculate dTau_pp_ijkp1half_n
        dTau_pp_ijkp1_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen+1]
          *((grid.dLocalGridOld[grid.nW][i][j][k+1]-grid.dLocalGridOld[grid.nW][i][j][k])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][nKCen+1])+(dU_ijkp1_nm1half-dU0_i_nm1half)
          /dR_i_n+dV_ijkp1_nm1half*grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]/dR_i_n
          -0.333333333333333*dDivU_ijkp1_n);
        
        //calculate dTau_pp_ijk_n
        dTau_pp_ijk_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen]
          *((grid.dLocalGridOld[grid.nW][i][j][k]-grid.dLocalGridOld[grid.nW][i][j][k-1])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][nKCen])+(dU_ijk_nm1half-dU0_i_nm1half)
          /dR_i_n+dV_ijk_nm1half*grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]/dR_i_n
          -0.333333333333333*dDivU_ijk_n);
        
        //calculate dTA1
        dTA1=(dTau_rp_ip1halfjkp1half_n-dTau_rp_im1halfjkp1half_n)*d1_rhoDM_ijkp1half_n;
        
        //calculate dTS1
        dTS1=3.0*dEddyVisc_ijkp1half_n*(dW_R_ip1halfjkp1half_n-dW_R_im1halfjkp1half_n)
          *d1_rhoDM_ijkp1half_n;
        
        //calculate dTA2
        dTA2=(dTau_tp_ijp1halfkp1half_n-dTau_tp_ijm1halfkp1half_n)/(dRRho_ijkp1half_n
          *grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //calculate dTS2
        dTS2=2.0*grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*(dW_SinTheta_ijp1halfkp1half_n
          -dW_SinTheta_ijm1halfkp1half_n)/(dR_i_n*grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //calculate dTA3
        dTA3=(dTau_pp_ijkp1_n-dTau_pp_ijk_n)/(dRRho_ijkp1half_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *dDPhi_kp1half);
        //dTA3=0.0;
        
        //calculate dTS3
        dTS3=(3.0*(dU_ijkp1_nm1half-dU_ijk_nm1half)+2.0
          *grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]*(dV_ijkp1_nm1half-dV_ijk_nm1half))/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*dDPhi_kp1half);
        
        dEddyViscosityTerms=-4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dTA1+dTS1)-dTA2-dTA3-dEddyVisc_ijkp1half_n/(dRho_ijkp1half_n*dR_i_n)*(dTS2+dTS3);
        //dEddyViscosityTerms=0.0;
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nW][i][j][k]=grid.dLocalGridOld[grid.nW][i][j][k]-time.dDeltat_n
          *(4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1)
          +dS1+dA2+dS2+dA3+dS3+dEddyViscosityTerms);
          
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
        ssName<<"W_A1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1));
        
        //add S1
        ssName.str("");
        ssName<<"W_S1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS1);
        
        //add A2
        ssName.str("");
        ssName<<"W_A2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dA2);
        
        //add S2
        ssName.str("");
        ssName<<"W_S2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS2);
        
        //add A3
        ssName.str("");
        ssName<<"W_A3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dA3);
        
        //add S3
        ssName.str("");
        ssName<<"W_S3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS3);
        
        //add EV
        ssName.str("");
        ssName<<"W_EV"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dEddyViscosityTerms);
        
        //add A3
        ssName.str("");
        ssName<<"W_DwDt"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,(grid.dLocalGridNew[grid.nW][i][j][k]-grid.dLocalGridOld[grid.nW][i][j][k])
          /time.dDeltat_n);
        #endif
      }
    }
  }
  
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nV][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nV][0][0];i++){
    
    //calculate j of interface quantities
    nIInt=i+grid.nCenIntOffset[0];
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
    dR3_ip1half_n=dRSq_ip1half_n*grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dR3_im1half_n=dRSq_im1half_n*grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dU0_i_nm1half=(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0])*0.5;
    dRhoAve_ip1half_n=(grid.dLocalGridOld[grid.nDenAve][i][0][0])*0.5;
    dRhoAve_im1half_n=(grid.dLocalGridOld[grid.nDenAve][i][0][0]
      +grid.dLocalGridOld[grid.nDenAve][i-1][0][0])*0.5;
    dDM_ip1half=(grid.dLocalGridOld[grid.nDM][i][0][0])*0.5;
    dDM_im1half=(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i-1][0][0])*0.5;
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nV][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nV][0][1];j++){
      
      //calculate j of centered quantities
      nJInt=j+grid.nCenIntOffset[1];
      dDTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
        +grid.dLocalGridOld[grid.nDTheta][0][j][0])*0.5;
      dDTheta_jm1half=(grid.dLocalGridOld[grid.nDTheta][0][j][0]
        +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*0.5;
      
      for(k=grid.nStartGhostUpdateExplicit[grid.nV][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nV][0][2];k++){
        
        //calculate k of interface quantities
        nKCen=k-grid.nCenIntOffset[2];
        
        //Calculate interpolated quantities
        dDPhi_kp1half=(grid.dLocalGridOld[grid.nDPhi][0][0][nKCen]
          +grid.dLocalGridOld[grid.nDPhi][0][0][nKCen+1])*0.5;
        dDPhi_km1half=(grid.dLocalGridOld[grid.nDPhi][0][0][nKCen]
          +grid.dLocalGridOld[grid.nDPhi][0][0][nKCen-1])*0.5;
        dU_ijkp1half_nm1half=(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen+1]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen+1]
          +grid.dLocalGridOld[grid.nU][nIInt][j][nKCen]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen])*0.25;
        dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen])*0.5;
        dU_ijkp1_nm1half=(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen+1]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen+1])*0.5;
        dV_ijk_nm1half=(grid.dLocalGridOld[grid.nV][i][nJInt][nKCen]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen])*0.5;
        dV_ijkp1_nm1half=(grid.dLocalGridOld[grid.nV][i][nJInt][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen+1])*0.5;
        dV_ijkp1half_nm1half=(grid.dLocalGridOld[grid.nV][i][nJInt][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt][nKCen]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen])*0.25;
        dV_ijm1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen])*0.5;
        dV_ijm1halfkm1half_nm1half=(grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen-1])*0.5;
        dW_ijkp1half_nm1half=grid.dLocalGridOld[grid.nW][i][j][k];
        dW_ijp1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nW][i][j+1][k]
          +grid.dLocalGridOld[grid.nW][i][j][k])*0.5;
        dW_ijm1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nW][i][j-1][k]
          +grid.dLocalGridOld[grid.nW][i][j][k])*0.5;
        dW_ip1halfjkp1half_nm1half=grid.dLocalGridOld[grid.nW][i][j][k];/**\BC assume theta and
          phi velocities are constant across surface*/
        dW_im1halfjkp1half_nm1half=(grid.dLocalGridOld[grid.nW][i-1][j][k]
          +grid.dLocalGridOld[grid.nW][i][j][k])*0.5;
        dW_ijkp1_nm1half=(grid.dLocalGridOld[grid.nW][i][j][k]
          +grid.dLocalGridOld[grid.nW][i][j][k+1])*0.5;
        dW_ijk_nm1half=(grid.dLocalGridOld[grid.nW][i][j][k]
          +grid.dLocalGridOld[grid.nW][i][j][k-1])*0.5;
        dRho_ijkp1half_n=(grid.dLocalGridOld[grid.nD][i][j][nKCen]
          +grid.dLocalGridOld[grid.nD][i][j][nKCen+1])*0.5;
        dP_ijkp1_n=grid.dLocalGridOld[grid.nP][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nQ0][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nQ1][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nQ2][i][j][nKCen+1];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][nKCen]+grid.dLocalGridOld[grid.nQ0][i][j][nKCen]
          +grid.dLocalGridOld[grid.nQ1][i][j][nKCen]+grid.dLocalGridOld[grid.nQ2][i][j][nKCen];
        dEddyVisc_ip1halfjkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen])*0.25;/** \BC assume eddy viscosity is 
          zero at surface*/
        dEddyVisc_im1halfjkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen]
          +grid.dLocalGridOld[grid.nEddyVisc][i-1][j][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i-1][j][nKCen])*0.25;
        dEddyVisc_ijp1halfkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i][j+1][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j+1][nKCen]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen])*0.25;
        dEddyVisc_ijm1halfkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j-1][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j-1][nKCen])*0.25;
        dEddyVisc_ijkp1half_n=(grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen])*0.5;
        
        //calculate derived quantities
        dUmU0_ijkp1half_nm1half=dU_ijkp1half_nm1half-dU0_i_nm1half;
        d1_rhoDM_ijkp1half_n=1.0/(dRho_ijkp1half_n*grid.dLocalGridOld[grid.nDM][i][0][0]);
        dRRho_ijkp1half_n=dR_i_n*dRho_ijkp1half_n;
        dRSq_UmU0_ip1halfjkp1_n=dRSq_ip1half_n*(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen+1]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0]);
        dRSq_UmU0_im1halfjkp1_n=dRSq_im1half_n*(grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen+1]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
        dRSq_UmU0_ip1halfjk_n=dRSq_ip1half_n*(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0]);
        dRSq_UmU0_im1halfjk_n=dRSq_im1half_n*(grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
        dV_SinTheta_ijp1halfkp1_n=grid.dLocalGridOld[grid.nV][i][nJInt][nKCen+1]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0];
        dV_SinTheta_ijm1halfkp1_n=grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen+1]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0];
        dV_SinTheta_ijp1halfk_n=grid.dLocalGridOld[grid.nV][i][nJInt][nKCen]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0];
        dV_SinTheta_ijm1halfk_n=grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0];
        dW_R_ip1jkp1half_n=grid.dLocalGridOld[grid.nW][i][j][k]/dR_ip1_n;
        dW_R_im1jkp1half_n=grid.dLocalGridOld[grid.nW][i-1][j][k]/dR_im1_n;
        dW_R_ijkp1half_n=grid.dLocalGridOld[grid.nW][i][j][k]/dR_i_n;
        dW_R_ip1halfjkp1half_n=dW_ip1halfjkp1half_nm1half/grid.dLocalGridOld[grid.nR][nIInt][0][0];
        dW_R_im1halfjkp1half_n=dW_im1halfjkp1half_nm1half
          /grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
        dW_SinTheta_ijp1kp1half_n=grid.dLocalGridOld[grid.nW][i][j+1][k]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][j+1][0];
        dW_SinTheta_ijm1kp1half_n=grid.dLocalGridOld[grid.nW][i][j-1][k]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][j-1][0];
        dW_SinTheta_ijkp1half_n=grid.dLocalGridOld[grid.nW][i][j][k]
          /grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0];
        dW_SinTheta_ijp1halfkp1half_n=dW_ijp1halfkp1half_nm1half
          /grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0];
        dW_SinTheta_ijm1halfkp1half_n=dW_ijm1halfkp1half_nm1half
          /grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0];
        
        //calculate A1
        dA1CenGrad=(dW_ip1halfjkp1half_nm1half-dW_im1halfjkp1half_nm1half)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        dA1UpWindGrad=0.0;
        
        if(dUmU0_ijkp1half_nm1half<0.0){//moving in a negative direction
          dA1UpWindGrad=dA1CenGrad;/**\BC assume upwind gradient is the same as centered gradient
            across surface*/
        }
        else{//moving in a positive direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j][k]
            -grid.dLocalGridOld[grid.nW][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=dUmU0_ijkp1half_nm1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate S1
        dS1=dU_ijkp1half_nm1half*dW_ijkp1half_nm1half/dR_i_n;
        
        //calculate dA2
        dA2CenGrad=(dW_ijp1halfkp1half_nm1half-dW_ijm1halfkp1half_nm1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        dA2UpWindGrad=0.0;
        if(dV_ijkp1half_nm1half<0.0){//moning in a negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j+1][k]
            -grid.dLocalGridOld[grid.nW][i][j][k])/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        }
        else{//moving in a positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j][k]
            -grid.dLocalGridOld[grid.nW][i][j-1][k])/(grid.dLocalGridOld[grid.nDTheta][0][j-1][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        }
        dA2=dV_ijkp1half_nm1half/dR_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA2CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
        
        //calculate S2
        dS2=dV_ijkp1half_nm1half*grid.dLocalGridOld[grid.nW][i][j][k]
          *grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]/dR_i_n;
        
        //calculate A3
        dA3CenGrad=(dW_ijkp1_nm1half-dW_ijk_nm1half)/dDPhi_kp1half;
        dA3UpWindGrad=0.0;
        if(dW_ijkp1half_nm1half<0.0){//moving in a negative direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j][k+1]
            -grid.dLocalGridOld[grid.nW][i][j][k])/grid.dLocalGridOld[grid.nDPhi][0][0][nKCen+1];
        }
        else{//moving in a positive direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j][k]
            -grid.dLocalGridOld[grid.nW][i][j][k-1])/grid.dLocalGridOld[grid.nDPhi][0][0][nKCen];
        }
        dA3=dW_ijkp1half_nm1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA3CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA3UpWindGrad)/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
        
        //calculate S3
        dS3=(dP_ijkp1_n-dP_ijk_n)/(dRho_ijkp1half_n*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*dDPhi_kp1half);
        
        //calculate dDivU_ijkp1_n
        dDivU_ijkp1_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dRSq_UmU0_ip1halfjkp1_n-dRSq_UmU0_im1halfjkp1_n)/grid.dLocalGridOld[grid.nDM][i][0][0]
          +(dV_SinTheta_ijp1halfkp1_n-dV_SinTheta_ijm1halfkp1_n)/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDTheta][0][j][0])
          +(grid.dLocalGridOld[grid.nW][i][j][k+1]-grid.dLocalGridOld[grid.nW][i][j][k])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][nKCen+1]);
        
        //calculate dDivU_ijk_n
        dDivU_ijk_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dRSq_UmU0_ip1halfjk_n-dRSq_UmU0_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0]
          +(dV_SinTheta_ijp1halfk_n-dV_SinTheta_ijm1halfk_n)/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDTheta][0][j][0])
          +(grid.dLocalGridOld[grid.nW][i][j][k]-grid.dLocalGridOld[grid.nW][i][j][k-1])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][nKCen]);
        
        //calculate dTau_rp_ip1halfjkp1half_n
        dTau_rp_ip1halfjkp1half_n=dEddyVisc_ip1halfjkp1half_n*(4.0*parameters.dPi*dR3_ip1half_n
          *dRhoAve_ip1half_n*(dW_R_ip1jkp1half_n-dW_R_ijkp1half_n)/dDM_ip1half
          +(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen+1]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0])-(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0])
          /(dDPhi_kp1half*grid.dLocalGridOld[grid.nR][nIInt][0][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]));
        
        //calculate dTau_rp_im1halfjkp1half_n
        dTau_rp_im1halfjkp1half_n=dEddyVisc_im1halfjkp1half_n*(4.0*parameters.dPi*dR3_im1half_n
          *dRhoAve_im1half_n*(dW_R_ijkp1half_n-dW_R_im1jkp1half_n)/dDM_im1half
          +((grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen+1]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0])
          -(grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen]
          -grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]))
          /(dDPhi_kp1half*grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]));
        
        //calculate dTau_tp_ijp1halfkp1half_n
        dTau_tp_ijp1halfkp1half_n=dEddyVisc_ijp1halfkp1half_n
          *(grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]*(dW_SinTheta_ijp1kp1half_n
          -dW_SinTheta_ijkp1half_n)/(dR_i_n*dDTheta_jp1half)
          +(grid.dLocalGridOld[grid.nV][i][nJInt][nKCen+1]
          -grid.dLocalGridOld[grid.nV][i][nJInt][nKCen])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]*dDPhi_kp1half));
        
        //calculate dTau_tp_ijm1halfkp1half_n
        dTau_tp_ijm1halfkp1half_n=dEddyVisc_ijm1halfkp1half_n
          *(grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]*(dW_SinTheta_ijkp1half_n
          -dW_SinTheta_ijm1kp1half_n)/(dR_i_n*dDTheta_jm1half)
          +(grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen+1]
          -grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]*dDPhi_kp1half));
        
        //calculate dTau_pp_ijkp1half_n
        dTau_pp_ijkp1_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen+1]
          *((grid.dLocalGridOld[grid.nW][i][j][k+1]-grid.dLocalGridOld[grid.nW][i][j][k])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][nKCen+1])+(dU_ijkp1_nm1half-dU0_i_nm1half)
          /dR_i_n+dV_ijkp1_nm1half*grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]/dR_i_n
          -0.333333333333333*dDivU_ijkp1_n);
        
        dTau_pp_ijk_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][i][j][nKCen]
          *((grid.dLocalGridOld[grid.nW][i][j][k]-grid.dLocalGridOld[grid.nW][i][j][k-1])/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][nKCen])
          +(dU_ijk_nm1half-dU0_i_nm1half)/dR_i_n+dV_ijk_nm1half
          *grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]/dR_i_n-0.333333333333333*dDivU_ijk_n);
        
        //calculate dTA1
        dTA1=(dTau_rp_ip1halfjkp1half_n-dTau_rp_im1halfjkp1half_n)*d1_rhoDM_ijkp1half_n;
        
        //calculate dTS1
        dTS1=3.0*dEddyVisc_ijkp1half_n*(dW_R_ip1halfjkp1half_n-dW_R_im1halfjkp1half_n)
          *d1_rhoDM_ijkp1half_n;
        
        //calculate dTA2
        dTA2=(dTau_tp_ijp1halfkp1half_n-dTau_tp_ijm1halfkp1half_n)/(dRRho_ijkp1half_n
          *grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //calculate dTS2
        dTS2=2.0*grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*(dW_SinTheta_ijp1halfkp1half_n
          -dW_SinTheta_ijm1halfkp1half_n)/(dR_i_n*grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //calculate dTA3
        dTA3=(dTau_pp_ijkp1_n-dTau_pp_ijk_n)/(dRRho_ijkp1half_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][nKCen]);
        //dTA3=0.0;
        
        //calculate dTS3
        dTS3=3.0*((dU_ijkp1_nm1half-dU0_i_nm1half)-(dU_ijk_nm1half-dU0_i_nm1half)+2.0
          *grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]*(dV_ijkp1_nm1half-dV_ijk_nm1half))/(dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][nKCen]);
        
        dEddyViscosityTerms=-4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dTA1+dTS1)-dTA2-dTA3-dEddyVisc_ijkp1half_n/(dRho_ijkp1half_n*dR_i_n)*(dTS2+dTS3);
        //dEddyViscosityTerms=0.0;
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nW][i][j][k]=grid.dLocalGridOld[grid.nW][i][j][k]-time.dDeltat_n
          *(4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1)
          +dS1+dA2+dS2+dA3+dS3+dEddyViscosityTerms);
        
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
        ssName<<"W_A1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1));
        
        //add S1
        ssName.str("");
        ssName<<"W_S1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS1);
        
        //add A2
        ssName.str("");
        ssName<<"W_A2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dA2);
        
        //add S2
        ssName.str("");
        ssName<<"W_S2"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS2);
        
        //add A3
        ssName.str("");
        ssName<<"W_A3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dA3);
        
        //add S3
        ssName.str("");
        ssName<<"W_S3"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dS3);
        
        //add EV
        ssName.str("");
        ssName<<"W_EV"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,-dEddyViscosityTerms);
        
        //add A3
        ssName.str("");
        ssName<<"W_DwDt"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
          ,(grid.dLocalGridNew[grid.nW][i][j][k]-grid.dLocalGridOld[grid.nW][i][j][k])
          /time.dDeltat_n);
        #endif
      }
    }
  }
}
