void calNewV_RT_LES(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop){
  int i;
  int j;
  int k;
  int nIInt;
  int nJCen;
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
  double dDTheta_jp1half;
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
  double dU_U0_Diff_ijp1halfk_nm1half;
  double dA1CenGrad;
  double dA1UpWindGrad;
  double dA1;
  double dS1;
  double dA2CenGrad;
  double dA2UpWindGrad;
  double dA2;
  double dS2;
  double dTau_rt_ip1halfjp1halfk_n;
  double dTau_rt_im1halfjp1halfk_n;
  double dDivU_ijp1k_n;
  double dDivU_ijk_n;
  double dTau_tt_ijp1k_n;
  double dTau_tt_ijk_n;
  double dTA1;
  double dTS1;
  double dTA2;
  double dTS2;
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
        dRho_ijp1halfk_n=(grid.dLocalGridOld[grid.nD][i][nJCen][k]
          +grid.dLocalGridOld[grid.nD][i][nJCen+1][k])*0.5;
        dP_ijp1k_n=grid.dLocalGridOld[grid.nP][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nQ0][i][nJCen+1][k]+grid.dLocalGridOld[grid.nQ1][i][nJCen+1][k];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][nJCen][k]+grid.dLocalGridOld[grid.nQ0][i][nJCen][k]
          +grid.dLocalGridOld[grid.nQ1][i][nJCen][k];
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
          +(grid.dLocalGridOld[grid.nV][i][j+1][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j+1][0]
          -grid.dLocalGridOld[grid.nV][i][j][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0])
          /(grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0]
          *dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen+1][0]);
        
        //calculate DivU_ijk_n
        dDivU_ijk_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dRSqUmU0_ip1halfjk_n-dRSqUmU0_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0]
          +(grid.dLocalGridOld[grid.nV][i][j][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]
          -grid.dLocalGridOld[grid.nV][i][j-1][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j-1][0])
          /(grid.dLocalGridOld[grid.nDTheta][0][nJCen][0]*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen][0]);
        
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
        
        //calculate TS4
        dTS4=2.0*grid.dLocalGridOld[grid.nV][i][j][k]
          *grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]
          *grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]/dR_i_n;
        
        dEddyViscosityTerms=-4.0*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dTA1+dTS1)
          -dTA2-dEddyVisc_ijp1halfk_n/(dRho_ijp1halfk_n*dR_i_n)*(dTS2-dTS4);
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nV][i][j][k]=grid.dLocalGridOld[grid.nV][i][j][k]
          -time.dDeltat_n*(4.0*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1)
          +dS1+dA2+dS2+dEddyViscosityTerms);
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
        dRho_ijp1halfk_n=(grid.dLocalGridOld[grid.nD][i][nJCen][k]
          +grid.dLocalGridOld[grid.nD][i][nJCen+1][k])*0.5;
        dP_ijp1k_n=grid.dLocalGridOld[grid.nP][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nQ0][i][nJCen+1][k]+grid.dLocalGridOld[grid.nQ1][i][nJCen+1][k];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][nJCen][k]+grid.dLocalGridOld[grid.nQ0][i][nJCen][k]
          +grid.dLocalGridOld[grid.nQ1][i][nJCen][k];
        dEddyVisc_ip1halfjp1halfk_n=(grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k])*0.25;/**\BC Assuming eddy viscosity is
          zero at surface.*/
        dEddyVisc_im1halfjp1halfk_n=(grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i-1][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i-1][nJCen+1][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k])*0.25;
        dEddyVisc_ijp1halfk_n=(grid.dLocalGridOld[grid.nEddyVisc][i][nJCen][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][nJCen+1][k])*0.5;
        
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
          +(grid.dLocalGridOld[grid.nV][i][j+1][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j+1][0]
          -grid.dLocalGridOld[grid.nV][i][j][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0])
          /(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen+1][0]
          *grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0]);
        
        //calculate DivU_ijk_n
        dDivU_ijk_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dRSqUmU0_ip1halfjk_n-dRSqUmU0_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0]
          +(grid.dLocalGridOld[grid.nV][i][j][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]
          -grid.dLocalGridOld[grid.nV][i][j-1][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j-1][0])
          /(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJK][0][nJCen][0]
          *grid.dLocalGridOld[grid.nDTheta][0][nJCen][0]);
        
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
        
        //calculate TS4
        dTS4=2.0*grid.dLocalGridOld[grid.nV][i][j][k]
          *grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]
          *grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]/dR_i_n;
        
        dEddyViscosityTerms=-4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dTA1+dTS1)-dTA2-dEddyVisc_ijp1halfk_n/(dRho_ijp1halfk_n*dR_i_n)*(dTS2-dTS4);
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nV][i][j][k]=grid.dLocalGridOld[grid.nV][i][j][k]
          -time.dDeltat_n*(4.0*parameters.dPi*dRSq_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dA1)+dS1+dA2+dS2+dEddyViscosityTerms);
      }
    }
  }
}
