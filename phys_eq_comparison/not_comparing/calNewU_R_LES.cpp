void calNewU_R_LES(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop){
  //calculate new u
  int i;
  int j;
  int k;
  int nICen;
  double dRho_ip1halfjk_n;
  double dP_ip1jk_n;
  double dP_ijk_n;
  double dA1;
  double dS1;
  double dS4;
  double dTA1;
  double dTS1;
  double dTS4;
  double dDivU_ijk_n;
  double dDivU_ip1jk_n;
  double dDivU_ip1halfjk_n;//used at surface boundary condition
  double dTau_rr_ip1jk_n;
  double dTau_rr_ip1halfjk_n;
  double dTau_rr_ijk_n;
  double dR_i_n;
  double dR_ip1_n;
  double dU_ip1jk_nm1half;
  double dU_ijk_nm1half;
  double dRSq_ip1half_n;
  double dRSq_ip1_n;
  double dRSq_i_n;
  double dRSq_ip3half_n;
  double dRSq_im1half_n;
  double dRSqU_ip3halfjk_n;
  double dRSqU_ip1halfjk_n;
  double dRSqU_im1halfjk_n;
  double dDM_ip1half;
  double dEddyVisc_ip1halfjk_n;
  double dA1CenGrad;
  double dA1UpWindGrad;
  double dU_U0_Diff;
  double dDonorFrac_ip1half;
  
  for(i=grid.nStartUpdateExplicit[grid.nU][0];i<grid.nEndUpdateExplicit[grid.nU][0];i++){
    
    //calculate i of centered quantities
    nICen=i-grid.nCenIntOffset[0];
    dR_ip1_n=(grid.dLocalGridOld[grid.nR][i+1][0][0]+grid.dLocalGridOld[grid.nR][i][0][0])*0.5;
    dR_i_n=(grid.dLocalGridOld[grid.nR][i][0][0]+grid.dLocalGridOld[grid.nR][i-1][0][0])*0.5;
    dRSq_ip1_n=dR_ip1_n*dR_ip1_n;
    dRSq_i_n=dR_i_n*dR_i_n;
    dRSq_ip1half_n=grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nR][i][0][0];
    dRSq_ip3half_n=grid.dLocalGridOld[grid.nR][i+1][0][0]*grid.dLocalGridOld[grid.nR][i+1][0][0];
    dRSq_im1half_n=grid.dLocalGridOld[grid.nR][i-1][0][0]*grid.dLocalGridOld[grid.nR][i-1][0][0];
    dDM_ip1half=(grid.dLocalGridOld[grid.nDM][nICen+1][0][0]
      +grid.dLocalGridOld[grid.nDM][nICen][0][0])*0.5;
    dDonorFrac_ip1half=(grid.dLocalGridOld[grid.nDonorCellFrac][nICen+1][0][0]
      +grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nU][1];j<grid.nEndUpdateExplicit[grid.nU][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nU][2];k<grid.nEndUpdateExplicit[grid.nU][2];k++){
        
        //calculate interpolated quantities
        dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
        dU_ip1jk_nm1half=(grid.dLocalGridOld[grid.nU][i+1][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
        dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][i+1][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
        dEddyVisc_ip1halfjk_n=(grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][nICen+1][j][k])*0.5;
        dP_ip1jk_n=grid.dLocalGridOld[grid.nP][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nQ0][nICen+1][j][k];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][nICen][j][k]
          +grid.dLocalGridOld[grid.nQ0][nICen][j][k];
        
        //calculate derived quantities
        dRSqU_ip3halfjk_n=dRSq_ip3half_n*grid.dLocalGridOld[grid.nU][i+1][j][k];
        dRSqU_ip1halfjk_n=dRSq_ip1half_n*grid.dLocalGridOld[grid.nU][i][j][k];
        dRSqU_im1halfjk_n=dRSq_im1half_n*grid.dLocalGridOld[grid.nU][i-1][j][k];
        
        //cal DivU_ip1jk_n
        dDivU_ip1jk_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][nICen+1][0][0]
          *(dRSqU_ip3halfjk_n-dRSqU_ip1halfjk_n)/grid.dLocalGridOld[grid.nDM][nICen+1][0][0];
        
        //cal DivU_ijk_n
        dDivU_ijk_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][nICen][0][0]
          *(dRSqU_ip1halfjk_n-dRSqU_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][nICen][0][0];
        
        //cal Tau_rr_ip1jk_n
        dTau_rr_ip1jk_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][nICen+1][j][k]*(4.0*parameters.dPi
          *dRSq_ip1_n*grid.dLocalGridOld[grid.nD][nICen+1][0][0]
          *(grid.dLocalGridOld[grid.nU][i+1][j][k]-grid.dLocalGridOld[grid.nU][i][j][k])
          /grid.dLocalGridOld[grid.nDM][nICen+1][0][0]-0.3333333333333333*dDivU_ip1jk_n);
        
        //cal Tau_rr_ijk_n
        dTau_rr_ijk_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]*(4.0*parameters.dPi
          *dRSq_i_n*grid.dLocalGridOld[grid.nD][nICen][0][0]*(grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU][i-1][j][k])/grid.dLocalGridOld[grid.nDM][nICen][0][0]
          -0.3333333333333333*dDivU_ijk_n);
        
        //cal dTA1
        dTA1=1.0/dRho_ip1halfjk_n*(dTau_rr_ip1jk_n-dTau_rr_ijk_n)/dDM_ip1half;
        
        //cal dTS1
        dTS1=dEddyVisc_ip1halfjk_n/(dRho_ip1halfjk_n*grid.dLocalGridOld[grid.nR][i][0][0])*(4.0
          *(dU_ip1jk_nm1half-dU_ijk_nm1half)/dDM_ip1half);
        
        //cal dTS4
        dTS4=4.0*grid.dLocalGridOld[grid.nU][i][j][k]/grid.dLocalGridOld[grid.nR][i][0][0];
        
        //Calculate dA1
        dA1CenGrad=(dU_ip1jk_nm1half-dU_ijk_nm1half)
          /(grid.dLocalGridOld[grid.nDM][nICen+1][0][0]
          +grid.dLocalGridOld[grid.nDM][nICen][0][0])*2.0;
        dA1UpWindGrad=0.0;
        dU_U0_Diff=grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0];
        if(dU_U0_Diff<0.0){//moving from outside in
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nU][i+1][j][k]
            -grid.dLocalGridOld[grid.nU][i][j][k])
            /grid.dLocalGridOld[grid.nDM][nICen+1][0][0];
        }
        else{//moving from inside out
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU][i-1][j][k])
            /grid.dLocalGridOld[grid.nDM][nICen][0][0];
        }
        dA1=dU_U0_Diff*((1.0-dDonorFrac_ip1half)*dA1CenGrad+dDonorFrac_ip1half*dA1UpWindGrad);
        
        //calculate dS1
        dS1=(dP_ip1jk_n-dP_ijk_n)/(grid.dLocalGridOld[grid.nDM][nICen+1][0][0]
          +grid.dLocalGridOld[grid.nDM][nICen][0][0])*2.0/dRho_ip1halfjk_n;
        
        //calculate dS4
        dS4=parameters.dG*grid.dLocalGridOld[grid.nM][i][0][0]/dRSq_ip1half_n;
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nU][i][j][k]=grid.dLocalGridOld[grid.nU][i][j][k]
          -time.dDeltat_n*(4.0*parameters.dPi*dRho_ip1halfjk_n*dRSq_ip1half_n*(dA1+dS1+dTA1+dTS1)
          +dS4+dEddyVisc_ip1halfjk_n/(dRho_ip1halfjk_n*grid.dLocalGridOld[grid.nR][i][0][0])
          *(dTS4));
      }
    }
  }
  
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nU][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nU][0][0];i++){
    
    //calculate i of centered quantities
    nICen=i-grid.nCenIntOffset[0];
    dR_i_n=(grid.dLocalGridOld[grid.nR][i][0][0]+grid.dLocalGridOld[grid.nR][i-1][0][0])*0.5;
    dRSq_i_n=dR_i_n*dR_i_n;
    dRSq_ip1half_n=grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nR][i][0][0];
    dRSq_im1half_n=grid.dLocalGridOld[grid.nR][i-1][0][0]*grid.dLocalGridOld[grid.nR][i-1][0][0];
    dDM_ip1half=(0.0+grid.dLocalGridOld[grid.nDM][nICen][0][0])*0.5;
    dDonorFrac_ip1half=grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0];
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nU][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nU][0][1];j++){
      for(k=grid.nStartGhostUpdateExplicit[grid.nU][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nU][0][2];k++){
        
        //calculate interpolated quantities
        dRho_ip1halfjk_n=(0.0+grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;/**\BC Missing 
          grid.dLocalGridOld[grid.nD][nICen+1][j][k] in calculation of \f$\rho_{i+1/2}\f$, setting
          it to 0.0*/
        dU_ip1jk_nm1half=grid.dLocalGridOld[grid.nU][i][j][k];/**\BC missing 
          grid.dLocalGridOld[grid.nU][i+1][j][k] using velocity at i*/
        dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i-1][j][k])*0.5;
        dEddyVisc_ip1halfjk_n=(grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k])*0.5;/**\BC Assuming 
          eddy viscosity outside model is zero.*/
        
        //calculate derived quantities
        /*dRSqU_ip3halfjk_n=dRSq_ip3half_n*grid.dLocalGridOld[grid.nU][i][j][k]; missing outside 
          velocity */
        dRSqU_ip1halfjk_n=dRSq_ip1half_n*grid.dLocalGridOld[grid.nU][i][j][k];
        dRSqU_im1halfjk_n=dRSq_im1half_n*grid.dLocalGridOld[grid.nU][i-1][j][k];
        
        //calculate advection term in x1-direction
        dP_ijk_n=grid.dLocalGridOld[grid.nP][nICen][j][k]
          +grid.dLocalGridOld[grid.nQ0][nICen][j][k];
        dP_ip1jk_n=-1.0*dP_ijk_n;/**\BC Missing grid.dLocalGridOld[grid.nP][nICen+1][j][k] in 
          calculation of \f$S_1\f$, setting it to -1.0*grid.dLocalGridOld[grid.nP][nICen][j][k].*/
        
        //Calculate dA1
        dA1CenGrad=(dU_ip1jk_nm1half-dU_ijk_nm1half)/grid.dLocalGridOld[grid.nDM][nICen][0][0]*2.0;
          /**\BC Missing grid.dLocalGridOld[grid.nDM][nICen+1][0][0] in calculation of centered 
          \f$A_1\f$ gradient, setting it to zero.*/
        dA1UpWindGrad=0.0;
        if(grid.dLocalGridOld[grid.nU][i][j][k]<0.0){//moving from outside in
          dA1UpWindGrad=dA1CenGrad;/**\BC Missing grid.dLocalGridOld[grid.nU][i+1][j][k] and 
            grid.dLocalGridOld[grid.nDM][nICen+1][0][0] in calculation of upwind gradient, when 
            moving inward. Using centered gradient instead.*/
        }
        else{//moving from inside out
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU][i-1][j][k])/grid.dLocalGridOld[grid.nDM][nICen][0][0];
        }
        dA1=(grid.dLocalGridOld[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU0][i][0][0])
          *((1.0-dDonorFrac_ip1half)*dA1CenGrad+dDonorFrac_ip1half*dA1UpWindGrad);
        
        //calculate source terms in x1-direction
        dS1=(dP_ip1jk_n-dP_ijk_n)/(grid.dLocalGridOld[grid.nDM][nICen][0][0]*(0.5+parameters.dAlpha
        +parameters.dAlphaExtra))
          /dRho_ip1halfjk_n;
        dS4=parameters.dG*grid.dLocalGridOld[grid.nM][i][0][0]/dRSq_ip1half_n;
        
        //cal DivU_ip1jk_n
        dDivU_ip1halfjk_n=4.0*parameters.dPi*dRho_ip1halfjk_n*(dRSqU_ip1halfjk_n-dRSqU_im1halfjk_n)
          /dDM_ip1half;
        
        //cal DivU_ijk_n
        dDivU_ijk_n=4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][nICen][0][0]
          *(dRSqU_ip1halfjk_n-dRSqU_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][nICen][0][0];
        
        //cal Tau_rr_ip1halfjk_n
        dTau_rr_ip1halfjk_n=2.0*dEddyVisc_ip1halfjk_n*(4.0*parameters.dPi*dRSq_ip1half_n
          *dRho_ip1halfjk_n*(grid.dLocalGridOld[grid.nU][i][j][k]-dU_ijk_nm1half)/dDM_ip1half
          -0.3333333333333333*dDivU_ip1halfjk_n);
        
        //cal dTS4
        dTS4=4.0*grid.dLocalGridOld[grid.nU][i][j][k]/grid.dLocalGridOld[grid.nR][i][0][0];
        
        //cal Tau_rr_ijk_n
        dTau_rr_ijk_n=2.0*grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]*(4.0*parameters.dPi
          *dRSq_i_n*grid.dLocalGridOld[grid.nD][nICen][0][0]*(grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU][i-1][j][k])/grid.dLocalGridOld[grid.nDM][nICen][0][0]
          -0.3333333333333333*dDivU_ijk_n);
        
        //cal dTA1
        dTA1=1.0/dRho_ip1halfjk_n*(dTau_rr_ip1halfjk_n-dTau_rr_ijk_n)/dDM_ip1half;
        
        //cal dTS1
        dTS1=dEddyVisc_ip1halfjk_n/(dRho_ip1halfjk_n*grid.dLocalGridOld[grid.nR][i][0][0])*(4.0
          *(grid.dLocalGridOld[grid.nU][i][j][k]-dU_ijk_nm1half)/dDM_ip1half);
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nU][i][j][k]=grid.dLocalGridOld[grid.nU][i][j][k]
          -time.dDeltat_n*(4.0*parameters.dPi*dRho_ip1halfjk_n*dRSq_ip1half_n*(dA1+dS1)+dS4);
      }
    }
  }
}
