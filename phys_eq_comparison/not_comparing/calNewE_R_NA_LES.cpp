void calNewE_R_NA_LES(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop){
  
  int i;
  int j;
  int k;
  int nIInt;
  double dU_ijk_np1half;
  double dU_ip1jk_np1half;
  double dU_im1jk_np1half;
  double dU0_i_np1half;
  double dE_ip1halfjk_n;
  double dE_im1halfjk_n;
  double dR_ip1_np1half;
  double dRSq_ip1_np1half;
  double dR_i_np1half;
  double dR_im1half_np1half;
  double dR_ip1half_np1half;
  double dRSq_i_np1half;
  double dR_im1_np1half;
  double dRSq_im1_np1half;
  double dRSq_ip1half_np1half;
  double dR4_ip1half_np1half;
  double dRSq_im1half_np1half;
  double dR4_im1half_np1half;
  double dRhoAve_i_n;
  double dRhoAve_ip1half_n;
  double dRhoAve_im1_n;
  double dRhoAve_im1half_n;
  double dRhoAve_ip1_n;
  double dRho_ijk_n;
  double dRho_ip1jk_n;
  double dRho_im1jk_n;
  double dRho_ip1halfjk_n;
  double dRho_im1halfjk_n;
  double dTSq_ip1jk_n;
  double dT4_ip1jk_n;
  double dTSq_ijk_n;
  double dT4_ijk_n;
  double dTSq_im1jk_n;
  double dT4_im1jk_n;
  double dKappa_ip1halfjk_n;
  double dKappa_im1halfjk_n;
  double dEddyVisc_ip1halfjk_np1half;
  double dEddyVisc_im1halfjk_np1half;
  double dDM_ip1half;
  double dDM_im1half;
  double dUR2_ip1jk_np1half;
  double dUR2_ijk_np1half;
  double dUR2_im1jk_np1half;
  double dUR2_im1halfjk_np1half;
  double dUR2_ip1halfjk_np1half;
  double dA1CenGrad;
  double dA1UpWindGrad;
  double dUmU0_ijk_np1half;
  double dA1;
  double dP_ijk_n;
  double dS1;
  double dTGrad_ip1half;
  double dTGrad_im1half;
  double dGrad_ip1half;
  double dGrad_im1half;
  double dS4;
  double dT1;
  double dEddyVisc_ip1halfjk_n;
  double dEddyVisc_im1halfjk_n;
  double dEGrad_ip1halfjk_np1half;
  double dEGrad_im1halfjk_np1half;
  double dPiSq=parameters.dPi*parameters.dPi;
  
  for(i=grid.nStartUpdateExplicit[grid.nE][0];i<grid.nEndUpdateExplicit[grid.nE][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_np1half=(grid.dLocalGridOld[grid.nR][nIInt][0][0]
      +grid.dLocalGridNew[grid.nR][nIInt][0][0])*0.5;
    dR_im1half_np1half=(grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      +grid.dLocalGridNew[grid.nR][nIInt-1][0][0])*0.5;
    dR_ip1_np1half=(grid.dLocalGridOld[grid.nR][nIInt+1][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridNew[grid.nR][nIInt+1][0][0]
      +grid.dLocalGridNew[grid.nR][nIInt][0][0])*0.25;
    dRSq_ip1_np1half=dR_ip1_np1half*dR_ip1_np1half;
    dR_im1_np1half=(grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt-2][0][0]+grid.dLocalGridNew[grid.nR][nIInt-1][0][0]
      +grid.dLocalGridNew[grid.nR][nIInt-2][0][0])*0.25;
    dRSq_im1_np1half=dR_im1_np1half*dR_im1_np1half;
    dR_i_np1half=(dR_ip1half_np1half+dR_im1half_np1half)*0.5;
    dRSq_i_np1half=dR_i_np1half*dR_i_np1half;
    dRSq_ip1half_np1half=dR_ip1half_np1half*dR_ip1half_np1half;
    dR4_ip1half_np1half=dRSq_ip1half_np1half*dRSq_ip1half_np1half;
    dRSq_im1half_np1half=dR_im1half_np1half*dR_im1half_np1half;
    dR4_im1half_np1half=dRSq_im1half_np1half*dRSq_im1half_np1half;
    dRhoAve_i_n=grid.dLocalGridOld[grid.nD][i][0][0];
    dRhoAve_ip1_n=grid.dLocalGridOld[grid.nD][i+1][0][0];
    dRhoAve_im1_n=grid.dLocalGridOld[grid.nD][i-1][0][0];
    dRhoAve_ip1half_n=(dRhoAve_i_n+dRhoAve_ip1_n)*0.5;
    dRhoAve_im1half_n=(dRhoAve_i_n+dRhoAve_im1_n)*0.5;
    dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
      +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
    dDM_ip1half=(grid.dLocalGridOld[grid.nDM][i+1][0][0]+grid.dLocalGridOld[grid.nDM][i][0][0])*0.5;
    dDM_im1half=(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i-1][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nE][1];j<grid.nEndUpdateExplicit[grid.nE][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nE][2];k<grid.nEndUpdateExplicit[grid.nE][2];k++){
        
        //Calculate interpolated quantities
        dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dU_ip1jk_np1half=(grid.dLocalGridNew[grid.nU][nIInt+1][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt][j][k])*0.5;
        dU_im1jk_np1half=(grid.dLocalGridNew[grid.nU][nIInt-1][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-2][j][k])*0.5;
        dE_ip1halfjk_n=(grid.dLocalGridOld[grid.nE][i+1][j][k]
          +grid.dLocalGridOld[grid.nE][i][j][k])*0.5;
        dE_im1halfjk_n=(grid.dLocalGridOld[grid.nE][i][j][k]
          +grid.dLocalGridOld[grid.nE][i-1][j][k])*0.5;
        dRho_ijk_n=grid.dLocalGridOld[grid.nD][i][j][k];
        dRho_ip1jk_n=grid.dLocalGridOld[grid.nD][i+1][j][k];
        dRho_im1jk_n=grid.dLocalGridOld[grid.nD][i-1][j][k];
        dRho_ip1halfjk_n=(dRho_ip1jk_n+dRho_ijk_n)*0.5;
        dRho_im1halfjk_n=(dRho_ijk_n+dRho_im1jk_n)*0.5;
        dTSq_ip1jk_n=grid.dLocalGridOld[grid.nT][i+1][j][k]
          *grid.dLocalGridOld[grid.nT][i+1][j][k];
        dT4_ip1jk_n=dTSq_ip1jk_n*dTSq_ip1jk_n;
        dTSq_ijk_n=grid.dLocalGridOld[grid.nT][i][j][k]
          *grid.dLocalGridOld[grid.nT][i][j][k];
        dT4_ijk_n=dTSq_ijk_n*dTSq_ijk_n;
        dTSq_im1jk_n=grid.dLocalGridOld[grid.nT][i-1][j][k]
          *grid.dLocalGridOld[grid.nT][i-1][j][k];
        dT4_im1jk_n=dTSq_im1jk_n*dTSq_im1jk_n;
        dKappa_ip1halfjk_n=(dT4_ip1jk_n+dT4_ijk_n)
          /(dT4_ijk_n/grid.dLocalGridOld[grid.nKappa][i][j][k]
          +dT4_ip1jk_n/grid.dLocalGridOld[grid.nKappa][i+1][j][k]);
        dKappa_im1halfjk_n=(dT4_im1jk_n+dT4_ijk_n)
          /(dT4_ijk_n/grid.dLocalGridOld[grid.nKappa][i][j][k]
          +dT4_im1jk_n/grid.dLocalGridOld[grid.nKappa][i-1][j][k]);
        dEddyVisc_ip1halfjk_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j][k]
          +grid.dLocalGridNew[grid.nEddyVisc][i+1][j][k])*0.5;
        dEddyVisc_im1halfjk_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j][k]
          +grid.dLocalGridNew[grid.nEddyVisc][i-1][j][k])*0.5;
        dEddyVisc_ip1halfjk_n=(grid.dLocalGridOld[grid.nEddyVisc][i+1][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][k])*0.5;
        dEddyVisc_im1halfjk_n=(grid.dLocalGridOld[grid.nEddyVisc][i-1][j][k]
          +grid.dLocalGridOld[grid.nEddyVisc][i][j][k])*0.5;
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP_ijk_n+=grid.dLocalGridOld[grid.nQ0][i][j][k];
        #endif
        
        //calculate derinved quantities
        dUR2_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dRSq_im1half_np1half;
        dUR2_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dRSq_ip1half_np1half;
        dUR2_ip1jk_np1half=dU_ip1jk_np1half*dRSq_ip1_np1half;
        dUR2_ijk_np1half=dU_ijk_np1half*dRSq_i_np1half;
        dUR2_im1jk_np1half=dU_im1jk_np1half*dRSq_im1_np1half;
        
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
        dA1=dUmU0_ijk_np1half*dRSq_i_np1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate dS1
        dS1=dP_ijk_n/grid.dLocalGridOld[grid.nD][i][j][k]
          *(dUR2_ip1halfjk_np1half-dUR2_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //Calculate dS4
        dTGrad_ip1half=(dT4_ip1jk_n-dT4_ijk_n)/(grid.dLocalGridOld[grid.nDM][i+1][0][0]
          +grid.dLocalGridOld[grid.nDM][i][0][0])*2.0;
        dTGrad_im1half=(dT4_ijk_n-dT4_im1jk_n)/(grid.dLocalGridOld[grid.nDM][i][0][0]
          +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        dGrad_ip1half=dRhoAve_ip1half_n*dR4_ip1half_np1half/(dKappa_ip1halfjk_n
          *dRho_ip1halfjk_n)*dTGrad_ip1half;
        dGrad_im1half=dRhoAve_im1half_n*dR4_im1half_np1half/(dKappa_im1halfjk_n
          *dRho_im1halfjk_n)*dTGrad_im1half;
        dS4=16.0*parameters.dPi*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]
          *(dGrad_ip1half-dGrad_im1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
                
        //calculate dT1
        dEGrad_ip1halfjk_np1half=dR4_ip1half_np1half*dEddyVisc_ip1halfjk_n*dRhoAve_ip1half_n
          *(grid.dLocalGridOld[grid.nE][i+1][j][k]-grid.dLocalGridOld[grid.nE][i][j][k])
          /(dRho_ip1halfjk_n*dDM_ip1half);
        dEGrad_im1halfjk_np1half=dR4_im1half_np1half*dEddyVisc_im1halfjk_n*dRhoAve_im1half_n
          *(grid.dLocalGridOld[grid.nE][i][j][k]-grid.dLocalGridOld[grid.nE][i-1][j][k])
          /(dRho_im1halfjk_n*dDM_im1half);
        dT1=16.0*dPiSq*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dEGrad_ip1halfjk_np1half
          -dEGrad_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //calculate new energy
        grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridOld[grid.nE][i][j][k]
          -time.dDeltat_np1half*(4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]*(dA1+dS1-dT1)
          -4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS4));
        
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
    dR_ip1half_np1half=(grid.dLocalGridOld[grid.nR][nIInt][0][0]
      +grid.dLocalGridNew[grid.nR][nIInt][0][0])*0.5;
    dR_im1half_np1half=(grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      +grid.dLocalGridNew[grid.nR][nIInt-1][0][0])*0.5;
    dR_im1_np1half=(grid.dLocalGridOld[grid.nR][nIInt-1][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt-2][0][0]+grid.dLocalGridNew[grid.nR][nIInt-1][0][0]
      +grid.dLocalGridNew[grid.nR][nIInt-2][0][0])*0.25;
    dRSq_im1_np1half=dR_im1_np1half*dR_im1_np1half;
    dRSq_im1half_np1half=dR_im1half_np1half*dR_im1half_np1half;
    dR4_im1half_np1half=dRSq_im1half_np1half*dRSq_im1half_np1half;
    dR_i_np1half=(dR_ip1half_np1half+dR_im1half_np1half)*0.5;
    dRSq_i_np1half=dR_i_np1half*dR_i_np1half;
    dRSq_ip1half_np1half=dR_ip1half_np1half*dR_ip1half_np1half;
    dR4_ip1half_np1half=dRSq_ip1half_np1half*dRSq_ip1half_np1half;
    dRSq_im1half_np1half=dR_im1half_np1half*dR_im1half_np1half;
    dR4_im1half_np1half=dRSq_im1half_np1half*dRSq_im1half_np1half;
    dRhoAve_i_n=grid.dLocalGridOld[grid.nD][i][0][0];
    dRhoAve_ip1_n=0.0;
    dRhoAve_im1_n=grid.dLocalGridOld[grid.nD][i-1][0][0];
    dRhoAve_ip1half_n=(dRhoAve_i_n+dRhoAve_ip1_n)*0.5;
    dRhoAve_im1half_n=(dRhoAve_i_n+dRhoAve_im1_n)*0.5;
    dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
      +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
    dDM_ip1half=(0.0+grid.dLocalGridOld[grid.nDM][i][0][0])*0.5;
    dDM_im1half=(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i-1][0][0])*0.5;
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nE][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nE][0][1];j++){
      for(k=grid.nStartGhostUpdateExplicit[grid.nE][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nE][0][2];k++){
        
        //Calculate interpolated quantities
        dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
          +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
        dU_im1jk_np1half=(grid.dLocalGridNew[grid.nU][nIInt-2][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dE_ip1halfjk_n=(grid.dLocalGridOld[grid.nE][i][j][k]);/**\BC Missing
          grid.dLocalGridOld[grid.nE][i+1][j][k] in calculation of \f$E_{i+1/2,j,k}\f$ 
          setting it equal to value at i.*/
        dE_im1halfjk_n=(grid.dLocalGridOld[grid.nE][i][j][k]
          +grid.dLocalGridOld[grid.nE][i-1][j][k])*0.5;
        dTSq_ijk_n=grid.dLocalGridOld[grid.nT][i][j][k]
          *grid.dLocalGridOld[grid.nT][i][j][k];
        dT4_ijk_n=dTSq_ijk_n*dTSq_ijk_n;
        dTSq_im1jk_n=grid.dLocalGridOld[grid.nT][i-1][j][k]
          *grid.dLocalGridOld[grid.nT][i-1][j][k];
        dT4_im1jk_n=dTSq_im1jk_n*dTSq_im1jk_n;
        dKappa_im1halfjk_n=(dT4_im1jk_n+dT4_ijk_n)
          /(dT4_ijk_n/grid.dLocalGridOld[grid.nKappa][i][j][k]
          +dT4_im1jk_n/grid.dLocalGridOld[grid.nKappa][i-1][j][k]);
        dRho_ijk_n=grid.dLocalGridOld[grid.nD][i][j][k];
        dRho_im1jk_n=grid.dLocalGridOld[grid.nD][i-1][j][k];
        dRho_im1halfjk_n=(dRho_ijk_n+dRho_im1jk_n)*0.5;
        
        dEddyVisc_ip1halfjk_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j][k]
          +grid.dLocalGridNew[grid.nEddyVisc][i+1][j][k])*0.5;
        dEddyVisc_im1halfjk_np1half=(grid.dLocalGridNew[grid.nEddyVisc][i][j][k]
          +grid.dLocalGridNew[grid.nEddyVisc][i-1][j][k])*0.5;
        
        //Calcuate dA1
        dA1CenGrad=(dE_ip1halfjk_n-dE_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0];
        dA1UpWindGrad=0.0;
        dUmU0_ijk_np1half=(dU_ijk_np1half-dU0_i_np1half);
        if(dUmU0_ijk_np1half<0.0){//moving in the negative direction
          dA1UpWindGrad=dA1CenGrad;/**\BC grid.dLocalGridOld[grid.nDM][i+1][0][0] and 
            grid.dLocalGridOld[grid.nE][i+1][j][k] missing in the calculation of upwind gradient in 
            dA1. Using the centered gradient instead.*/
        }
        else{//moving in the postive direction*/
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k]
            -grid.dLocalGridOld[grid.nE][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=dUmU0_ijk_np1half*dRSq_i_np1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate dUR2_im1half_np1half
        dUR2_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dRSq_im1half_np1half;
        
        //calculate dR2_ip1half_np1half
        dUR2_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dRSq_ip1half_np1half;
        
        //calculate dURsq_ijk_np1half
        dUR2_ijk_np1half=dU_ijk_np1half*dRSq_i_np1half;
        
        //calculate dURsq_im1jk_np1half
        dUR2_im1jk_np1half=dU_im1jk_np1half*dRSq_im1_np1half;
        
        //calculate dS1
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP_ijk_n+=grid.dLocalGridOld[grid.nQ0][i][j][k];
        #endif
        dS1=dP_ijk_n/grid.dLocalGridOld[grid.nD][i][j][k]
          *(dUR2_ip1halfjk_np1half-dUR2_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //Calculate dS4
        dTGrad_im1half=(dT4_ijk_n-dT4_im1jk_n)/(grid.dLocalGridOld[grid.nDM][i][0][0]
          +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        dGrad_ip1half=-3.0*dRSq_ip1half_np1half*dT4_ijk_n/(8.0*parameters.dPi);/**\BC
          Missing grid.dLocalGridOld[grid.nT][i+1][0][0]*/
        dGrad_im1half=dRhoAve_im1half_n*dR4_im1half_np1half/(dKappa_im1halfjk_n
          *dRho_im1halfjk_n)*dTGrad_im1half;
        dS4=16.0*parameters.dPi*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]
          *(dGrad_ip1half-dGrad_im1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //calculate dT1
        dEGrad_ip1halfjk_np1half=0.0;/**\BC missing energy outside the model, assuming it is the 
          same as that in the last zone. That causes this term to be zero.*/
        dEGrad_im1halfjk_np1half=dR4_im1half_np1half*dEddyVisc_im1halfjk_n*dRhoAve_im1half_n
          *(grid.dLocalGridOld[grid.nE][i][j][k]-grid.dLocalGridOld[grid.nE][i-1][j][k])
          /(dRho_im1halfjk_n*dDM_im1half);
        dT1=16.0*dPiSq*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dEGrad_ip1halfjk_np1half
          -dEGrad_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //calculate new energy
        grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridOld[grid.nE][i][j][k]
          -time.dDeltat_np1half*(4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]*(dA1+dS1-dT1)
          -4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS4));
        
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
