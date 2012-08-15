void calNewE_RT_NA(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop){
  int i;
  int j;
  int k;
  int nIInt;
  int nJInt;
  double dU0_i_np1half;
  double dU_ijk_np1half;
  double dE_ip1halfjk_n;
  double dE_im1halfjk_n;
  double dR_i_n;
  double dR_im1half_n;
  double dR_ip1half_n;
  double dRSq_i_n;
  double dRSq_ip1half;
  double dR4_ip1half;
  double dR_im1half_sq;
  double dR_im1half_4;
  double dV_ijk_np1half;
  double dE_ijp1halfk_n;
  double dE_ijm1halfk_n;
  double dVSinTheta_ijp1halfk_np1half;
  double dVSinTheta_ijm1halfk_np1half;
  double dRhoAve_ip1half;
  double dRhoAve_im1half;
  double dRho_im1halfjk;
  double dRho_ip1halfjk;
  double dRho_ijp1halfk;
  double dRho_ijm1halfk;
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
  double dKappa_ip1halfjk_n;
  double dKappa_im1halfjk_n;
  double dKappa_ijp1halfk_n;
  double dKappa_ijm1halfk_n;
  double dA1UpWindGrad;
  double dA1CenGrad;
  double dU_U0_Diff;
  double dA1;
  double dUR2_im1half_np1half;
  double dUR2_ip1half_np1half;
  double dP_ijk_n;
  double dS1;
  double dA2CenGrad;
  double dA2UpWindGrad;
  double dA2;
  double dS2;
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
  
  for(i=grid.nStartUpdateExplicit[grid.nE][0];i<grid.nEndUpdateExplicit[grid.nE][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
      +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
    dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])
      *0.5;
    dR_im1half_n=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_ip1half_n=grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dRSq_i_n=dR_i_n*dR_i_n;
    dRSq_ip1half=dR_ip1half_n*dR_ip1half_n;
    dR4_ip1half=dRSq_ip1half*dRSq_ip1half;
    dR_im1half_sq=dR_im1half_n*dR_im1half_n;
    dR_im1half_4=dR_im1half_sq*dR_im1half_sq;
    dRhoAve_ip1half=(grid.dLocalGridOld[grid.nDenAve][i+1][0][0]
      +grid.dLocalGridOld[grid.nDenAve][i][0][0])*0.5;
    dRhoAve_im1half=(grid.dLocalGridOld[grid.nDenAve][i][0][0]
      +grid.dLocalGridOld[grid.nDenAve][i-1][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nE][1];j<grid.nEndUpdateExplicit[grid.nE][1];j++){
      
      //calculate i for interface centered quantities
      nJInt=j+grid.nCenIntOffset[1];
      
      for(k=grid.nStartUpdateExplicit[grid.nE][2];k<grid.nEndUpdateExplicit[grid.nE][2];k++){
        
        //Calculate interpolated quantities
        dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dE_ip1halfjk_n=(grid.dLocalGridOld[grid.nE][i+1][j][k]+grid.dLocalGridOld[grid.nE][i][j][k])
          *0.5;
        dE_im1halfjk_n=(grid.dLocalGridOld[grid.nE][i][j][k]+grid.dLocalGridOld[grid.nE][i-1][j][k])
          *0.5;
        dV_ijk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k])*0.5;
        dE_ijp1halfk_n=(grid.dLocalGridOld[grid.nE][i][j+1][k]+grid.dLocalGridOld[grid.nE][i][j][k])
          *0.5;
        dE_ijm1halfk_n=(grid.dLocalGridOld[grid.nE][i][j][k]+grid.dLocalGridOld[grid.nE][i][j-1][k])
          *0.5;
        dVSinTheta_ijp1halfk_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          *grid.dLocalGridNew[grid.nV][i][nJInt][k];
        dVSinTheta_ijm1halfk_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
          *grid.dLocalGridNew[grid.nV][i][nJInt-1][k];
        dRho_ip1halfjk=(grid.dLocalGridOld[grid.nD][i+1][j][k]+grid.dLocalGridOld[grid.nD][i][j][k])
          *0.5;
        dRho_im1halfjk=(grid.dLocalGridOld[grid.nD][i][j][k]+grid.dLocalGridOld[grid.nD][i-1][j][k])
          *0.5;
        dRho_ijp1halfk=(grid.dLocalGridOld[grid.nD][i][j+1][k]+grid.dLocalGridOld[grid.nD][i][j][k])
          *0.5;
        dRho_ijm1halfk=(grid.dLocalGridOld[grid.nD][i][j][k]+grid.dLocalGridOld[grid.nD][i][j-1][k])
          *0.5;
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
        
        //Calcuate dA1
        dA1CenGrad=(dE_ip1halfjk_n-dE_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0];
        dU_U0_Diff=(dU_ijk_np1half-dU0_i_np1half);
        if(dU_U0_Diff<0.0){//moving in the negative direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nE][i+1][j][k]
            -grid.dLocalGridOld[grid.nE][i][j][k])/(grid.dLocalGridOld[grid.nDM][i+1][0][0]
            +grid.dLocalGridOld[grid.nDM][i][0][0])*2.0;
        }
        else{//moving in the postive direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k]
            -grid.dLocalGridOld[grid.nE][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=dU_U0_Diff*dRSq_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA1CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate dS1
        dUR2_im1half_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dR_im1half_n
          *dR_im1half_n;
        dUR2_ip1half_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dR_ip1half_n
          *dR_ip1half_n;
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP_ijk_n+=grid.dLocalGridOld[grid.nQ0][i][j][k]+grid.dLocalGridOld[grid.nQ1][i][j][k];
        #endif
        dS1=dP_ijk_n/grid.dLocalGridOld[grid.nD][i][j][k]
          *(dUR2_ip1half_np1half-dUR2_im1half_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
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
        dA2=dV_ijk_np1half/dR_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA2CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
          
        //Calcualte dS2
        dS2=dP_ijk_n/(grid.dLocalGridOld[grid.nD][i][j][k]*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDTheta][0][j][0])
          *(dVSinTheta_ijp1halfk_np1half-dVSinTheta_ijm1halfk_np1half);
        
        //Calculate dS4
        dTGrad_ip1half=(dT4_ip1jk_n-dT4_ijk_n)/(grid.dLocalGridOld[grid.nDM][i+1][0][0]
          +grid.dLocalGridOld[grid.nDM][i][0][0])*2.0;
        dTGrad_im1half=(dT4_ijk_n-dT4_im1jk_n)/(grid.dLocalGridOld[grid.nDM][i][0][0]
          +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        dGrad_ip1half=dRhoAve_ip1half*dR4_ip1half/(dKappa_ip1halfjk_n*dRho_ip1halfjk)
          *dTGrad_ip1half;
        dGrad_im1half=dRhoAve_im1half*dR_im1half_4/(dKappa_im1halfjk_n*dRho_im1halfjk)
          *dTGrad_im1half;
        dS4=16.0*parameters.dPi*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dGrad_ip1half-dGrad_im1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //Calculate dS5
        dTGrad_jp1half=(dT4_ijp1k_n-dT4_ijk_n)/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
          +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        dTGrad_jm1half=(dT4_ijk_n-dT4_ijm1k_n)/(grid.dLocalGridOld[grid.nDTheta][0][j][0]
          +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;;
        dGrad_jp1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          /(dKappa_ijp1halfk_n*dRho_ijp1halfk*dR_i_n)*dTGrad_jp1half;
        dGrad_jm1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
          /(dKappa_ijm1halfk_n*dRho_ijm1halfk*dR_i_n)*dTGrad_jm1half;;
        dS5=(dGrad_jp1half-dGrad_jm1half)/(grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *dR_i_n*grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //calculate new energy
        grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridOld[grid.nE][i][j][k]
          -time.dDeltat_np1half*(4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1+dS1)
          +dA2+dS2-4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS4+dS5));
        
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
    dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt-1][0][0])*0.5;
    dR_im1half_n=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_ip1half_n=grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dRSq_i_n=dR_i_n*dR_i_n;
    dRSq_ip1half=dR_ip1half_n*dR_ip1half_n;
    dR_im1half_sq=dR_im1half_n*dR_im1half_n;
    dR_im1half_4=dR_im1half_sq*dR_im1half_sq;
    dRhoAve_im1half=(grid.dLocalGridOld[grid.nDenAve][i][0][0]
      +grid.dLocalGridOld[grid.nDenAve][i-1][0][0])*0.5;
    
    for(int j=grid.nStartGhostUpdateExplicit[grid.nE][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nE][0][1];j++){
      
      //calculate i for interface centered quantities
      nJInt=j+grid.nCenIntOffset[1];
      
      for(k=grid.nStartGhostUpdateExplicit[grid.nE][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nE][0][2];k++){
        
        //Calculate interpolated quantities
        dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
          +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
        dE_ip1halfjk_n=(grid.dLocalGridOld[grid.nE][i][j][k]);/**\BC Missing
          grid.dLocalGridOld[grid.nE][i+1][j][k] in calculation of \f$E_{i+1/2,j,k}\f$ 
          setting it equal to value at i.*/
        dE_im1halfjk_n=(grid.dLocalGridOld[grid.nE][i][j][k]
          +grid.dLocalGridOld[grid.nE][i-1][j][k])*0.5;
        dV_ijk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k])*0.5;
        dE_ijp1halfk_n=(grid.dLocalGridOld[grid.nE][i][j+1][k]
          +grid.dLocalGridOld[grid.nE][i][j][k])*0.5;
        dE_ijm1halfk_n=(grid.dLocalGridOld[grid.nE][i][j][k]
          +grid.dLocalGridOld[grid.nE][i][j-1][k])*0.5;
        dVSinTheta_ijp1halfk_np1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          *grid.dLocalGridNew[grid.nV][i][nJInt][k];
        dVSinTheta_ijm1halfk_np1half=
          grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
          *grid.dLocalGridNew[grid.nV][i][nJInt-1][k];
        dRho_im1halfjk=(grid.dLocalGridOld[grid.nD][i][j][k]
          +grid.dLocalGridOld[grid.nD][i-1][j][k])*0.5;
        dRho_ijp1halfk=(grid.dLocalGridOld[grid.nD][i][j+1][k]
          +grid.dLocalGridOld[grid.nD][i][j][k])*0.5;
        dRho_ijm1halfk=(grid.dLocalGridOld[grid.nD][i][j][k]
          +grid.dLocalGridOld[grid.nD][i][j-1][k])*0.5;
        dTSq_ijk_n=grid.dLocalGridOld[grid.nT][i][j][k]
          *grid.dLocalGridOld[grid.nT][i][j][k];
        dT4_ijk_n=dTSq_ijk_n*dTSq_ijk_n;
        dTSq_im1jk_n=grid.dLocalGridOld[grid.nT][i-1][j][k]
          *grid.dLocalGridOld[grid.nT][i-1][j][k];
        dT4_im1jk_n=dTSq_im1jk_n*dTSq_im1jk_n;
        dTSq_ijp1k_n=grid.dLocalGridOld[grid.nT][i][j+1][k]
          *grid.dLocalGridOld[grid.nT][i][j+1][k];
        dT4_ijp1k_n=dTSq_ijp1k_n*dTSq_ijp1k_n;
        dTSq_ijm1k_n=grid.dLocalGridOld[grid.nT][i][j-1][k]
          *grid.dLocalGridOld[grid.nT][i][j-1][k];
        dT4_ijm1k_n=dTSq_ijm1k_n*dTSq_ijm1k_n;
        dKappa_im1halfjk_n=(dT4_im1jk_n+dT4_ijk_n)/(dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]+dT4_im1jk_n
          /grid.dLocalGridOld[grid.nKappa][i-1][j][k]);
        dKappa_ijp1halfk_n=(dT4_ijp1k_n+dT4_ijk_n)/(dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]+dT4_ijp1k_n
          /grid.dLocalGridOld[grid.nKappa][i][j+1][k]);
        dKappa_ijm1halfk_n=(dT4_ijm1k_n+dT4_ijk_n)/(dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]+dT4_ijm1k_n
          /grid.dLocalGridOld[grid.nKappa][i][j-1][k]);
        
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
        dA1=dU_U0_Diff*dRSq_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA1CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate dS1
        dUR2_im1half_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dR_im1half_n
          *dR_im1half_n;
        dUR2_ip1half_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dR_ip1half_n
          *dR_ip1half_n;
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP_ijk_n+=grid.dLocalGridOld[grid.nQ0][i][j][k];
        #endif
        dS1=dP_ijk_n/grid.dLocalGridOld[grid.nD][i][j][k]
          *(dUR2_ip1half_np1half-dUR2_im1half_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //Calcualte dA2
        dA2CenGrad=(dE_ijp1halfk_n-dE_ijm1halfk_n)/grid.dLocalGridOld[grid.nDTheta][0][j][0];
        dA2UpWindGrad=0.0;
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
        dA2=dV_ijk_np1half/dR_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA2CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
          
        //Calcualte dS2
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP_ijk_n+=grid.dLocalGridOld[grid.nQ1][i][j][k];
        #endif
        dS2=dP_ijk_n/(grid.dLocalGridOld[grid.nD][i][j][k]*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDTheta][0][j][0])
          *(dVSinTheta_ijp1halfk_np1half-dVSinTheta_ijm1halfk_np1half);
        
        //Calculate dS4
        dTGrad_im1half=(dT4_ijk_n-dT4_im1jk_n)/(grid.dLocalGridOld[grid.nDM][i][0][0]
          +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        dGrad_ip1half=-3.0*dRSq_ip1half*dT4_ijk_n/(8.0*parameters.dPi);/**\BC
          Missing grid.dLocalGridOld[grid.nT][i+1][0][0] using flux equals \f$2\sigma T^4\f$ at 
          surface.*/
        dGrad_im1half=dRhoAve_im1half*dR_im1half_4/(dKappa_im1halfjk_n*dRho_im1halfjk)
          *dTGrad_im1half;
        dS4=16.0*parameters.dPi*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *(dGrad_ip1half-dGrad_im1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //Calculate dS5
        dTGrad_jp1half=(dT4_ijp1k_n-dT4_ijk_n)/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
          +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        dTGrad_jm1half=(dT4_ijk_n-dT4_ijm1k_n)/(grid.dLocalGridOld[grid.nDTheta][0][j][0]
          +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;;
        dGrad_jp1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0]
          /(dKappa_ijp1halfk_n*dRho_ijp1halfk*dR_i_n)*dTGrad_jp1half;
        dGrad_jm1half=grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0]
          /(dKappa_ijm1halfk_n*dRho_ijm1halfk*dR_i_n)*dTGrad_jm1half;;
        dS5=(dGrad_jp1half-dGrad_jm1half)/(grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]
          *dR_i_n*grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //calculate new energy
        grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridOld[grid.nE][i][j][k]
          -time.dDeltat_np1half*(4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1+dS1)
          +dA2+dS2-4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS4+dS5));
        
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
