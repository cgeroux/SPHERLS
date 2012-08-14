void calNewE_RT_AD(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop){
  int i;
  int j;
  int k;
  int nIInt;
  int nJInt;
  double dU_ijk_np1half;
  double dU0_i_np1half;
  double dE_ip1halfjk_n;
  double dE_im1halfjk_n;
  double dR_i_n;
  double dR_im1half_n;
  double dR_ip1half_n;
  double dRSq_i_n;
  double dV_ijk_np1half;
  double dE_ijp1halfk_n;
  double dE_ijm1halfk_n;
  double dVSinTheta_ijp1halfk_np1half;
  double dVSinTheta_ijm1halfk_np1half;
  double dA1CenGrad;
  double dA1UpWindGrad;
  double dU_U0_Diff;
  double dA1;
  double dUR2_im1half_np1half;
  double dUR2_ip1half_np1half;
  double dP;
  double dS1;
  double dA2CenGrad;
  double dA2UpWindGrad;
  double dA2;
  double dS2;
  
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
        dUR2_im1half_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dR_im1half_n*dR_im1half_n;
        dUR2_ip1half_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dR_ip1half_n*dR_ip1half_n;
        dP=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP+=grid.dLocalGridOld[grid.nQ0][i][j][k];
        #endif
        dS1=dP/grid.dLocalGridOld[grid.nD][i][j][k]*(dUR2_ip1half_np1half-dUR2_im1half_np1half)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        
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
        dP=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP+=grid.dLocalGridOld[grid.nQ1][i][j][k];
        #endif
        dS2=dP/(grid.dLocalGridOld[grid.nD][i][j][k]*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDTheta][0][j][0])
          *(dVSinTheta_ijp1halfk_np1half-dVSinTheta_ijm1halfk_np1half);
        
        //calculate new energy
        grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridOld[grid.nE][i][j][k]
          -time.dDeltat_np1half*(4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1+dS1)
          +dA2+dS2);
        
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
    dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])
      *0.5;
    dR_im1half_n=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_ip1half_n=grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dRSq_i_n=dR_i_n*dR_i_n;
    dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
      +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])
      *0.5;
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nE][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nE][0][1];j++){
      
      //calculate i for interface centered quantities
      nJInt=j+grid.nCenIntOffset[1];
      
      for(k=grid.nStartGhostUpdateExplicit[grid.nE][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nE][0][2];k++){
        
        
        //Calculate interpolated quantities
        dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dE_ip1halfjk_n=0.0;/**\BC grid.dLocalGridOld[grid.nE][i+1][j][k] is missing*/
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
        
        //Calcuate dA1
        dA1CenGrad=(dE_ip1halfjk_n-dE_im1halfjk_n)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        dU_U0_Diff=(dU_ijk_np1half-dU0_i_np1half);
        /**\BC grid.dLocalGridOld[grid.nDM][i+1][0][0] and grid.dLocalGridOld[grid.nE][i+1][j][k]
        missing using inner gradient for both*/
        dA1UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k]-grid.dLocalGridOld[grid.nE][i-1][j][k])
          /(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        dA1=dU_U0_Diff*dRSq_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA1CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate dS1
        dUR2_im1half_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dR_im1half_n*dR_im1half_n;
        dUR2_ip1half_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dR_ip1half_n*dR_ip1half_n;
        dP=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP+=grid.dLocalGridOld[grid.nQ0][i][j][k];
        #endif
        dS1=dP/grid.dLocalGridOld[grid.nD][i][j][k]*(dUR2_ip1half_np1half-dUR2_im1half_np1half)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        
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
        dP=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP+=grid.dLocalGridOld[grid.nQ1][i][j][k];
        #endif
        dS2=dP/(grid.dLocalGridOld[grid.nD][i][j][k]*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDTheta][0][j][0])
          *(dVSinTheta_ijp1halfk_np1half-dVSinTheta_ijm1halfk_np1half);
        
        //calculate new energy
        grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridOld[grid.nE][i][j][k]
          -time.dDeltat_np1half*(4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1+dS1)
          +dA2+dS2);
        
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
  
  #if SEDOV==1 //use zero P, E, and rho gradients
    for(i=grid.nStartGhostUpdateExplicit[grid.nE][1][0];
      i<grid.nEndGhostUpdateExplicit[grid.nE][1][0];i++){
      for(j=grid.nStartGhostUpdateExplicit[grid.nE][1][1];
        j<grid.nEndGhostUpdateExplicit[grid.nE][1][1];j++){
        for(k=grid.nStartGhostUpdateExplicit[grid.nE][1][2];
          k<grid.nEndGhostUpdateExplicit[grid.nE][1][2];k++){
          grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridNew[grid.nE][i+1][j][k];
        }
      }
    }
  #endif
}
