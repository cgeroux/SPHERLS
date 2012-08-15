void calNewE_R_AD(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop){
  int i;
  int j;
  int k;
  double dU_ijk_np1half;
  double dU0_i_np1half;
  int nIInt;
  double dE_ip1halfjk_n;
  double dE_im1halfjk_n;
  double dR_i_n;
  double dR_im1half_np1half;
  double dR_ip1half_np1half;
  double dRSq_i_n;
  double dA1;
  double dUR2_im1halfjk_n;
  double dUR2_ip1halfjk_n;
  double dP_ijk_n;
  double dS1;
  double dA1CenGrad;
  double dA1UpWindGrad;
  double dUmU0_ijk_Diff;
  double dR_im1half_n;
  double dR_ip1half_n;
  double dRSq_im1half_n;
  double dRSq_ip1half_n;
  
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
    dRSq_im1half_n=dR_im1half_n*dR_im1half_n;
    dRSq_ip1half_n=dR_ip1half_n*dR_ip1half_n;
    
    for(j=grid.nStartUpdateExplicit[grid.nE][1];j<grid.nEndUpdateExplicit[grid.nE][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nE][2];k<grid.nEndUpdateExplicit[grid.nE][2];k++){
        
        //calculate interpolated quantities
        dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dE_ip1halfjk_n=(grid.dLocalGridOld[grid.nE][i+1][j][k]+grid.dLocalGridOld[grid.nE][i][j][k])
          *0.5;
        dE_im1halfjk_n=(grid.dLocalGridOld[grid.nE][i][j][k]+grid.dLocalGridOld[grid.nE][i-1][j][k])
          *0.5;
        
        //Calcuate dA1
        dA1CenGrad=(dE_ip1halfjk_n-dE_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0];
        dUmU0_ijk_Diff=(dU_ijk_np1half-dU0_i_np1half);
        if(dUmU0_ijk_Diff<0.0){//moving in the negative direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nE][i+1][j][k]
            -grid.dLocalGridOld[grid.nE][i][j][k])/(grid.dLocalGridOld[grid.nDM][i+1][0][0]
            +grid.dLocalGridOld[grid.nDM][i][0][0])*2.0;
        }
        else{//moving in the postive direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k]
            -grid.dLocalGridOld[grid.nE][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=dUmU0_ijk_Diff*dRSq_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //source term in x-direction
        dUR2_im1halfjk_n=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dRSq_im1half_n;
        dUR2_ip1halfjk_n=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dRSq_ip1half_n;
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP_ijk_n+=grid.dLocalGridOld[grid.nQ0][i][j][k];
        #endif
        
        dS1=dP_ijk_n/grid.dLocalGridOld[grid.nD][i][j][k]*(dUR2_ip1halfjk_n-dUR2_im1halfjk_n)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //calculate new energy
        grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridOld[grid.nE][i][j][k]
        -time.dDeltat_np1half*4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][i][j][k]*(dA1+dS1);
        
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
  for(i=grid.nStartGhostUpdateExplicit[grid.nE][0][0];i<grid.nEndGhostUpdateExplicit[grid.nE][0][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
      +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
    dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])
      *0.5;
    dR_im1half_n=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_ip1half_n=grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dRSq_i_n=dR_i_n*dR_i_n;
    dRSq_im1half_n=dR_im1half_n*dR_im1half_n;
    dRSq_ip1half_n=dR_ip1half_n*dR_ip1half_n;
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nE][0][1];j<grid.nEndGhostUpdateExplicit[grid.nE][0][1];j++){
      for(k=grid.nStartGhostUpdateExplicit[grid.nE][0][2];k<grid.nEndGhostUpdateExplicit[grid.nE][0][2];k++){
        
        
        //calculate interpolated quantities
        dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dE_ip1halfjk_n=grid.dLocalGridOld[grid.nE][i][j][k]*0.5;
        dE_im1halfjk_n=(grid.dLocalGridOld[grid.nE][i][j][k]+grid.dLocalGridOld[grid.nE][i-1][j][k])
          *0.5;
        
        //Calcuate dA1
        dA1CenGrad=(dE_ip1halfjk_n-dE_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0];
        dUmU0_ijk_Diff=(dU_ijk_np1half-dU0_i_np1half);
        if(dUmU0_ijk_Diff<0.0){//moving in the negative direction
          dA1UpWindGrad=dA1CenGrad;
        }
        else{//moving in the postive direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k]
            -grid.dLocalGridOld[grid.nE][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=dUmU0_ijk_Diff*dRSq_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //source term in x-direction
        dUR2_im1halfjk_n=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dRSq_im1half_n;
        dUR2_ip1halfjk_n=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dRSq_ip1half_n;
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP_ijk_n+=grid.dLocalGridOld[grid.nQ0][i][j][k];
        #endif
        
        dS1=dP_ijk_n/grid.dLocalGridOld[grid.nD][i][j][k]*(dUR2_ip1halfjk_n-dUR2_im1halfjk_n)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //calculate new energy
        grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridOld[grid.nE][i][j][k]
        -time.dDeltat_np1half*4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][i][j][k]*(dA1+dS1);
        
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
    for(i=grid.nStartGhostUpdateExplicit[grid.nE][1][0];i<grid.nEndGhostUpdateExplicit[grid.nE][1][0];i++){
      
      for(j=grid.nStartGhostUpdateExplicit[grid.nE][1][1];j<grid.nEndGhostUpdateExplicit[grid.nE][1][1];j++){
        for(k=grid.nStartGhostUpdateExplicit[grid.nE][1][2];k<grid.nEndGhostUpdateExplicit[grid.nE][1][2];k++){
          grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridNew[grid.nE][i+1][j][k];
        }
      }
    }
  #endif
}
