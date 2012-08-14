void calNewV_RT(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop){
  
  int i;
  int j;
  int k;
  int nIInt;
  int nJCen;
  double dR_i_n;
  double dU0i_n;
  double dU_ijp1halfk_n;
  double dV_ip1halfjp1halfk_n;
  double dV_im1halfjp1halfk_n;
  double dV_ijp1halfk_n;
  double dV_ijp1k_n;
  double dV_ijk_n;
  double dDeltaTheta_jp1half;
  double dRho_ijp1halfk_n;
  double dA1CenGrad;
  double dA1UpWindGrad;
  double dU_U0_Diff;
  double dA1;
  double dS1;
  double dA2CenGrad;
  double dA2UpWindGrad;
  double dA2;
  double dP_ijp1k_n;
  double dP_ijk_n;
  double dS2;
  
  //calculate new v
  for(i=grid.nStartUpdateExplicit[grid.nV][0];i<grid.nEndUpdateExplicit[grid.nV][0];i++){
    
    //calculate j of interface quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt-1][0][0])*0.5;
    dU0i_n=0.5*(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
    
    for(j=grid.nStartUpdateExplicit[grid.nV][1];j<grid.nEndUpdateExplicit[grid.nV][1];j++){
      
      //calculate j of centered quantities
      nJCen=j-grid.nCenIntOffset[1];
      
      for(k=grid.nStartUpdateExplicit[grid.nV][2];k<grid.nEndUpdateExplicit[grid.nV][2];k++){
        
        //Calculate interpolated quantities
        dU_ijp1halfk_n=0.25*(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]);
        dV_ip1halfjp1halfk_n=0.5*(grid.dLocalGridOld[grid.nV][i+1][j][k]
          +grid.dLocalGridOld[grid.nV][i][j][k]);
        dV_im1halfjp1halfk_n=0.5*(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i-1][j][k]);
        dV_ijp1halfk_n=grid.dLocalGridOld[grid.nV][i][j][k];
        dV_ijp1k_n=(grid.dLocalGridOld[grid.nV][i][j+1][k]
          +grid.dLocalGridOld[grid.nV][i][j][k])*0.5;
        dV_ijk_n=(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i][j-1][k])*0.5;
        dDeltaTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0]
          +grid.dLocalGridOld[grid.nDTheta][0][nJCen][0])*0.5;
        dRho_ijp1halfk_n=(grid.dLocalGridOld[grid.nD][i][nJCen][k]
          +grid.dLocalGridOld[grid.nD][i][nJCen+1][k])*0.5;
        dP_ijp1k_n=grid.dLocalGridOld[grid.nP][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nQ0][i][nJCen+1][k]+grid.dLocalGridOld[grid.nQ1][i][nJCen+1][k];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][nJCen][k]+grid.dLocalGridOld[grid.nQ0][i][nJCen][k]
          +grid.dLocalGridOld[grid.nQ1][i][nJCen][k];
        
        //calculate derived quantities
        dU_U0_Diff=dU_ijp1halfk_n-dU0i_n;
        
        //calculate A1
        dA1CenGrad=(dV_ip1halfjp1halfk_n-dV_im1halfjp1halfk_n)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        dA1UpWindGrad=0.0;
        if(dU_U0_Diff<0.0){//moving in a negative direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nV][i+1][j][k]
            -grid.dLocalGridOld[grid.nV][i][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i+1][0][0])*2.0;
        }
        else{//moving in a positive direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
            -grid.dLocalGridOld[grid.nV][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=4.0*parameters.dPi*dR_i_n*dR_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *dU_U0_Diff*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA1CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate S1
        dS1=dU_ijp1halfk_n*dV_ijp1halfk_n/dR_i_n;
        
        //calculate dA2
        dA2CenGrad=(dV_ijp1k_n-dV_ijk_n)/dDeltaTheta_jp1half;
        dA2UpWindGrad=0.0;
        if(dV_ijp1halfk_n<0.0){//moning in a negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j+1][k]
            -grid.dLocalGridOld[grid.nV][i][j][k])/grid.dLocalGridOld[grid.nDTheta][0][j+1][0];
        }
        else{//moving in a positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
            -grid.dLocalGridOld[grid.nV][i][j-1][k])/grid.dLocalGridOld[grid.nDTheta][0][j][0];
        }
        dA2=dV_ijp1halfk_n/dR_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA2CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
        
        //calculate S2
        dS2=(dP_ijp1k_n-dP_ijk_n)/dDeltaTheta_jp1half/dRho_ijp1halfk_n/dR_i_n;
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nV][i][j][k]=grid.dLocalGridOld[grid.nV][i][j][k]
          -time.dDeltat_n*(dA1+dS1+dA2+dS2);
      }
    }
  }
  
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nV][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nV][0][0];i++){
    
    //calculate j of interface quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt-1][0][0])*0.5;
    dU0i_n=0.5*(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nV][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nV][0][1];j++){
      
      //calculate j of centered quantities
      nJCen=j-grid.nCenIntOffset[1];
      dDeltaTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0]
        +grid.dLocalGridOld[grid.nDTheta][0][nJCen][0])*0.5;
      
      for(k=grid.nStartGhostUpdateExplicit[grid.nV][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nV][0][2];k++){
        
        //Calculate interpolated quantities
        dU_ijp1halfk_n=0.25*(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]);
        dV_ip1halfjp1halfk_n=grid.dLocalGridOld[grid.nV][i][j][k];/**\BC 
          grid.dLocalGridOld[grid.nV][i+1][j+1][k] is missing*/
        dV_im1halfjp1halfk_n=0.5*(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i-1][j][k]);
        dV_ijp1halfk_n=grid.dLocalGridOld[grid.nV][i][j][k];
        dV_ijp1k_n=(grid.dLocalGridOld[grid.nV][i][j+1][k]
          +grid.dLocalGridOld[grid.nV][i][j][k])*0.5;
         dV_ijk_n=(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i][j-1][k])*0.5;
        dRho_ijp1halfk_n=(grid.dLocalGridOld[grid.nD][i][nJCen][k]
          +grid.dLocalGridOld[grid.nD][i][nJCen+1][k])*0.5;
        dP_ijp1k_n=grid.dLocalGridOld[grid.nP][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nQ0][i][nJCen+1][k]+grid.dLocalGridOld[grid.nQ1][i][nJCen+1][k];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][nJCen][k]+grid.dLocalGridOld[grid.nQ0][i][nJCen][k]
          +grid.dLocalGridOld[grid.nQ1][i][nJCen][k];
        
        //calculate derived quantities
        dU_U0_Diff=dU_ijp1halfk_n-dU0i_n;
        
        //calculate A1
        dA1CenGrad=(dV_ip1halfjp1halfk_n-dV_im1halfjp1halfk_n)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        dA1UpWindGrad=0.0;
        
        if(dU_U0_Diff<0.0){//moving in a negative direction
          dA1UpWindGrad=dA1CenGrad;/**\BC missing upwind gradient, using centred gradient instead*/
        }
        else{//moving in a positive direction*/
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
            -grid.dLocalGridOld[grid.nV][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=4.0*parameters.dPi*dR_i_n*dR_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]*dU_U0_Diff
          *((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA1CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate S1
        dS1=dU_ijp1halfk_n*dV_ijp1halfk_n/dR_i_n;
        
        //calculate dA2
        dA2CenGrad=(dV_ijp1k_n-dV_ijk_n)/dDeltaTheta_jp1half;
        dA2UpWindGrad=0.0;
        if(dV_ijp1halfk_n<0.0){//moning in a negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j+1][k]
            -grid.dLocalGridOld[grid.nV][i][j][k])/grid.dLocalGridOld[grid.nDTheta][0][j+1][0];
        }
        else{//moving in a positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
            -grid.dLocalGridOld[grid.nV][i][j-1][k])/grid.dLocalGridOld[grid.nDTheta][0][j][0];
        }
        dA2=dV_ijp1halfk_n/dR_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA2CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
        
        //calculate S2
        dS2=(dP_ijp1k_n-dP_ijk_n)/dDeltaTheta_jp1half/dRho_ijp1halfk_n/dR_i_n;
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nV][i][j][k]=grid.dLocalGridOld[grid.nV][i][j][k]
          -time.dDeltat_n*(dA1+dS1+dA2+dS2);
      }
    }
  }
#if SEDOV==1
    
    //calculate gost region 1 inner most ghost region in x1 direction
    for(i=grid.nStartGhostUpdateExplicit[grid.nV][1][0];
      i<grid.nEndGhostUpdateExplicit[grid.nV][1][0];i++){//nU0 needs to be 1D
      
      //calculate j of interface quantities
      nIInt=i+grid.nCenIntOffset[0];
      dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])
        *0.5;
      dU0i_n=0.5*(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
        +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
      
      for(j=grid.nStartGhostUpdateExplicit[grid.nV][1][1];
        j<grid.nEndGhostUpdateExplicit[grid.nV][1][1];j++){
        
        //calculate j of centered quantities
        nJCen=j-grid.nCenIntOffset[1];
        dDeltaTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0]
          +grid.dLocalGridOld[grid.nDTheta][0][nJCen][0])*0.5;
            
        for(k=grid.nStartGhostUpdateExplicit[grid.nV][1][2];
          k<grid.nEndGhostUpdateExplicit[grid.nV][1][2];k++){
          
          //Calculate interpolated quantities
          dU_ijp1halfk_n=0.25*(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
            +grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
            +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
            +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]);
          dV_ip1halfjp1halfk_n=0.5*(grid.dLocalGridOld[grid.nV][i+1][j][k]
            +grid.dLocalGridOld[grid.nV][i][j][k]);
          dV_im1halfjp1halfk_n=0.5*(grid.dLocalGridOld[grid.nV][i][j][k]
            +grid.dLocalGridOld[grid.nV][i-1][j][k]);
          dV_ijp1halfk_n=grid.dLocalGridOld[grid.nV][i][j][k];
          dV_ijp1k_n=(grid.dLocalGridOld[grid.nV][i][j+1][k]+grid.dLocalGridOld[grid.nV][i][j][k])
            *0.5;
          dV_ijk_n=(grid.dLocalGridOld[grid.nV][i][j][k]+grid.dLocalGridOld[grid.nV][i][j-1][k])
            *0.5;
          dRho_ijp1halfk_n=(grid.dLocalGridOld[grid.nD][i][nJCen][k]
            +grid.dLocalGridOld[grid.nD][i][nJCen+1][k])*0.5;
          dP_ijp1k_n=grid.dLocalGridOld[grid.nP][i][nJCen+1][k]
            +grid.dLocalGridOld[grid.nQ0][i][nJCen+1][k]
            +grid.dLocalGridOld[grid.nQ1][i][nJCen+1][k];
          dP_ijk_n=grid.dLocalGridOld[grid.nP][i][nJCen][k]
            +grid.dLocalGridOld[grid.nQ0][i][nJCen][k]
            +grid.dLocalGridOld[grid.nQ1][i][nJCen][k];
          
          //calculate derived quantities
          dU_U0_Diff=dU_ijp1halfk_n-dU0i_n;
          
          //calculate A1
          dA1CenGrad=(dV_ip1halfjp1halfk_n-dV_im1halfjp1halfk_n)
            /grid.dLocalGridOld[grid.nDM][i][0][0];
          dA1UpWindGrad=0.0;
          if(dU_U0_Diff<0.0){//moving in a negative direction
            dA1UpWindGrad=(grid.dLocalGridOld[grid.nV][i+1][j][k]
              -grid.dLocalGridOld[grid.nV][i][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
              +grid.dLocalGridOld[grid.nDM][i+1][0][0])*2.0;
          }
          else{//moving in a positive direction
            dA1UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
              -grid.dLocalGridOld[grid.nV][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
              +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
          }
          dA1=4.0*parameters.dPi*dR_i_n*dR_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
            *dU_U0_Diff*((1.0-parameters.dDonorFrac)*dA1CenGrad
            +parameters.dDonorFrac*dA1UpWindGrad);
          
          //calculate S1
          dS1=dU_ijp1halfk_n*dV_ijp1halfk_n/dR_i_n;
          
          //calculate dA2
          dA2CenGrad=(dV_ijp1k_n-dV_ijk_n)/dDeltaTheta_jp1half;
          dA2UpWindGrad=0.0;
          if(dV_ijp1halfk_n<0.0){//moning in a negative direction
            dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j+1][k]
              -grid.dLocalGridOld[grid.nV][i][j][k])
              /grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0];
          }
          else{//moving in a positive direction
            dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
              -grid.dLocalGridOld[grid.nV][i][j-1][k])
              /grid.dLocalGridOld[grid.nDTheta][0][nJCen][0];
          }
          dA2=dV_ijp1halfk_n/dR_i_n*((1.0-parameters.dDonorFrac)*dA2CenGrad
            +parameters.dDonorFrac*dA2UpWindGrad);
          
          //calculate S2
          dS2=(dP_ijp1k_n-dP_ijk_n)/dDeltaTheta_jp1half/dRho_ijp1halfk_n
            /dR_i_n;
          
          //calculate new velocity
          grid.dLocalGridNew[grid.nV][i][j][k]=grid.dLocalGridOld[grid.nV][i][j][k]
            -time.dDeltat_n*(dA1+dS1+dA2+dS2);
        }
      }
    }
#endif
}
