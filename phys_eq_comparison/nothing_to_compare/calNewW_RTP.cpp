void calNewW_RTP(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop){
  int i;
  int j;
  int k;
  int nIInt;
  int nJInt;
  int nKCen;
  double dU_ijkp1half_nm1half;
  double dV_ijkp1half_nm1half;
  double dU0i_nm1half;
  double dR_i_n;
  double dW_ijkp1half_nm1half;
  double dW_ijp1halfkp1half_nm1half;
  double dW_ijm1halfkp1half_nm1half;
  double dW_ip1halfjkp1half_nm1half;
  double dW_im1halfjkp1half_nm1half;
  double dW_ijkp1_nm1half;
  double dW_ijk_nm1half;
  double dDeltaPhi_kp1half;
  double dRho_ijkp1half_n;
  double dA1CenGrad;
  double dA1UpWindGrad;
  double dUmU0_ijkp1half_nm1half;
  double dA1;
  double dS1;
  double dA2CenGrad;
  double dA2UpWindGrad;
  double dA2;
  double dS2;
  double dA3CenGrad;
  double dA3UpWindGrad;
  double dA3;
  double dP_ijkp1_n;
  double dP_ijk_n;
  double dS3;
  
  //calculate new w
  for(i=grid.nStartUpdateExplicit[grid.nW][0];i<grid.nEndUpdateExplicit[grid.nW][0];i++){
    
    //calculate j of interface quantities
    nIInt=i+grid.nCenIntOffset[0];
    dU0i_nm1half=(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0])*0.5;
    dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])
      *0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nW][1];j<grid.nEndUpdateExplicit[grid.nW][1];j++){
      
      //calculate j of centered quantities
      nJInt=j+grid.nCenIntOffset[1];
      
      for(k=grid.nStartUpdateExplicit[grid.nW][2];k<grid.nEndUpdateExplicit[grid.nW][2];k++){
        
        //calculate k of interface quantities
        nKCen=k-grid.nCenIntOffset[2];
        
        //Calculate interpolated quantities
        dU_ijkp1half_nm1half=(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen+1]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen+1]
          +grid.dLocalGridOld[grid.nU][nIInt][j][nKCen]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen])*0.25;
        dV_ijkp1half_nm1half=(grid.dLocalGridOld[grid.nV][i][nJInt][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt][nKCen]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen])*0.25;
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
        dDeltaPhi_kp1half=(grid.dLocalGridOld[grid.nDPhi][0][0][nKCen]
          +grid.dLocalGridOld[grid.nDPhi][0][0][nKCen+1])*0.5;
        dRho_ijkp1half_n=(grid.dLocalGridOld[grid.nD][i][j][nKCen]
          +grid.dLocalGridOld[grid.nD][i][j][nKCen+1])*0.5;
        dP_ijkp1_n=grid.dLocalGridOld[grid.nP][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nQ0][i][j][nKCen+1]+grid.dLocalGridOld[grid.nQ1][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nQ2][i][j][nKCen+1];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][nKCen]+grid.dLocalGridOld[grid.nQ0][i][j][nKCen]
          +grid.dLocalGridOld[grid.nQ1][i][j][nKCen]+grid.dLocalGridOld[grid.nQ2][i][j][nKCen];
        
        //calculate A1
        dA1CenGrad=(dW_ip1halfjkp1half_nm1half-dW_im1halfjkp1half_nm1half)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        dA1UpWindGrad=0.0;
        dUmU0_ijkp1half_nm1half=dU_ijkp1half_nm1half-dU0i_nm1half;
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
        dA1=4.0*parameters.dPi*dR_i_n*dR_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *dUmU0_ijkp1half_nm1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate S1
        dS1=dU_ijkp1half_nm1half*dW_ijkp1half_nm1half/dR_i_n;
        
        //calculate dA2
        dA2CenGrad=(dW_ijp1halfkp1half_nm1half-dW_ijm1halfkp1half_nm1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        dA2UpWindGrad=0.0;
        if(dV_ijkp1half_nm1half<0.0){//moning in a negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j+1][nKCen]
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
        dA3CenGrad=(dW_ijkp1_nm1half-dW_ijk_nm1half)
          /dDeltaPhi_kp1half;
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
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA3UpWindGrad)
          /(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
        
        //calculate S3
        dS3=(dP_ijkp1_n-dP_ijk_n)/(dRho_ijkp1half_n*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*dDeltaPhi_kp1half);
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nW][i][j][k]=grid.dLocalGridOld[grid.nW][i][j][k]
          -time.dDeltat_n*(dA1+dS1+dA2+dS2+dA3+dS3);
      }
    }
  }
  
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nV][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nV][0][0];i++){
    
    //calculate j of interface quantities
    nIInt=i+grid.nCenIntOffset[0];
    dU0i_nm1half=(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0])*0.5;
    dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])
      *0.5;
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nV][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nV][0][1];j++){
      
      //calculate j of centered quantities
      nJInt=j+grid.nCenIntOffset[1];
      
      for(k=grid.nStartGhostUpdateExplicit[grid.nV][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nV][0][2];k++){
        
        //calculate k of interface quantities
        nKCen=k-grid.nCenIntOffset[2];
        
        //Calculate interpolated quantities
        dU_ijkp1half_nm1half=(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen+1]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen+1]
          +grid.dLocalGridOld[grid.nU][nIInt][j][nKCen]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen])*0.25;
        dV_ijkp1half_nm1half=(grid.dLocalGridOld[grid.nV][i][nJInt][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt][nKCen]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen])*0.25;
        dW_ijkp1half_nm1half=grid.dLocalGridOld[grid.nW][i][j][k];
        dW_ijp1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nW][i][j+1][k]
          +grid.dLocalGridOld[grid.nW][i][j][k])*0.5;
        dW_ijm1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nW][i][j-1][k]
          +grid.dLocalGridOld[grid.nW][i][j][k])*0.5;
        dW_ip1halfjkp1half_nm1half=grid.dLocalGridOld[grid.nW][i][j][k];/**\BC missing 
          grid.dLocalGridOld[grid.nW][i+1][j][k] assuming that the phi velocity at the outter most 
          interface is the same as the phi velocity in the center of the zone.*/
        dW_im1halfjkp1half_nm1half=(grid.dLocalGridOld[grid.nW][i-1][j][k]
          +grid.dLocalGridOld[grid.nW][i][j][k])*0.5;
        dW_ijkp1_nm1half=(grid.dLocalGridOld[grid.nW][i][j][k]
          +grid.dLocalGridOld[grid.nW][i][j][k+1])*0.5;
        dW_ijk_nm1half=(grid.dLocalGridOld[grid.nW][i][j][k]
          +grid.dLocalGridOld[grid.nW][i][j][k-1])*0.5;
        dDeltaPhi_kp1half=(grid.dLocalGridOld[grid.nDPhi][0][0][nKCen]
          +grid.dLocalGridOld[grid.nDPhi][0][0][nKCen+1])*0.5;
        dRho_ijkp1half_n=(grid.dLocalGridOld[grid.nD][i][j][nKCen]
          +grid.dLocalGridOld[grid.nD][i][j][nKCen+1])*0.5;
        dP_ijkp1_n=grid.dLocalGridOld[grid.nP][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nQ0][i][j][nKCen+1]+grid.dLocalGridOld[grid.nQ1][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nQ2][i][j][nKCen+1];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][nKCen]+grid.dLocalGridOld[grid.nQ0][i][j][nKCen]
          +grid.dLocalGridOld[grid.nQ1][i][j][nKCen]+grid.dLocalGridOld[grid.nQ2][i][j][nKCen];
          
        //calculate A1
        dA1CenGrad=(dW_ip1halfjkp1half_nm1half-dW_im1halfjkp1half_nm1half)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        dA1UpWindGrad=0.0;
        dUmU0_ijkp1half_nm1half=dU_ijkp1half_nm1half-dU0i_nm1half;
        if(dUmU0_ijkp1half_nm1half<0.0){//moving in a negative direction
          dA1UpWindGrad=dA1CenGrad;/**\BC missing grid.dLocalGridOld[grid.nW][i+1][j][k] in outter 
            most zone. This is needed to  calculate the upwind gradient for donnor cell. The 
            centered gradient is used instead when moving in the negative direction.*/
        }
        else{//moving in a positive direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j][k]
            -grid.dLocalGridOld[grid.nW][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=4.0*parameters.dPi*dR_i_n*dR_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *dUmU0_ijkp1half_nm1half*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
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
        dA3CenGrad=(dW_ijkp1_nm1half-dW_ijk_nm1half)/dDeltaPhi_kp1half;
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
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA3UpWindGrad)
          /(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
        
        //calculate S3
        dS3=(dP_ijkp1_n-dP_ijk_n)/(dRho_ijkp1half_n*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*dDeltaPhi_kp1half);
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nW][i][j][k]=grid.dLocalGridOld[grid.nW][i][j][k]
          -time.dDeltat_n*(dA1+dS1+dA2+dS2+dA3+dS3);
      }
    }
  }
  
  #if SEDOV==1
  
    //ghost region 1, inner ghost region in x1 direction
    for(i=grid.nStartGhostUpdateExplicit[grid.nV][1][0];
      i<grid.nEndGhostUpdateExplicit[grid.nV][1][0];i++){
      
      //calculate j of interface quantities
      nIInt=i+grid.nCenIntOffset[0];
      dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])
        *0.5;
      dU0i_nm1half=(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
        +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0])*0.5;
      for(j=grid.nStartGhostUpdateExplicit[grid.nV][1][1];
        j<grid.nEndGhostUpdateExplicit[grid.nV][1][1];j++){
        
        //calculate j of centered quantities
        nJInt=j+grid.nCenIntOffset[1];
        
        for(k=grid.nStartGhostUpdateExplicit[grid.nV][1][2];
          k<grid.nEndGhostUpdateExplicit[grid.nV][1][2];k++){
        
        //calculate k of interface quantities
        nKCen=k-grid.nCenIntOffset[2];
        
        //Calculate interpolated quantities
        dU_ijkp1half_nm1half=(grid.dLocalGridOld[grid.nU][nIInt][j][nKCen+1]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen+1]
          +grid.dLocalGridOld[grid.nU][nIInt][j][nKCen]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][nKCen])*0.25;
        ddV_ijkp1half_nm1half=(grid.dLocalGridOld[grid.nV][i][nJInt][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt][nKCen]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][nKCen])*0.25;
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
        dDeltaPhi_kp1half=(grid.dLocalGridOld[grid.nDPhi][0][0][nKCen]
          +grid.dLocalGridOld[grid.nDPhi][0][0][nKCen+1])*0.5;
        dRho_ijkp1half_n=(grid.dLocalGridOld[grid.nD][i][j][nKCen]
          +grid.dLocalGridOld[grid.nD][i][j][nKCen+1])*0.5;
        dP_ijkp1_n=grid.dLocalGridOld[grid.nP][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nQ0][i][j][nKCen+1]+grid.dLocalGridOld[grid.nQ1][i][j][nKCen+1]
          +grid.dLocalGridOld[grid.nQ2][i][j][nKCen+1];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][nKCen]+grid.dLocalGridOld[grid.nQ0][i][j][nKCen]
          +grid.dLocalGridOld[grid.nQ1][i][j][nKCen]+grid.dLocalGridOld[grid.nQ2][i][j][nKCen];
        
        //calculate A1
        dA1CenGrad=(dW_ip1halfjkp1half_nm1half-dW_im1halfjkp1half_nm1half)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        dA1UpWindGrad=0.0;
        dUmU0_ijkp1half_nm1half=dU_ijkp1half_nm1half-dU0i_nm1half;
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
        dA1=4.0*parameters.dPi*dR_i_n*dR_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *dUmU0_ijkp1half_nm1half*((1.0-parameters.dDonorFrac)*dA1CenGrad+parameters.dDonorFrac
          *dA1UpWindGrad);
        
        //calculate S1
        dS1=dU_ijkp1half_nm1half*dW_ijkp1half_nm1half/dR_i_n;
        
        //calculate dA2
        dA2CenGrad=(dW_ijp1halfkp1half_nm1half-dW_ijm1halfkp1half_nm1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        dA2UpWindGrad=0.0;
        if(dV_ijkp1half_nm1half<0.0){//moning in a negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j+1][nKCen]
            -grid.dLocalGridOld[grid.nW][i][j][k])/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        }
        else{//moving in a positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j][k]
            -grid.dLocalGridOld[grid.nW][i][j-1][k])/(grid.dLocalGridOld[grid.nDTheta][0][j-1][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        }
        dA2=dV_ijkp1half_nm1half/dR_i_n*((1.0-parameters.dDonorFrac)*dA2CenGrad
          +parameters.dDonorFrac*dA2UpWindGrad);
        
        //calculate S2
        dS2=dV_ijkp1half_nm1half*grid.dLocalGridOld[grid.nW][i][j][k]
          *grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]/dR_i_n;
        
        //calculate A3
        dA3CenGrad=(dW_ijkp1_nm1half-dW_ijk_nm1half)/dDeltaPhi_kp1half;
        dA3UpWindGrad=0.0;
        if(dW_ijkp1half_nm1half<0.0){//moving in a negative direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j][k+1]-grid.dLocalGridOld[grid.nW][i][j][k])
            /grid.dLocalGridOld[grid.nDPhi][0][0][nKCen+1];
        }
        else{//moving in a positive direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nW][i][j][k]-grid.dLocalGridOld[grid.nW][i][j][k-1])
            /grid.dLocalGridOld[grid.nDPhi][0][0][nKCen];
        }
        dA3=dW_ijkp1half_nm1half*((1.0-parameters.dDonorFrac)*dA3CenGrad+parameters.dDonorFrac
          *dA3UpWindGrad)/(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
        
        //calculate S3
        dS3=(dP_ijkp1_n-dP_ijk_n)/(dRho_ijkp1half_n*dR_i_n
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*dDeltaPhi_kp1half);
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nW][i][j][k]=grid.dLocalGridOld[grid.nW][i][j][k]
          -time.dDeltat_n*(dA1+dS1+dA2+dS2+dA3+dS3);
        }
      }
    }
  
  #endif
}
