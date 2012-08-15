void calNewV_RTP(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop){
  
  int i;
  int j;
  int k;
  int nIInt;
  int nJCen;
  int nKInt;
  double dR_i_n;
  double dU_ijp1halfk_nm1half;
  double dU0i_nm1half;
  double dV_ip1halfjp1halfk_nm1half;
  double dV_im1halfjp1halfk_nm1half;
  double dV_ijp1halfk_nm1half;
  double dV_ijp1k_nm1half;
  double dV_ijk_nm1half;
  double dDeltaTheta_jp1half;
  double dRho_ijp1halfk_n;
  double dV_ijp1halfkm1half_nm1half;
  double dV_ijp1halfkp1half_nm1half;
  double dW_ijp1halfk_nm1half;
  double dP_ijp1k_n;
  double dP_ijk_n;
  double dA1CenGrad;
  double dA1UpWindGrad;
  double dU_U0_Diff;
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
  
  //calculate new v
  for(i=grid.nStartUpdateExplicit[grid.nV][0];i<grid.nEndUpdateExplicit[grid.nV][0];i++){
    
    //calculate j of interface quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])
      *0.5;
    dU0i_nm1half=0.5*(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
    
    for(j=grid.nStartUpdateExplicit[grid.nV][1];j<grid.nEndUpdateExplicit[grid.nV][1];j++){
      
      //calculate j of centered quantities
      nJCen=j-grid.nCenIntOffset[1];
      dDeltaTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0]
        +grid.dLocalGridOld[grid.nDTheta][0][nJCen][0])*0.5;
      
      for(k=grid.nStartUpdateExplicit[grid.nV][2];k<grid.nEndUpdateExplicit[grid.nV][2];k++){
        
        //calculate k of interface quantities
        nKInt=k+grid.nCenIntOffset[2];
        
        //Calculate interpolated quantities
        dU_ijp1halfk_nm1half=0.25*(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]);
        dV_ip1halfjp1halfk_nm1half=0.5*(grid.dLocalGridOld[grid.nV][i+1][j][k]
          +grid.dLocalGridOld[grid.nV][i][j][k]);
        dV_im1halfjp1halfk_nm1half=0.5*(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i-1][j][k]);
        dV_ijp1halfk_nm1half=grid.dLocalGridOld[grid.nV][i][j][k];
        dV_ijp1k_nm1half=(grid.dLocalGridOld[grid.nV][i][j+1][k]
          +grid.dLocalGridOld[grid.nV][i][j][k])*0.5;
        dV_ijk_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k]+grid.dLocalGridOld[grid.nV][i][j-1][k])
          *0.5;
        dRho_ijp1halfk_n=(grid.dLocalGridOld[grid.nD][i][nJCen][k]
          +grid.dLocalGridOld[grid.nD][i][nJCen+1][k])*0.5;
        dV_ijp1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k+1]
          +grid.dLocalGridOld[grid.nV][i][j][k])*0.5;
        dV_ijp1halfkm1half_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i][j][k-1])*0.5;
        dW_ijp1halfk_nm1half=0.25*(grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt]
          +grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt-1]
          +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt]
          +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt-1]);
        dP_ijp1k_n=grid.dLocalGridOld[grid.nP][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nQ0][i][nJCen+1][k]+grid.dLocalGridOld[grid.nQ1][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nQ2][i][nJCen+1][k];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][nJCen][k]+grid.dLocalGridOld[grid.nQ0][i][nJCen][k]
          +grid.dLocalGridOld[grid.nQ1][i][nJCen][k]+grid.dLocalGridOld[grid.nQ2][i][nJCen][k];
        
        //calculate dreived quantities
        dU_U0_Diff=dU_ijp1halfk_nm1half-dU0i_nm1half;
        
        //calculate A1
        dA1CenGrad=(dV_ip1halfjp1halfk_nm1half-dV_im1halfjp1halfk_nm1half)
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
        dS1=dU_ijp1halfk_nm1half*dV_ijp1halfk_nm1half/dR_i_n;
        
        //calculate dA2
        dA2CenGrad=(dV_ijp1k_nm1half-dV_ijk_nm1half)/dDeltaTheta_jp1half;
        dA2UpWindGrad=0.0;
        if(dV_ijp1halfk_nm1half<0.0){//moning in a negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j+1][k]
            -grid.dLocalGridOld[grid.nV][i][j][k])/grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0];
        }
        else{//moving in a positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
            -grid.dLocalGridOld[grid.nV][i][j-1][k])/grid.dLocalGridOld[grid.nDTheta][0][nJCen][0];
        }
        dA2=dV_ijp1halfk_nm1half/dR_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA2CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
        
        //calculate S2
        dS2=(dP_ijp1k_n-dP_ijk_n)/(dDeltaTheta_jp1half*dRho_ijp1halfk_n*dR_i_n);
        
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
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA3UpWindGrad)
          /(dR_i_n*grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]);
        
        //calculate S3
        dS3=-1.0*dW_ijp1halfk_nm1half*dW_ijp1halfk_nm1half
          *grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]/dR_i_n;
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nV][i][j][k]=grid.dLocalGridOld[grid.nV][i][j][k]-time.dDeltat_n
          *(dA1+dS1+dA2+dS2+dA3+dS3);
      }
    }
  }
  
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nV][0][0];i<grid.nEndGhostUpdateExplicit[grid.nV][0][0];
    i++){
    
    //calculate j of interface quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])
      *0.5;
    dU0i_nm1half=0.5*(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
          
    for(j=grid.nStartGhostUpdateExplicit[grid.nV][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nV][0][1];j++){
      
      //calculate j of centered quantities
      nJCen=j-grid.nCenIntOffset[1];
      dDeltaTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0]
        +grid.dLocalGridOld[grid.nDTheta][0][nJCen][0])*0.5;
        
      for(k=grid.nStartGhostUpdateExplicit[grid.nV][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nV][0][2];k++){
        
        //calculate k of interface quantities
        nKInt=k+grid.nCenIntOffset[2];
        
        //Calculate interpolated quantities
        dU_ijp1halfk_nm1half=0.25*(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]);
        dV_ip1halfjp1halfk_nm1half=grid.dLocalGridOld[grid.nV][i][j][k];/**\BC Assuming theta
          and phi velocities are the same at the surface of the star as just inside the star.*/
        dV_im1halfjp1halfk_nm1half=0.5*(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i-1][j][k]);
        dV_ijp1halfk_nm1half=grid.dLocalGridOld[grid.nV][i][j][k];
        dV_ijp1k_nm1half=(grid.dLocalGridOld[grid.nV][i][j+1][k]
          +grid.dLocalGridOld[grid.nV][i][j][k])*0.5;
        dV_ijk_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k]+grid.dLocalGridOld[grid.nV][i][j-1][k])
          *0.5;
        dRho_ijp1halfk_n=(grid.dLocalGridOld[grid.nD][i][nJCen][k]
          +grid.dLocalGridOld[grid.nD][i][nJCen+1][k])*0.5;
        dP_ijp1k_n=grid.dLocalGridOld[grid.nP][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nQ0][i][nJCen+1][k]+grid.dLocalGridOld[grid.nQ1][i][nJCen+1][k]
          +grid.dLocalGridOld[grid.nQ2][i][nJCen+1][k];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][nJCen][k]+grid.dLocalGridOld[grid.nQ0][i][nJCen][k]
          +grid.dLocalGridOld[grid.nQ1][i][nJCen][k]+grid.dLocalGridOld[grid.nQ2][i][nJCen][k];
        dV_ijp1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k+1]
          +grid.dLocalGridOld[grid.nV][i][j][k])*0.5;
        dV_ijp1halfkm1half_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k]
          +grid.dLocalGridOld[grid.nV][i][j][k-1])*0.5;
        dW_ijp1halfk_nm1half=0.25*(grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt]
          +grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt-1]
          +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt]
          +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt-1]);
        
        //calculate derived quantities
        dU_U0_Diff=dU_ijp1halfk_nm1half-dU0i_nm1half;
        
        //calculate A1
        dA1CenGrad=(dV_ip1halfjp1halfk_nm1half-dV_im1halfjp1halfk_nm1half)
          /grid.dLocalGridOld[grid.nDM][i][0][0];
        dA1UpWindGrad=0.0;
        if(dU_U0_Diff<0.0){//moving in a negative direction
          dA1UpWindGrad=dA1CenGrad;/**\BC ussing cetnered gradient for upwind gradient outside star
            at surface.*/
        }
        else{//moving in a positive direction*/
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
            -grid.dLocalGridOld[grid.nV][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=4.0*parameters.dPi*dR_i_n*dR_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]
          *dU_U0_Diff*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*dA1CenGrad
          +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate S1
        dS1=dU_ijp1halfk_nm1half*dV_ijp1halfk_nm1half/dR_i_n;
        
        //calculate dA2
        dA2CenGrad=(dV_ijp1k_nm1half-dV_ijk_nm1half)/dDeltaTheta_jp1half;
        dA2UpWindGrad=0.0;
        if(dV_ijp1halfk_nm1half<0.0){//moning in a negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j+1][k]
            -grid.dLocalGridOld[grid.nV][i][j][k])/grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0];
        }
        else{//moving in a positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
            -grid.dLocalGridOld[grid.nV][i][j-1][k])/grid.dLocalGridOld[grid.nDTheta][0][nJCen][0];
        }
        dA2=dV_ijp1halfk_nm1half/dR_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA2CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA2UpWindGrad);
        
        //calculate S2
        dS2=(dP_ijp1k_n-dP_ijk_n)/dDeltaTheta_jp1half/dRho_ijp1halfk_n/dR_i_n;
        
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
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nV][i][j][k]=grid.dLocalGridOld[grid.nV][i][j][k]
          -time.dDeltat_n*(dA1+dS1+dA2+dS2+dA3+dS3);
      }
    }
  }
  
  #if SEDOV==1
    //ghost region 1, innner most ghost region in x1 direction
    for(i=grid.nStartGhostUpdateExplicit[grid.nV][1][0];
      i<grid.nEndGhostUpdateExplicit[grid.nV][1][0];i++){
      
      //calculate j of interface quantities
      nIInt=i+grid.nCenIntOffset[0];
      dR_i_n=(grid.dLocalGridOld[grid.nR][nIInt][0][0]+grid.dLocalGridOld[grid.nR][nIInt-1][0][0])
        *0.5;
          
      for(j=grid.nStartGhostUpdateExplicit[grid.nV][1][1];
        j<grid.nEndGhostUpdateExplicit[grid.nV][1][1];j++){
        
        //calculate j of centered quantities
        nJCen=j-grid.nCenIntOffset[1];
        
        for(k=grid.nStartGhostUpdateExplicit[grid.nV][1][2];
          k<grid.nEndGhostUpdateExplicit[grid.nV][1][2];k++){
          
          //calculate k of interface quantities
          nKInt=k+grid.nCenIntOffset[2];
          
          //Calculate interpolated quantities
          dU_ijp1halfk_nm1half=0.25*(grid.dLocalGridOld[grid.nU][nIInt][nJCen][k]
            +grid.dLocalGridOld[grid.nU][nIInt][nJCen+1][k]
            +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen][k]
            +grid.dLocalGridOld[grid.nU][nIInt-1][nJCen+1][k]);
          dU0i_nm1half=0.5*(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
            +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0]);
          dV_ip1halfjp1halfk_nm1half=0.5*(grid.dLocalGridOld[grid.nV][i+1][j][k]
            +grid.dLocalGridOld[grid.nV][i][j][k]);
          dV_im1halfjp1halfk_nm1half=0.5*(grid.dLocalGridOld[grid.nV][i][j][k]
            +grid.dLocalGridOld[grid.nV][i-1][j][k]);
          dV_ijp1halfk_nm1half=grid.dLocalGridOld[grid.nV][i][j][k];
          dV_ijp1k_nm1half=(grid.dLocalGridOld[grid.nV][i][j+1][k]
            +grid.dLocalGridOld[grid.nV][i][j][k])*0.5;
          dV_ijk_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k]
            +grid.dLocalGridOld[grid.nV][i][j-1][k])*0.5;
          dDeltaTheta_jp1half=(grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0]
            +grid.dLocalGridOld[grid.nDTheta][0][nJCen][0])*0.5;
          dRho_ijp1halfk_n=(grid.dLocalGridOld[grid.nD][i][nJCen][k]
            +grid.dLocalGridOld[grid.nD][i][nJCen+1][k])*0.5;
          dV_ijp1halfkp1half_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k+1]
            +grid.dLocalGridOld[grid.nV][i][j][k])*0.5;
          dV_ijp1halfkm1half_nm1half=(grid.dLocalGridOld[grid.nV][i][j][k]
            +grid.dLocalGridOld[grid.nV][i][j][k-1])*0.5;
          dW_ijp1halfk_nm1half=0.25*(grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt]
            +grid.dLocalGridOld[grid.nW][i][nJCen+1][nKInt-1]
            +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt]
            +grid.dLocalGridOld[grid.nW][i][nJCen][nKInt-1]);
          dP_ijp1k_n=grid.dLocalGridOld[grid.nP][i][nJCen+1][k]
            +grid.dLocalGridOld[grid.nQ1][i][nJCen+1][k];
          dP_ijk_n=grid.dLocalGridOld[grid.nP][i][nJCen][k]
            +grid.dLocalGridOld[grid.nQ1][i][nJCen][k];
          
          //calculate derived quantities
          dU_U0_Diff=dU_ijp1halfk_nm1half-dU0i_nm1half;
          
          //calculate A1
          dA1CenGrad=(dV_ip1halfjp1halfk_nm1half-dV_im1halfjp1halfk_nm1half)
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
          dA1=4.0*parameters.dPi*dR_i_n*dR_i_n*grid.dLocalGridOld[grid.nDenAve][i][0][0]*dU_U0_Diff
            *((1.0-parameters.dDonorFrac)*dA1CenGrad+parameters.dDonorFrac*dA1UpWindGrad);
          
          //calculate S1
          dS1=dU_ijp1halfk_nm1half*dV_ijp1halfk_nm1half/dR_i_n;
          
          //calculate dA2
          dA2CenGrad=(dV_ijp1k_nm1half-dV_ijk_nm1half)/dDeltaTheta_jp1half;
          dA2UpWindGrad=0.0;
          if(dV_ijp1halfk_nm1half<0.0){//moning in a negative direction
            dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j+1][k]
              -grid.dLocalGridOld[grid.nV][i][j][k])
              /grid.dLocalGridOld[grid.nDTheta][0][nJCen+1][0];
          }
          else{//moving in a positive direction
            dA2UpWindGrad=(grid.dLocalGridOld[grid.nV][i][j][k]
              -grid.dLocalGridOld[grid.nV][i][j-1][k])
              /grid.dLocalGridOld[grid.nDTheta][0][nJCen][0];
          }
          dA2=dV_ijp1halfk_nm1half/dR_i_n*((1.0-parameters.dDonorFrac)*dA2CenGrad
            +parameters.dDonorFrac*dA2UpWindGrad);
          
          //calculate S2
          dS2=(dP_ijp1k_n-dP_ijk_n)/(dDeltaTheta_jp1half*dRho_ijp1halfk_n*dR_i_n);
          
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
          dA3=dW_ijp1halfk_nm1half*((1.0-parameters.dDonorFrac)*dA3CenGrad
            +parameters.dDonorFrac*dA3UpWindGrad)/(dR_i_n
            *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0]);
          
          //calculate S3
          dS3=-1.0*dW_ijp1halfk_nm1half*dW_ijp1halfk_nm1half
            *grid.dLocalGridOld[grid.nCotThetaIJp1halfK][0][j][0]/dR_i_n;
          
          //calculate new velocity
          grid.dLocalGridNew[grid.nV][i][j][k]=grid.dLocalGridOld[grid.nV][i][j][k]
            -time.dDeltat_n*(dA1+dS1+dA2+dS2+dA3+dS3);
        }
      }
    }
  #endif
}
