void calNewU_RTP(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop){
  int i;
  int j;
  int k;
  int nICen;
  int nJInt;
  int nKInt;
  double dRSq_ip1half_n;
  double dU_ip1jk_nm1half;
  double dU_ijk_nm1half;
  double dU_ip1halfjp1halfk_nm1half;
  double dU_ip1halfjm1halfk_nm1half;
  double dU_ip1halfjkp1half_nm1half;
  double dU_ip1halfjkm1half_nm1half;
  double dV_ip1halfjk_nm1half;
  double dW_ip1halfjk_nm1half;
  double dRho_ip1halfjk_n;
  double dRhoAve_ip1halfjk_n;
  double dP_ip1jk_n;
  double dP_ijk_n;
  double dUmU0_ijk_nm1half;
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
  double dS4;
  double dDonorFrac_ip1half;
  
  //calculate new u
  for(i=grid.nStartUpdateExplicit[grid.nU][0];i<grid.nEndUpdateExplicit[grid.nU][0];i++){
    
    //calculate i of centered quantities
    nICen=i-grid.nCenIntOffset[0];
    dRSq_ip1half_n=grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nR][i][0][0];
    dRhoAve_ip1halfjk_n=(grid.dLocalGridOld[grid.nDenAve][nICen+1][0][0]
      +grid.dLocalGridOld[grid.nDenAve][nICen][0][0])*0.5;
    dDonorFrac_ip1half=(grid.dLocalGridOld[grid.nDonorCellFrac][nICen+1][0][0]
          +grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nU][1];j<grid.nEndUpdateExplicit[grid.nU][1];j++){
      
      //calculate j of interface quantities
      nJInt=j+grid.nCenIntOffset[1];
      
      for(k=grid.nStartUpdateExplicit[grid.nU][2];k<grid.nEndUpdateExplicit[grid.nU][2];k++){
        
        //calculate k of interface quantities
        nKInt=k+grid.nCenIntOffset[2];
        
        //CALCULATE INTERPOLATED QUANTITIES
        dU_ip1jk_nm1half=(grid.dLocalGridOld[grid.nU][i+1][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
        dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i-1][j][k])*0.5;
        dU_ip1halfjp1halfk_nm1half=(grid.dLocalGridOld[grid.nU][i][j+1][k]
          +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
        dU_ip1halfjm1halfk_nm1half=0.5*(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i][j-1][k]);
        dU_ip1halfjkp1half_nm1half=0.5*(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k+1]);
        dU_ip1halfjkm1half_nm1half=0.5*(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k-1]);
        dV_ip1halfjk_nm1half=(grid.dLocalGridOld[grid.nV][nICen+1][nJInt][k]
          +grid.dLocalGridOld[grid.nV][nICen+1][nJInt-1][k]
          +grid.dLocalGridOld[grid.nV][nICen][nJInt][k]
          +grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k])*0.25;
        dW_ip1halfjk_nm1half=(grid.dLocalGridOld[grid.nW][nICen+1][j][nKInt]
          +grid.dLocalGridOld[grid.nW][nICen+1][j][nKInt-1]
          +grid.dLocalGridOld[grid.nW][nICen][j][nKInt]
          +grid.dLocalGridOld[grid.nW][nICen][j][nKInt-1])*0.25;
        dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
        dP_ip1jk_n=grid.dLocalGridOld[grid.nP][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nQ0][nICen+1][j][k];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][nICen][j][k]+grid.dLocalGridOld[grid.nQ0][nICen][j][k];
        
        //Calculate derived quantities
        dUmU0_ijk_nm1half=grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0];
        
        //Calculate dA1
        dA1CenGrad=(dU_ip1jk_nm1half-dU_ijk_nm1half)/(grid.dLocalGridOld[grid.nDM][nICen+1][0][0]
          +grid.dLocalGridOld[grid.nDM][nICen][0][0])*2.0;
        dA1UpWindGrad=0.0;
        if(dUmU0_ijk_nm1half<0.0){//moving from outside in
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nU][i+1][j][k]
            -grid.dLocalGridOld[grid.nU][i][j][k])/grid.dLocalGridOld[grid.nDM][nICen+1][0][0];
        }
        else{//moving from inside out
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU][i-1][j][k])/grid.dLocalGridOld[grid.nDM][nICen][0][0];
        }
        dA1=dUmU0_ijk_nm1half*((1.0-dDonorFrac_ip1half)*dA1CenGrad+dDonorFrac_ip1half
          *dA1UpWindGrad);
        
        //calculate dS1
        dS1=(dP_ip1jk_n-dP_ijk_n)/((grid.dLocalGridOld[grid.nDM][nICen+1][0][0]
          +grid.dLocalGridOld[grid.nDM][nICen][0][0])*dRho_ip1halfjk_n)*2.0;
        
        //Calculate dS4
        dS4=parameters.dG*grid.dLocalGridOld[grid.nM][i][0][0]/dRSq_ip1half_n;
        
        //Calculate dA2
        dA2CenGrad=(dU_ip1halfjp1halfk_nm1half-dU_ip1halfjm1halfk_nm1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        dA2UpWindGrad=0.0;
        if(dV_ip1halfjk_nm1half>0.0){//moving in positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU][i][j-1][k])/(grid.dLocalGridOld[grid.nDTheta][0][j][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;
        }
        else{//moving in negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j+1][k]
            -grid.dLocalGridOld[grid.nU][i][j][k])/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        }
        dA2=dV_ip1halfjk_nm1half*((1.0-dDonorFrac_ip1half)*dA2CenGrad
          +dDonorFrac_ip1half*dA2UpWindGrad)/grid.dLocalGridOld[grid.nR][i][0][0];
        
        //Calculate dS2
        dS2=-1.0*dV_ip1halfjk_nm1half*dV_ip1halfjk_nm1half/grid.dLocalGridOld[grid.nR][i][0][0];
        
        //Calculate dA3
        dA3CenGrad=(dU_ip1halfjkp1half_nm1half-dU_ip1halfjkm1half_nm1half)
          /grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dA3UpWindGrad=0.0;
        if(dW_ip1halfjk_nm1half>0.0){//moving in positive direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU][i][j][k-1])
            /(grid.dLocalGridOld[grid.nDPhi][0][0][k]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*2.0;
        }
        else{//moving in negative direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k+1]
            -grid.dLocalGridOld[grid.nU][i][j][k])
            /(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k])*2.0;
        }
        dA3=dW_ip1halfjk_nm1half*((1.0-dDonorFrac_ip1half)*dA3CenGrad
          +dDonorFrac_ip1half*dA3UpWindGrad)/(grid.dLocalGridOld[grid.nR][i][0][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
        
        //Calculate dS3
        dS3=-1.0*dW_ip1halfjk_nm1half*dW_ip1halfjk_nm1half
          /grid.dLocalGridOld[grid.nR][i][0][0];
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nU][i][j][k]=grid.dLocalGridOld[grid.nU][i][j][k]
          -time.dDeltat_n*(4.0*parameters.dPi*dRhoAve_ip1halfjk_n*dRSq_ip1half_n*(dA1+dS1)+dA2+dS2
          +dA3+dS3+dS4);
      }
    }
  }
  
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nU][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nU][0][0];i++){
    
    //calculate i of centered quantities
    nICen=i-grid.nCenIntOffset[0];
    dRSq_ip1half_n=grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nR][i][0][0];
    dRhoAve_ip1halfjk_n=(grid.dLocalGridOld[grid.nDenAve][nICen][0][0])*0.5;
    dDonorFrac_ip1half=grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0];
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nU][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nU][0][1];j++){
      
      //calculate j of interface quantities
      nJInt=j+grid.nCenIntOffset[1];
      
      for(k=grid.nStartGhostUpdateExplicit[grid.nU][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nU][0][2];k++){
        
        //calculate k of interface quantities
        nKInt=k+grid.nCenIntOffset[2];
        
        //CALCULATE INTERPOLATED QUANTITIES
        dU_ip1jk_nm1half=grid.dLocalGridOld[grid.nU][i][j][k];/**\BC
          missing grid.dLocalGridOld[grid.nU][i+1][j][k] in calculation of \f$u_{i+1,j,k}\f$ setting
          \f$u_{i+1,j,k}\f$=\f$u_{i+1/2,j,k}\f$.*/
        dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]+grid.dLocalGridOld[grid.nU][i-1][j][k])
          *0.5;/**\BC Missing grid.dLocalGridOld[grid.nD][i+1][j][k] in calculation of 
          \f$\rho_{i+1/2,j,k}\f$, setting it to zero.*/
        dU_ip1halfjp1halfk_nm1half=0.5*(grid.dLocalGridOld[grid.nU][i][j+1][k]
          +grid.dLocalGridOld[grid.nU][i][j][k]);
        dU_ip1halfjm1halfk_nm1half=0.5*(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i][j-1][k]);
        dU_ip1halfjkp1half_nm1half=0.5*(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k+1]);
        dU_ip1halfjkm1half_nm1half=0.5*(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k-1]);
        dV_ip1halfjk_nm1half=0.5*(grid.dLocalGridOld[grid.nV][nICen][nJInt][k]
          +grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k]);/**\BC assuming theta velocity is 
          constant across the surface.*/
        dW_ip1halfjk_nm1half=(grid.dLocalGridOld[grid.nW][nICen][j][nKInt]
          +grid.dLocalGridOld[grid.nW][nICen][j][nKInt-1])*0.5;/**\BC assuming phi velocity is
          constant across the surface.*/
        dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;/**\BC Missing 
          grid.dLocalGridOld[grid.nDenAve][nICen+1][0][0] in calculation of 
          \f$\langle\rho\rangle_{i+1/2}\f$ setting it to zero.*/
        dP_ijk_n=grid.dLocalGridOld[grid.nP][nICen][j][k]+grid.dLocalGridOld[grid.nQ0][nICen][j][k];
        dP_ip1jk_n=-1.0*dP_ijk_n;
        
        //Calculate dA1
        dA1CenGrad=(dU_ip1jk_nm1half-dU_ijk_nm1half)/((0.5+parameters.dAlpha
          +parameters.dAlphaExtra)*grid.dLocalGridOld[grid.nDM][nICen][0][0]);/**\BC Missing
          grid.dLocalGridOld[grid.nDM][nICen+1][0][0] in calculation of centered \f$A_1\f$
          gradient, setting it equal to \ref Parameters.dAlpha
          *grid.dLocalGridOld[grid.nDM][nICen][0][0].*/
        dA1UpWindGrad=0.0;
        if(grid.dLocalGridOld[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU0][i][0][0]<0.0){//moving from outside in
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
        
        //calculate dS1
        dS1=(dP_ip1jk_n-dP_ijk_n)/(grid.dLocalGridOld[grid.nDM][nICen][0][0]*(0.5+parameters.dAlpha
          +parameters.dAlphaExtra))/dRho_ip1halfjk_n;/**\BC Missing 
          grid.dLocalGridOld[grid.nDM][i+1][0][0] in calculation of \f$S_1\f$ using
          \ref Parameters.dAlpha *grid.dLocalGridOld[grid.nDM][nICen][0][0] instead.*/
        
        //Calculate dS4
        dS4=parameters.dG*grid.dLocalGridOld[grid.nM][i][0][0]/dRSq_ip1half_n;
        
        //Calculate dA2
        dA2CenGrad=(dU_ip1halfjp1halfk_nm1half-dU_ip1halfjm1halfk_nm1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        dA2UpWindGrad=0.0;
        if(dV_ip1halfjk_nm1half>0.0){//moving in positive direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU][i][j-1][k])/(grid.dLocalGridOld[grid.nDTheta][0][j][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;
        }
        else{//moving in negative direction
          dA2UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j+1][k]
            -grid.dLocalGridOld[grid.nU][i][j][k])/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
            +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
        }
        dA2=dV_ip1halfjk_nm1half*((1.0-dDonorFrac_ip1half)*dA2CenGrad+dDonorFrac_ip1half
          *dA2UpWindGrad)/grid.dLocalGridOld[grid.nR][i][0][0];
        
        //Calculate dS2
        dS2=-1.0*dV_ip1halfjk_nm1half*dV_ip1halfjk_nm1half/grid.dLocalGridOld[grid.nR][i][0][0];
        
        //Calculate dA3
        dA3CenGrad=(dU_ip1halfjkp1half_nm1half-dU_ip1halfjkm1half_nm1half)
          /grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dA3UpWindGrad=0.0;
        if(dW_ip1halfjk_nm1half>0.0){//moving in positive direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU][i][j][k-1])/(grid.dLocalGridOld[grid.nDPhi][0][0][k]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*2.0;
        }
        else{//moving in negative direction
          dA3UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k+1]
            -grid.dLocalGridOld[grid.nU][i][j][k])/(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
            +grid.dLocalGridOld[grid.nDPhi][0][0][k])*2.0;
        }
        dA3=dW_ip1halfjk_nm1half*((1.0-dDonorFrac_ip1half)*dA3CenGrad+dDonorFrac_ip1half
          *dA3UpWindGrad)/(grid.dLocalGridOld[grid.nR][i][0][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
        
        //Calculate dS3
        dS3=-1.0*dW_ip1halfjk_nm1half*dW_ip1halfjk_nm1half/grid.dLocalGridOld[grid.nR][i][0][0];
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nU][i][j][k]=grid.dLocalGridOld[grid.nU][i][j][k]
          -time.dDeltat_n*(4.0*parameters.dPi*dRhoAve_ip1halfjk_n*dRSq_ip1half_n*(dA1+dS1)+dA2+dS2
          +dA3+dS3+dS4);
      }
    }
  }

  #if SEDOV==1
    //calculate gost region 1 inner most ghost region in x1 direction
    for(i=grid.nStartGhostUpdateExplicit[grid.nU][1][0];
      i<grid.nEndGhostUpdateExplicit[grid.nU][1][0];i++){
      
      //calculate i of centered quantities
      nICen=i-grid.nCenIntOffset[0];
      dRSq_ip1half_n=grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nR][i][0][0];
      dRhoAve_ip1halfjk_n=(grid.dLocalGridOld[grid.nDenAve][nICen+1][0][0]
        +grid.dLocalGridOld[grid.nDenAve][nICen][0][0])*0.5;
      dDonorFrac_ip1half=(grid.dLocalGridOld[grid.nDonorCellFrac][nICen+1][0][0]
            +grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])*0.5;
      
      for(j=grid.nStartGhostUpdateExplicit[grid.nU][1][1];
        j<grid.nEndGhostUpdateExplicit[grid.nU][1][1];j++){
        
        //calculate j of interface quantities
        nJInt=j+grid.nCenIntOffset[1];
        
        for(k=grid.nStartGhostUpdateExplicit[grid.nU][1][2];
          k<grid.nEndGhostUpdateExplicit[grid.nU][1][2];k++){
            
          //calculate k of interface quantities
          nKInt=k+grid.nCenIntOffset[2];
          
          //CALCULATE INTERPOLATED QUANTITIES
          dU_ip1jk_nm1half=(grid.dLocalGridOld[grid.nU][i+1][j][k]
            +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
          dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
            +grid.dLocalGridOld[grid.nU][i-1][j][k])*0.5;
          dU_ip1halfjp1halfk_nm1half=(grid.dLocalGridOld[grid.nU][i][j+1][k]
            +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
          dU_ip1halfjm1halfk_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
            +grid.dLocalGridOld[grid.nU][i][j-1][k])*0.5;
          dU_ip1halfjkp1half_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
            +grid.dLocalGridOld[grid.nU][i][j][k+1])*0.5;
          dU_ip1halfjkm1half_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
            +grid.dLocalGridOld[grid.nU][i][j][k-1])*0.5;
          dV_ip1halfjk_nm1half=(grid.dLocalGridOld[grid.nV][nICen+1][nJInt][k]
            +grid.dLocalGridOld[grid.nV][nICen+1][nJInt-1][k]
            +grid.dLocalGridOld[grid.nV][nICen][nJInt][k]
            +grid.dLocalGridOld[grid.nV][nICen][nJInt-1][k])*0.25;
          dW_ip1halfjk_nm1half=(grid.dLocalGridOld[grid.nW][nICen+1][j][nKInt]
            +grid.dLocalGridOld[grid.nW][nICen+1][j][nKInt-1]
            +grid.dLocalGridOld[grid.nW][nICen][j][nKInt]
            +grid.dLocalGridOld[grid.nW][nICen][j][nKInt-1])*0.25;
          dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][nICen+1][j][k]
            +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
          dP_ip1jk_n=grid.dLocalGridOld[grid.nP][nICen+1][j][k]
            +grid.dLocalGridOld[grid.nQ0][nICen+1][j][k];
          dP_ijk_n=grid.dLocalGridOld[grid.nP][nICen][j][k]
            +grid.dLocalGridOld[grid.nQ0][nICen][j][k];
          
          //Calculate dA1
          dA1CenGrad=(dU_ip1jk_nm1half-dU_ijk_nm1half)/(grid.dLocalGridOld[grid.nDM][nICen+1][0][0]
            +grid.dLocalGridOld[grid.nDM][nICen][0][0])*2.0;
          dA1UpWindGrad=0.0;
          dUmU0_ijk_nm1half=grid.dLocalGridOld[grid.nU][i][j][k]
            -grid.dLocalGridOld[grid.nU0][i][0][0];
          if(dUmU0_ijk_nm1half<0.0){//moving from outside in
            dA1UpWindGrad=(grid.dLocalGridOld[grid.nU][i+1][j][k]
              -grid.dLocalGridOld[grid.nU][i][j][k])
              /grid.dLocalGridOld[grid.nDM][nICen+1][0][0];
          }
          else{//moving from inside out
            dA1UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
              -grid.dLocalGridOld[grid.nU][i-1][j][k])
              /grid.dLocalGridOld[grid.nDM][nICen][0][0];
          }
          dA1=dUmU0_ijk_nm1half*((1.0-parameters.dDonorFrac)*dA1CenGrad+parameters.dDonorFrac
            *dA1UpWindGrad);
          
          //calculate dS1
          dS1=(dP_ip1jk_n-dP_ijk_n)/((grid.dLocalGridOld[grid.nDM][nICen+1][0][0]
            +grid.dLocalGridOld[grid.nDM][nICen][0][0])*dRho_ip1halfjk_n)*2.0;
          
          //Calculate dS4
          dS4=parameters.dG*grid.dLocalGridOld[grid.nM][i][0][0]/dRSq_ip1half_n;
          
          //Calculate dA2
          dA2CenGrad=(dU_ip1halfjp1halfk_nm1half-dU_ip1halfjm1halfk_nm1half)
            /grid.dLocalGridOld[grid.nDTheta][0][j][0];
          dA2UpWindGrad=0.0;
          if(dV_ip1halfjk_nm1half>0.0){//moving in positive direction
            dA2UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
              -grid.dLocalGridOld[grid.nU][i][j-1][k])/(grid.dLocalGridOld[grid.nDTheta][0][j][0]
              +grid.dLocalGridOld[grid.nDTheta][0][j-1][0])*2.0;
          }
          else{//moving in negative direction
            dA2UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j+1][k]
              -grid.dLocalGridOld[grid.nU][i][j][k])/(grid.dLocalGridOld[grid.nDTheta][0][j+1][0]
              +grid.dLocalGridOld[grid.nDTheta][0][j][0])*2.0;
          }
          dA2=dV_ip1halfjk_nm1half*((1.0-parameters.dDonorFrac)*dA2CenGrad+parameters.dDonorFrac
            *dA2UpWindGrad)/grid.dLocalGridOld[grid.nR][i][0][0];
          
          //Calculate dS2
          dS2=-1.0*dV_ip1halfjk_nm1half*dV_ip1halfjk_nm1half/grid.dLocalGridOld[grid.nR][i][0][0];
          
          //Calculate dA3
          dA3CenGrad=(dU_ip1halfjkp1half_nm1half-dU_ip1halfjkm1half_nm1half)
            /grid.dLocalGridOld[grid.nDPhi][0][0][k];
          dA3UpWindGrad=0.0;
          if(dW_ip1halfjk_nm1half>0.0){//moving in positive direction
            dA3UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k]
              -grid.dLocalGridOld[grid.nU][i][j][k-1])/(grid.dLocalGridOld[grid.nDPhi][0][0][k]
              +grid.dLocalGridOld[grid.nDPhi][0][0][k-1])*2.0;
          }
          else{//moving in negative direction
            dA3UpWindGrad=(grid.dLocalGridOld[grid.nU][i][j][k+1]
              -grid.dLocalGridOld[grid.nU][i][j][k])/(grid.dLocalGridOld[grid.nDPhi][0][0][k+1]
              +grid.dLocalGridOld[grid.nDPhi][0][0][k])*2.0;
          }
          dA3=dW_ip1halfjk_nm1half*((1.0-parameters.dDonorFrac)*dA3CenGrad
            +parameters.dDonorFrac*dA3UpWindGrad)/(grid.dLocalGridOld[grid.nR][i][0][0]
            *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]);
          
          //Calculate dS3
          dS3=-1.0*dW_ip1halfjk_nm1half*dW_ip1halfjk_nm1half
            /grid.dLocalGridOld[grid.nR][i][0][0];
          
          //calculate new velocity
          grid.dLocalGridNew[grid.nU][i][j][k]=grid.dLocalGridOld[grid.nU][i][j][k]
            -time.dDeltat_n*(4.0*parameters.dPi*dRhoAve_ip1halfjk_n*dRSq_ip1half_n*(dA1+dS1)+dA2
             +dS2+dA3+dS3+dS4);
        }
      }
    }
  #endif

}
