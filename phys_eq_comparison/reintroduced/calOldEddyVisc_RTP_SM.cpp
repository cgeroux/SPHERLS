void calOldEddyVisc_RTP_SM(Grid &grid, Parameters &parameters){
  int i;
  int j;
  int k;
  int nIInt;//i+1/2 index
  int nJInt;//j+1/2 index
  int nKInt;//k+1/2 index
  double d1;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  double d7;
  double d8;
  double d9;
  double d10;
  double d11;
  double d12;
  double d13;
  double d14;
  double dA;
  double dB;
  double dC;
  double dD;
  double dE;
  double dF;
  double dG;
  double dH;
  double dI;
  double dTerms;
  double dR_ip1half_np1half;
  double dR_im1half_np1half;
  double dR_i_np1half;
  double dU0_ip1half_np1half;
  double dU0_im1half_np1half;
  double dU0_i_np1half;
  double dU_ip1halfjk_np1half;
  double dU_im1halfjk_np1half;
  double dU_ijp1halfk_np1half;
  double dU_ijm1halfk_np1half;
  double dU_ijk_np1half;
  double dU_ijkp1half_np1half;
  double dU_ijkm1half_np1half;
  double dV_ijk_np1half;
  double dV_ip1halfjk_np1half;
  double dV_im1halfjk_np1half;
  double dV_ijp1halfk_np1half;
  double dV_ijm1halfk_np1half;
  double dV_ijkp1half_np1half;
  double dV_ijkm1half_np1half;
  double dW_ijk_np1half;
  double dW_ip1halfjk_np1half;
  double dW_im1halfjk_np1half;
  double dW_ijp1halfk_np1half;
  double dW_ijm1halfk_np1half;
  double dW_ijkp1half_np1half;
  double dW_ijkm1half_np1half;
  double dDelR_i_np1half;
  double dR_i_np1half_sq;
  double dConstantSq=parameters.dEddyViscosity*parameters.dEddyViscosity/pow(2.0,0.5);
  double dLengthScaleSq;
  
  //main grid explicit
  for(i=grid.nStartUpdateExplicit[grid.nEddyVisc][0];
    i<grid.nEndUpdateExplicit[grid.nEddyVisc][0];i++){
      
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_np1half=grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dR_im1half_np1half=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_i_np1half=(dR_ip1half_np1half+dR_im1half_np1half)*0.5;
    dR_i_np1half_sq=dR_i_np1half*dR_i_np1half;
    dDelR_i_np1half=dR_ip1half_np1half-dR_im1half_np1half;
    dU0_ip1half_np1half=grid.dLocalGridOld[grid.nU0][nIInt][0][0];
    dU0_im1half_np1half=grid.dLocalGridOld[grid.nU0][nIInt-1][0][0];
    dU0_i_np1half=(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nEddyVisc][1];
      j<grid.nEndUpdateExplicit[grid.nEddyVisc][1];j++){
      
      nJInt=j+grid.nCenIntOffset[1];
      
      for(k=grid.nStartUpdateExplicit[grid.nEddyVisc][2];
        k<grid.nEndUpdateExplicit[grid.nEddyVisc][2];k++){
        
        nKInt=k+grid.nCenIntOffset[2];
        
        dLengthScaleSq=dR_i_np1half_sq*dDelR_i_np1half*grid.dLocalGridOld[grid.nDTheta][0][j][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dLengthScaleSq=pow(dLengthScaleSq,0.666666666666666);
        
        //interpolate
        dU_ip1halfjk_np1half=grid.dLocalGridOld[grid.nU][nIInt][j][k];
        dU_im1halfjk_np1half=grid.dLocalGridOld[grid.nU][nIInt-1][j][k];
        dU_ijk_np1half=(grid.dLocalGridOld[grid.nU][nIInt][j][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][k])*0.5;
        dU_ijkp1half_np1half=(grid.dLocalGridOld[grid.nU][nIInt][j][k+1]
          +grid.dLocalGridOld[grid.nU][nIInt][j][k]+grid.dLocalGridOld[grid.nU][nIInt-1][j][k+1]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][k])*0.25;
        dU_ijkm1half_np1half=(grid.dLocalGridOld[grid.nU][nIInt][j][k]
          +grid.dLocalGridOld[grid.nU][nIInt][j][k-1]+grid.dLocalGridOld[grid.nU][nIInt-1][j][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][k-1])*0.25;
        dU_ijp1halfk_np1half=(grid.dLocalGridOld[grid.nU][nIInt][j][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][k]+grid.dLocalGridOld[grid.nU][nIInt][j+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j+1][k])*0.25;
        dU_ijm1halfk_np1half=(grid.dLocalGridOld[grid.nU][nIInt][j][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][k]+grid.dLocalGridOld[grid.nU][nIInt][j-1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j-1][k])*0.25;
        dV_ijk_np1half=(grid.dLocalGridOld[grid.nV][i][nJInt][k]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k])*0.5;
        dV_ip1halfjk_np1half=(grid.dLocalGridOld[grid.nV][i][nJInt][k]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k]+grid.dLocalGridOld[grid.nV][i+1][nJInt][k]
          +grid.dLocalGridOld[grid.nV][i+1][nJInt-1][k])*0.25;
        dV_im1halfjk_np1half=(grid.dLocalGridOld[grid.nV][i][nJInt][k]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k]+grid.dLocalGridOld[grid.nV][i-1][nJInt][k]
          +grid.dLocalGridOld[grid.nV][i-1][nJInt-1][k])*0.25;
        dV_ijp1halfk_np1half=grid.dLocalGridOld[grid.nV][i][nJInt][k];
        dV_ijm1halfk_np1half=grid.dLocalGridOld[grid.nV][i][nJInt-1][k];
        dV_ijkp1half_np1half=(grid.dLocalGridOld[grid.nV][i][nJInt][k+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt][k]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k])*0.25;
        dV_ijkm1half_np1half=(grid.dLocalGridOld[grid.nV][i][nJInt][k]
          +grid.dLocalGridOld[grid.nV][i][nJInt][k-1]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k-1])*0.25;
        dW_ijk_np1half=(grid.dLocalGridOld[grid.nW][i][j][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt-1])*0.5;
        dW_ip1halfjk_np1half=(grid.dLocalGridOld[grid.nW][i+1][j][nKInt]
          +grid.dLocalGridOld[grid.nW][i+1][j][nKInt-1]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt-1])*0.25;
        dW_im1halfjk_np1half=(grid.dLocalGridOld[grid.nW][i-1][j][nKInt]
          +grid.dLocalGridOld[grid.nW][i-1][j][nKInt-1]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt-1])*0.25;
        dW_ijp1halfk_np1half=(grid.dLocalGridOld[grid.nW][i][j+1][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j+1][nKInt-1]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt-1])*0.25;
        dW_ijm1halfk_np1half=(grid.dLocalGridOld[grid.nW][i][j-1][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j-1][nKInt-1]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt-1])*0.25;
        dW_ijkp1half_np1half=grid.dLocalGridOld[grid.nW][i][j][nKInt];
        dW_ijkm1half_np1half=grid.dLocalGridOld[grid.nW][i][j][nKInt-1];
        
        //term 1
        d1=((dU_ip1halfjk_np1half-dU0_ip1half_np1half)-(dU_im1halfjk_np1half-dU0_im1half_np1half))
          /(dR_ip1half_np1half-dR_im1half_np1half);
        
        //term 2
        d2=1.0/dR_i_np1half*(dU_ijp1halfk_np1half-dU_ijm1halfk_np1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        
        //term 3
        d3=dV_ijk_np1half/dR_i_np1half;
        
        //term 4
        d4=(dU_ijkp1half_np1half-dU_ijkm1half_np1half)/(dR_i_np1half
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //term 5
        d5=(dW_ip1halfjk_np1half-dW_im1halfjk_np1half)/(dR_ip1half_np1half-dR_im1half_np1half);
        
        //term 6
        d6=dW_ijk_np1half/dR_i_np1half;
        
        //term 7
        d7=(dV_ip1halfjk_np1half-dV_im1halfjk_np1half)/dDelR_i_np1half;
        
        //term 8
        d8=1.0/dR_i_np1half*(dV_ijp1halfk_np1half-dV_ijm1halfk_np1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        
        //term 9
        d9=(dU_ijk_np1half-dU0_i_np1half)/dR_i_np1half;
        
        //term 10
        d10=(dW_ijp1halfk_np1half-dW_ijm1halfk_np1half)/(dR_i_np1half
          *grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //term 11
        d11=(dV_ijkp1half_np1half-dV_ijkm1half_np1half)/(dR_i_np1half
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //term 12
        d12=dW_ijk_np1half*grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]/dR_i_np1half;
        
        //term 13
        d13=(dW_ijkp1half_np1half-dW_ijkm1half_np1half)/(dR_i_np1half
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //term 14
        d14=dV_ijk_np1half*grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]/dR_i_np1half;
        
        //term A
        dA=2.0*d1*d1;
        
        //term B
        dB=(d2+d1-d3)*(d2-d3);
        
        //term C
        dC=(d4+d5-d6)*(d4-d6);
        
        //term D
        dD=d7*(d2+d7-d3);
        
        //term E
        dE=d8+d9;
        dE=2.0*dE*dE;
        
        //term F
        dF=(d10+d11-d12)*(d11-d12);
        
        //term dG
        dG=d5*(d4+d5-d6);
        
        //term dH
        dH=d10*(d10+d11-d12);
        
        //term dI
        dI=d13+d14+d9;
        dI=2.0*dI*dI;
        
        dTerms=dA+dB+dC+dD+dE+dF+dG+dH+dI;
        grid.dLocalGridOld[grid.nEddyVisc][i][j][k]=dConstantSq*dLengthScaleSq
          *grid.dLocalGridOld[grid.nD][i][j][k]*pow(dTerms,0.5);
      }
    }
  }
  
  //inner radial ghost cells, set to zero
  for(i=0;i<grid.nNumGhostCells;i++){
    for(j=grid.nStartUpdateExplicit[grid.nEddyVisc][1];
      j<grid.nEndUpdateExplicit[grid.nEddyVisc][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nEddyVisc][2];
        k<grid.nEndUpdateExplicit[grid.nEddyVisc][2];k++){
        grid.dLocalGridOld[grid.nEddyVisc][i][j][k]=0.0;
      }
    }
  }
  
  //outter radial ghost cells,explicit
  for(i=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_np1half=grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dR_im1half_np1half=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_i_np1half=(dR_ip1half_np1half+dR_im1half_np1half)*0.5;
    dR_i_np1half_sq=dR_i_np1half*dR_i_np1half;
    dDelR_i_np1half=dR_ip1half_np1half-dR_im1half_np1half;
    dU0_ip1half_np1half=grid.dLocalGridOld[grid.nU0][nIInt][0][0];
    dU0_im1half_np1half=grid.dLocalGridOld[grid.nU0][nIInt-1][0][0];
    dU0_i_np1half=(grid.dLocalGridOld[grid.nU0][nIInt][0][0]
      +grid.dLocalGridOld[grid.nU0][nIInt-1][0][0])*0.5;
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][1];j++){
      
      nJInt=j+grid.nCenIntOffset[1];
      
      for(k=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][2];k++){
        
        nKInt=k+grid.nCenIntOffset[2];
        dLengthScaleSq=dR_i_np1half_sq*dDelR_i_np1half*grid.dLocalGridOld[grid.nDTheta][0][j][0]
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dLengthScaleSq=pow(dLengthScaleSq,0.666666666666666);
        
        //interpolate
        dU_ip1halfjk_np1half=grid.dLocalGridOld[grid.nU][nIInt][j][k];
        dU_im1halfjk_np1half=grid.dLocalGridOld[grid.nU][nIInt-1][j][k];
        dU_ijk_np1half=(grid.dLocalGridOld[grid.nU][nIInt][j][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][k])*0.5;
        dU_ijkp1half_np1half=(grid.dLocalGridOld[grid.nU][nIInt][j][k+1]
          +grid.dLocalGridOld[grid.nU][nIInt][j][k]+grid.dLocalGridOld[grid.nU][nIInt-1][j][k+1]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][k])*0.25;
        dU_ijkm1half_np1half=(grid.dLocalGridOld[grid.nU][nIInt][j][k]
          +grid.dLocalGridOld[grid.nU][nIInt][j][k-1]+grid.dLocalGridOld[grid.nU][nIInt-1][j][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][k-1])*0.25;
        dU_ijp1halfk_np1half=(grid.dLocalGridOld[grid.nU][nIInt][j][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][k]+grid.dLocalGridOld[grid.nU][nIInt][j+1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j+1][k])*0.25;
        dU_ijm1halfk_np1half=(grid.dLocalGridOld[grid.nU][nIInt][j][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j][k]+grid.dLocalGridOld[grid.nU][nIInt][j-1][k]
          +grid.dLocalGridOld[grid.nU][nIInt-1][j-1][k])*0.25;
        dV_ijk_np1half=(grid.dLocalGridOld[grid.nV][i][nJInt][k]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k])*0.5;
        dV_ip1halfjk_np1half=(grid.dLocalGridOld[grid.nV][i][nJInt][k]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k])*0.5;/**\BC assuming that theta velocity
          is constant across surface*/
        dV_im1halfjk_np1half=(grid.dLocalGridOld[grid.nV][i][nJInt][k]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k]+grid.dLocalGridOld[grid.nV][i-1][nJInt][k]
          +grid.dLocalGridOld[grid.nV][i-1][nJInt-1][k])*0.25;
        dV_ijp1halfk_np1half=grid.dLocalGridOld[grid.nV][i][nJInt][k];
        dV_ijm1halfk_np1half=grid.dLocalGridOld[grid.nV][i][nJInt-1][k];
        dV_ijkp1half_np1half=(grid.dLocalGridOld[grid.nV][i][nJInt][k+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt][k]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k+1]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k])*0.25;
        dV_ijkm1half_np1half=(grid.dLocalGridOld[grid.nV][i][nJInt][k]
          +grid.dLocalGridOld[grid.nV][i][nJInt][k-1]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k]
          +grid.dLocalGridOld[grid.nV][i][nJInt-1][k-1])*0.25;
        dW_ijk_np1half=(grid.dLocalGridOld[grid.nW][i][j][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt-1])*0.5;
        dW_ip1halfjk_np1half=(grid.dLocalGridOld[grid.nW][i][j][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt-1])*0.5;/**\BC assume phi velocity is constant
          across surface*/
        dW_im1halfjk_np1half=(grid.dLocalGridOld[grid.nW][i-1][j][nKInt]
          +grid.dLocalGridOld[grid.nW][i-1][j][nKInt-1]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt-1])*0.25;
        dW_ijp1halfk_np1half=(grid.dLocalGridOld[grid.nW][i][j+1][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j+1][nKInt-1]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt-1])*0.25;
        dW_ijm1halfk_np1half=(grid.dLocalGridOld[grid.nW][i][j-1][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j-1][nKInt-1]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt]
          +grid.dLocalGridOld[grid.nW][i][j][nKInt-1])*0.25;
        dW_ijkp1half_np1half=grid.dLocalGridOld[grid.nW][i][j][nKInt];
        dW_ijkm1half_np1half=grid.dLocalGridOld[grid.nW][i][j][nKInt-1];
        
        //term 1
        d1=((dU_ip1halfjk_np1half-dU0_ip1half_np1half)-(dU_im1halfjk_np1half-dU0_im1half_np1half))
          /(dR_ip1half_np1half-dR_im1half_np1half);
        
        //term 2
        d2=1.0/dR_i_np1half*(dU_ijp1halfk_np1half-dU_ijm1halfk_np1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        
        //term 3
        d3=dV_ijk_np1half/dR_i_np1half;
        
        //term 4
        d4=(dU_ijkp1half_np1half-dU_ijkm1half_np1half)/(dR_i_np1half
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //term 5
        d5=(dW_ip1halfjk_np1half-dW_im1halfjk_np1half)/(dR_ip1half_np1half-dR_im1half_np1half);
        
        //term 6
        d6=dW_ijk_np1half/dR_i_np1half;
        
        //term 7
        d7=(dV_ip1halfjk_np1half-dV_im1halfjk_np1half)/dDelR_i_np1half;
        
        //term 8
        d8=1.0/dR_i_np1half*(dV_ijp1halfk_np1half-dV_ijm1halfk_np1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        
        //term 9
        d9=(dU_ijk_np1half-dU0_i_np1half)/dR_i_np1half;
        
        //term 10
        d10=(dW_ijp1halfk_np1half-dW_ijm1halfk_np1half)/(dR_i_np1half
          *grid.dLocalGridOld[grid.nDTheta][0][j][0]);
        
        //term 11
        d11=(dV_ijkp1half_np1half-dV_ijkm1half_np1half)/(dR_i_np1half
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //term 12
        d12=dW_ijk_np1half*grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]/dR_i_np1half;
        
        //term 13
        d13=(dW_ijkp1half_np1half-dW_ijkm1half_np1half)/(dR_i_np1half
          *grid.dLocalGridOld[grid.nSinThetaIJK][0][j][0]*grid.dLocalGridOld[grid.nDPhi][0][0][k]);
        
        //term 14
        d14=dV_ijk_np1half*grid.dLocalGridOld[grid.nCotThetaIJK][0][j][0]/dR_i_np1half;
        
        //term A
        dA=2.0*d1*d1;
        
        //term B
        dB=(d2+d1-d3)*(d2-d3);
        
        //term C
        dC=(d4+d5-d6)*(d4-d6);
        
        //term D
        dD=d7*(d2+d7-d3);
        
        //term E
        dE=d8+d9;
        dE=2.0*dE*dE;
        
        //term dF
        dF=(d10+d11-d12)*(d11-d12);
        
        //term dG
        dG=d5*(d4+d5-d6);
        
        //term dH
        dH=d10*(d10+d11-d12);
        
        //term dI
        dI=d13+d14+d9;
        dI=2.0*dI*dI;
        
        dTerms=dA+dB+dC+dD+dE+dF+dG+dH+dI;
        grid.dLocalGridOld[grid.nEddyVisc][i][j][k]=dConstantSq*dLengthScaleSq
          *grid.dLocalGridOld[grid.nD][i][j][k]*pow(dTerms,0.5);
      }
    }
  }
}
