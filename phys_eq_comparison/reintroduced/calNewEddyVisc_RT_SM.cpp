void calNewEddyVisc_RT_SM(Grid &grid, Parameters &parameters){
  
  int i;
  int j;
  int k;
  int nIInt;//i+1/2 index
  int nJInt;//j+1/2 index
  double d1;
  double d2;
  double d3;
  double d7;
  double d8;
  double d9;
  double dA;
  double dB;
  double dD;
  double dE;
  double dTerms;
  double dR_ip1half_np1half;
  double dR_im1half_np1half;
  double dR_i_np1half;
  double dRSq_i_np1half;
  double dU0_i_np1half;
  double dU_ip1halfjk_np1half;
  double dU_im1halfjk_np1half;
  double dU_ijp1halfk_np1half;
  double dU_ijm1halfk_np1half;
  double dU_ijk_np1half;
  double dV_ijk_np1half;
  double dV_ip1halfjk_np1half;
  double dV_im1halfjk_np1half;
  double dV_ijp1halfk_np1half;
  double dV_ijm1halfk_np1half;
  double dDelR_i_np1half;
  double dD_ijk_np1half;
  double dConstantSq=parameters.dEddyViscosity*parameters.dEddyViscosity/pow(2.0,0.5);
  double dLengthScaleSq;
  
  //main grid explicit
  for(i=grid.nStartUpdateExplicit[grid.nEddyVisc][0];
    i<grid.nEndUpdateExplicit[grid.nEddyVisc][0];i++){
      
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_np1half=(grid.dLocalGridNew[grid.nR][nIInt][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt][0][0])*0.5;
    dR_im1half_np1half=(grid.dLocalGridNew[grid.nR][nIInt-1][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt-1][0][0])*0.5;
    dR_i_np1half=(dR_ip1half_np1half+dR_im1half_np1half)*0.5;
    dRSq_i_np1half=dR_i_np1half*dR_i_np1half;
    dDelR_i_np1half=dR_ip1half_np1half-dR_im1half_np1half;
    dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
      +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nEddyVisc][1];
      j<grid.nEndUpdateExplicit[grid.nEddyVisc][1];j++){
      
      nJInt=j+grid.nCenIntOffset[1];
      
      for(k=grid.nStartUpdateExplicit[grid.nEddyVisc][2];
        k<grid.nEndUpdateExplicit[grid.nEddyVisc][2];k++){
        
        
        dLengthScaleSq=dR_i_np1half*dDelR_i_np1half*grid.dLocalGridOld[grid.nDTheta][0][j][0];
        
        //interpolate
        dU_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k];
        dU_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k];
        dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dU_ijp1halfk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k]+grid.dLocalGridNew[grid.nU][nIInt][j+1][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j+1][k])*0.25;
        dU_ijm1halfk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k]+grid.dLocalGridNew[grid.nU][nIInt][j-1][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j-1][k])*0.25;
        dV_ijk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k])*0.5;
        dV_ip1halfjk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k]+grid.dLocalGridNew[grid.nV][i+1][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i+1][nJInt-1][k])*0.25;
        dV_im1halfjk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k]+grid.dLocalGridNew[grid.nV][i-1][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i-1][nJInt-1][k])*0.25;
        dV_ijp1halfk_np1half=grid.dLocalGridNew[grid.nV][i][nJInt][k];
        dV_ijm1halfk_np1half=grid.dLocalGridNew[grid.nV][i][nJInt-1][k];
        
        //term 1
        d1=((dU_ip1halfjk_np1half-grid.dLocalGridNew[grid.nU0][nIInt][0][0])
          -(dU_im1halfjk_np1half-grid.dLocalGridNew[grid.nU0][nIInt-1][0][0]))
          /(dR_ip1half_np1half-dR_im1half_np1half);
        
        //term 2
        d2=1.0/dR_i_np1half*(dU_ijp1halfk_np1half-dU_ijm1halfk_np1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        
        //term 3
        d3=dV_ijk_np1half/dR_i_np1half;
        
        //term 7
        d7=(dV_ip1halfjk_np1half-dV_im1halfjk_np1half)/dDelR_i_np1half;
        
        //term 8
        d8=1.0/dR_i_np1half*(dV_ijp1halfk_np1half-dV_ijm1halfk_np1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        
        //term 9
        d9=(dU_ijk_np1half-dU0_i_np1half)/dR_i_np1half;
        
        //term A
        dA=2.0*d1*d1;
        
        //term B
        dB=(d2+d1-d3)*(d2-d3);
        
        //term D
        dD=d7*(d2+d7-d3);
        
        //term E
        dE=d8+d9;
        dE=2.0*dE*dE;
        
        dTerms=dA+dB+dD+dE;
        grid.dLocalGridNew[grid.nEddyVisc][i][j][k]=dConstantSq*dLengthScaleSq
          *grid.dLocalGridNew[grid.nD][i][j][k]*pow(dTerms,0.5);
      }
    }
  }
  
  //outter radial ghost cells,explicit
  for(i=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_np1half=(grid.dLocalGridNew[grid.nR][nIInt][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt][0][0])*0.5;
    dR_im1half_np1half=(grid.dLocalGridNew[grid.nR][nIInt-1][0][0]
      +grid.dLocalGridOld[grid.nR][nIInt-1][0][0])*0.5;
    dR_i_np1half=(dR_ip1half_np1half+dR_im1half_np1half)*0.5;
    dRSq_i_np1half=dR_i_np1half*dR_i_np1half;
    dDelR_i_np1half=dR_ip1half_np1half-dR_im1half_np1half;
    dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
      +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][1];j++){
      
      nJInt=j+grid.nCenIntOffset[1];
      
      for(k=grid.nStartGhostUpdateExplicit[grid.nEddyVisc][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nEddyVisc][0][2];k++){
        
        dLengthScaleSq=dR_i_np1half*dDelR_i_np1half*grid.dLocalGridOld[grid.nDTheta][0][j][0];
        
        //interpolate
        dU_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k];
        dU_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k];
        dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dU_ijp1halfk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k]+grid.dLocalGridNew[grid.nU][nIInt][j+1][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j+1][k])*0.25;
        dU_ijm1halfk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k]+grid.dLocalGridNew[grid.nU][nIInt][j-1][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j-1][k])*0.25;
        dV_ijk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k])*0.5;
        dV_ip1halfjk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k])*0.5;/**\BC assuming that theta velocity
          is constant across surface*/
        dV_im1halfjk_np1half=(grid.dLocalGridNew[grid.nV][i][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i][nJInt-1][k]+grid.dLocalGridNew[grid.nV][i-1][nJInt][k]
          +grid.dLocalGridNew[grid.nV][i-1][nJInt-1][k])*0.25;
        dV_ijp1halfk_np1half=grid.dLocalGridNew[grid.nV][i][nJInt][k];
        dV_ijm1halfk_np1half=grid.dLocalGridNew[grid.nV][i][nJInt-1][k];
        
        //term 1
        d1=((dU_ip1halfjk_np1half-grid.dLocalGridNew[grid.nU0][nIInt][0][0])
          -(dU_im1halfjk_np1half-grid.dLocalGridNew[grid.nU0][nIInt-1][0][0]))
          /(dR_ip1half_np1half-dR_im1half_np1half);
        
        //term 2
        d2=1.0/dR_i_np1half*(dU_ijp1halfk_np1half-dU_ijm1halfk_np1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        
        //term 3
        d3=dV_ijk_np1half/dR_i_np1half;
        
        //term 7
        d7=(dV_ip1halfjk_np1half-dV_im1halfjk_np1half)/dDelR_i_np1half;
        
        //term 8
        d8=1.0/dR_i_np1half*(dV_ijp1halfk_np1half-dV_ijm1halfk_np1half)
          /grid.dLocalGridOld[grid.nDTheta][0][j][0];
        
        //term 9
        d9=(dU_ijk_np1half-dU0_i_np1half)/dR_i_np1half;
        
        //term A
        dA=2.0*d1*d1;
        
        //term B
        dB=(d2+d1-d3)*(d2-d3);
        
        //term D
        dD=d7*(d2+d7-d3);
        
        //term E
        dE=d8+d9;
        dE=2.0*dE*dE;
        
        dTerms=dA+dB+dD+dE;
        grid.dLocalGridNew[grid.nEddyVisc][i][j][k]=dConstantSq*dLengthScaleSq
          *grid.dLocalGridNew[grid.nD][i][j][k]*pow(dTerms,0.5);
      }
    }
  }
}
