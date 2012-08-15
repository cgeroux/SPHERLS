void calNewU0_RTP(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop,MessPass &messPass){
  
  /**\todo
    At some point I will likely want to make this funciton compatiable with a 3D domain 
    decomposition instead of a purely radial domain decomposition.
  */
  
  //post a blocking recieve from inner radial neighbour
  int i;
  for(i=0;i<procTop.nNumRadialNeighbors;i++){
    if(procTop.nCoords[procTop.nRank][0]>procTop.nCoords[procTop.nRadialNeighborRanks[i]][0]){/*
      if current processor has a radial neighbor inside post a recieve*/
      MPI::COMM_WORLD.Recv(grid.dLocalGridNew,1
        ,messPass.typeRecvNewVar[procTop.nRadialNeighborNeighborIDs[i]][grid.nU0]
        ,procTop.nRadialNeighborRanks[i],2);
    }
  }
  
  #if SEDOV==1
    //calculate grid velocities for local grid at inner most ghost region
    for(i=grid.nStartGhostUpdateExplicit[grid.nU0][1][0];
      i<grid.nEndGhostUpdateExplicit[grid.nU0][1][0];i++){//nU0 needs to be 1D
      grid.dLocalGridNew[grid.nU0][i][0][0]=grid.dLocalGridNew[grid.nU][i][0][0];
    }
  #endif
  
  //calculate grid velocities for local grid
  double dCSum;
  double dARhoSum;
  int nICen;
  int j;
  int k;
  int nJInt;
  int nKInt;
  double dR_im1half_np1half;
  double dR_ip1half_np1half;
  double dRSq_im1half_np1half;
  double dRSq_ip1half_np1half;
  double d1half_RSq_ip1half_m_RSq_im1half;
  double dTemp;
  double dA_ip1halfjk_np1half;
  double dA_im1halfjk_np1half;
  double dA_ijm1halfk_np1half;
  double dA_ijp1halfk_np1half;
  double dA_ijkm1half_np1half;
  double dA_ijkp1half_np1half;
  double dRho_im1halfjk_np1half;
  double dRho_ip1halfjk_np1half;
  double dRho_ijm1halfk_np1half;
  double dRho_ijp1halfk_np1half;
  double dRho_ijkm1half_np1half;
  double dRho_ijkp1half_np1half;
  double dUmU0_ip1halfjk_nm1half;
  double dUmU0_im1halfjk_np1half;
  double dRho_cen_im1half;
  double dRho_upwind_im1half;
  double dRho_cen_ip1half;
  double dRho_upwind_ip1half;
  double dRho_cen_jm1half;
  double dRho_upwind_jm1half;
  double dRho_cen_jp1half;
  double dRho_upwind_jp1half;
  double dRho_cen_km1half;
  double dRho_upwind_km1half;
  double dRho_cen_kp1half;
  double dRho_upwind_kp1half;
  double dDonorFrac_ip1half;
  double dDonorFrac_im1half;
  
  for(i=grid.nStartUpdateExplicit[grid.nU0][0];i<grid.nEndUpdateExplicit[grid.nU0][0];i++){/*nU0
    needs to be 1D*/
    dCSum=0.0;
    dARhoSum=0.0;
    
    //calculate i of centered quantities
    nICen=i-grid.nCenIntOffset[0];
    
    dR_im1half_np1half=grid.dLocalGridOld[grid.nR][i-1][0][0];
    dR_ip1half_np1half=grid.dLocalGridOld[grid.nR][i][0][0];
    dRSq_im1half_np1half=dR_im1half_np1half*dR_im1half_np1half;
    dRSq_ip1half_np1half=dR_ip1half_np1half*dR_ip1half_np1half;
    d1half_RSq_ip1half_m_RSq_im1half=0.5*(dRSq_ip1half_np1half-dRSq_im1half_np1half);
    dDonorFrac_ip1half=(grid.dLocalGridOld[grid.nDonorCellFrac][nICen+1][0][0]
      +grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])*0.5;
    dDonorFrac_im1half=(grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0]
      +grid.dLocalGridOld[grid.nDonorCellFrac][nICen-1][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nU][1];j<grid.nEndUpdateExplicit[grid.nU][1];j++){
      
      //calculate j of interface quantities
      nJInt=j+grid.nCenIntOffset[1];
      
      for(k=grid.nStartUpdateExplicit[grid.nU][2];k<grid.nEndUpdateExplicit[grid.nU][2];k++){
        
        nKInt=k+grid.nCenIntOffset[2];
        
        dTemp=grid.dLocalGridOld[grid.nDCosThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dA_im1halfjk_np1half=dRSq_im1half_np1half*dTemp;
        dA_ip1halfjk_np1half=dRSq_ip1half_np1half*dTemp;
        dA_ijm1halfk_np1half=d1half_RSq_ip1half_m_RSq_im1half
          *grid.dLocalGridOld[grid.nDPhi][0][0][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt-1][0];
         dA_ijp1halfk_np1half=d1half_RSq_ip1half_m_RSq_im1half
          *grid.dLocalGridOld[grid.nDPhi][0][0][k]
          *grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][nJInt][0];
        dA_ijkm1half_np1half=d1half_RSq_ip1half_m_RSq_im1half
          *grid.dLocalGridOld[grid.nDTheta][0][j][0];
        dA_ijkp1half_np1half=dA_ijkm1half_np1half;
        
        //calculate difference between U and U0
        dUmU0_ip1halfjk_nm1half=grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0];
        dUmU0_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][i-1][j][k]
          -grid.dLocalGridNew[grid.nU0][i-1][0][0];
        
        //calculate rho at i-1/2, not time centered
        dRho_cen_im1half=(grid.dLocalGridOld[grid.nD][nICen][j][k]
          +grid.dLocalGridOld[grid.nD][nICen-1][j][k])*0.5;
        if(dUmU0_im1halfjk_np1half<0.0){//moving from outside in
          dRho_upwind_im1half=grid.dLocalGridOld[grid.nD][nICen][j][k];
        }
        else{//moving from inside out
          dRho_upwind_im1half=grid.dLocalGridOld[grid.nD][nICen-1][j][k];
        }
        dRho_im1halfjk_np1half=((1.0-dDonorFrac_im1half)*dRho_cen_im1half+dDonorFrac_im1half
          *dRho_upwind_im1half);
        
        //calculate rho at i+1/2, not time centered
        dRho_cen_ip1half=(grid.dLocalGridOld[grid.nD][nICen][j][k]
          +grid.dLocalGridOld[grid.nD][nICen+1][j][k])*0.5;
        if(dUmU0_ip1halfjk_nm1half<0.0){//moving from outside in
          dRho_upwind_ip1half=grid.dLocalGridOld[grid.nD][nICen+1][j][k];
        }
        else{//moving from inside out
          dRho_upwind_ip1half=grid.dLocalGridOld[grid.nD][nICen][j][k];
        }
        dRho_ip1halfjk_np1half=((1.0-dDonorFrac_ip1half)*dRho_cen_ip1half+dDonorFrac_ip1half
          *dRho_upwind_ip1half);
        
        //calculte rho at j-1/2
        dRho_cen_jm1half=(grid.dLocalGridOld[grid.nD][nICen][j-1][k]
          +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
        if(grid.dLocalGridNew[grid.nV][nICen][nJInt-1][k]<0.0){
          dRho_upwind_jm1half=grid.dLocalGridOld[grid.nD][nICen][j][k];
        }
        else{
          dRho_upwind_jm1half=grid.dLocalGridOld[grid.nD][nICen][j-1][k];
        }
        dRho_ijm1halfk_np1half=((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])
          *dRho_cen_jm1half+grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0]
          *dRho_upwind_jm1half);
        
        //calculte rho at j+1/2
        dRho_cen_jp1half=(grid.dLocalGridOld[grid.nD][nICen][j+1][k]
          +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
        if(grid.dLocalGridNew[grid.nV][nICen][nJInt][k]<0.0){
          dRho_upwind_jp1half=grid.dLocalGridOld[grid.nD][nICen][j+1][k];
        }
        else{
          dRho_upwind_jp1half=grid.dLocalGridOld[grid.nD][nICen][j][k];
        }
        dRho_ijp1halfk_np1half=((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])
          *dRho_cen_jp1half+grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0]
          *dRho_upwind_jp1half);
        
        //calculte rho at k-1/2
        dRho_cen_km1half=(grid.dLocalGridOld[grid.nD][nICen][j][k-1]
          +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
        if(grid.dLocalGridNew[grid.nW][nICen][j][nKInt-1]<0.0){
          dRho_upwind_km1half=grid.dLocalGridOld[grid.nD][nICen][j][k];
        }
        else{
          dRho_upwind_km1half=grid.dLocalGridOld[grid.nD][nICen][j][k-1];
        }
        dRho_ijkm1half_np1half=((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])
          *dRho_cen_km1half+grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0]
          *dRho_upwind_km1half);
        
        //calculte rho at k+1/2
        dRho_cen_kp1half=(grid.dLocalGridOld[grid.nD][nICen][j][k+1]
          +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
        if(grid.dLocalGridNew[grid.nW][nICen][j][nKInt]<0.0){
          dRho_upwind_kp1half=grid.dLocalGridOld[grid.nD][nICen][j][k+1];
        }
        else{
          dRho_upwind_kp1half=grid.dLocalGridOld[grid.nD][nICen][j][k];
        }
        dRho_ijkp1half_np1half=((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])
          *dRho_cen_kp1half+grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0]
          *dRho_upwind_kp1half);
        
        dCSum+=(grid.dLocalGridNew[grid.nU][i-1][j][k]-grid.dLocalGridNew[grid.nU0][i-1][0][0])
          *dA_im1halfjk_np1half*dRho_im1halfjk_np1half
          -grid.dLocalGridNew[grid.nU][i][j][k]*dA_ip1halfjk_np1half*dRho_ip1halfjk_np1half/*
          +grid.dLocalGridNew[grid.nV][nICen][nJInt-1][k]*dA_ijm1halfk_np1half
          *dRho_ijm1halfk_np1half
          -grid.dLocalGridNew[grid.nV][nICen][nJInt][k]*dA_ijp1halfk_np1half
          *dRho_ijp1halfk_np1half
          +grid.dLocalGridNew[grid.nW][nICen][j][k-1]*dA_ijkm1half_np1half
          *dRho_ijkm1half_np1half-grid.dLocalGridNew[grid.nW][nICen][j][k]*dA_ijkp1half_np1half
          *dRho_ijkp1half_np1half*/;
        dARhoSum+=dA_ip1halfjk_np1half*dRho_ip1halfjk_np1half;
      }
    }
    grid.dLocalGridNew[grid.nU0][i][0][0]=-1.0*dCSum/dARhoSum;
  }
  
  //post a blocking send to outer radial neighbour
  int nNumRecieves=0;
  for(i=0;i<procTop.nNumRadialNeighbors;i++){
    if(procTop.nCoords[procTop.nRank][0]<procTop.nCoords[procTop.nRadialNeighborRanks[i]][0]){/*
      if current processor has a radial neighbor at inside post a recieve*/
      MPI::COMM_WORLD.Send(grid.dLocalGridNew,1
        ,messPass.typeSendNewVar[procTop.nRadialNeighborNeighborIDs[i]][grid.nU0]
        ,procTop.nRadialNeighborRanks[i],2);
      nNumRecieves++;
    }
  }
  
  //post a non-blocking recieve for outer radial neighbour
  MPI::Request *requestTempRecv=new MPI::Request[nNumRecieves];
  int nCount=0;
  for(i=0;i<procTop.nNumRadialNeighbors;i++){
    if(procTop.nCoords[procTop.nRank][0]<procTop.nCoords[procTop.nRadialNeighborRanks[i]][0]){/*
      if current processor has a radial neighbor at inside post a recieve*/
      requestTempRecv[nCount]=MPI::COMM_WORLD.Irecv(grid.dLocalGridNew,1
        ,messPass.typeRecvNewVar[procTop.nRadialNeighborNeighborIDs[i]][grid.nU0]
        ,procTop.nRadialNeighborRanks[i],2);
      nCount++;
    }
  }
  
  //post a blocking recieve for inner radial neighbour
  for(i=0;i<procTop.nNumRadialNeighbors;i++){
    if(procTop.nCoords[procTop.nRank][0]>procTop.nCoords[procTop.nRadialNeighborRanks[i]][0]){/*
      if current processor has a radial neighbor at inside post a recieve*/
      MPI::COMM_WORLD.Send(grid.dLocalGridNew,1
        ,messPass.typeSendNewVar[procTop.nRadialNeighborNeighborIDs[i]][grid.nU0]
        ,procTop.nRadialNeighborRanks[i],2);
    }
  }
  
  //calculate outermost grid velocity
  for(i=grid.nStartGhostUpdateExplicit[grid.nU0][0][0];i<grid.nEndGhostUpdateExplicit[grid.nU0][0][0]
    ;i++){
    dCSum=0.0;
    dARhoSum=0.0;
    
    //calculate i of centered quantities
    nICen=i-grid.nCenIntOffset[0];
    dR_im1half_np1half=grid.dLocalGridOld[grid.nR][i-1][0][0];
    dR_ip1half_np1half=grid.dLocalGridOld[grid.nR][i][0][0];
    dRSq_im1half_np1half=dR_im1half_np1half*dR_im1half_np1half;
    dRSq_ip1half_np1half=dR_ip1half_np1half*dR_ip1half_np1half;
    dDonorFrac_ip1half=grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0];
    dDonorFrac_im1half=(grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0]
      +grid.dLocalGridOld[grid.nDonorCellFrac][nICen-1][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nU][1];j<grid.nEndUpdateExplicit[grid.nU][1];j++){
      
      //calculate j of interface quantities
      nJInt=j+grid.nCenIntOffset[1];
      
      for(k=grid.nStartUpdateExplicit[grid.nU][2];k<grid.nEndUpdateExplicit[grid.nU][2];k++){
      
        nKInt=k+grid.nCenIntOffset[2];
        dTemp=grid.dLocalGridOld[grid.nDCosThetaIJK][0][j][0]
          *grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dA_im1halfjk_np1half=dRSq_im1half_np1half*dTemp;
        dA_ip1halfjk_np1half=dRSq_ip1half_np1half*dTemp;
        dTemp=0.5*(dRSq_ip1half_np1half-dRSq_im1half_np1half)
          *grid.dLocalGridOld[grid.nDPhi][0][0][k];
        dA_ijm1halfk_np1half=dTemp*grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j-1][0];
        dA_ijp1halfk_np1half=dTemp*grid.dLocalGridOld[grid.nSinThetaIJp1halfK][0][j][0];
        dA_ijkm1half_np1half=0.5*(dRSq_ip1half_np1half-dRSq_im1half_np1half)
          *grid.dLocalGridOld[grid.nDTheta][0][j][0];
        dA_ijkp1half_np1half=dA_ijkm1half_np1half;
        
        //calculate difference between U and U0
        dUmU0_ip1halfjk_nm1half=grid.dLocalGridOld[grid.nU][i][j][k]
          -grid.dLocalGridOld[grid.nU0][i][0][0];
        dUmU0_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][i-1][j][k]
          -grid.dLocalGridNew[grid.nU0][i-1][0][0];
        
        //calculate rho at i-1/2, not time centered
        dRho_cen_im1half=(grid.dLocalGridOld[grid.nD][nICen][j][k]
          +grid.dLocalGridOld[grid.nD][nICen-1][j][k])*0.5;
        if(dUmU0_im1halfjk_np1half<0.0){//moving from outside in
          dRho_upwind_im1half=grid.dLocalGridOld[grid.nD][nICen][j][k];
        }
        else{//moving from inside out
          dRho_upwind_im1half=grid.dLocalGridOld[grid.nD][nICen-1][j][k];
        }
        dRho_im1halfjk_np1half=((1.0-dDonorFrac_im1half)*dRho_cen_im1half+dDonorFrac_im1half
          *dRho_upwind_im1half);
        
        //calculate rho at i+1/2, not time centered
        dRho_cen_ip1half=(grid.dLocalGridOld[grid.nD][nICen][j][k]+0.0)*0.5;
        if(dUmU0_ip1halfjk_nm1half<0.0){//moving from outside in
          dRho_upwind_ip1half=0.0;
        }
        else{//moving from inside out
          dRho_upwind_ip1half=grid.dLocalGridOld[grid.nD][nICen][j][k];
        }
        dRho_ip1halfjk_np1half=((1.0-dDonorFrac_ip1half)*dRho_cen_ip1half+dDonorFrac_ip1half
          *dRho_upwind_ip1half);
        
        //calculte rho at j-1/2
        dRho_cen_jm1half=(grid.dLocalGridOld[grid.nD][nICen][j-1][k]
          +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
        if(grid.dLocalGridNew[grid.nV][nICen][nJInt-1][k]<0.0){
          dRho_upwind_jm1half=grid.dLocalGridOld[grid.nD][nICen][j][k];
        }
        else{
          dRho_upwind_jm1half=grid.dLocalGridOld[grid.nD][nICen][j-1][k];
        }
        dRho_ijm1halfk_np1half=((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])
          *dRho_cen_jm1half+grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0]
          *dRho_upwind_jm1half);
        
        //calculte rho at j+1/2
        dRho_cen_jp1half=(grid.dLocalGridOld[grid.nD][nICen][j+1][k]
          +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
        if(grid.dLocalGridNew[grid.nV][nICen][nJInt][k]<0.0){
          dRho_upwind_jp1half=grid.dLocalGridOld[grid.nD][nICen][j+1][k];
        }
        else{
          dRho_upwind_jp1half=grid.dLocalGridOld[grid.nD][nICen][j][k];
        }
        dRho_ijp1halfk_np1half=((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])
          *dRho_cen_jp1half+grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0]
          *dRho_upwind_jp1half);
        
        //calculte rho at k-1/2
        dRho_cen_km1half=(grid.dLocalGridOld[grid.nD][nICen][j][k-1]
          +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
        if(grid.dLocalGridNew[grid.nW][nICen][j][nKInt-1]<0.0){
          dRho_upwind_km1half=grid.dLocalGridOld[grid.nD][nICen][j][k];
        }
        else{
          dRho_upwind_km1half=grid.dLocalGridOld[grid.nD][nICen][j][k-1];
        }
        dRho_ijkm1half_np1half=((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])
          *dRho_cen_km1half+grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0]
          *dRho_upwind_km1half);
        
        //calculte rho at j+1/2
        dRho_cen_kp1half=(grid.dLocalGridOld[grid.nD][nICen][j][k+1]
          +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
        if(grid.dLocalGridNew[grid.nW][nICen][j][nKInt]<0.0){
          dRho_upwind_kp1half=grid.dLocalGridOld[grid.nD][nICen][j][k+1];
        }
        else{
          dRho_upwind_kp1half=grid.dLocalGridOld[grid.nD][nICen][j][k];
        }
        dRho_ijkp1half_np1half=((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])
          *dRho_cen_kp1half+grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0]
          *dRho_upwind_kp1half);
        
        dCSum+=(grid.dLocalGridNew[grid.nU][i-1][j][k]-grid.dLocalGridNew[grid.nU0][i-1][0][0])
          *dA_im1halfjk_np1half*dRho_im1halfjk_np1half
          -grid.dLocalGridNew[grid.nU][i][j][k]*dA_ip1halfjk_np1half*dRho_ip1halfjk_np1half/*
          +grid.dLocalGridNew[grid.nV][nICen][nJInt-1][k]*dA_ijm1halfk_np1half
          *dRho_ijm1halfk_np1half
          -grid.dLocalGridNew[grid.nV][nICen][nJInt][k]*dA_ijp1halfk_np1half
          *dRho_ijp1halfk_np1half
          +grid.dLocalGridNew[grid.nW][nICen][j][k-1]*dA_ijkm1half_np1half
          *dRho_ijkm1half_np1half
          -grid.dLocalGridNew[grid.nW][nICen][j][k]*dA_ijkp1half_np1half
          *dRho_ijkp1half_np1half*/;
        dARhoSum+=dA_ip1halfjk_np1half*dRho_ip1halfjk_np1half;
      }
    }
    grid.dLocalGridNew[grid.nU0][i][0][0]=-1.0*dCSum/dARhoSum;
    
    //set U equal to U0
    for(j=0;j<grid.nLocalGridDims[procTop.nRank][grid.nU][1]+2*grid.nNumGhostCells;j++){
      for(k=0;k<grid.nLocalGridDims[procTop.nRank][grid.nU][2]+2*grid.nNumGhostCells;k++){
        grid.dLocalGridNew[grid.nU][i][j][k]=grid.dLocalGridNew[grid.nU0][i][0][0];
      }
    }
  }
  
  //wait for all recieves to complete
  MPI::Status *statusTempRecv=new MPI::Status[nNumRecieves];
  MPI::Request::Waitall(nNumRecieves,requestTempRecv,statusTempRecv);
  delete [] requestTempRecv;
  delete [] statusTempRecv;
}
