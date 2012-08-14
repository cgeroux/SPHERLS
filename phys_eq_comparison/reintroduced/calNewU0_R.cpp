void calNewU0_R(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop,MessPass &messPass){
  
  /**\todo
    At some point I will likely want to make this funciton compatiable with a 3D domain 
    decomposition instead of a purely radial domain decomposition.
  */
  
  //post a blocking recieve from inner radial neighbour
  int i;
  for(i=0;i<procTop.nNumRadialNeighbors;i++){
    if(procTop.nCoords[procTop.nRank][0]>procTop.nCoords[procTop.nRadialNeighborRanks[i]][0]){/*if
      current processor has a radial neighbor at inside, post a recieve*/
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
  int nICen;
  double dARatio;
  double dRho_im1half;
  double dRho_upwind_im1half;
  double dRho_cen_im1half;
  double dRho_ip1half;
  double dRho_upwind_ip1half;
  double dRho_cen_ip1half;
  double dUmU0_ip1halfjk_nm1half;
  double dUmU0_im1halfjk_np1half;
  double dDonorFrac_ip1half;
  double dDonorFrac_im1half;
  
  for(i=grid.nStartUpdateExplicit[grid.nU0][0];i<grid.nEndUpdateExplicit[grid.nU0][0];i++){/*nU0
    needs to be 1D*/
    
    //calculate new grid velocity
    grid.dLocalGridNew[grid.nU0][i][0][0]=grid.dLocalGridNew[grid.nU][i][0][0];
  }
  
  //post a blocking send to outer radial neighbour
  int nNumRecieves=0;
  for(i=0;i<procTop.nNumRadialNeighbors;i++){
    if(procTop.nCoords[procTop.nRank][0]<procTop.nCoords[procTop.nRadialNeighborRanks[i]][0]){/*if
      current processor has a radial neighbor at inside post a recieve*/
      
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
    if(procTop.nCoords[procTop.nRank][0]<procTop.nCoords[procTop.nRadialNeighborRanks[i]][0]){/*if
      current processor has a radial neighbor at inside post a recieve*/
      
      requestTempRecv[nCount]=MPI::COMM_WORLD.Irecv(grid.dLocalGridNew,1
        ,messPass.typeRecvNewVar[procTop.nRadialNeighborNeighborIDs[i]][grid.nU0]
        ,procTop.nRadialNeighborRanks[i],2);
      nCount++;
    }
  }
  
  //post a blocking recieve for inner radial neighbour
  for(i=0;i<procTop.nNumRadialNeighbors;i++){
    if(procTop.nCoords[procTop.nRank][0]>procTop.nCoords[procTop.nRadialNeighborRanks[i]][0]){/*if
      current processor has a radial neighbor at inside post a recieve*/
      
      MPI::COMM_WORLD.Send(grid.dLocalGridNew,1
        ,messPass.typeSendNewVar[procTop.nRadialNeighborNeighborIDs[i]][grid.nU0]
        ,procTop.nRadialNeighborRanks[i],2);
    }
  }
  
  //calculate outermost grid velocity
  for(i=grid.nStartGhostUpdateExplicit[grid.nU][0][0];i<grid.nEndGhostUpdateExplicit[grid.nU][0][0]
    ;i++){
    
    grid.dLocalGridNew[grid.nU0][i][0][0]=grid.dLocalGridNew[grid.nU][i][0][0];
    
  }
  
  //wait for all recieves to complete
  MPI::Status *statusTempRecv=new MPI::Status[nNumRecieves];
  MPI::Request::Waitall(nNumRecieves,requestTempRecv,statusTempRecv);
  delete [] requestTempRecv;
  delete [] statusTempRecv;
}
