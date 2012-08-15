void initDonorFracAndMaxConVel_R_GL(Grid &grid, Parameters &parameters){
  
  int nEndCalc=std::max(grid.nEndGhostUpdateExplicit[grid.nD][0][0]
    ,grid.nEndUpdateExplicit[grid.nD][0]);
  int i;
  int j;
  int k;
  int nIInt;
  double dC;
  double dTest_ConVelOverSoundSpeed_R;
  double dTest_ConVelOverSoundSpeed=0.0;
  double dTest_ConVel_R;
  double dTest_ConVel=0.0;
  
  for(i=grid.nStartUpdateExplicit[grid.nD][0];i<nEndCalc;i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dTest_ConVelOverSoundSpeed=0.0;
    
    for(j=grid.nStartUpdateExplicit[grid.nD][1];j<grid.nEndUpdateExplicit[grid.nD][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nD][2];k<grid.nEndUpdateExplicit[grid.nD][2];k++){
        
        dC=sqrt(parameters.dGamma*(grid.dLocalGridOld[grid.nP][i][j][k]
          +grid.dLocalGridOld[grid.nQ0][i][j][k])/grid.dLocalGridOld[grid.nD][i][j][k]);
        dTest_ConVelOverSoundSpeed_R=fabs(grid.dLocalGridOld[grid.nU][nIInt][j][k]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0])/dC;
        dTest_ConVel_R=fabs(grid.dLocalGridOld[grid.nU][nIInt][j][k]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0]);
        
        //keep largest convective velocity over sound speed
        if(dTest_ConVelOverSoundSpeed_R>dTest_ConVelOverSoundSpeed){
          dTest_ConVelOverSoundSpeed=dTest_ConVelOverSoundSpeed_R;
        }
        
        //keep largest convective velocity
        if(dTest_ConVel_R>dTest_ConVel){
          dTest_ConVel=dTest_ConVel_R;
        }
      }
    }
    
    //set donnor fraction
    double dTest_ConVelOverSoundSpeed2=parameters.dDonorCellMultiplier*dTest_ConVelOverSoundSpeed;
    if(dTest_ConVelOverSoundSpeed2>1.0){
      grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]=1.0;
    }
    else if(dTest_ConVelOverSoundSpeed2<parameters.dDonorCellMin){
      grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]=parameters.dDonorCellMin;
    }
    else{
      grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]=dTest_ConVelOverSoundSpeed2;
    }
  }
  
  //keep largest convective velocity
  double dTest_ConVel2;
  MPI::COMM_WORLD.Allreduce(&dTest_ConVel,&dTest_ConVel2,1,MPI::DOUBLE,MPI_MAX);
  parameters.dMaxConvectiveVelocity=dTest_ConVel2;
}
