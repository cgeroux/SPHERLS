void calNewVelocities_RTP(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop){
  
  //calcualte new radial velocities
  calNewU_RTP(grid,parameters,time,procTop);
  
  //calculate new theta velocities
  calNewV_RTP(grid,parameters,time,procTop);
  
  //calculate new phi velocities
  calNewW_RTP(grid,parameters,time,procTop);
}
