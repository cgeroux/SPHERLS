void calNewVelocities_RT(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop){
  
  //calcualte new radial velocities
  calNewU_RT(grid,parameters,time,procTop);
  
  //calculate new theta velocities
  calNewV_RT(grid,parameters,time,procTop);
}
