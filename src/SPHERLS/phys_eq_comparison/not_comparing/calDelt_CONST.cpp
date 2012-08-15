void calDelt_CONST(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop){
  time.dDeltat_nm1half=time.dConstTimeStep;// time between t^n and t^{n+1}
  time.dDeltat_np1half=time.dConstTimeStep;
  time.dDeltat_n=time.dConstTimeStep;//time between t^{n-1/2} and t^{n+1/2}
  time.dt+=time.dConstTimeStep;//increase time by time step
  time.nTimeStepIndex++;//increase time step index
  
  //not yet updated to properly calculate donor cell fraction
  std::stringstream ssTemp;
  ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
    <<": Constant timestep function is not yet properly implemented to handel calculation of donor"
    <<" cell fraction! Stopping. \n";
  throw exception2(ssTemp.str(),INPUT);
}
