double dEOS_GL(double dRho, double dE, Parameters parameters){
  return dRho*(parameters.dGamma-1.0)*dE;
}
