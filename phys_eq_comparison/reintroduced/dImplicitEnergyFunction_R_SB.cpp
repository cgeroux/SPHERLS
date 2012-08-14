double dImplicitEnergyFunction_R_SB(Grid &grid,Parameters &parameters,Time &time
  ,double dTemps[],int i,int j,int k){
  
  double dT_ijk_np1=dTemps[0];
  double dT_im1jk_np1=dTemps[1];
  int nIInt=i+grid.nCenIntOffset[0];
  
  //Calculate interpolated quantities
  double dR_ip1half_n=grid.dLocalGridOld[grid.nR][nIInt][0][0];
  double dR_im1half_n=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
  double dR_i_n=(dR_ip1half_n+dR_im1half_n)*0.5;
  double dRSq_i_n=dR_i_n*dR_i_n;
  double dRSq_ip1half_n=dR_ip1half_n*dR_ip1half_n;
  double dRSq_im1half_n=dR_im1half_n*dR_im1half_n;
  double dR4_ip1half_n=dRSq_ip1half_n*dRSq_ip1half_n;
  double dR4_im1half_n=dRSq_im1half_n*dRSq_im1half_n;
  double dRhoAve_ip1half_n=(grid.dLocalGridOld[grid.nDenAve][i][0][0])*0.5;/**\BC missing density
    outside model assuming it is zero*/
  double dRhoAve_im1half_n=(grid.dLocalGridOld[grid.nD][i][0][0]
    +grid.dLocalGridOld[grid.nD][i-1][0][0])*0.5;
  double dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][i][j][k])*0.5;/**\BC missing desnity outside
    model assuming it is zero*/
  double dRho_im1halfjk_n=(grid.dLocalGridOld[grid.nD][i][j][k]
    +grid.dLocalGridOld[grid.nD][i-1][j][k])*0.5;
  double dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
    +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
  double dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
    +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
  
  double dT_ijk_np1half=(dT_ijk_np1+grid.dLocalGridOld[grid.nT][i][j][k])*0.5;
  double dTSq_ijk_np1half=dT_ijk_np1half*dT_ijk_np1half;
  double dT4_ijk_np1half=dTSq_ijk_np1half*dTSq_ijk_np1half;
  
  double dT_im1jk_np1half=(dT_im1jk_np1+grid.dLocalGridOld[grid.nT][i-1][j][k])*0.5;
  double dTSq_im1jk_np1half=dT_im1jk_np1half*dT_im1jk_np1half;
  double dT4_im1jk_np1half=dTSq_im1jk_np1half*dTSq_im1jk_np1half;
  
  double dE_ijk_np1=parameters.eosTable.dGetEnergy(dT_ijk_np1
    ,grid.dLocalGridNew[grid.nD][i][j][k]);
  double dE_ijk_np1half=parameters.eosTable.dGetEnergy(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  double dE_im1jk_np1half=parameters.eosTable.dGetEnergy(dT_im1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i-1][j][k]);
  double dE_ip1halfjk_np1half=dE_ijk_np1half;/**\BC Assuming energy outside model is the same as
    the energy in the last zone inside the model.*/
  double dE_im1halfjk_np1half=(dE_im1jk_np1half+dE_ijk_np1half)*0.5;
  
  double dP_ijk_np1half=parameters.eosTable.dGetPressure(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  #if VISCOUS_ENERGY_EQ==1
  dP_ijk_np1half=dP_ijk_np1half+grid.dLocalGridOld[grid.nQ0][i][j][k];
  #endif
  
  double dKappa_ijk_np1half=parameters.eosTable.dGetOpacity(dT_ijk_np1half
    ,grid.dLocalGridOld[grid.nD][i][j][k]);
  double dKappa_im1jk_np1half=parameters.eosTable.dGetOpacity(dT_im1jk_np1half
    ,grid.dLocalGridOld[grid.nD][i-1][j][k]);
  double dKappa_im1halfjk_np1half=(dT4_im1jk_np1half+dT4_ijk_np1half)/(dT4_ijk_np1half
    /dKappa_ijk_np1half+dT4_im1jk_np1half/dKappa_im1jk_np1half);
  
  //Calcuate dA1
  double dA1CenGrad=(dE_ip1halfjk_np1half-dE_im1halfjk_np1half)
    /grid.dLocalGridOld[grid.nDM][i][0][0];
  double dA1UpWindGrad=0.0;
  double dU_U0_Diff=(dU_ijk_np1half-dU0_i_np1half);
  if(dU_U0_Diff<0.0){//moving in the negative radial direction
    dA1UpWindGrad=0.0;/**\BC A1 upwind set to zero as no material is flowing into the star*/
  }
  else{//moving in the postive radial direction
    dA1UpWindGrad=(dE_ijk_np1half-dE_im1jk_np1half)/(grid.dLocalGridOld[grid.nDM][i][0][0]
      +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
  }
  
  double dDEDM=((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
    *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
  
  //apply DEDM clamp if set, and above the required mass
  if(parameters.bDEDMClamp){
    if(parameters.dDEDMClampMr!=-1.0){//clamp has been set
      if(grid.dLocalGridOld[grid.nM][nIInt][0][0]>=parameters.dDEDMClampMr){
        dDEDM=parameters.dDEDMClampValue;
      }
    }
    else{//if the clamp hasn't been set lets see if we can set it
      
      //this should only be set on the first, static and spherically symetric model
      if(grid.dLocalGridOld[grid.nT][i][0][0]<=parameters.dEDMClampTemperature){
        parameters.dDEDMClampMr=grid.dLocalGridOld[grid.nM][nIInt][0][0];
        parameters.dDEDMClampValue=dDEDM;
      }
    }
  }
  
  double dA1=dU_U0_Diff*dRSq_i_n*dDEDM;
  
  //calculate dS1
  double dUR2_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dRSq_im1half_n;
  double dUR2_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dRSq_ip1half_n;
  double dS1=dP_ijk_np1half/grid.dLocalGridOld[grid.nD][i][j][k]
    *(dUR2_ip1halfjk_np1half-dUR2_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
  //Calculate dS4
  double dTGrad_im1half_np1half=(dT4_ijk_np1half-dT4_im1jk_np1half)
    /(grid.dLocalGridOld[grid.nDM][i][0][0]+grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
  double dGrad_ip1half_np1half=-3.0*dRSq_ip1half_n*dT4_ijk_np1half/(8.0*parameters.dPi);/**\BC 
    Missing grid.dLocalGridOld[grid.nT][i+1][0][0] using flux equals \f$2\sigma T^4\f$ at surface.*/
  double dGrad_im1half_np1half=dRhoAve_im1half_n*dR4_im1half_n/(dKappa_im1halfjk_np1half
    *dRho_im1halfjk_n)*dTGrad_im1half_np1half;
  double dS4=16.0*parameters.dPi*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]
    *(dGrad_ip1half_np1half-dGrad_im1half_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
  
  
  #if DEBUG_EQUATIONS==1
  if(parameters.bSetThisCall){
    
    //if we don't want zone by zone, set ssEnd.str("")
    std::stringstream ssName;
    std::stringstream ssEnd;
    if(parameters.bEveryJK){
      ssEnd<<"_"<<j<<"_"<<k;
    }
    else{
      ssEnd.str("");
    }
    
    //add E
    ssName.str("");
    ssName<<"E"<<ssEnd.str();
    parameters.profileDataDebug.setMaxAbs(ssName.str()
      ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
      ,grid.dLocalGridOld[grid.nE][i][j][k]);
    
    //add A1
    ssName.str("");
    ssName<<"E_A1"<<ssEnd.str();
    parameters.profileDataDebug.setMaxAbs(ssName.str()
      ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
      ,-4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][j][k]*(dA1));
    
    //add S1
    ssName.str("");
    ssName<<"E_S1"<<ssEnd.str();
    parameters.profileDataDebug.setMaxAbs(ssName.str()
      ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
      ,-4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][j][k]*(dS1));
    
    //add S4
    ssName.str("");
    ssName<<"E_S4"<<ssEnd.str();
    parameters.profileDataDebug.setMaxAbs(ssName.str()
      ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
      ,4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS4));
    
    //add E_DEDt
    ssName.str("");
    ssName<<"E_DEDt"<<ssEnd.str();
    parameters.profileDataDebug.setMaxAbs(ssName.str()
      ,i+grid.nGlobalGridPositionLocalGrid[0]-grid.nNumGhostCells
      ,(dE_ijk_np1-grid.dLocalGridOld[grid.nE][i][j][k])
      /time.dDeltat_np1half);
  }
  #endif
  
  //calculate energy equation discrepancy
  return (dE_ijk_np1-grid.dLocalGridOld[grid.nE][i][j][k])/time.dDeltat_np1half
    +4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]*(dA1+dS1)
    -4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*dS4;
  
}
