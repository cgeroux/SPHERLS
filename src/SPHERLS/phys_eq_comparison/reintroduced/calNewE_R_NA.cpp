void calNewE_R_NA(Grid &grid, Parameters &parameters, Time &time, ProcTop &procTop){
  int i;
  int j;
  int k;
  int nIInt;
  double dU_ijk_np1half;
  double dU0_i_np1half;
  double dE_ip1halfjk_n;
  double dE_im1halfjk_n;
  double dR_i_n;
  double dR_im1half_n;
  double dR_ip1half_n;
  double dRSq_i_n;
  double dRSq_ip1half_n;
  double dR4_ip1half_n;
  double dRSq_im1half_n;
  double dR4_im1half_n;
  double dUmU0_ijk_np1half;
  double dUR2_im1halfjk_np1half;
  double dUR2_ip1halfjk_np1half;
  double dRhoAve_ip1half_n;
  double dRhoAve_im1half_n;
  double dRho_ip1halfjk_n;
  double dRho_im1halfjk_n;
  double dTSq_ip1jk_n;
  double dTSq_ijk_n;
  double dTSq_im1jk_n;
  double dT4_ip1jk_n;
  double dT4_ijk_n;
  double dT4_im1jk_n;
  double dKappa_ip1halfjk_n;
  double dKappa_im1halfjk_n;
  double dA1CenGrad;
  double dA1UpWindGrad;
  double dA1;
  double dP_ijk_n;
  double dS1;
  double dTGrad_ip1half;
  double dTGrad_im1half;
  double dGrad_ip1half;
  double dGrad_im1half;
  double dS4;
  double dPiSq=parameters.dPi*parameters.dPi;
  for(i=grid.nStartUpdateExplicit[grid.nE][0];i<grid.nEndUpdateExplicit[grid.nE][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_n=grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dR_im1half_n=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_i_n=(dR_ip1half_n+dR_im1half_n)*0.5;
    dRSq_i_n=dR_i_n*dR_i_n;
    dRSq_ip1half_n=dR_ip1half_n*dR_ip1half_n;
    dR4_ip1half_n=dRSq_ip1half_n*dRSq_ip1half_n;
    dRSq_im1half_n=dR_im1half_n*dR_im1half_n;
    dR4_im1half_n=dRSq_im1half_n*dRSq_im1half_n;
    dRhoAve_ip1half_n=(grid.dLocalGridOld[grid.nD][i][0][0]
      +grid.dLocalGridOld[grid.nD][i+1][0][0])*0.5;
    dRhoAve_im1half_n=(grid.dLocalGridOld[grid.nD][i][0][0]
      +grid.dLocalGridOld[grid.nD][i-1][0][0])*0.5;
    dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
      +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nE][1];j<grid.nEndUpdateExplicit[grid.nE][1];j++){
      
      for(k=grid.nStartUpdateExplicit[grid.nE][2];k<grid.nEndUpdateExplicit[grid.nE][2];k++){
        
        //Calculate interpolated quantities
        dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dE_ip1halfjk_n=(grid.dLocalGridOld[grid.nE][i+1][j][k]
          +grid.dLocalGridOld[grid.nE][i][j][k])*0.5;
        dE_im1halfjk_n=(grid.dLocalGridOld[grid.nE][i][j][k]
          +grid.dLocalGridOld[grid.nE][i-1][j][k])*0.5;
        dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][i+1][j][k]
          +grid.dLocalGridOld[grid.nD][i][j][k])*0.5;
        dRho_im1halfjk_n=(grid.dLocalGridOld[grid.nD][i][j][k]
          +grid.dLocalGridOld[grid.nD][i-1][j][k])*0.5;

        //calculate derived quantities
        dUR2_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dRSq_im1half_n;
        dUR2_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dRSq_ip1half_n;
        dTSq_ip1jk_n=grid.dLocalGridOld[grid.nT][i+1][j][k]*grid.dLocalGridOld[grid.nT][i+1][j][k];
        dT4_ip1jk_n=dTSq_ip1jk_n*dTSq_ip1jk_n;
        dTSq_ijk_n=grid.dLocalGridOld[grid.nT][i][j][k]*grid.dLocalGridOld[grid.nT][i][j][k];
        dT4_ijk_n=dTSq_ijk_n*dTSq_ijk_n;
        dTSq_im1jk_n=grid.dLocalGridOld[grid.nT][i-1][j][k]*grid.dLocalGridOld[grid.nT][i-1][j][k];
        dT4_im1jk_n=dTSq_im1jk_n*dTSq_im1jk_n;
        dKappa_ip1halfjk_n=(dT4_ip1jk_n+dT4_ijk_n)/(dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]+dT4_ip1jk_n
          /grid.dLocalGridOld[grid.nKappa][i+1][j][k]);
        dKappa_im1halfjk_n=(dT4_im1jk_n+dT4_ijk_n)/(dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]+dT4_im1jk_n
          /grid.dLocalGridOld[grid.nKappa][i-1][j][k]);
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
          dP_ijk_n=dP_ijk_n+grid.dLocalGridOld[grid.nQ0][i][j][k];
        #endif
        
        //Calcuate dA1
        dA1CenGrad=(dE_ip1halfjk_n-dE_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0];
        dUmU0_ijk_np1half=(dU_ijk_np1half-dU0_i_np1half);
        if(dUmU0_ijk_np1half<0.0){//moving in the negative direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nE][i+1][j][k]
            -grid.dLocalGridOld[grid.nE][i][j][k])/(grid.dLocalGridOld[grid.nDM][i+1][0][0]
            +grid.dLocalGridOld[grid.nDM][i][0][0])*2.0;
        }
        else{//moving in the postive direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k]
            -grid.dLocalGridOld[grid.nE][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=dUmU0_ijk_np1half*dRSq_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate dS1
        dS1=dP_ijk_n/grid.dLocalGridOld[grid.nD][i][j][k]
          *(dUR2_ip1halfjk_np1half-dUR2_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //Calculate dS4
        dTGrad_ip1half=(dT4_ip1jk_n-dT4_ijk_n)/(grid.dLocalGridOld[grid.nDM][i+1][0][0]
          +grid.dLocalGridOld[grid.nDM][i][0][0])*2.0;
        dTGrad_im1half=(dT4_ijk_n-dT4_im1jk_n)/(grid.dLocalGridOld[grid.nDM][i][0][0]
          +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        dGrad_ip1half=dRhoAve_ip1half_n*dR4_ip1half_n/(dKappa_ip1halfjk_n
          *dRho_ip1halfjk_n)*dTGrad_ip1half;
        dGrad_im1half=dRhoAve_im1half_n*dR4_im1half_n/(dKappa_im1halfjk_n
          *dRho_im1halfjk_n)*dTGrad_im1half;
        dS4=16.0*dPiSq*grid.dLocalGridOld[grid.nD][i][0][0]
          *(dGrad_ip1half-dGrad_im1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //calculate new energy
        grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridOld[grid.nE][i][j][k]
          -time.dDeltat_np1half*(4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]
          *(dA1+dS1)-4.0*parameters.dSigma/(3.0
          *grid.dLocalGridOld[grid.nD][i][j][k])*(dS4));
        
        #if DEBUG_EQUATIONS==1
        
        int nGhostCells=1;
        if(procTop.nRank==0){
          nGhostCells=0;
        }
        
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
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,grid.dLocalGridOld[grid.nE][i][j][k]);
        
        //add A1
        ssName.str("");
        ssName<<"E_A1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,-4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1));
        
        //add S1
        ssName.str("");
        ssName<<"E_S1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,-4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]*(dS1));
        
        //add S4
        ssName.str("");
        ssName<<"E_S4"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS4));
        
        //add DEDt
        ssName.str("");
        ssName<<"E_DEDt"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,(grid.dLocalGridNew[grid.nE][i][j][k]-grid.dLocalGridOld[grid.nE][i][j][k])
          /time.dDeltat_np1half);
        #endif
        
        if(grid.dLocalGridNew[grid.nE][i][j][k]<0.0){
          
          #if SIGNEGENG==1
          raise(SIGINT);
          #endif
          
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
            <<": negative energy calculated in , ("<<i<<","<<j<<","<<k<<")\n";
          throw exception2(ssTemp.str(),CALCULATION);
          
        }
      }
    }
  }
  
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nE][0][0];i<grid.nEndGhostUpdateExplicit[grid.nE][0][0]
    ;i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dR_ip1half_n=grid.dLocalGridOld[grid.nR][nIInt][0][0];
    dR_im1half_n=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];
    dR_i_n=(dR_ip1half_n+dR_im1half_n)*0.5;
    dRSq_i_n=dR_i_n*dR_i_n;
    dRSq_ip1half_n=dR_ip1half_n*dR_ip1half_n;
    dRSq_im1half_n=dR_im1half_n*dR_im1half_n;
    dR4_im1half_n=dRSq_im1half_n*dRSq_im1half_n;
    dRhoAve_im1half_n=(grid.dLocalGridOld[grid.nD][i][0][0]
      +grid.dLocalGridOld[grid.nD][i-1][0][0])*0.5;
    dU0_i_np1half=(grid.dLocalGridNew[grid.nU0][nIInt][0][0]
      +grid.dLocalGridNew[grid.nU0][nIInt-1][0][0])*0.5;
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nE][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nE][0][1];j++){
      
      for(k=grid.nStartGhostUpdateExplicit[grid.nE][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nE][0][2];k++){
        
        //Calculate interpolated quantities
        dU_ijk_np1half=(grid.dLocalGridNew[grid.nU][nIInt][j][k]
          +grid.dLocalGridNew[grid.nU][nIInt-1][j][k])*0.5;
        dE_ip1halfjk_n=grid.dLocalGridOld[grid.nE][i][j][k];/**\BC Setting energy at surface equal 
          to energy in last zone.*/
        dE_im1halfjk_n=(grid.dLocalGridOld[grid.nE][i][j][k]+grid.dLocalGridOld[grid.nE][i-1][j][k])
          *0.5;
        dRho_im1halfjk_n=(grid.dLocalGridOld[grid.nD][i][j][k]
          +grid.dLocalGridOld[grid.nD][i-1][j][k])*0.5;
        
        //calculate derived quantities
        dUR2_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]*dRSq_im1half_n;
        dUR2_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]*dRSq_ip1half_n;
        dTSq_ijk_n=grid.dLocalGridOld[grid.nT][i][j][k]*grid.dLocalGridOld[grid.nT][i][j][k];
        dT4_ijk_n=dTSq_ijk_n*dTSq_ijk_n;
        dTSq_im1jk_n=grid.dLocalGridOld[grid.nT][i-1][j][k]*grid.dLocalGridOld[grid.nT][i-1][j][k];
        dT4_im1jk_n=dTSq_im1jk_n*dTSq_im1jk_n;
        dKappa_im1halfjk_n=(dT4_im1jk_n+dT4_ijk_n)/(dT4_ijk_n
          /grid.dLocalGridOld[grid.nKappa][i][j][k]+dT4_im1jk_n
          /grid.dLocalGridOld[grid.nKappa][i-1][j][k]);
        dP_ijk_n=grid.dLocalGridOld[grid.nP][i][j][k];
        #if VISCOUS_ENERGY_EQ==1
        dP_ijk_n+=grid.dLocalGridOld[grid.nQ0][i][j][k];
        #endif
        
        //Calcuate dA1
        dA1CenGrad=(dE_ip1halfjk_n-dE_im1halfjk_n)/grid.dLocalGridOld[grid.nDM][i][0][0];
        dUmU0_ijk_np1half=(dU_ijk_np1half-dU0_i_np1half);
        if(dUmU0_ijk_np1half<0.0){//moving in the negative direction
          dA1UpWindGrad=0.0;/**\BC Upwind gradient in dA1 term should be zero as there is no flow
           into the star.*/
        }
        else{//moving in the postive direction
          dA1UpWindGrad=(grid.dLocalGridOld[grid.nE][i][j][k]
            -grid.dLocalGridOld[grid.nE][i-1][j][k])/(grid.dLocalGridOld[grid.nDM][i][0][0]
            +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        }
        dA1=dUmU0_ijk_np1half*dRSq_i_n*((1.0-grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])
          *dA1CenGrad+grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]*dA1UpWindGrad);
        
        //calculate dS1
        dS1=dP_ijk_n/grid.dLocalGridOld[grid.nD][i][j][k]
          *(dUR2_ip1halfjk_np1half-dUR2_im1halfjk_np1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //Calculate dS4
        dTGrad_im1half=(dT4_ijk_n-dT4_im1jk_n)/(grid.dLocalGridOld[grid.nDM][i][0][0]
          +grid.dLocalGridOld[grid.nDM][i-1][0][0])*2.0;
        dGrad_ip1half=-3.0*dRSq_ip1half_n*dT4_ijk_n/(8.0*parameters.dPi);/**\BC
          Missing grid.dLocalGridOld[grid.nT][i+1][0][0] using flux equals \f$2\sigma T^4\f$ at 
          surface.*/
        dGrad_im1half=dRhoAve_im1half_n*dR4_im1half_n/(dKappa_im1halfjk_n*dRho_im1halfjk_n)
          *dTGrad_im1half;
        dS4=16.0*parameters.dPi*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]
          *(dGrad_ip1half-dGrad_im1half)/grid.dLocalGridOld[grid.nDM][i][0][0];
        
        //calculate new energy
        grid.dLocalGridNew[grid.nE][i][j][k]=grid.dLocalGridOld[grid.nE][i][j][k]
          -time.dDeltat_np1half*(4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]*(dA1
          +dS1)-4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])
          *(dS4));
        
        
        #if DEBUG_EQUATIONS==1
        
        int nGhostCells=1;
        if(procTop.nRank==0){
          nGhostCells=0;
        }
        
        //add E
        parameters.profileDataDebug.setMaxAbs("E"
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,grid.dLocalGridOld[grid.nE][i][j][k]);
        
        //add A1
        parameters.profileDataDebug.setMaxAbs("E_A1"
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,-4.0*parameters.dPi*grid.dLocalGridOld[grid.nDenAve][i][0][0]*(dA1));
        
        //add S1
        parameters.profileDataDebug.setMaxAbs("E_S1"
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,-4.0*parameters.dPi*grid.dLocalGridOld[grid.nD][i][0][0]*(dS1));
        
        //add S4
        parameters.profileDataDebug.setMaxAbs("E_S4"
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,4.0*parameters.dSigma/(3.0*grid.dLocalGridOld[grid.nD][i][j][k])*(dS4));
        
        //add DEDt
        parameters.profileDataDebug.setMaxAbs("E_DEDt"
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,(grid.dLocalGridNew[grid.nE][i][j][k]-grid.dLocalGridOld[grid.nE][i][j][k])
          /time.dDeltat_np1half);
        #endif
        
        if(grid.dLocalGridNew[grid.nE][i][j][k]<0.0){
          
          #if SIGNEGENG==1
          raise(SIGINT);
          #endif
          
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
            <<": negative energy calculated in , ("<<i<<","<<j<<","<<k<<")\n";
          throw exception2(ssTemp.str(),CALCULATION);
          
        }
      }
    }
  }
}
