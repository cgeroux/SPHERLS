void calNewU_R(Grid &grid,Parameters &parameters,Time &time,ProcTop &procTop){
  int i;
  int j;
  int k;
  int nICen;
  double dRho_ip1halfjk_n;
  double dP_ip1jk_n;
  double dP_ijk_n;
  double dA1CenGrad;
  double dA1UpWindGrad;
  double dA1;
  double dS1;
  double dS4;
  double dU_U0_Diff;
  double dU_ip1jk_nm1half;
  double dU_ijk_nm1half;
  double dU_ip1halfjk_nm1half;
  double dU0_ip1half_nm1half;
  double dRSq_ip1half_n;
  double dDM_ip1half;
  double dDonorFrac_ip1half;
  for(i=grid.nStartUpdateExplicit[grid.nU][0];i<grid.nEndUpdateExplicit[grid.nU][0];i++){
    
    //calculate i of centered quantities
    nICen=i-grid.nCenIntOffset[0];
    dU0_ip1half_nm1half=grid.dLocalGridOld[grid.nU0][i][0][0];
    dRSq_ip1half_n=grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nR][i][0][0];
    dDM_ip1half=(grid.dLocalGridOld[grid.nDM][nICen+1][0][0]
      +grid.dLocalGridOld[grid.nDM][nICen][0][0])*0.5;
    dDonorFrac_ip1half=(grid.dLocalGridOld[grid.nDonorCellFrac][nICen+1][0][0]
      +grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nU][1];j<grid.nEndUpdateExplicit[grid.nU][1];j++){
      for(k=grid.nStartUpdateExplicit[grid.nU][2];k<grid.nEndUpdateExplicit[grid.nU][2];k++){
        
        //CALCULATE INTERPOLATED QUANTITIES
        dU_ip1jk_nm1half=(grid.dLocalGridOld[grid.nU][i+1][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
        dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i-1][j][k])*0.5;
        dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
        dP_ip1jk_n=grid.dLocalGridOld[grid.nP][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nQ0][nICen+1][j][k];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][nICen][j][k]
          +grid.dLocalGridOld[grid.nQ0][nICen][j][k];
        
        //calculate dS1
        dS1=(dP_ip1jk_n-dP_ijk_n)/(dDM_ip1half*dRho_ip1halfjk_n);
        
        //Calculate dS4
        dS4=parameters.dG*grid.dLocalGridOld[grid.nM][i][0][0]/dRSq_ip1half_n;
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nU][i][j][k]=grid.dLocalGridOld[grid.nU][i][j][k]
          -time.dDeltat_n*(4.0*parameters.dPi*dRho_ip1halfjk_n*dRSq_ip1half_n*(dS1)+dS4);
          
        #if DEBUG_EQUATIONS==1
        
        //if we don't want zone by zone, set ssEnd.str("")
        std::stringstream ssName;
        std::stringstream ssEnd;
        if(parameters.bEveryJK){
          ssEnd<<"_"<<j<<"_"<<k;
        }
        else{
          ssEnd.str("");
        }
        
        //add M_r
        ssName.str("");
        ssName<<"M_r";
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,grid.dLocalGridOld[grid.nM][i][0][0]);
          
        //add S1
        ssName.str("");
        ssName<<"U_S1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-4.0*parameters.dPi*dRho_ip1halfjk_n*dRSq_ip1half_n*(dS1));
        
        //add S4
        ssName.str("");
        ssName<<"U_S4"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-1.0*dS4);
        
        //add DuDt
        ssName.str("");
        ssName<<"U_DuDt"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,(grid.dLocalGridNew[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU][i][j][k])
          /time.dDeltat_n);
        #endif
      }
    }
  }
  
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nU][0][0];
    i<grid.nEndGhostUpdateExplicit[grid.nU][0][0];i++){
    
    //calculate i of centered quantities
    nICen=i-grid.nCenIntOffset[0];
    dRSq_ip1half_n=grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nR][i][0][0];
    dDM_ip1half=(grid.dLocalGridOld[grid.nDM][nICen][0][0])*(0.5+parameters.dAlpha
      +parameters.dAlphaExtra);/**\BC Missing grid.dLocalGridOld[grid.nDM][i+1][0][0] in calculation
      of \f$S_1\f$ using \ref Parameters.dAlpha *grid.dLocalGridOld[grid.nDM][nICen][0][0] instead.
      */
    dDonorFrac_ip1half=grid.dLocalGridOld[grid.nDonorCellFrac][nICen][0][0];
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nU][0][1];
      j<grid.nEndGhostUpdateExplicit[grid.nU][0][1];j++){
      for(k=grid.nStartGhostUpdateExplicit[grid.nU][0][2];
        k<grid.nEndGhostUpdateExplicit[grid.nU][0][2];k++){
        
        //CALCULATE INTERPOLATED QUANTITIES
        dU_ip1jk_nm1half=grid.dLocalGridOld[grid.nU][i][j][k];
        dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][i][j][k]
          +grid.dLocalGridOld[grid.nU][i-1][j][k])*0.5;
        dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;/**\BC Missing density 
          outside model, setting it to zero. */
        dP_ijk_n=grid.dLocalGridOld[grid.nP][nICen][j][k]
          +grid.dLocalGridOld[grid.nQ0][nICen][j][k];
        dP_ip1jk_n=-1.0*dP_ijk_n;/**\BC Missing pressure outside surface setting it equal to 
          negative pressure in the center of the first cell so that it will be zero at surface.*/
        dEddyVisc_ip1halfjk_n=grid.dLocalGridOld[grid.nEddyVisc][nICen][j][k]*0.5;/**\BC assume 
          viscosity is zero outside the star.*/
        
        //calculate dS1
        dS1=(dP_ip1jk_n-dP_ijk_n)/(dDM_ip1half*dRho_ip1halfjk_n);
        
        //Calculate dS4
        dS4=parameters.dG*grid.dLocalGridOld[grid.nM][i][0][0]/dRSq_ip1half_n;
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nU][i][j][k]=grid.dLocalGridOld[grid.nU][i][j][k]
          -time.dDeltat_n*(4.0*parameters.dPi*dRho_ip1halfjk_n*dRSq_ip1half_n*(dS1)+dS4);
        
        #if DEBUG_EQUATIONS==1
        
        //if we don't want zone by zone, set ssEnd.str("")
        std::stringstream ssName;
        std::stringstream ssEnd;
        if(parameters.bEveryJK){
          ssEnd<<"_"<<j<<"_"<<k;
        }
        else{
          ssEnd.str("");
        }
        
        //add M_r
        ssName.str("");
        ssName<<"M_r";
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,grid.dLocalGridOld[grid.nM][i][0][0]);
          
        //add S1
        ssName.str("");
        ssName<<"U_S1"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-4.0*parameters.dPi*dRho_ip1halfjk_n*dRSq_ip1half_n*(dS1));
        
        //add S4
        ssName.str("");
        ssName<<"U_S4"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
        ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,-1.0*dS4);
        
        //add DuDt
        ssName.str("");
        ssName<<"U_DuDt"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]+grid.nCenIntOffset[0]-1
          ,(grid.dLocalGridNew[grid.nU][i][j][k]-grid.dLocalGridOld[grid.nU][i][j][k])
          /time.dDeltat_n);
        #endif
      }
    }
  }
  
  #if SEDOV==1
  //calculate gost region 1 inner most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nU][1][0];
    i<grid.nEndGhostUpdateExplicit[grid.nU][1][0];i++){
    
    //calculate i of centered quantities
    dU0_ip1half_nm1half=grid.dLocalGridOld[grid.nU0][i][0][0];
    nICen=i-grid.nCenIntOffset[0];
    
    for(j=grid.nStartGhostUpdateExplicit[grid.nU][1][1];
      j<grid.nEndGhostUpdateExplicit[grid.nU][1][1];j++){
      for(k=grid.nStartGhostUpdateExplicit[grid.nU][1][2];
        k<grid.nEndGhostUpdateExplicit[grid.nU][1][2];k++){
        
        //calculate interpolated quantities
        dRho_ip1halfjk_n=(grid.dLocalGridOld[grid.nD][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nD][nICen][j][k])*0.5;
        dU_ip1jk_nm1half=(grid.dLocalGridOld[grid.nU][i+1][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
        dU_ijk_nm1half=(grid.dLocalGridOld[grid.nU][i+1][j][k]
          +grid.dLocalGridOld[grid.nU][i][j][k])*0.5;
        dU_ip1halfjk_nm1half=grid.dLocalGridOld[grid.nU][i][j][k];
        
        //calculate advection term in x1-direction
        dP_ip1jk_n=grid.dLocalGridOld[grid.nP][nICen+1][j][k]
          +grid.dLocalGridOld[grid.nQ0][nICen+1][j][k];
        dP_ijk_n=grid.dLocalGridOld[grid.nP][nICen][j][k]
          +grid.dLocalGridOld[grid.nQ0][nICen][j][k];
        
        //calculate source terms in x1-direction
        dS1=(dP_ip1jk_n-dP_ijk_n)/(grid.dLocalGridOld[grid.nDM][nICen+1][0][0]
          +grid.dLocalGridOld[grid.nDM][nICen][0][0])*2.0/dRho_ip1halfjk_n;
        dRSq_ip1half_n=grid.dLocalGridOld[grid.nR][i][0][0]*grid.dLocalGridOld[grid.nR][i][0][0];
        dS4=parameters.dG*grid.dLocalGridOld[grid.nM][i][0][0]/dRSq_ip1half_n;
        
        //calculate new velocity
        grid.dLocalGridNew[grid.nU][i][j][k]=grid.dLocalGridOld[grid.nU][i][j][k]
          -time.dDeltat_n*(4.0*parameters.dPi*dRho_ip1halfjk_n*dRSq_ip1half_n*(dS1)+dS4);
      }
    }
  }
  #endif
}
