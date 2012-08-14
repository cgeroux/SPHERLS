void calNewD_R(Grid &grid, Parameters &parameters, Time &time,ProcTop &procTop){
  int i;
  int j;
  int k;
  int nIInt;
  int nJInt;
  int nKInt;
  double dDelRCu_i_n;
  double dDelRCu_i_np1;
  double dVRatio;
  double dR_ip1half_np1half;
  double dR_im1half_np1half;
  double dRSq_ip1half_np1half;
  double dRSq_im1half_np1half;
  double dDelRSq_i_np1half;
  double dV_np1;
  double dA_im1half;
  double dA_ip1half;
  double dRho_im1half;
  double dRho_cen_im1half;
  double dRho_upwind_im1half;
  double dRho_ip1half;
  double dRho_cen_ip1half;
  double dRho_upwind_ip1half;
  double dDeltaRhoR;
  double dVr_np1;
  double dUmU0_ip1halfjk_np1half;
  double dUmU0_ip1halfjk_nm1half;
  double dUmU0_im1halfjk_np1half;
  double d1Thrid=0.333333333333333333333333333333;
  double dDonorFrac_ip1half;
  double dDonorFrac_im1half;
  for(i=grid.nStartUpdateExplicit[grid.nD][0];i<grid.nEndUpdateExplicit[grid.nD][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dDelRCu_i_n=(pow(grid.dLocalGridOld[grid.nR][nIInt][0][0],3.0)
          -pow(grid.dLocalGridOld[grid.nR][nIInt-1][0][0],3.0));
    dDelRCu_i_np1=(pow(grid.dLocalGridNew[grid.nR][nIInt][0][0],3.0)
          -pow(grid.dLocalGridNew[grid.nR][nIInt-1][0][0],3.0));
    dR_ip1half_np1half=grid.dLocalGridOld[grid.nR][nIInt][0][0];//could be time centered
    dR_im1half_np1half=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];//could be time centered
    dRSq_ip1half_np1half=dR_ip1half_np1half*dR_ip1half_np1half;
    dRSq_im1half_np1half=dR_im1half_np1half*dR_im1half_np1half;
    dDelRSq_i_np1half=dRSq_ip1half_np1half-dRSq_im1half_np1half;
    dVRatio=dDelRCu_i_n/dDelRCu_i_np1;//calculate ratio of volume at n to volume at n+1
    dDonorFrac_ip1half=(grid.dLocalGridOld[grid.nDonorCellFrac][i+1][0][0]
      +grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0])*0.5;
    dDonorFrac_im1half=(grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]
      +grid.dLocalGridOld[grid.nDonorCellFrac][i-1][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nD][1];j<grid.nEndUpdateExplicit[grid.nD][1];j++){
      
      for(k=grid.nStartUpdateExplicit[grid.nD][2];k<grid.nEndUpdateExplicit[grid.nD][2];k++){
        
        dV_np1=d1Thrid*dDelRCu_i_np1;
          
        //CALCULATE RATE OF CHANGE IN RHO IN RADIAL DIRECTION
        
        //calculate area at i-1/2
        dA_im1half=dRSq_im1half_np1half;
        
        //calculate area at i+1/2
        dA_ip1half=dRSq_ip1half_np1half;
        
        //calculate difference between U and U0
        dUmU0_ip1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt][j][k]
          -grid.dLocalGridNew[grid.nU0][nIInt][0][0];
        
        dUmU0_ip1halfjk_nm1half=grid.dLocalGridOld[grid.nU][nIInt][j][k]
          -grid.dLocalGridOld[grid.nU0][nIInt][0][0];
        
        dUmU0_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]
          -grid.dLocalGridNew[grid.nU0][nIInt-1][0][0];
        
        //calculate rho at i-1/2, not time centered
        dRho_cen_im1half=(grid.dLocalGridOld[grid.nD][i][j][k]
          +grid.dLocalGridOld[grid.nD][i-1][j][k])*0.5;
        if(dUmU0_im1halfjk_np1half<0.0){//moving from outside in
          dRho_upwind_im1half=grid.dLocalGridOld[grid.nD][i][j][k];
        }
        else{//moving from inside out
          dRho_upwind_im1half=grid.dLocalGridOld[grid.nD][i-1][j][k];
        }
        dRho_im1half=((1.0-dDonorFrac_im1half)*dRho_cen_im1half+dDonorFrac_im1half
          *dRho_upwind_im1half);
        
        //calculate rho at i+1/2, not time centered
        dRho_cen_ip1half=(grid.dLocalGridOld[grid.nD][i][j][k]
          +grid.dLocalGridOld[grid.nD][i+1][j][k])*0.5;
        if(dUmU0_ip1halfjk_nm1half<0.0){//moving from outside in
          dRho_upwind_ip1half=grid.dLocalGridOld[grid.nD][i+1][j][k];
        }
        else{//moving from inside out
          dRho_upwind_ip1half=grid.dLocalGridOld[grid.nD][i][j][k];
        }
        dRho_ip1half=((1.0-dDonorFrac_ip1half)*dRho_cen_ip1half+dDonorFrac_ip1half
          *dRho_upwind_ip1half);
        
        //calculate radial term
        dDeltaRhoR=dUmU0_im1halfjk_np1half*dRho_im1half*dA_im1half
          -dUmU0_ip1halfjk_np1half*dRho_ip1half*dA_ip1half;
        
        //calculate new density
        grid.dLocalGridNew[grid.nD][i][j][k]=dVRatio*grid.dLocalGridOld[grid.nD][i][j][k]
          +time.dDeltat_np1half*(dDeltaRhoR)/dV_np1;
        
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
        
        //add rho
        ssName.str("");
        ssName<<"rho"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,grid.dLocalGridOld[grid.nD][i][j][k]);
        
        //add DeltaRhoDt_R
        ssName.str("");
        ssName<<"DeltaRhoDt_R"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,time.dDeltat_np1half*(dDeltaRhoR)/dV_np1);
        
        //add DeltaRhoDt_Vol
        ssName.str("");
        ssName<<"DeltaRhoDt_Vol"<<ssEnd.str();
        parameters.profileDataDebug.setMaxAbs(ssName.str()
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,dVRatio*grid.dLocalGridOld[grid.nD][i][j][k]);
        #endif
        
        if(grid.dLocalGridNew[grid.nD][i][j][k]<0.0){
          
          #if SIGNEGDEN==1
          raise(SIGINT);
          #endif
          
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
            <<": negative density calculated in , ("<<i<<","<<j<<","<<k<<")\n";
          throw exception2(ssTemp.str(),CALCULATION);
        }
      }
    }
  }
  
  //ghost region 0, outter most ghost region in x1 direction
  for(i=grid.nStartGhostUpdateExplicit[grid.nD][0][0]
    ;i<grid.nEndGhostUpdateExplicit[grid.nD][0][0];i++){
    
    //calculate i for interface centered quantities
    nIInt=i+grid.nCenIntOffset[0];
    dDelRCu_i_n=(pow(grid.dLocalGridOld[grid.nR][nIInt][0][0],3.0)
          -pow(grid.dLocalGridOld[grid.nR][nIInt-1][0][0],3.0));
    dDelRCu_i_np1=(pow(grid.dLocalGridNew[grid.nR][nIInt][0][0],3.0)
          -pow(grid.dLocalGridNew[grid.nR][nIInt-1][0][0],3.0));
    dR_ip1half_np1half=grid.dLocalGridOld[grid.nR][nIInt][0][0];//could be time centered
    dR_im1half_np1half=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];//could be time centered
    dRSq_ip1half_np1half=dR_ip1half_np1half*dR_ip1half_np1half;
    dRSq_im1half_np1half=dR_im1half_np1half*dR_im1half_np1half;
    dDelRSq_i_np1half=dRSq_ip1half_np1half-dRSq_im1half_np1half;
    dVRatio=dDelRCu_i_n/dDelRCu_i_np1;//calculate ratio of volume at n to volume at n+1
    dDonorFrac_im1half=(grid.dLocalGridOld[grid.nDonorCellFrac][i][0][0]
      +grid.dLocalGridOld[grid.nDonorCellFrac][i-1][0][0])*0.5;
    
    for(j=grid.nStartUpdateExplicit[grid.nD][1];j<grid.nEndUpdateExplicit[grid.nD][1];j++){
      
      for(k=grid.nStartUpdateExplicit[grid.nD][2];k<grid.nEndUpdateExplicit[grid.nD][2];k++){
        
        dV_np1=d1Thrid*dDelRCu_i_np1;
          
        //CALCULATE RATE OF CHANGE IN RHO IN RADIAL DIRECTION
        
        //calculate area at i-1/2
        dA_im1half=dRSq_im1half_np1half;
        
        //calculate difference between U and U0
        dUmU0_im1halfjk_np1half=grid.dLocalGridNew[grid.nU][nIInt-1][j][k]
          -grid.dLocalGridNew[grid.nU0][nIInt-1][0][0];
        
        //calculate rho at i-1/2, not time centered
        dRho_cen_im1half=(grid.dLocalGridOld[grid.nD][i][j][k]
          +grid.dLocalGridOld[grid.nD][i-1][j][k])*0.5;
        if(dUmU0_im1halfjk_np1half<0.0){//moving from outside in
          dRho_upwind_im1half=grid.dLocalGridOld[grid.nD][i][j][k];
        }
        else{//moving from inside out
          dRho_upwind_im1half=grid.dLocalGridOld[grid.nD][i-1][j][k];
        }
        dRho_im1half=((1.0-dDonorFrac_im1half)*dRho_cen_im1half+dDonorFrac_im1half
          *dRho_upwind_im1half);
        
        //calculate radial term
        dDeltaRhoR=dUmU0_im1halfjk_np1half*dRho_im1half*dA_im1half;
        
        //calculate new density
        /**\BC doesn't allow mass flux through outter interface*/
        grid.dLocalGridNew[grid.nD][i][j][k]=dVRatio*grid.dLocalGridOld[grid.nD][i][j][k]
          +time.dDeltat_np1half*(dDeltaRhoR)/dV_np1;
        
        #if DEBUG_EQUATIONS==1
        
        int nGhostCells=1;
        if(procTop.nRank==0){
          nGhostCells=0;
        }
        
        //add rho
        parameters.profileDataDebug.setMaxAbs("rho"
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,grid.dLocalGridOld[grid.nD][i][j][k]);
        
        //add DeltaRhoDt_R
        parameters.profileDataDebug.setMaxAbs("DeltaRhoDt_R"
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,time.dDeltat_np1half*(dDeltaRhoR)/dV_np1);
        
        //add DeltaRhoDt_Vol
        parameters.profileDataDebug.setMaxAbs("DeltaRhoDt_Vol"
          ,i+grid.nGlobalGridPositionLocalGrid[0]-nGhostCells*grid.nNumGhostCells
          ,dVRatio*grid.dLocalGridOld[grid.nD][i][j][k]);
        #endif
        
        if(grid.dLocalGridNew[grid.nD][i][j][k]<0.0){
          
          #if SIGNEGDEN==1
          raise(SIGINT);
          #endif
          
          std::stringstream ssTemp;
          ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
            <<": negative density calculated in , ("<<i<<","<<j<<","<<k<<")\n";
          throw exception2(ssTemp.str(),CALCULATION);
        }
      }
    }
  }

  #if SEDOV==1
    
    //ghost region 1, inner most ghost region in x1 direction
    for(i=grid.nStartGhostUpdateExplicit[grid.nD][1][0];
      i<grid.nEndGhostUpdateExplicit[grid.nD][1][0];i++){
      
      //calculate i for interface centered quantities
      nIInt=i+grid.nCenIntOffset[0];
      dDelRCu_i_n=(pow(grid.dLocalGridOld[grid.nR][nIInt][0][0],3.0)
            -pow(grid.dLocalGridOld[grid.nR][nIInt-1][0][0],3.0));
      dDelRCu_i_np1=(pow(grid.dLocalGridNew[grid.nR][nIInt][0][0],3.0)
            -pow(grid.dLocalGridNew[grid.nR][nIInt-1][0][0],3.0));
      dR_ip1half_np1half=grid.dLocalGridOld[grid.nR][nIInt][0][0];//could be time centered
      dR_im1half_np1half=grid.dLocalGridOld[grid.nR][nIInt-1][0][0];//could be time centered
      dRSq_ip1half_np1half=dR_ip1half_np1half*dR_ip1half_np1half;
      dRSq_im1half_np1half=dR_im1half_np1half*dR_im1half_np1half;
      dDelRSq_i_np1half=dRSq_ip1half_np1half-dRSq_im1half_np1half;
      dVRatio=dDelRCu_i_n/dDelRCu_i_np1;//calculate ratio of volume at n to volume at n+1
      
      for(j=grid.nStartGhostUpdateExplicit[grid.nD][1][1];
        j<grid.nEndGhostUpdateExplicit[grid.nD][1][1];j++){
        
        for(k=grid.nStartGhostUpdateExplicit[grid.nD][1][2];
          k<grid.nEndGhostUpdateExplicit[grid.nD][1][2];k++){
          
          dV_np1=d1Thrid*dDelRCu_i_np1;
            
          //CALCULATE RATE OF CHANGE IN RHO IN RADIAL DIRECTION
          
          //calculate area at i+1/2
          dA_ip1half=dRSq_ip1half_np1half;
          
          //calculate rho at i+1/2, not time centered
          dRho_ip1half=(grid.dLocalGridOld[grid.nD][i][j][k]+grid.dLocalGridOld[grid.nD][i+1][j][k])
            *0.5;
          
          //calculate radial term
          dDeltaRhoR=-(grid.dLocalGridNew[grid.nU][nIInt][j][k]
            -grid.dLocalGridNew[grid.nU0][nIInt][0][0])*dRho_ip1half*dA_ip1half;
          
          //calculate new density
          grid.dLocalGridNew[grid.nD][i][j][k]=dVRatio*grid.dLocalGridOld[grid.nD][i][j][k]
            +time.dDeltat_np1half*(dDeltaRhoR)/dV_np1;
          
          if(grid.dLocalGridNew[grid.nD][i][j][k]<0.0){
            
            #if SIGNEGDEN==1
            raise(SIGINT);
            #endif
            
            std::stringstream ssTemp;
            ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
              <<": negative density calculated in , ("<<i<<","<<j<<","<<k<<")\n";
            throw exception2(ssTemp.str(),CALCULATION);
          }
        }
      }
    }
  #endif
}
