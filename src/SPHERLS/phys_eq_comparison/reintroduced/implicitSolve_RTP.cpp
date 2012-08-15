void implicitSolve_RTP(Grid &grid,Implicit &implicit,Parameters &parameters,Time &time
  ,ProcTop &procTop,MessPass &messPass,Functions &functions){
  
  double *dValuesRHS=new double[implicit.nNumRowsALocal+implicit.nNumRowsALocalSB];
  int *nIndicesRHS=new int[implicit.nNumRowsALocal+implicit.nNumRowsALocalSB];
  
  //loop until corrections are small enough
  double dRelTError=std::numeric_limits<double>::max();
  int nNumIterations=0;
  int nI;
  int nJ;
  int nK;
  double dTemps[7];
  double dF_ijk_Tijk;
  double *dValues;
  double dF_ijk_Tijk1;
  double dF_ijk_Tip1;
  double dF_ijk_Tim1;
  double dF_ijk_Tjp1;
  double dF_ijk_Tjm1;
  double dF_ijk_Tkp1;
  double dF_ijk_Tkm1;
  double dF_ijk_Ti1;
  double dTemp2;
  double dRelTErrorLocal;
  while(dRelTError>implicit.dTolerance
    &&nNumIterations<implicit.nMaxNumIterations){
    //CALCULATE COEFFECIENT MATRIX AND RHS
    
    //calculate on inner grid
    for(int i=0;i<implicit.nNumRowsALocal;i++){//for each row
      nI=implicit.nLocFun[i][0];
      nJ=implicit.nLocFun[i][1];
      nK=implicit.nLocFun[i][2];
      
      dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
      dTemps[1]=grid.dLocalGridNew[grid.nT][nI+1][nJ][nK];
      dTemps[2]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
      dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
      dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
      dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
      dTemps[6]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
      
      #if DEBUG_EQUATIONS==1
      if(dRelTError<implicit.dTolerance*2.0e6){
        parameters.bSetThisCall=true;
      }
      #endif
      dF_ijk_Tijk=functions.fpImplicitEnergyFunction(grid,parameters,time,dTemps,nI,nJ,nK);
      #if DEBUG_EQUATIONS==1
      parameters.bSetThisCall=false;
      #endif
      
      dValuesRHS[i]=-1.0*dF_ijk_Tijk;
      nIndicesRHS[i]=implicit.nLocDer[i][0][0];
      dValues=new double[implicit.nNumDerPerRow[i]];
      for(int j=0;j<implicit.nNumDerPerRow[i];j++){//for each derivative
        
        switch(implicit.nTypeDer[i][j]){
          case 0 :{//calculate derivative of energy equation wrt. T at i
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI+1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[6]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tijk1=functions.fpImplicitEnergyFunction(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tijk1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ][nK]);
            break;
          }
          case 1 :{//calculate derivative of energy equation wrt. T at i+1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI+1][nJ][nK]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[6]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tip1=functions.fpImplicitEnergyFunction(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tip1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI+1][nJ][nK]);
            break;
          }
          case 2 :{//calculate derivative of energy equation wrt. T at i-1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI+1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[6]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tim1=functions.fpImplicitEnergyFunction(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tim1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI-1][nJ][nK]);
            break;
          }
          case 3 :{//calculate derivative of energy equation wrt. T at j+1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI+1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[6]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tjp1=functions.fpImplicitEnergyFunction(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tjp1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ+1][nK]);
            break;
          }
          case 4 :{//calculate derivative of energy equation wrt. T at j-1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI+1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[6]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tjm1=functions.fpImplicitEnergyFunction(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tjm1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ-1][nK]);
            break;
          }
          case 34 :{//calculate derivative of energy equation wrt. T at j+1 and j-1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI+1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[6]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tjp1=functions.fpImplicitEnergyFunction(grid,parameters,time,dTemps,nI,nJ,nK);
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK]*(1.0+implicit.dDerivativeStepFraction);
            dF_ijk_Tjm1=functions.fpImplicitEnergyFunction(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tjp1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ+1][nK])
              +(dF_ijk_Tjm1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ-1][nK]);
            break;
          }
          case 5 :{//calculate derivative of energy equation wrt. T at k+1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI+1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[6]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];

            dF_ijk_Tkp1=functions.fpImplicitEnergyFunction(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tkp1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ][nK+1]);
            break;
          }
          case 6 :{//calculate derivative of energy equation wrt. T at k-1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI+1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[6]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1]*(1.0+implicit.dDerivativeStepFraction);
            dF_ijk_Tkm1=functions.fpImplicitEnergyFunction(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tkm1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ][nK-1]);
            break;
          }
          case 56 :{//calculate derivative of energy equation wrt. T at k+1 and k-1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI+1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[6]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tkp1=functions.fpImplicitEnergyFunction(grid,parameters,time,dTemps,nI,nJ,nK);
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[6]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1]*(1.0+implicit.dDerivativeStepFraction);
            dF_ijk_Tkm1=functions.fpImplicitEnergyFunction(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tkp1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ][nK+1])
              +(dF_ijk_Tkm1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ][nK-1]);
            break;
          }
        }
      }
      MatSetValues(
        implicit.matCoeff,//matrix to set values in
        1,//number or rows
        &implicit.nLocDer[i][0][0],//global index of rows
        implicit.nNumDerPerRow[i],//number of columns
        implicit.nLocDer[i][1],//global index of column
        dValues,//logically two-dimensional array of values
        INSERT_VALUES);
      delete [] dValues;
    }
    
    //calculate at surface
    for(int i=implicit.nNumRowsALocal;i<implicit.nNumRowsALocal+implicit.nNumRowsALocalSB;i++){//for each row
      nI=implicit.nLocFun[i][0];
      nJ=implicit.nLocFun[i][1];
      nK=implicit.nLocFun[i][2];
      
      dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
      dTemps[1]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
      dTemps[2]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
      dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
      dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
      dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
      
      #if DEBUG_EQUATIONS==1
      if(dRelTError<implicit.dTolerance*2.0e6){
        parameters.bSetThisCall=true;
      }
      #endif
      dF_ijk_Tijk=functions.fpImplicitEnergyFunction_SB(grid,parameters,time,dTemps
        ,nI,nJ,nK);
      #if DEBUG_EQUATIONS==1
      parameters.bSetThisCall=false;
      #endif
      dValuesRHS[i]=-1.0*dF_ijk_Tijk;
      nIndicesRHS[i]=implicit.nLocDer[i][0][0];
      dValues=new double[implicit.nNumDerPerRow[i]];
      for(int j=0;j<implicit.nNumDerPerRow[i];j++){//for each derivative
        
        switch(implicit.nTypeDer[i][j]){
          case 0 :{//calculate derivative of energy equation wrt. T at i
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tijk1=functions.fpImplicitEnergyFunction_SB(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tijk1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ][nK]);
            break;
          }
          //case 1: no i+1 at surface
          case 2 :{//calculate derivative of energy equation wrt. T at i-1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tim1=functions.fpImplicitEnergyFunction_SB(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tim1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI-1][nJ][nK]);
            break;
          }
          case 3 :{//calculate derivative of energy equation wrt. T at j+1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tjp1=functions.fpImplicitEnergyFunction_SB(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tjp1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ+1][nK]);
            break;
          }
          case 4 :{//calculate derivative of energy equation wrt. T at j-1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tjm1=functions.fpImplicitEnergyFunction_SB(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tjm1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ-1][nK]);
            break;
          }
          case 34 :{//calculate derivative of energy equation wrt. T at j+1 and j-1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tjp1=functions.fpImplicitEnergyFunction_SB(grid,parameters,time,dTemps,nI,nJ,nK);
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK]*(1.0+implicit.dDerivativeStepFraction);
            dF_ijk_Tjm1=functions.fpImplicitEnergyFunction_SB(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tjp1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ+1][nK])
              +(dF_ijk_Tjm1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ-1][nK]);
            break;
          }
          case 5 :{//calculate derivative of energy equation wrt. T at k+1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tkp1=functions.fpImplicitEnergyFunction_SB(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tkp1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ][nK+1]);
            break;
          }
          case 6 :{//calculate derivative of energy equation wrt. T at k-1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1]*(1.0+implicit.dDerivativeStepFraction);
            dF_ijk_Tkm1=functions.fpImplicitEnergyFunction_SB(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tkm1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ][nK-1]);
            break;
          }
          case 56 :{//calculate derivative of energy equation wrt. T at k+1 and k-1
            dTemps[0]=grid.dLocalGridNew[grid.nT][nI][nJ][nK];
            dTemps[1]=grid.dLocalGridNew[grid.nT][nI-1][nJ][nK];
            dTemps[2]=grid.dLocalGridNew[grid.nT][nI][nJ+1][nK];
            dTemps[3]=grid.dLocalGridNew[grid.nT][nI][nJ-1][nK];
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1]*(1.0+implicit.dDerivativeStepFraction);
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1];
            dF_ijk_Tkp1=functions.fpImplicitEnergyFunction_SB(grid,parameters,time,dTemps,nI,nJ,nK);
            dTemps[4]=grid.dLocalGridNew[grid.nT][nI][nJ][nK+1];
            dTemps[5]=grid.dLocalGridNew[grid.nT][nI][nJ][nK-1]*(1.0+implicit.dDerivativeStepFraction);
            dF_ijk_Tkm1=functions.fpImplicitEnergyFunction_SB(grid,parameters,time,dTemps,nI,nJ,nK);
            dValues[j]=(dF_ijk_Tkp1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ][nK+1])
              +(dF_ijk_Tkm1-dF_ijk_Tijk)
              /(implicit.dDerivativeStepFraction*grid.dLocalGridNew[grid.nT][nI][nJ][nK-1]);
            break;
          }
        }
      }
      MatSetValues(
        implicit.matCoeff,//matrix to set values in
        1,//number or rows
        &implicit.nLocDer[i][0][0],//global index of rows
        implicit.nNumDerPerRow[i],//number of columns
        implicit.nLocDer[i][1],//global index of column
        dValues,//logically two-dimensional array of values
        INSERT_VALUES);
      delete [] dValues;
    }
    
    //assemble coeffecient matrix
    MatAssemblyBegin(implicit.matCoeff,MAT_FINAL_ASSEMBLY);
    
    //set values of the RHS
    VecSetValues(implicit.vecRHS
      ,implicit.nNumRowsALocal+implicit.nNumRowsALocalSB
      ,nIndicesRHS
      ,dValuesRHS
      ,INSERT_VALUES);
    
    VecAssemblyBegin(implicit.vecRHS);
    VecAssemblyEnd(implicit.vecRHS);
    MatAssemblyEnd(implicit.matCoeff,MAT_FINAL_ASSEMBLY);
    
    //solve system
    KSPSetOperators(implicit.kspContext,implicit.matCoeff,implicit.matCoeff,SAME_NONZERO_PATTERN);
    KSPSolve(implicit.kspContext,implicit.vecRHS,implicit.vecTCorrections);
    
    //get distributed corrections
    VecScatterBegin(implicit.vecscatTCorrections,implicit.vecTCorrections
      ,implicit.vecTCorrectionsLocal,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(implicit.vecscatTCorrections,implicit.vecTCorrections
      ,implicit.vecTCorrectionsLocal,INSERT_VALUES,SCATTER_FORWARD);
    VecGetArray(implicit.vecTCorrectionsLocal,&dValues);
    
    //apply corrections
    dRelTErrorLocal=0.0;
    for(int i=0;i<implicit.nNumRowsALocal+implicit.nNumRowsALocalSB;i++){
      nI=implicit.nLocFun[i][0];
      nJ=implicit.nLocFun[i][1];
      nK=implicit.nLocFun[i][2];
      grid.dLocalGridNew[grid.nT][nI][nJ][nK]+=dValues[i];
      if(grid.dLocalGridNew[grid.nT][nI][nJ][nK]<0.0){
        
        #if SIGNEGTEMP==1
        raise(SIGINT);
        #endif
        
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
          <<": negative temperature calculated in , ("<<nI<<","<<nJ<<","<<nK<<") on iteration "
          <<nNumIterations<<"\n";
        throw exception2(ssTemp.str(),CALCULATION);
        
      }
      dTemp2=fabs(dValues[i]/grid.dLocalGridNew[grid.nT][nI][nJ][nK]);
      if(dRelTErrorLocal<dTemp2){
        dRelTErrorLocal=dTemp2;
      }
    }
    
    updateLocalBoundariesNewGrid(grid.nT,procTop,messPass,grid);
    
    MPI::COMM_WORLD.Allreduce(&dRelTErrorLocal,&dRelTError,1,MPI::DOUBLE,MPI_MAX);
    
    VecRestoreArray(implicit.vecTCorrectionsLocal,&dValues);
    nNumIterations++;
  }
  
  #if TRACKMAXSOLVERERROR==1
    
    /* Calculate absolute error in solver*/
    Vec vecCalRHS1;
    Vec vecCalRHS2;
    VecDuplicate(implicit.vecRHS,&vecCalRHS1);
    VecDuplicate(implicit.vecRHS,&vecCalRHS2);
    MatMult(implicit.matCoeff,implicit.vecTCorrections,vecCalRHS1);
    VecCopy(vecCalRHS1,vecCalRHS2);
    
    VecAXPY(vecCalRHS1,-1.0,implicit.vecRHS);
    
    //get maximum absolute error, and average value of the RHS
    int nIndexLargestError;
    double dMaxError;
    VecMax(vecCalRHS1,&nIndexLargestError,&dMaxError);
    dMaxError=fabs(dMaxError);
    if(dMaxError>implicit.dMaxErrorInRHS){
      implicit.dMaxErrorInRHS=dMaxError;
      double dSumRHS=0.0;
      VecAbs(vecCalRHS2);
      VecSum(vecCalRHS2,&dSumRHS);
      int nSizeVecCalRHS;
      VecGetSize(vecCalRHS2,&nSizeVecCalRHS);
      implicit.dAverageRHS=dSumRHS/double(nSizeVecCalRHS);
    }
    
    //keep track of largest number of interations the solver took
    int nNumIterationsSolve;
    KSPGetIterationNumber(implicit.kspContext,&nNumIterationsSolve);
    if(nNumIterationsSolve>implicit.nMaxNumSolverIterations){
      implicit.nMaxNumSolverIterations=nNumIterationsSolve;
    }
  #endif
  
  delete [] dValuesRHS;
  delete [] nIndicesRHS;
  if(dRelTError>implicit.dCurrentRelTError){
    implicit.dCurrentRelTError=dRelTError;
  }
  if(nNumIterations>implicit.nCurrentNumIterations){
    implicit.nCurrentNumIterations=nNumIterations;
  }
  
  if(procTop.nRank==0){
    if(nNumIterations>=implicit.nMaxNumIterations){
      std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<procTop.nRank
      <<": The maximum number of iteration for implicit solution ("<<implicit.nMaxNumIterations
      <<") has be exceeded in current time step ("<<time.nTimeStepIndex
      <<") with a maximum relative error in the implicit calculation of temperature of "
      <<implicit.dCurrentRelTError<<std::endl;
    }
  }
  
  //calculate E, Kappa, P form new temperature
  calNewPEKappaGamma_TEOS(grid,parameters);
}
