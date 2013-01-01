/**
  @file
  
  This file contains the main function which is the driver for SPHERLS.
  
*/

#include <mpi.h>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <csignal>
#include <fenv.h>
#include "main.h"
#include "global.h"
#include "watchzone.h"
#include "exception2.h"
#include "xmlParser.h"
#include "xmlFunctions.h"
#include "dataManipulation.h"
#include "dataMonitoring.h"
#include "physEquations.h"

int main(int argc, char* argv[]){
  
  Global global;
  
  //initialize MPI
  MPI::Init(argc,argv);
  
  //set handler for Floatpoint Exceptions
  signal(SIGFPE, signalHandler);
  
  try{
    
    //Initialize program, read in starting model
    init(global.procTop,global.grid,global.output,global.time,global.parameters
      ,global.messPass,global.performance,global.implicit,argc,argv);
    
    //set function pointers to be used for calculations
    setMainFunctions(global.functions,global.procTop,global.parameters,global.grid,global.time
      ,global.implicit);
    
    //update new grid with old grid after read
    updateNewGridWithOld(global.grid,global.procTop);
    
    //update boundaries, needed here to make sure that the old grid has all internal variables that 
    //are updated in time, initialized in the ghost zones
    updateLocalBoundaries(global.procTop,global.messPass,global.grid);
    
    bool bFirstIterationDump=true;
    bool bFirstIterationPrint=true;
    if(global.output.nPrintMode==1&&global.procTop.nRank==0){//print out header if print time step info
      std::cout<<"Time_Step_Index"
        <<" "<<"Time"
        <<" "<<"Delta_Time"
        <<" "<<"max(Del_Rho/Rho)"
        <<" "<<"max(Del_T/T)"
        <<" "<<"max(Del_UmU0/UmU0)"
        <<" "<<"max(Del_V/V)"<<std::endl;
    }
    while(global.time.dt<=global.time.dEndTime&&global.time.nTimeStepIndex<=global.time.nEndTimeStep){
      
      //if bDump is true write out grid
      if(global.output.bDump){
        
        global.output.nNumTimeStepsSinceLastDump++;
        
        //decide if dumping this time step
        bool bDump=false;
        if(global.output.nDumpFrequencyStep!=0){
          if(global.time.nTimeStepIndex%global.output.nDumpFrequencyStep==0){
            bDump=true;
          }
        }
        if(global.output.dDumpFrequencyTime!=0.0){
          if(global.time.dt>=(global.output.dDumpFrequencyTime+global.output.dTimeLastDump)){
            bDump=true;
            global.output.dTimeLastDump=global.time.dt;
          }
        }
        
        //if dumping this time step, then dump
        if(bDump||bFirstIterationDump){
          std::stringstream ssFileNameOut;
          
          ssFileNameOut<<global.output.sBaseOutputFileName<<"_t"<<std::setfill('0')<<std::setw(8)
            <<global.time.nTimeStepIndex;
            
          if(global.procTop.nRank==0){
            std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<global.procTop.nRank<<":"
              <<std::endl<<"  Dumping model to file: "<<ssFileNameOut.str()<<std::endl;
          }
          
          global.output.nNumTimeStepsSinceLastDump=0;
          global.functions.fpModelWrite(ssFileNameOut.str(), global.procTop,global.grid,global.time
            ,global.parameters);
          
          #if DEBUG_EQUATIONS==1
          if(!bFirstIterationDump){//nothing to print on the first iteration
            std::stringstream ssFileNameProOut;
            ssFileNameProOut<<global.parameters.sDebugProfileOutput<<"_t"<<std::setfill('0')
              <<std::setw(8)<<global.time.nTimeStepIndex<<"_pro.txt";
            global.parameters.profileDataDebug.toFile(ssFileNameProOut.str(),global.time
              ,global.procTop);
            global.parameters.profileDataDebug.clear();
          }
          #endif
          
          bFirstIterationDump=false;
        }
      }
      
      //Print status
      if(global.output.bPrint){
        
        global.output.nNumTimeStepsSinceLastPrint++;
        
        //decide if printing this time step
        bool bPrint=false;
        if(global.output.nPrintFrequencyStep!=0){
          if(global.time.nTimeStepIndex%global.output.nPrintFrequencyStep==0){
            bPrint=true;
          }
        }
        if(global.output.dPrintFrequencyTime!=0.0){
          if(global.time.dt>=(global.output.dPrintFrequencyTime+global.output.dTimeLastPrint)){
            bPrint=true;
            global.output.dTimeLastPrint=global.time.dt;
          }
        }
        
        //if print this time step, then print
        if(bPrint||bFirstIterationPrint){
          bFirstIterationPrint=false;
          if(global.procTop.nRank==0){
            std::cout.setf(std::ios::scientific);
            std::cout.precision(14);
            if(global.output.nPrintMode==0){
              std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<global.procTop.nRank<<":"
                <<std::endl
                <<"  For current time step "<<global.time.nTimeStepIndex<<" at time t="
                <<global.time.dt<<" [s] with time step dt=";
              std::cout.precision(5);
              std::cout<<global.time.dDeltat_np1half<<" [s]"<<std::endl;
              std::cout
                <<"    * max( (Del Rho)_t/Rho )   ="<<global.time.dDelRho_t_Rho_max<<std::endl
                <<"    * max( (Del T)_t/T )       ="<<global.time.dDelT_t_T_max<<std::endl
                <<"    * max( (Del UmU0)_t/UmU0 ) ="<<global.time.dDelUmU0_t_UmU0_max<<std::endl;
              if(global.grid.nNumDims>1){//if 2D or more
                std::cout<<"    * max( (Del V)_t/V )       ="<<global.time.dDelV_t_V_max<<std::endl;
              }
              if(global.grid.nNumDims>2){//if 3D or more
                std::cout<<"    * max( (Del W)_t/W )       ="<<global.time.dDelW_t_W_max<<std::endl;
              }
              if(global.grid.nNumDims>1){
                std::cout<<"    * max(convective velocity) ="
                  <<global.parameters.dMaxConvectiveVelocity/1.0e5<<" [km/s]"<<std::endl;
              }
              if(global.implicit.nNumImplicitZones>0){
                std::cout
                  <<"  : Over last "<<global.output.nNumTimeStepsSinceLastPrint<<" time steps :\n"
                  <<"    * Largest number of iterations for temperature calculation was "
                  <<global.implicit.nCurrentNumIterations<<std::endl
                  <<"    * Largest relative error in calculation of temperature was     "
                  <<global.implicit.dCurrentRelTError<<std::endl;
                #if TRACKMAXSOLVERERROR==1
                  std::cout
                    <<"    * Largest number of iterations for solver was                  "
                    <<global.implicit.nMaxNumSolverIterations<<std::endl
                    <<"    * Largest absolute error in RHS of system was                  "
                    <<global.implicit.dMaxErrorInRHS<<std::endl
                    <<"    * with average RHS magnitude of                                "
                    <<global.implicit.dAverageRHS<<std::endl;
                #endif
              }
            }
            else if(global.output.nPrintMode==1){//print out time step diagnostic info
              std::cout<<global.time.nTimeStepIndex
                <<" "<<global.time.dt
                <<" "<<global.time.dDeltat_np1half
                <<" "<<global.time.dDelRho_t_Rho_max
                <<" "<<global.time.dDelT_t_T_max
                <<" "<<global.time.dDelUmU0_t_UmU0_max
                <<" "<<global.time.dDelV_t_V_max<<std::endl;
            }
          }
          
          #if TRACKMAXSOLVERERROR==1
          global.implicit.dMaxErrorInRHS=0.0;
          #endif
          global.output.nNumTimeStepsSinceLastPrint=0;
          global.implicit.dCurrentRelTError=0.0;
          global.implicit.nCurrentNumIterations=0;
        }
      }
      
      //output watch zone info
      global.functions.fpWriteWatchZones(global.output,global.grid,global.parameters, global.time
        , global.procTop);
      
      //calculate new velocities and update boundaries
      global.functions.fpCalculateNewVelocities(global.grid,global.parameters,global.time
        ,global.procTop);
      global.functions.fpUpdateLocalBoundaryVelocitiesNewGrid(global.procTop,global.messPass
        ,global.grid);
      
      //calculate new grid velocity and update boundaries
      global.functions.fpCalculateNewGridVelocities(global.grid,global.parameters,global.time
        ,global.procTop,global.messPass);
      updateLocalBoundariesNewGrid(global.grid.nU0,global.procTop,global.messPass,global.grid);
      
      //calculate new radius and update boundaries
      global.functions.fpCalculateNewRadii(global.grid,global.time);
      updateLocalBoundariesNewGrid(global.grid.nR,global.procTop,global.messPass,global.grid);
      
      //calculate new densities, and update boundaries
      global.functions.fpCalculateNewDensities(global.grid,global.parameters, global.time
        ,global.procTop);
      updateLocalBoundariesNewGrid(global.grid.nD,global.procTop,global.messPass,global.grid);
      
      //calculate horizontally averaged density
      global.functions.fpCalculateAveDensities(global.grid);
      updateLocalBoundariesNewGrid(global.grid.nDenAve,global.procTop,global.messPass,global.grid);
      
      //calculate new eddy viscosity
      global.functions.fpCalculateNewEddyVisc(global.grid,global.parameters);
      updateLocalBoundariesNewGrid(global.grid.nEddyVisc,global.procTop,global.messPass
        ,global.grid);
      
      //calculate new energies in explicit region
      global.functions.fpCalculateNewEnergies(global.grid,global.parameters, global.time
        ,global.procTop);
      
      //calculate new variables (T,Kappa,P, gamma) via equation of state in explicit region
      global.functions.fpCalculateNewEOSVars(global.grid,global.parameters);
      updateLocalBoundariesNewGrid(global.grid.nP,global.procTop,global.messPass,global.grid);
      updateLocalBoundariesNewGrid(global.grid.nGamma,global.procTop,global.messPass,global.grid);
      
      //update temperature at boundaries, need new temperature for implicit solution
      updateLocalBoundariesNewGrid(global.grid.nT,global.procTop,global.messPass,global.grid);
      
      //implicityly solve for T, and update (E,Kappa,P,Gamma) via equation of state in implicit region
      global.functions.fpImplicitSolve(global.grid,global.implicit,global.parameters,global.time
        ,global.procTop,global.messPass,global.functions);
      
      //calculate new artificial viscosity
      global.functions.fpCalculateNewAV(global.grid,global.parameters);
      
      //calculate timestep
      global.functions.fpCalculateDeltat(global.grid,global.parameters, global.time,global.procTop);
      
      //update boundaries remaining boundaries to old grid and copy new grid to old grid
      updateLocalBoundaries(global.procTop,global.messPass,global.grid);
      
    }
    
    global.output.nNumTimeStepsSinceLastDump++;
    
    //output watch zone info, for last timestep
    global.functions.fpWriteWatchZones(global.output,global.grid,global.parameters, global.time
      , global.procTop);
    
    //finish program by deleting dynamic memory
    fin(true,global.time, global.output, global.procTop,global.grid,global.parameters
      ,global.functions, global.performance, global.implicit);
  }
  
  //error handeling
  catch(exception2& eTemp){
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<global.procTop.nRank<<":"
      <<eTemp.getMsg();
    MPI::COMM_WORLD.Abort(1);
  }
  catch(std::exception& eTemp){
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<global.procTop.nRank<<":"
      <<"Standard exception:"<<eTemp.what()<<std::endl;
    MPI::COMM_WORLD.Abort(1);
  }
  catch(...){
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<global.procTop.nRank<<":"
      <<"main: unknown error\n";
    MPI::COMM_WORLD.Abort(1);
  }
  
  //finalize mpi
  MPI::Finalize();

  return 0;
}
void signalHandler(int nSig){
  if(nSig==SIGFPE){
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Floating point signal "<<nSig
      <<" detected.\n";
    MPI::COMM_WORLD.Abort(1);
    exit(EXIT_FAILURE);
    return;
  }
  else if(SIGINT){
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Interupt signal "<<nSig
      <<" detected.\n";
    MPI::COMM_WORLD.Abort(1);
    exit(EXIT_FAILURE);
    return;
  }
  else{
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": unknown signal "<<nSig
      <<" detected.\n";
    MPI::COMM_WORLD.Abort(1);
    exit(EXIT_FAILURE);
    return;
  }
}
