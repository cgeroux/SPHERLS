/**
  @file
  
  Implementation file for the Time class
  
*/

#include "time.h"
#include <limits>

Time::Time(){
  dDeltat_np1half=0.0;
  dDeltat_n=0.0;
  dt=0.0;
  dEndTime=0.0;
  dTimeStepFactor=0.0;
  nTimeStepIndex=0;
  dPerChange=1e-2;
  dDelRho_t_Rho_max=0.0;
  dDelT_t_T_max=0.0;
  dDelE_t_E_max=0.0;
  dDelUmU0_t_UmU0_max=0.0;
  dDelV_t_V_max=0.0;
  dDelW_t_W_max=0.0;
  nEndTimeStep=std::numeric_limits<int>::max();
}
