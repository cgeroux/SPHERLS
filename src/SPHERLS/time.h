/**
  @file
  
  Header file for the ProcTop class
  
*/

#ifndef TIME_H
#define TIME_H

class Time{
  public:
    double dDeltat_np1half; /**<
      The time step centered at \f$n+1/2\f$ in seconds. It is used for
      calculating new variables defined at time step \f$n\f$, e.g. the density \ref Grid::nD.
      */
    double dDeltat_nm1half;/**<
      The previously used timestep centered at \f$n-1/2\f$ in seconds. It 
      is used for calculating \ref dDeltat_n the \f$n\f$ centered time step.
      */
    double dDeltat_n; /**<
      The time step centered at \f$n\f$ in seconds. It is used for calculating
      new variables defined at time step \f$n+1/2\f$, e.g. the radial velocity \ref Grid::nU. This value is
      determined by averaging the current \ref Time::dDeltat_np1half, and the last 
      \ref Time::dDeltat_np1half.
      */
    double dt; /**<
      The current time of the simulation in seconds.
      */
    double dEndTime; /**<
      The end time of the current calculation in seconds.
      */
    int nEndTimeStep;/**<
      The last time step to calculate, will stop if the current time step is larger than this. The 
      default value is the largest integer of the system.
      */
    double dTimeStepFactor; /**<
      Used for determining the time step. It is the factor which the
      courrant time step is multiplied by in order to determine \ref Time::dDeltat_np1half.
      */
    int nTimeStepIndex; /**<
      An index indecating the current time step. An index of zero corresponds
      to a \ref Time::dt=0. \todo should probably make this an unsigned variable, and perhaps also
      use the keyword long to help ensure there are enough values. Often need 7 decimal places.
      */
    bool bVariableTimeStep;/**<
      If true a variable time step is used as specified by the Courant condition, times the \ref dTimeStepFactor.
      */
    double dConstTimeStep;/**<
      If set to a value other than 0, will use that constant time step in place of the courant time
      step.
      */
    double dPerChange;/**< 
      A percentage amount to allow the maximum horizontal temperture variation and radial, theta 
      and phi convective velocities to change by from one time step to the next. The time step is 
      reduced accordingly to keep this precent change intact.
    */
    double dDelRho_t_Rho_max;/**<
      Keeps track of the maximum relative change in density from one time step to the next.*/
    double dDelT_t_T_max;/**<
      Keeps track of the maximum relative change in temperature from one time step to the 
      next. This quantity is only tracked if the calculation is non-adiabatic, else the energy is
      tracked instead, see \ref Time::dDelE_t_E_max*/
    double dDelE_t_E_max;/**<
      Keeps track of the maximum relative change in energy from one time step to the 
      next. This quantity is only tracked if the calculation is adiabatic, else the temperature is
      tracked instead, see \ref Time::dDelT_t_T_max*/
    double dDelUmU0_t_UmU0_max;
    double dDelV_t_V_max;
    double dDelW_t_W_max;
    Time(); /**<
      Constructor for the class \ref Time.
      */
};/**@class Time
  This class manages information which pertains to time variables.
*/
#endif
