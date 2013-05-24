#ifndef EOS_H
#define EOS_H

/** @file
  
  Header file for \ref eos.cpp
*/

#include <string>
#include "exception2.h"
class eos{
  public:
    
    //member variables
    int nNumRho;/**<
      Number of densities in the equation of state table
      */
    int nNumT;/**<
      Number of temperatures in the equation of state table
      */
    double dXMassFrac;/**<
      Hydrogen mass fraction of the composition used to generate the equation of state table.
      */
    double dYMassFrac;/**<
      Helium mass fraction of the composition used to generate the equation of state table.
      */
    double dLogRhoMin;/**<
      Minimum density of the table in log10.
      */
    double dLogRhoDelta;/**<
      Increment of the density between table entries in log10.
      */
    double dLogTMin;/**<
      Minimum temperature of the table in log10.
      */
    double dLogTDelta;/**<
      Increment of the temperature between table entries in log10.
      */
    double **dLogP;/**<
      2D array of log10 pressures. dLogP[i][j] gives the log10 pressure at
      log10 density of \ref eos::dLogRhoDelta*i+\ref eos::dLogRhoMin, and at log10 temperature of
      \ref eos::dLogTDelta*j+\ref eos::dLogTMin.
      */
    double **dLogE;/**<
      2D array of log10 energies. dLogE[i][j] gives the log10 energy at
      log10 density of \ref eos::dLogRhoDelta*i+\ref eos::dLogRhoMin, and at log10 temperature of
      \ref eos::dLogTDelta*j+\ref eos::dLogTMin.
      */
    double **dLogKappa;/**<
      2D array of log10 opacities. dLogKappa[i][j] gives 
      the log10 opacity at log10 density of \ref eos::dLogRhoDelta*i+\ref eos::dLogRhoMin, 
      and at log10 temperature of \ref eos::dLogTDelta*j+\ref eos::dLogTMin.
      */
    
    //member functions
    eos();/**<
      Constructor, doesn't really do anything
      */
    eos(int nNumT,int nNumRho);/**<
      Constructor, allocates memory for the 2D arrays
      
      @param[in] nNumT number of temperatures in the equaiton of state table
      @param[in] nNumRho number of densities in the equaiton of state table
      
      Note: I don't think this version of the constructor is implemented,
      I should probably get rid of this definition
      */
    eos(const eos &ref);/**<
      Copy constructor, simply constructs a new eos object from another eos object
      */
    ~eos();/**<
      Destructor, delets dynamic arrays
      */
    eos& operator=(const eos & eosRightSide);/**<
      Assignment operator, assigns one eos object to another.
      */
    void readAscii(std::string sFileName)throw(exception2);/**<
      This fuction reads in an ascii file and stores it in the current object.
      
      @param[in] sFileName name of the equation of state file to read from.
      */
    void readBobsAscii(std::string sFileName)throw(exception2);/**<
      This fuction reads in an ascii file and stores it in the current object. The ascii file is in 
      Bob's format.
      
      @param[in] sFileName name of the equation of state file to read from.
      */
    void writeAscii(std::string sFileName)throw(exception2);/**<
      This fuction writes the equation of state stored in the current object to an ascii file.
      @param[in] sFileName name of the file to write the equation of state to.
      */
    void readBin(std::string sFileName)throw(exception2);/**<
      This fuction reads in a binary file and stores it in the current object.
      @param[in] sFileName name of the equation of state file to read from.
      */
    void writeBin(std::string sFileName)throw(exception2);/**<
      This fuction writes the equation of state stored in the current object to a binary file.
      @param[in] sFileName name of the file to write the equaiton of state to.
      */
    double dGetPressure(double dT, double dRho)throw(exception2);/**<
      This function linearly interpolates the pressure to a given temperature and density. 
      Note that both \c dT and \c dRho are not in log space.
      
      @param[in] dT temperature to interpolate to.
      @param[in] dRho density to interpolate to.
      @return the interpolated pressure.
      */
    double dGetEnergy(double dT, double dRho)throw(exception2);/**<
      This function linearly interpolates the energy to a given temperature and and density. 
      Note that both \c dT and \c dRho are not in log space.
      
      @param[in] dT temperature to interpolate to.
      @param[in] dRho density to interpolate to.
      @return the interpolated energy.
      */
    double dGetOpacity(double dT, double dRho)throw(exception2);/**<
      This function linearly interpolates the opacity to a given temperature and and density. 
      Note that both \c dT and \c dRho are not in log space.
      
      @param [in] dT temperature to interpolate to.
      @param [in] dRho density to interpolate to.
      @return the interpolated opacity.
      */
    double dDRhoDP(double dT,double dRho)throw(exception2);/**<
      This function calculates the partial derivative of density w.r.t. pressure
      @param [in] dT temperature at which the derivative is to be computed
      @param [in] dRho density at which the derivative is to be computed
      @return the partial derivative of density w.r.t. pressure.
      */
    double dSoundSpeed(double dT,double dRho)throw(exception2);/**<
      This function calculates the adiabatic sound speed
      @param [in] dT temperature at which the derivative is to be computed
      @param [in] dRho density at which the derivative is to be computed
      @return the sound speed.
      */
    void getEKappa(double dT, double dRho, double &dE, double &dKappa)throw(exception2);/**<
      This function linearly interpolates the three dependent quantities (Pressure, Energy
      , Opacity) to a given temperature and density. Note that both \c dT and \c dRho are 
      not in log space.
      
      @param[in] dT temperature to interpolate to.
      @param[in] dRho density to interpolate to.
      @param[out] dE energy at dT and dRho.
      @param[out] dKappa opacity at dT and dRho.
      */
    void getPEKappa(double dT, double dRho, double &dP,double &dE, double &dKappa)throw(exception2);/**<
      This function linearly interpolates the three dependent quantities (Pressure, Energy
      , Opacity) to a given temperature and density. Note that both \c dT and \c dRho are 
      not in log space.
      
      @param[in] dT temperature to interpolate to.
      @param[in] dRho density to interpolate to.
      @param[out] dP pressure at dT and dRho.
      @param[out] dE energy at dT and dRho.
      @param[out] dKappa opacity at dT and dRho.
      */
    void getPEKappaGamma(double dT,double dRho,double &dP,double &dE,double &dKappa
      ,double &dGamma)throw(exception2);/**<
      This function linearly interpolates the energy and opacity to a given temperature and 
      density. Note that both \c dT and \c dRho are not in log space.
      
      @param[in] dT temperature to interpolate to.
      @param[in] dRho density to interpolate to.
      @param[out] dP pressure at dT and dRho.
      @param[out] dE energy at dT and dRho.
      @param[out] dKappa opacity at dT and dRho.
      @param[out] dGamma adiabatic index at dT and dRho.
      */
    void getPEKappaGammaCp(double dT,double dRho,double &dP,double &dE,double &dKappa
      ,double &dGamma,double &dCp)throw(exception2);/**<
      This function linearly interpolates the energy and opacity to a given temperature and 
      density. Note that both \c dT and \c dRho are not in log space.
      
      @param[in] dT temperature to interpolate to.
      @param[in] dRho density to interpolate to.
      @param[out] dP pressure at dT and dRho.
      @param[out] dE energy at dT and dRho.
      @param[out] dKappa opacity at dT and dRho.
      @param[out] dGamma adiabatic index at dT and dRho.
      @param[out] dCp specific heat at constant pressure at dT and dRho.
      */
    void getPKappaGamma(double dT, double dRho, double &dP, double &dKappa,double &dGamma)throw(exception2);/**<
      This function linearly interpolates the energy and opacity to a given temperature and 
      density. Note that both \c dT and \c dRho are not in log space.
      
      @param[in] dT temperature to interpolate to.
      @param[in] dRho density to interpolate to.
      @param[out] dP pressure at dT and dRho.
      @param[out] dKappa opacity at dT and dRho.
      @param[out] dGamma adiabatic index at dT and dRho.
      */
    void gamma1DelAdC_v(double dT,double dRho,double &dGamma1, double &dDelAd, double &dC_v)throw(exception2);/**<
      This function calculates gamma1 and the adiabatic gradient
      
      @param [in] dT temperature at which the derivative is to be computed
      @param [in] dRho density at which the derivative is to be computed
      @param [out] dGamma1 gamma1
      @param [out] dDelAd adiabatic gradient
      @param [out] dC_v specific heat at constant volume
      */
    void getPAndDRhoDP(double dT,double dRho,double &dP, double &dDRhoDP)throw(exception2);/**<
      This function calculates the partial derivative of density w.r.t. pressure
      and the pressure
      @param [in] dT temperature at which the derivative is to be computed
      @param [in] dRho density at which the derivative is to be computed
      @param [out] dP pressure at dT and dRho
      @param [out] dDRhoDP derivative of density w.r.t. pressure at conatant temperature
      */
    void getEAndDTDE(double dT,double dRho,double &dE, double & dDTDE)throw(exception2);/**<
      This function calculates the partial derivative of temperature w.r.t. energy
      and the energy
      @param [in] dT temperature at which the derivative is to be computed
      @param [in] dRho density at which the derivative is to be computed
      @param [out] dE energy at dT and dRho
      @param [out] dDTDE derivative of temperature w.r.t. energy at constant density
      */
    void getDlnPDlnTDlnPDlnPDEDT(double dT, double dRho, double &dDlnPDlnT, double &dDlnPDlnRho,
      double &dDEDT)throw(exception2);/**<
        This function calculates various partial derivatives
        @param [in] dT temperature at which the derivative is to be computed
        @param [in] dRho density at which the derivative is to be computed
        @param [out] dDlnPDlnT derivative of ln(P) w.r.t. ln(T)
        @param [out] dDlnPDlnRho derivative of ln(P) w.r.t. ln(Rho)
        @param [out] dDEDT derivative of temperature w.r.t. energy at constant density
      */
};/**@class eos
  This class holds an equation of state as well as many functions useful for manipulating it
  */

#endif
