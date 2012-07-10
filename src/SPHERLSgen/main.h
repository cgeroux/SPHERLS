#ifndef MAIN_H
#define MAIN_H

/** @file
  
  Header file for \ref main.cpp.
  
*/

#include <vector>
#include "eos.h"

//variables
unsigned int nNumR;/**<
  Number of radial zones. When generating a stellar model this value is dependent on the depth to 
  which the model is integrated to (\ref dRStop), the initial mass spacing (\ref dMDelta) and the 
  rate at which the mass spacing is changed (\ref dMDeltaDelta). If generating a sedov model for 
  testing this value is explicitly set by the user.
  */
unsigned int nNumTheta;/**<
  Number of theta zones. Set in the configuration file under the 
  "<dimensions>" node using the "<num-theta>" tag.
  */
unsigned int nNumPhi;/**<
  Number of phi zones. Set in the configuration file under the 
  "<dimensions>" node using the "<num-phi>" tag.
  */
unsigned int nNumDims;/**<
  number of dimensions for the output model.
  */
double dMSun;/**<
  Value to use for a solar mass [g]. Set in the configuration file with the 
  "<M-sun>" tag.
  */
double dRSun;/**<
  Value to use for a solar radius [cm]. Set in the configuration file with the
  "<R-sun>" tag.
  */
double dLSun;/**<
  Value to use for a solar luminosity [ergs s^-1]. Set in the configuration file 
  with the "<L-sun>" tag.
  */
double dSigma;/**<
  Stefan-Boltzmann constant [ergs s^-1 cm^-2 K^-4]. Set in the configuration file
  with the "<sigma>" tag.
  */
double dRSurf;/**<
  surface radius of the model to be computed [solar radii].
  */
double dG;/**<
  Gravitational constant. Set in the configuration file with the "<G>" tag.
  */
double dPi=3.1415926535897932384626433832795;/**<
  Pi to a rediculous number of decimals
  */
double dL;/**<
  Luminosity of the model to be computed [L_sun]. Set in the configuration file with
  the "<L>" tag.
  */
double dTeff;/**<
  Effective temperature of the model to be computed [K]. Set in the configuration
  file with the "<T-eff>" tag.
  */
double dMTotal;/**<
  Total mass of the model to be computed [solar masses]. Set in the configuration
  file with the "<M-total>" tag.
  */
double dMDelta;/**<
  Initial mass step size [solar masses]. Set in the configuration file with the
  "<M-delta>" tag.
  */
double dRDelta;/**<
  Radial spacing when generating a sedov test model, otherwise it is undefined.
  */
double dRMin;/**<
  Minimum radius when generating a sedov test model, otherwise it is undefined.
  */
double dMDeltaDelta;/**<
  Fraction to increase \ref dMDelta by each shell. Set in the configuration 
  file with the "<M-delta-delta>" tag.
  */
double dDeltaTheta;/**<
  The spacing to use between theta's. It is set in the configuration file
  under the "<dimensions>" node with the "<delta-theta>" tag. \todo At some point in the future
  it may be desirable to vary this value across theta.
  */
double dDeltaPhi;/**<
  The spacing to use between phi's. It is set in the configuration file
  under the "<dimensions>" node with the "<delta-phi>" tag.
  */
int nNumCellsCent;/**<
  The number of radial cells to include in the higher energy center of a sedov model.
  */
double dEngCent;/**<
  The central energy of a sedov model, This value is units of [ergs] and will be divided by the mass
  in the central region to produce the specific energy of the central region.
  */
double dEng;/**<
  The energy in the remained of a sedov model exluding the central region, which is calculated from
  \ref dEngCent.
  */
double dRho;/**<
  The density of the sedov model, set to be uniform throughout.
  */
double dTolerance;/**<
  Allowed tolerance when converging a model. This value is compared to the 
  relative error/change in a quantity while converging it e.x. \f$\delta\rho/\rho\f$ and when this 
  erro is smaller than dTolerance intartions cease. A good value is usually around 5e-15
  (about machine precision).
  */
std::string sUDistType="POLY";/**<
  Specifies the method used to generate the radial velocity
  profile. It is set in the configuration file as the "type" attribute in the "<velocityDist>" tag.
  Accepted values are "POLY" for a polynomial profile and "PRO" to use a tabulated
  velocity profile.
  
  \see term for more details about "PLOY" velocity profiles, and dUPro, dUProR, nNumUProPoints, and
  dUSurf for more details about "PRO" velocity profiles.
  */
unsigned int nNumUProPoints=0;/**<
  Number of points in the velocity profile arrays \ref dUPro and 
  \ref dUProR. This value is set in \ref readUProfile by reading it from the velocity profile file.
  */
double *dUPro;/**<
  Radial velocity profile, velocity array of size \ref nNumUProPoints. It is used
  when generating a velocity profile if \ref sUDistType="PRO". The array values are set in 
  \ref readUProfile by reading them from the velocity profile file.
  */
double *dUProR;/**<
  Radial vleocity profile, radii array of size \ref nNumUProPoints. It is used
  when generating a velocity profile if \ref sUDistType="PRO". The array values are set in 
  \ref readUProfile by reading them from the velocity profile file.
  */
double dUSurf=1.0;/**<
  Scaling for radial velocity profile when generating a velocity profile of 
  type \ref sUDistTypes="PRO". This will be the value of the radial velocity at the surface of the
  model, and sets the scaling of the rest of the velocities. The value of the variable is set in the
  configuration file under the "<velocityDist>" node with the "<uSurf>" tag.
  */
std::string sOutPutfile;/**<
  Name of the file to write the calculated model to. It is set in the 
  configuration file with the "<fileName>" tag.
  */
unsigned int nNumZones1D;/**<
  The number of zones to be included in the 1D region. If it is -1 all
  zones will be output in 1D. It is set in the configuration file with under the "<dimensions>" node
  with the "<num-1d>" tag.
  */
unsigned int nNumGhostCells;/**<
  Number of ghost cells to use while outputing model. It is set in the 
  configuration file with the "<num-ghost-cells>" tag.
  */
unsigned int nPrecision=14;/**<\
  Number of decimal places to include in the output model. It is set in the 
  configuration file with the "<precision>" tag.
  */
unsigned int nPeriodic[3]={0,0,0};/**<
  Indicates which directions should have periodic boundary conditions
  in the evolved model. The three array values are set under the "<periodic>" node with tags
  "<x0>", "<x1>" and "<x2>" for the three directions, where x0 is the radial direciton, x1 is the
  theta direction, and x2 is the phi direction.
  */
double dTimeStepFactor=1.0;/**<
  Fraction to multiply the courant timestep by when calculating the
  time step. This value should be the same at that used when following the hydrodynamics of the 
  static model with SPHERLS. The time step is included in the initial static model to make the
  model format compatible with that required for resuming calculations with SPHERLS. It is set
  in the configuration file with the "<timeStepFactor>" tag.
  */
double dAlpha=0.2;/**<
  Parameter used to indicate mass above outter most zone. It multiplies 
  \ref dMDelta to calucluate the amount of mass lying above the outter interface. Appropriate values
  for this paramter are between 0 and 0.5. This value is set in the configuration file with the
  "<alpha>" tag.
  */
std::vector<double> vecdM;/**<
  Holds independent variable, interior mass. They are interface
  centered.
  */
std::vector<double> vecdMDel;/**<
  Holds the mass of the shells. They are zone centered
  */
std::vector<double> vecdP;/**<
  Holds the pressures of the radial shells. They are zone centered.
  */
std::vector<double> vecdE;/**<
  Holds the internal energy of the radial shells. They are zone 
  centered.
  */
std::vector<double> vecdRho;/**<
  Holds the density of the radial shells. They are zone centered.
  */
std::vector<double> vecdR;/**<
  Holds the radius of the radial shells. They are interface centered.
  */
std::vector<double> vecdT;/**<
  Holds the temperature of the radial shells. They are zone centered.
  */
std::vector<double> vecdKappa;/**<
  Holds the opcacity of the radial shells.  They are zone centered.
  */
double*** dU;/**<
  Holds the radial velocity. They are interface centered.
  */
double* dU0;/**<
  Holds the radial grid velocity. It is the same as the unpreturbed radial velocity \ref dU.
  **/
double*** dV;/**<
  Holds the theta velocity, assumed to be zero and are centered on 
  theta interfaces.
  */
double*** dW;/**<
  Holds the phi velocity, assumed to be zero and are centered on
  phi interfaces.
  */
struct term{
  double dCoeff;/**< Coeffeicent of the term in the polynomial.
    */
  double dPower;/**< Power of the term in the polynomial.
    */
};/**@struct term
  Structured variable holding a coeffecient and a power. Multiple terms combined together
  constructing a polynomial used for approximating the initial velocity profile. 
  \see vectVelDist, sUDistType
  */
std::vector<term> vectVelDist;/**<
  Holds mulitple \ref term's used to make up a polynomail to 
  calcuate the initial radial velocity profile of the model.
  \see term, sUDistType
  */
struct MDeltaDelta{
  std::string sStopType;
  double dMDeltaDelta;
  double dStopValue;
};
std::vector<MDeltaDelta> vecMDeltaDeltaList;
eos eosTable;/**<
  It is of type eos and holds the equation of state and opacity information and 
  functions used to provide a tabulated equation of state.
  \see SPHERLS reference manual for a more in depth description of the eos class.
  */
bool bGammaLawEOS=true;/**<
  if true will cause model to be generated using a gamma law gas
  */
double dGamma;/**<
  adiabatic gamma, 5/3=1.66 is stable to convection, will want to vary with 
  depth in future versions
  */
int nNumEProPoints=0;/**<
  number of points in the energy profile
  */
double *dEPro;/**<
  energy profile for interpolation [erg]
  */
double *dEProM;/**<
  mass points along energy profile [g]
  */
std::string sEOSFile;/**<
  filename for the equation of state table file
  */
bool bBinaryOutput;/**<
  If true the output file is in binary format. If false it outputs an ascii file.
  */
bool bWriteToScreen;/**<
  If true the model is written to the screen, else model is only printed to a file.
  */
bool bIsSedov;/**<
  If true it will generate a starting model to preform a Sedov (blast wave) test of the code.
  */
bool bAutoDeltaM;/**<
  If true it will use an algorithm to choose the mass spacing.
  */
//functions
void readConfig(std::string sConfigFileName,std::string sStartNode);/**<
  Reads in an xml configuration file and sets the values of many global variables. 
  \see main.h for variables that can be set in the configuration file and what the tag names specify
  those variables.
  
  @param [in] sConfigFileName filename and path to configuration file. Often "SPHERLSgen.xml"
  @param [in] sStartNode name of the starting node in the configuration file. Often "data"
  */
void readUProfile(std::string sProfileFileName);/**<
  This function reads in the radial velocity profile. The radial velocities are stored in \ref dUPro
  , the radii of those points are stored in \ref dUProR, and the size of these arrays is given by
  \ref nNumUProPoints. These radial velocities are used to interpolate a radial velocity profile for
  the output model.
  
  @param [in] sProfileFileName the name of the file containing the radial veolcity profile.
  */
void generateModel_SEDOV();/**<
  This fuction drives model generation for a sedov spherical blast wave model. It does this by 
  calling \ref calculateFirstShell_SEDOV to calculate the first shell of the model, then calling
  \ref calculateShell_SEDOV reapeatadly until all the zones are set.
  */
void generateModel_TEOS();/**<
  This function drives model generateration for a spherical static stellar model. It 
  does this by calling \ref calculateFirstShell_TEOS to calculate the first shell of the model, 
  then calling \ref calculateShell_TEOS repeatedly to generate the rest of the radial shells until
  the desired depth is reached (indicated by \ref dRStop). Finally a radial velocity profile is 
  generated by calling \ref makeVelocityDist.
  
  \see vecdM, vecdMDel, vecdP,vecdE,vecdRho,vecdR,vecdT,vecdKappa,vecdU,vecdV,vecdW for a 
  discription of the quantities calculated in the model.
  */
void generateModel_GL();/**<
  This function drives model generateration the spherical static stellar model. It 
  does this by calling \ref calculateFirstShell to calculate the first shell of the model, 
  then calling \ref calculateShell repeatedly to generate the rest of the radial shells until the 
  desired depth is reached (indicated by \ref dRStop). Finally a radial velocity profile is 
  generated by calling \ref makeVelocityDist.
  
  \see vecdM, vecdMDel, vecdP,vecdE,vecdRho,vecdR,vecdT,vecdKappa,vecdU,vecdV,vecdW for a 
  discription of the quantities calculated in the model.
  */
void calculateShell_SEDOV(int nShell);/**<
  
  */
void calculateShell_TEOS(unsigned int nShell);/**<
  
  */
void calculateFirstShell_SEDOV();/**<
  Calcualtes the first shell of a spherical blast wave model, or in otherwords a sedov test model.
  */
void calculateFirstShell_TEOS();/**<
  Calculates the quantites required in the first shell of the model. This is done by 
  <ol>
    <li>setting the mass at the surface equal to \ref dMTotal</li>
    
    <li>setting the temperatuer at the surface equal to the surface temperature given by 
      \f$T_s=\sqrt[4]{2}T_{eff}\f$ where \f$T_{eff}\f$ is \ref dTeff.</li>
      
    <li>Calculating the outter radius from \f$R_s=\sqrt{\frac{L}{4\pi\sigma T_s^4}}\f$ where 
      \f$L\f$ is \ref dL*\ref dLSun, \f$\pi\f$ is \ref dPi, \f$\sigma\f$ is \ref dSigma.</li>
      
    <li>Calculate the mass step. In the first shell this is simply -1.0*\ref dMDelta*\ref dMSun</li>
    
    <li>Calculate the pressure at zone center by solving the static momentum conservatin equation
      in 1D for the pressure at the center of the first zone, \n
      \f$P=\frac{-GM_r}{8\pi r^4}\left(0.5+\alpha\right)\Delta M\f$\n
      where \f$\alpha\f$ is \ref dAlpha, \f$\Delta M\f$ is \ref vecdMDel[0], \f$r\f$ is 
      \ref vecdR[0], \f$M_r\f$ is \ref vecdM[0] and \f$G\f$ is \ref dG. The pressure was directly
      solved for by applying the boundary condition that \f$P=0\f$ at the outer most interface. This
      lead to the pressure at \ref vecdM[0] - \ref dAlpha *\ref dMDelta equal to the negative of the
      pressure at \ref vecdM[0] + 0.5*\ref vecdMDel[0], allowing for direct solution of the
      pressure. Note that \ref vecdMDel[0] is negative to indicate moving into the star.
      </li>
    
    <li>Find the density such that the calculated pressure is recovered at the surface temperature
      \f$T_s\f$. This is done by calculating the derivative \f$\frac{\partial \rho}{\partial P}\f$
      from the equation of state and using it to calculate a linear correction to the density. This
      correction is repeatedly applied to the density until the correction is smaller than 
      \ref dTolerance.
      </li>
    
    <li>Once the density is found, the energy (\ref vecdE[0]) and the opacity (\ref vecdKappa[0])
      are easily found from the equation of state table.</li>
    <li>Finally the mass at the inner interface, and the radius at the inner interface are 
      calculated. The radius is simply that required to produce a volume to give the density
      of the zone at the mass of the shell.</li>
  </ol>
  */
double dPartialDerivativeVar1(double dVar1, double dVar2,int nShell
  ,double (*dFunction)(double dVar1, double dVar2, int nShell),double dH, double &dError);/**<
  Computes the partial deriavative of a function of two variables with respect to the first
  variable.
  
  @param[in] dVar1 location at which derivative should be evaluted, and also the variable that the
  partial derivative is take with respect to.
  @param[in] dVar2 locaiton at which derivative should be evaluated.
  @param[in] dFunction pointer to the double precision function of two variables, dVar1 and dVar2.
  
  @return returns the partial derivative with respect to dVar1 of the given function
  */
double dPartialDerivativeVar2(double dVar1, double dVar2,int nShell
  ,double (*dFunction)(double dVar1, double dVar2, int nShell),double dH, double &dError);/**<
  Computes the partial deriavative of a function of two variables with respect to the second
  variable.
  
  @param[in] dVar1 location at which derivative should be evaluted, and also the variable that the
  partial derivative is take with respect to.
  @param[in] dVar2 locaiton at which derivative should be evaluated.
  @param[in] dFunction pointer to the double precision function of two variables, dVar1 and dVar2.
  
  @return returns the partial derivative with respect to dVar1 of the given function
  */
void writeModel_R_TEOS();/**<
  Writes out the model generated using a tabulated equation of state in 1D to an ascii file.
    
  It writes the model out in the format accepted by SPHERLS. It writes the model the file file
  \ref sOutPutfile
  */
void writeModel_Bin_R_TEOS();/**<
  Writes out the model generated using a tabulated equation of state in 1D to a binary file.
    
  It writes the model out in the format accepted by SPHERLS. It writes the model the file file
  \ref sOutPutfile
  */
void writeModel_RT_TEOS();/**<
  Writes out the model generated using a tabulated equation of state in 2D in radius and theta to an
  ascii file.
    
  It writes the model out in the format accepted by SPHERLS. It writes the model the file file 
  \ref sOutPutfile
  */
void writeModel_Bin_RT_TEOS();/**<
  Writes out the model generated using a tabulated equation of state in 2D in radius and theta to a
  binary file.
    
  It writes the model out in the format accepted by SPHERLS. It writes the model the file file 
  \ref sOutPutfile
  */
void writeModel_RTP_TEOS();/**<
  Writes out the model generated using a tabulated equation of state in 3D in radius theta and phi
  to an ascii file.
    
  It writes the model out in the format accepted by SPHERLS. It writes the model the file file 
  \ref sOutPutfile
  */
void writeModel_Bin_RTP_TEOS();/**<
  Writes out the model generated using a tabulated equation of state in 3D in radius theta and phi
  to a binary file.
    
  It writes the model out in the format accepted by SPHERLS. It writes the model the file file 
  \ref sOutPutfile
  */
void makeVelocityDist();/**<
  */
double dTimeStep_TEOS();/**<
  */
double dMomentumCons(double dT,double dRho,int nShell);/**<
  */
double dEnergyCons(double dT,double dRho,int nShell);/**<
  */
void writeModel_R_GL();/**<
  Writes out the model generated using a gamma law gas in 1D in radius to an ascii file.
    
  It writes the model out in the format accepted by SPHERLS. It writes the model the file file 
  \ref sOutPutfile
  */
void writeModel_Bin_R_GL();/**<
  Writes out the model generated using a gamma law gas in 1D in radius to a binary file.
    
  It writes the model out in the format accepted by SPHERLS. It writes the model the file file 
  \ref sOutPutfile
  */
void writeModel_RT_GL();/**<
  Writes out the model generated using a gamma law gas in 2D in radius and theta to an ascii file.
    
  It writes the model out in the format accepted by SPHERLS. It writes the model the file file 
  \ref sOutPutfile
  */
void writeModel_Bin_RT_GL();/**<
  Writes out the model generated using a gamma law gas in 2D in radius and theta to a binary file.
    
  It writes the model out in the format accepted by SPHERLS. It writes the model the file file 
  \ref sOutPutfile
  */
void writeModel_RTP_GL();/**<
  Writes out the model generated using a gamma law gas in 3D in radius, theta and phi to an ascii 
  file.
    
  It writes the model out in the format accepted by SPHERLS. It writes the model the file file 
  \ref sOutPutfile
  */
void writeModel_Bin_RTP_GL();/**<
  Writes out the model generated using a gamma law gas in 3D in radius, theta and phi to a binary
  file.
    
  It writes the model out in the format accepted by SPHERLS. It writes the model the file file 
  \ref sOutPutfile
  */
double EOS_GL(double dP, double dE,int i);/**<
  */
double EOS_GL(double dP, double dE);/**<
  */
double interpolateE_GL(double dIntVar);/**<
  */
void writeModelToScreen_GL();/**<
  writes the calculated model to standard output
  */
void writeModelToScreen_TEOS();/**<
  writes the calculated model to standard output
  */
void readEnergyProfile_GL(std::string sProfileFileName);/**<
  */
void deleteEnergyProfile_GL();/**<
  */
void calculateFirstShell_GL();/**<
  */
void calculateShell_GL(unsigned int nShell);/**<
  */
double dTimeStep_GL();/**<
  */
void writeModel_GL();/**<
  */
double EOS_GL(double dP, double dE,int i);/**,
  */
double EOS_GL(double dP, double dE);/**<
  */
double interpolateE_GL(double dIntVar);/**<
  */
void generateModel_GL();/**<
  */
int main();/**<
  Main driving funciton for SPHERLSgen
  */
void makeVelocityDist_SEDOV();
void pretubeVelocityTorus(XMLNode xPerturb);

#endif
