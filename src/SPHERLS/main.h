#ifndef MAIN_H
#define MAIN_H

/** @file
  
  Header file for \ref main.cpp
  */

void signalHandler(int nSig);/**<
  Used for catching signals. 
  */
int main(int argc, char* argv[]);/**<
  Main driving function of SPHERLS.
  
  @param[in] argc number of arguments passed from the command line
  @param[in] argv array of character strings of size argc containing
  the arguments from the command line.
  
  The flow of this function is as follows:
    - Initilize program by calling \ref init()
    - Set function pointers by calling \ref setMainFunctions()
    - Update new grid with old grid by calling \ref updateNewGridWithOld()
    - Update boundaries of local grids
    - Calculate the first time step by calling \ref Functions::fpCalculateDeltat()
    - Enter while loop until end time (\ref Time::dEndTime) is reached, and for each interation of
      the loop:
      - Test to see if a model dump is needed (by checking Output::bDump and 
        \ref Output::nDumpFrequency), if so dump one by calling \ref modelWrite()
      - Write out information for any watchzones present by calling
        \ref writeWatchZones()
      - Write out information for peak kinetic energy per period by calling
        \ref writePeakKE()
      - calculate time step by calling function pointed to by \ref Functions::fpCalculateDeltat
      
      - Calculate new velocities by calling the function pointed to by
        \ref Functions::fpCalculateNewVelocities()
      - Update velocities on new grid boundaries between processors by calling
        \ref updateLocalBoundariesNewGrid() three times indicating the \f$r\f$-velocity (\ref U),
        \f$\theta\f$-velocity (\ref V) and the \f$\phi\f$-velocity (\ref W).
      - Calculate new grid velocities with \ref Functions::fpCalculateNewGridVelocities().
      - Calculate new radii with \ref Functions::fpCalculateNewRadii().
      - Update radii on new grid boundaries between processors by calling 
        \ref updateLocalBoundariesNewGrid() indicating radius is to be
        updated (\ref R).
      - Calculate new densities with \ref Functions::fpCalculateNewDensities()
      - Calculate new energies with \ref Functions::fpCalculateNewEnergies()
      - Update the old grid boundaries and centeres by calling 
        \ref updateLocalBoundaries()
      - Calculating the next time step with \ref Functions::fpCalculateDeltat()
    - Finish by dumping the last model computed
*/
#endif
