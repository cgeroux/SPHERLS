#!/usr/bin/env python
import dump
import optparse as op
import numpy
import sys

def diffDumps(fileName0,fileName1,allowedRelativeDifference,thresholdValue,out):
  """Checks to see if the two dump files differ
  
  """
  
  dumpFile0=dump.Dump(fileName0)
  dumpFile1=dump.Dump(fileName1)
  
  #check headers
  [headersMatch,fatal]=diffHeader(dumpFile0,dumpFile1,allowedRelativeDifference
    ,thresholdValue,out)
  if fatal:
    out.write("Differences do not allow comparison of variables, stopping!\n")
    return False
  
  #check variables
  varsMatch=True
  for i in range(dumpFile0.numVars):
    tmp=diffVar(dumpFile0,dumpFile1,i,allowedRelativeDifference,thresholdValue
      ,out)
    if not tmp:
      varsMatch=False
  
  if not varsMatch or not headersMatch:
    out.write("files are not numerically identical\n")
    return False
  else:
    out.write("files are numerically identical\n")
    return True
def diffHeader(dumpFile0,dumpFile1,allowedRelativeDifference,threshold,out):
  """Checks to see if the two file headers differ
  
  If the file headers do differ it prints to out where they differ
  
  Arguments:
  dumpFile0: one of the two dump files to be compared
  dumpFile1: the second of the two dump files to be compared
  allowedRelativeDifference: the allowed relative difference between two 
                             numerical values being compared between the two
                             dump files. If the relative difference is less,
                             no difference is reported.
  threshold: how large the variable must be before the difference is considered.
             Usefull when comparing velocities for example.
  out: an output object which supports the a write() allowing output of the
       comparison to be written out.
  """
  
  headersMatch=True
  fatal=False
  
  #check type
  if dumpFile0.type!=dumpFile1.type:
    out.write("file \""+dumpFile0.fileName+"\" has type \""+str(dumpFile0.type)
      +"\" while file \""+dumpFile1.fileName+"\" has type \""
      +str(dumpFile1.type)+"\"\n")
    headersMatch=False
  
  #check versions
  if dumpFile0.version!=dumpFile1.version:
    out.write("file \""+dumpFile0.fileName+"\" has version \""
      +str(dumpFile0.version)+"\" while file \""+dumpFile1.fileName
      +"\" has version \""+str(dumpFile1.version)+"\"\n")
    headersMatch=False
  
  #check time
  denomitator=abs(dumpFile0.time+dumpFile1.time)*0.5
  if denomitator>threshold:
    relativeDifference=abs(dumpFile0.time-dumpFile1.time)/denomitator
    if relativeDifference>allowedRelativeDifference:
      out.write("file \""+dumpFile0.fileName+"\" has time \""+str(dumpFile0.time)
        +"\" while file \""+dumpFile1.fileName+"\" has time \""
        +str(dumpFile1.time)+"\"\n")
      headersMatch=False
  
  #check timestep index
  if dumpFile0.timeStepIndex!=dumpFile1.timeStepIndex:
    out.write("file \""+dumpFile0.fileName+"\" has timeStepIndex \""
      +str(dumpFile0.timeStepIndex)+"\" while file \""+dumpFile1.fileName
      +"\" has timeStepIndex \""+str(dumpFile1.timeStepIndex)+"\"\n")
    headersMatch=False
  
  #check delta_t_nm1half
  denomitator=abs(dumpFile0.delta_t_nm1half+dumpFile1.delta_t_nm1half)*0.5
  if denomitator>threshold:
    relativeDifference=abs(dumpFile0.delta_t_nm1half-dumpFile1.delta_t_nm1half)
      /denomitator
    if relativeDifference>allowedRelativeDifference:
      out.write("file \""+dumpFile0.fileName+"\" has delta_t_nm1half \""
        +str(dumpFile0.delta_t_nm1half)+"\" while file \""+dumpFile1.fileName
        +"\" has delta_t_nm1half \""+str(dumpFile1.delta_t_nm1half)+"\"\n")
      headersMatch=False
  
  #check delta_t_np1half
  denomitator=abs(dumpFile0.delta_t_np1half+dumpFile1.delta_t_np1half)*0.5
  if denomitator>threshold:
    relativeDifference=abs(dumpFile0.delta_t_np1half-dumpFile1.delta_t_np1half)
      /denomitator
    if relativeDifference>allowedRelativeDifference:
      out.write("file \""+dumpFile0.fileName+"\" has delta_t_np1half \""
        +str(dumpFile0.delta_t_np1half)+"\" while file \""+dumpFile1.fileName
        +"\" has delta_t_np1half \""+str(dumpFile1.delta_t_np1half)+"\"\n")
      headersMatch=False
  
  #check alpha
  denomitator=abs(dumpFile0.alpha+dumpFile1.alpha)*0.5
  if denomitator>threshold:
    relativeDifference=abs(dumpFile0.alpha-dumpFile1.alpha)
      /denomitator
    if relativeDifference>allowedRelativeDifference:
      out.write("file \""+dumpFile0.fileName+"\" has alpha \""
        +str(dumpFile0.alpha)+"\" while file \""+dumpFile1.fileName
        +"\" has alpha \""+str(dumpFile1.alpha)+"\"\n")
      headersMatch=False
  
  #check eosStringLen
  if dumpFile0.eosStringLen!=dumpFile1.eosStringLen:
    out.write("file \""+dumpFile0.fileName+"\" has eosStringLen \""
      +str(dumpFile0.eosStringLen)+"\" while file \""+dumpFile1.fileName
      +"\" has eosStringLen \""+str(dumpFile1.eosStringLen)+"\"\n")
    #headersMatch=False
    #these are allowed to be different, but are still reported
  
  #check eosString, either gamma value, or a file path
  if dumpFile0.eosStringLen>0:
    if dumpFile0.eosString!=dumpFile1.eosString:
      out.write("file \""+dumpFile0.fileName+"\" has eosString \""
        +str(dumpFile0.eosString)+"\" while file \""+dumpFile1.fileName
        +"\" has eosString \""+str(dumpFile1.eosString)+"\"\n")
      #headersMatch=False
      #file paths are allowed to be different, but are still reported
  else:
    if dumpFile0.gamma!=dumpFile1.gamma:
      out.write("file \""+dumpFile0.fileName+"\" has gamma \""
        +str(dumpFile0.gamma)+"\" while file \""+dumpFile1.fileName
        +"\" has gamma \""+str(dumpFile1.gamma)+"\"\n")
      headersMatch=False
  
  #check av (artificial viscosity)
  denomitator=abs(dumpFile0.av+dumpFile1.av)*0.5
  if denomitator>threshold:
    relativeDifference=abs(dumpFile0.av-dumpFile1.av)
      /denomitator
    if relativeDifference>allowedRelativeDifference:
      out.write("file \""+dumpFile0.fileName+"\" has av \""+str(dumpFile0.av)\
        +"\" while file \""+dumpFile1.fileName+"\" has av \""\
        +str(dumpFile1.av)+"\"\n")
      headersMatch=False
  
  #check avthreshold
  denomitator=abs(dumpFile0.avthreshold+dumpFile1.avthreshold)*0.5
  if denomitator>threshold:
    relativeDifference=abs(dumpFile0.avthreshold-dumpFile1.avthreshold)
      /denomitator
    if relativeDifference>allowedRelativeDifference:
      out.write("file \""+dumpFile0.fileName+"\" has avthreshold \""
        +str(dumpFile0.avthreshold)+"\" while file \""+dumpFile1.fileName
        +"\" has avthreshold \""+str(dumpFile1.avthreshold)+"\"\n")
      headersMatch=False
  
  #check global grid dimensions
  if dumpFile0.globalDims!=dumpFile1.globalDims:
    out.write("file \""+dumpFile0.fileName+"\" has globalDims \""
      +str(dumpFile0.globalDims)+"\" while file \""+dumpFile1.fileName
      +"\" has globalDims \""+str(dumpFile1.globalDims)+"\"\n")
    headersMatch=False
    fatal=True
  
  #check boundary conditions
  if dumpFile0.boundaryConditions!=dumpFile1.boundaryConditions:
    out.write("file \""+dumpFile0.fileName+"\" has boundaryConditions \""
      +str(dumpFile0.boundaryConditions)+"\" while file \""+dumpFile1.fileName
      +"\" has boundaryConditions \""+str(dumpFile1.boundaryConditions)+"\"\n")
    headersMatch=False
    fatal=True
  
  #check number of 1D zones
  if dumpFile0.num1DZones!=dumpFile1.num1DZones:
    out.write("file \""+dumpFile0.fileName+"\" has num1DZones \""
      +str(dumpFile0.num1DZones)+"\" while file \""+dumpFile1.fileName
      +"\" has num1DZones \""+str(dumpFile1.num1DZones)+"\"\n")
    headersMatch=False
    fatal=True
  
  #check number of ghost cells
  if dumpFile0.numGhostCells!=dumpFile1.numGhostCells:
    out.write("file \""+dumpFile0.fileName+"\" has numGhostCells \""
      +str(dumpFile0.numGhostCells)+"\" while file \""+dumpFile1.fileName
      +"\" has numGhostCells \""+str(dumpFile1.numGhostCells)+"\"\n")
    headersMatch=False
    fatal=True
  
  #check number of variables
  if dumpFile0.numVars!=dumpFile1.numVars:
    out.write("file \""+dumpFile0.fileName+"\" has numVars \""
      +str(dumpFile0.numVars)+"\" while file \""+dumpFile1.fileName
      +"\" has numVars \""+str(dumpFile1.numVars)+"\"\n")
    headersMatch=False
    fatal=True
  
  #check information about all variables
  if dumpFile0.varInfo!=dumpFile1.varInfo:
    out.write("file \""+dumpFile0.fileName+"\" has varInfo \""
      +str(dumpFile0.varInfo)+"\" while file \""+dumpFile1.fileName
      +"\" has varInfo \""+str(dumpFile1.varInfo)+"\"\n")
    headersMatch=False
    fatal=True
  
  #check 
  if dumpFile0.varSize!=dumpFile1.varSize:
    out.write("file \""+dumpFile0.fileName+"\" has varSize \""\
      +str(dumpFile0.varSize)\
      +"\" while file \""+dumpFile1.fileName+"\" has varSize \""\
      +str(dumpFile1.varSize)+"\"\n")
    headersMatch=False
    fatal=True
  
  #check number of dimensions
  if dumpFile0.numDims!=dumpFile1.numDims:
    out.write("file \""+dumpFile0.fileName+"\" has numDims \""\
      +str(dumpFile0.numDims)\
      +"\" while file \""+dumpFile1.fileName+"\" has numDims \""\
      +str(dumpFile1.numDims)+"\"\n")
    headersMatch=False
    fatal=True
  
  return [headersMatch,fatal]
def diffVar(dumpFile0,dumpFile1,var,allowedDifferance,threshold,out):
  """Checks if the variable, var, in the two dump files differ
  
  If the variables differ the difference is written to out
  
  Arguments:
  dumpFile0: one of the two dump files to be compared
  dumpFile1: the second of the two dump files to be compared
  var: the variable to compare
  allowedRelativeDifference: the allowed relative difference between two 
                             numerical values being compared between the two
                             dump files. If the relative difference is less,
                             no difference is reported.
  threshold: how large the variable must be before the difference is considered.
             Usefull when comparing velocities for example.
  out: an output object which supports the a write() allowing output of the
       comparison to be written out.
  
  """
  
  match=True
  
  #set ghost cells based on which directions the variable is defined in
  ghostCellsInX0=1
  if dumpFile0.varInfo[var][0]==-1:
    ghostCellsInX0=0
  ghostCellsInX1=1
  if dumpFile0.varInfo[var][1]==-1:
    ghostCellsInX1=0
  ghostCellsInX2=1
  if dumpFile0.varInfo[var][2]==-1:
    ghostCellsInX2=0
  
  #read 1D part
  sizeX01=ghostCellsInX0*(dumpFile0.num1DZones+dumpFile0.numGhostCells)
  if dumpFile0.varInfo[var][0]==1 and dumpFile0.boundaryConditions[0]==0:
    sizeX01=ghostCellsInX0*(dumpFile0.num1DZones+dumpFile0.numGhostCells+1)
  sizeX1=1
  sizeX2=1
  for i in range(sizeX01):
    for j in range(sizeX1):
      for k in range(sizeX2):
        denomitator=abs(dumpFile0.vars[var][i][j][k]+dumpFile1.vars[var][i][j][k])*0.5
        if denomitator>threshold:
          numerator=abs(dumpFile0.vars[var][i][j][k]-dumpFile1.vars[var][i][j][k])
          relativeDiffernace=numerator/denomitator
          if relativeDiffernace>allowedDifferance:
            lineDiff="var "+str(dumpFile0.varNames[var])+" in zone ("+str(i)+","+str(j)+","+str(k)+\
              ") has relative difference of {0:.2e} and magnitudes {1:.15e} and {2:.15e}.\n"\
              .format(relativeDiffernace,dumpFile0.vars[var][i][j][k],dumpFile1.vars[var][i][j][k])
            out.write(lineDiff)
            match=False
  
  #read multi-D part
  sizeX02=dumpFile0.varSize[var][0]+ghostCellsInX0*2*dumpFile0.numGhostCells
  sizeX1=dumpFile0.varSize[var][1]+ghostCellsInX1*2*dumpFile0.numGhostCells
  sizeX2=dumpFile0.varSize[var][2]+ghostCellsInX2*2*dumpFile0.numGhostCells
  for i in range(sizeX01,sizeX02):
    for j in range(sizeX1):
      for k in range(sizeX2):
        denomitator=abs(dumpFile0.vars[var][i][j][k]+dumpFile1.vars[var][i][j][k])*0.5
        if denomitator>threshold:
          numerator=abs(dumpFile0.vars[var][i][j][k]-dumpFile1.vars[var][i][j][k])
          relativeDiffernace=numerator/denomitator
          if relativeDiffernace>allowedDifferance:
            lineDiff="var "+str(dumpFile0.varNames[var])+" in zone ("+str(i)+","+str(j)+","+str(k)+\
              ") has relative difference of {0:.2e} and magnitudes {1:.15e} and {2:.15e}.\n"\
              .format(relativeDiffernace,dumpFile0.vars[var][i][j][k],dumpFile1.vars[var][i][j][k])
            out.write(lineDiff)
            match=False
  return match
def main():
  parser=op.OptionParser(usage="Usage: %prog [options]"
    ,version="%prog 1.0"
    ,description="Tests to dump files to see if they differ.")
    
  #parse command line options
  (options,args)=parser.parse_args()
  
  if len(args)<2:
    print "must supply two files names"
    return False
  
  #test out some simple uses for this class
  diffDumps(args[0],args[1],5e-14,1e-17,sys.stdout)
if __name__ == "__main__":
  main()