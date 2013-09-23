#!/usr/bin/env python
import os
import getopt
import sys
import glob
import optparse as op
import combine_bins
import paths
import disect_filename
from myExceptions import *

def addParserOptions(parser):
  """Adds command line options for this script to parser
  
  """
  
  parser.add_option("--remake-profiles",action="store_true"
    ,dest="remakeProfiles"
    ,help="Will remake profiles even if they already exist. [not default]."
    ,default=False)
  parser.add_option("-e",action="store", dest="eosFile",type="string"
    ,help="Can be used to over ride the equation of state file found in the "
    +"model dumps. [not default].",default=None)
  parser.add_option("-v",action="store_true", dest="extraProfileInfo"
    ,help="Will include (dlnP/dlnT)_rho, (dlnP/dlnRho)_T, and (dE/dT)_rho in "
    +"radial profile. These are useful for calculating adiabatic gradient."
    ,default=False)
  
  #add parser options form combine_bins
  combine_bins.addParserOptions(parser)
def main():
  #make command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END]"
    ,version="%prog 1.0",description="Creates radial profiles from a range of "
    +"combined binary files starting with BASEFILENAME with time step indices "
    +"including START up to but not including END. Ex. there are a number of "
    +"files in the current directory from run1_t00000000, to run1_t00012300, "
    +"one could combine all these files with the command %prog run1_t[0-12300] "
    +"or %prog run1_t[0-*]. The \"*\" wild character is used to include all "
    +"files up to the maximum integer starting from START.")
  addParserOptions(parser)
  
  #parse command line options
  (options,args)=parser.parse_args()
    
  #should have at least one argument
  if len(args)!=1 :
    parser.error(" need one and one argument only.")
    
  #create profile files, and save list of files
  failedFiles=make_fileSet(args[0],options)
  
  if __name__=="__main__":#keeps from redundantly reporting errors
    #report profiles failed files
    for file in failedFiles:
      print file
def make_profile(fileName,options,nCount,nNumFiles):
  """Makes a radial profile form combined binary file fileName using SPHERLSanal
  
  uses:
  options.eosFile: if not None this equation of state file is used over the one
    indicated in the model dump.
  options.extraProfileInfo: add extra equation of state information to the 
    radial profile such as various gradients.
  options.remakeProfiles: if true it will remake the radial profile even if it
    already exists.
  """
  
  #if profile not already made for this combined binary file
  if not (os.path.exists(fileName+"_pro.txt")) or options.remakeProfiles:
    
    #make profile
    print (__name__+":"+make_profile.__name__+": creating profile from \""
      +fileName+"\" "+str(nCount)+"/"+str(nNumFiles)+" ...")
    cmd=paths.SPHERLSanalPath
    if options.extraProfileInfo:
      cmd+=" -v "
    if options.eosFile!=None:
      cmd+=" -e "+options.eosFile
    cmd+=' -a cb '+fileName
    success=os.system(cmd)
    if success!=0:#creating profile failed
      raise FileCreateFailed("error making profile "+fileName+"_pro.txt")
  else:
    print __name__+":"+make_profile.__name__+": profile \""+fileName\
      +"_pro.txt\" already exists, not remaking"
  
  return None#conversion successful
def make_fileSet(fileName,options,makeFileFunction=make_profile,frequency=1):
  """Makes a set of files from SPHERLS output using fileMakingFunction
  
  fileName: is expected to be the name of the inputfile
  options: command line options
  makeFileFunction: function to be used to create the a new file in the set.
  
  returns a list of files that failed being made
  """
  
  #check that we have a file name
  if len(fileName)<1:
    raise Exception("Requires an input file")
  
  #get base file name
  [start,end,baseFileName]=disect_filename.disectFileName(fileName)
  
  
  #make sure that all the distributed binary files in this range have been 
  #converted to combined binary files
  failedFiles=combine_bins.combine_bin_files(options.keep,fileName
    ,options.remakeBins)
  
  #check for combined binary files in iteration range starting with baseFileName
  filesExistCombBin=glob.glob(baseFileName
    +"[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]")
  files=[]
  
  for file in filesExistCombBin:
    intOfFile=int(file[len(baseFileName):len(file)])
    if intOfFile>=start and intOfFile<end:
      files.append(file)
  
  files.sort()
  if len(files)==0:
    raise NoFilesFound("no combined binary files found in range ["+str(start)
      +"-"+str(end)+"]")
    
  
  nNumFiles=len(files)
  nCount=1
  for i in range(0,len(files),frequency):
    try:
      makeFileFunction(files[i],options,nCount,nNumFiles)
    except FileCreateFailed as e:
      failedFiles.append(e.message)
    nCount+=1
  return failedFiles
if __name__ == "__main__":
  main()
  