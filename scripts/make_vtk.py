#!/usr/bin/env python
import os
import optparse as op
import sys
import combine_bins
from make_profiles import make_fileSet
import dump

def addParserOptions(parser):
  parser.add_option("--remake-vtk",action="store_true"
    ,dest="remakeVTK"
    ,help="Will remake vtk files even if they already exist. [not default]."
    ,default=False)
  parser.add_option("-e",action="store", dest="eosFile",type="string"
    ,help="Can be used to over ride the equation of state file found in the "
    +"model dumps. [not default].",default=None)
  parser.add_option("--no-curvature",action="store_false",dest="curvature"
    ,help="Will not calculate spherical curvature [not default].",default=True)
  
  #add options for combine_output parser
  combine_bins.addParserOptions(parser)
def makeVTKFile(fileName,options,nCount,nNumFiles):
  """Just a wrapper to writeRadialProfile member function of snapshot
  
  This allows different fileMakingFunctions to be specified. For example to
  create VTK files.
  """
    
  #if profile not already made for this combined binary file
  if not (os.path.exists(fileName+".vts")) or options.remakeVTK:
    
    #make vts file
    print (__name__+":"+makeVTKFile.__name__+": creating vts file from \""
      +fileName+"\" "+str(nCount)+"/"+str(nNumFiles)+" ...")
    
    #create object to hold model
    model=dump.Dump()
    
    #read in model
    model.read(fileName,eosFile=options.eosFile)
    
    #write out a vtk file from the model
    model.writeVTKFile(fileName,curvature=options.curvature
      ,eosFile=options.eosFile,includeScalars=["T","e","rho","kappa","c","cv","vr_con_cen"])
  else:
    print __name__+":"+makeVTKFile.__name__+": vts file \""+fileName\
      +".vts\" already exists, not remaking"
def main():
  #make command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END]"
    ,version="%prog 1.0",description="Creates vtk files (*.vts) from a range "
    +"of combined binary files starting with BASEFILENAME with time step "
    +"indices including START up to but not including END. Ex. there are a "
    +"number of files in the current directory from run1_t00000000, to "
    +"run1_t00012300, one could combine all these files with the command "
    +"%prog run1_t[0-12300] or %prog run1_t[0-*]. The \"*\" wild character is "
    +"used to include all files up to the maximum integer starting from START.")
  addParserOptions(parser)
  
  #parse command line options
  (options,args)=parser.parse_args()
    
  #should have at least one argument
  if len(args)!=1 :
    parser.error(" need one and one argument only.")
    
  #create profile files, and save list of files
  failedFiles=make_fileSet(args[0],options,makeFileFunction=makeVTKFile)
  
  if __name__=="__main__":#keeps from redundantly reporting errors
    #report profiles failed files
    for file in failedFiles:
      print file
if __name__ == "__main__":
  main()