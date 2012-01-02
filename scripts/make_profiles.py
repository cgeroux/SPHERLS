#!/usr/bin/env python
import os
import getopt
import sys
import glob
import optparse as op
import combine_bins
import SPHERLSanal_path
import disect_filename

def main():
  #make command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END]",version="%prog 1.0"
    ,description="Creates radial profiles from a range of combined binary files starting with "
    +"BASEFILENAME with time step indices including START upto but not including END. Ex. there "
    +"are a number of files in the current directory from run1_t00000000, to run1_t00012300, one "
    +"could combine all these files with the command %prog run1_t[0-12300] or %prog run1_t[0-*]. "
    +"The \"*\" wild character is used to include all files upto the maximum integer starting from "
    +"START.")
  parser.add_option("-k","--keep",action="store_true",dest="keep"
    ,help="Keeps distributed binary files [default].",default=True)
  parser.add_option("-r","--remove",action="store_false",dest="keep"
    ,help="Removes distributed binary files [not default]")
  parser.add_option("-m","--remake",action="store_true",dest="remake"
    ,help="Will remake profiles even if they already exist. [not default].",default=False)
  parser.add_option("--remake-bins",action="store_true",dest="remakeBins"
    ,help="Will remake binaries even if they already exist. [not default].",default=False)
  
  #parse command line options
  (options,args)=parser.parse_args()
    
  #should have at least one argument
  if len(args)!=1 :
    parser.error(" need one and one argument only.")
    
  #create profile files, and save list of files
  make_profiles(options.keep,args[0],options.remake,options.remakeBins)
def make_profiles(keep,fileName,remake,remakeBins):
  
  #get base file name
  [start,end,baseFileName]=disect_filename.disectFileName(fileName)
  
  #make sure that all the distributed binary files in this range have been converted to combined binary files
  failedFiles=combine_bins.combine_bin_files(keep,fileName,remakeBins)
  
  #check for combined binary files in interation range starting with baseFileName
  filesExistCombBin=glob.glob(baseFileName+"[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]")
  files=[]
  for file in filesExistCombBin:
    intOfFile=int(file[len(baseFileName):len(file)])
    if intOfFile>=start and intOfFile<end:
      files.append(file)
  
  files.sort()
  for file in files:
      
      #if profile not already made for this combined binary file
      if not (os.path.exists(file+"_pro.txt")) or remake:
        
        #make profile
        print __name__+":"+make_profiles.__name__+": creating profile from \""+file+"\" ..."
        success=os.system(SPHERLSanal_path.SPHERLSanalPath+' -a cb '+file)
        if success==0:
          pass
        else :
          #say there was an error and quit
          failedFiles.append(__name__+":"+make_profiles.__name__+": error making profile "+file+"_pro.txt")
      else:
        print __name__+":"+make_profiles.__name__+": profile \""+file\
          +"_pro.txt\" already exists, not remaking"
  
  if __name__=="__main__":#keeps from redundently reporting errors
    #report profiles failed files
    for file in failedFiles:
      print file
  
  return failedFiles
if __name__ == "__main__":
  main()
  