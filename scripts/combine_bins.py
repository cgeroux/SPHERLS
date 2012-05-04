#!/usr/bin/env python
import os
import getopt
import sys
import glob
import optparse as op
import paths
import disect_filename

def main():
  #make command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END]",version="%prog 1.0"
    ,description="Combines a range of binary files starting with BASEFILENAME with time step "
    +"indices including START upto but not including END. Ex. there are a number of files in the "
    +"current directory from run1_t00000000, to run1_t00012300, one could combine all these files "
    +"with the command %prog run1_t[0-12300] or %prog run1_t[0-*]. The \"*\" wild character is "
    +"used to include all files upto the maximum integer starting from START.")
  
  parser.add_option("-k","--keep",action="store_true",dest="keep"
    ,help="Keeps distributed binary files [default].",default=True)
  parser.add_option("-r","--remove",action="store_false",dest="keep"
    ,help="Removes distributed binary files")
  parser.add_option("-m","--remake",action="store_true",dest="remakeBins"
    ,help="Will remake combined binary files from distributed binary files even if combined "
    +"binary file already exists. [not default].",default=False)
    
  #parse command line options
  (options,args)=parser.parse_args()
    
  #should have at least one argument
  if len(args)!=1 :
    parser.error(" need one and one argument only.")
    
  #create profile files, and save list of files
  combine_bin_files(options.keep,args[0],options.remakeBins)
def combine_bin_files(keep,fileName,remakeBins):
  
  [start,end,baseFileName]=disect_filename.disectFileName(fileName)
  
  #if only doning one file we can make this faster
  filesExist=[]
  if end-start==1:
    filesExist.append(baseFileName+str(start).zfill(8)+"-0")
  else:
    #check for distributed binary files in interation range starting with baseFileName
    filesExist=glob.glob(baseFileName+"[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]-0")
  
  files=[]
  for file in filesExist:
    intOfFile=int(file[len(baseFileName):len(file)-2])
    if intOfFile>=start and intOfFile<end:
      files.append(file)
  
  if(len(files)==0):
    return [] #nothing to do, return an empty list of failedFiles
  
  files.sort()
  failedFiles=[]
  for i in range(len(files)):
    
    #if combined binary file doesn't already exist with this file name
    if not (os.path.exists(files[i][:len(files[i])-2])) or remakeBins:
      
      #make combined binary files
      if not keep:
        print __name__+":"+combine_bin_files.__name__+": combining \""+files[i][:len(files[i])-2]\
          +"\" and removing distributed binaries ..."
      else:
        print __name__+":"+combine_bin_files.__name__+": combining \""+files[i][:len(files[i])-2]+"\" ..."
      
      success=os.system(paths.SPHERLSanalPath+' -c dbcb '+files[i][:len(files[i])-2])
      if success==0:
        if not keep:
          #remove distributed binary files
          os.system('rm -f '+files[i][:len(files[i])-2]+'-*')
      else :
        #say there was an error and quit
        failedFiles.append(__name__+":"+combine_bin_files.__name__+": error combining binary file "
          +files[i][:len(files[i])-2])
    else:
      if not keep:
        print __name__+":"+combine_bin_files.__name__+": combined binary \""\
          +files[i][:len(files[i])-2]+"\" already exists, removing distributed binaries ..."
        #remove distributed binary files
        os.system('rm -f '+files[i][:len(files[i])-2]+'-*')
      else:
        print __name__+":"+combine_bin_files.__name__+": combined binary \""\
          +files[i][:len(files[i])-2]+"\" already exists"
  
  if __name__=="__main__":#keeps from redundently reporting errors
    #report problem files
    for failedFile in failedFiles:
      print failedFile
  
  return failedFiles
if __name__ == "__main__":
  main()
  