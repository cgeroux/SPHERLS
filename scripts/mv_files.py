#!/usr/bin/env python
import os
import getopt
import sys
import glob
import optparse as op
import disect_filename

def main():
  #make command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END] NEWBASEFILENAME"
    ,version="%prog 1.0"
    ,description="Moves a range of files starting with BASEFILENAME with time step indices starting"
    +"at START upto and including END to NEWBASFILENAME+INDEX. Ex. there are a number of files in "
    +"the current directory from run1_t00000000, to run1_t00012300, one could move all these files "
    +"to the new locations run2_t00000000 through to run2_t00012300 with the command %prog "
    +"run1_t[0-12300] run2 or %prog run1_t[0-*] run2. The \"*\" wild character is used to include "
    +"all files upto the maximum integer starting from START. This command also moves files with "
    +"additional suffixes, such as the distributed files, profiles, 2D slices etc., while "
    +"maintianing the correct suffixes and indices.")
  parser.add_option("--cb-only",action="store_true",dest="cbOnly"
    ,help="If set it will move only combined binary files matching BASEFILENAME [default]."
    ,default=False)
  
  #parse command line options
  (options,args)=parser.parse_args()
    
  #should have at least two arguments
  if len(args)!=2 :
    parser.error(" need two and two arguments only.")
    
  #create profile files, and save list of files
  move_files(args[0],args[1],options)
def move_files(fileName,newBaseFileName,options):
  
  #get base file name
  [start,end,baseFileName]=disect_filename.disectFileName(fileName)
  
  #check for distributed binary files in interation range starting with baseFileName
  filesExist=glob.glob(baseFileName+"[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]*")
  files=[]
  for file in filesExist:
    
    intOfFile=-1
    
    #if base part of file name matches the base file name given
    if file[0:len(baseFileName)]==baseFileName:
      testStr=file[len(baseFileName):len(baseFileName)+8]
      endStr=file[len(baseFileName)+8:len(file)]
      if testStr.isdigit():
        intOfFile=int(testStr)
      
      if intOfFile>=start and intOfFile<end:
        if options.cbOnly and endStr=='' or not options.cbOnly:
          files.append(file)
      
  if(len(files)==0):
    print __name__+":"+move_files.__name__+": no files found in range "+baseFileName+"["+parts2[0]+"-"+parts2[2]+"]"
    return False #nothing to do
  
  #put them in order, might be nice to have them in order? Certainly not required
  files.sort()
  for file in files:
    parts=file.partition('_t')
    
    newFileName=newBaseFileName+parts[1]+parts[2]
    print __name__+":"+move_files.__name__+": moving file \""+file+"\" to file "+newFileName
    
    #remove the file
    cmd='mv '+file+" "+newFileName
    os.system(cmd)
    
  return True
if __name__ == "__main__":
  main()
  