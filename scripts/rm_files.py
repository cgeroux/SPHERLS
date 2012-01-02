#!/usr/bin/env python
import os
import getopt
import sys
import glob
import optparse as op

def main():
  #make command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END]",version="%prog 1.0"
    ,description="Removes a range of file matching BASEFILENAME_tXXXXXXXX with XXXXXXXX being a "
    +"zero padded number starting with START and continuing to but not including END. Ex. there are"
    +" a number of files in the current directory from run1_t00000000, to run1_t00012300, one could"
    +" remove all these files with the command %prog run1_t[0-12300] or %prog run1_t[0-*]. The "
    +"\"*\" wild character is used to include all files upto the maximum integer starting from "
    +"START.")
  
  #parse command line options
  (options,args)=parser.parse_args()
  
  #should have at least one argument
  if len(args)!=1 :
    parser.error(" need one and one argument only.")
    
  #create profile files, and save list of files
  remove_files(args[0],options)
    
def remove_files(fileName,options):
  
  #get base file name
  parts0=fileName.partition('[')
  baseFileName=parts0[0]
  
  #get interation range
  parts1=parts0[2].partition(']')
  parts2=parts1[0].partition('-')
  start=int(parts2[0])
  if parts2[2]=="*":
    end=sys.maxint
  else:
    end=int(parts2[2])
  
  #check for distributed binary files in interation range starting with baseFileName
  filesExist=glob.glob(baseFileName+"[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]*")
  
  #may want to see if this also gets combined binary files, and if not may want to add it
  files=[]
  for file in filesExist:
    
    intOfFile=-1
    
    
    intOfFile=-1
    
    #if base part of file name matches the base file name given
    if file[0:len(baseFileName)]==baseFileName:
      testStr=file[len(baseFileName):len(baseFileName)+8]
      if testStr.isdigit():
        intOfFile=int(testStr)
      
      if intOfFile>=start and intOfFile<end:
        files.append(file)
  
  if(len(files)==0):
    print __name__+":"+remove_files.__name__+": no files found in range "+baseFileName+"["+parts2[0]+"-"+parts2[2]+"]"
    return False #nothing to do
  
  for file in files:
    
    print __name__+":"+remove_files.__name__+": removing file \""+file+"\""
    
    #remove the file
    os.system('rm -f '+file)
    
  return True
if __name__ == "__main__":
  main()
  