#!/usr/bin/env python
import os
import datafile
import optparse as op
import make_profiles
import glob
import math
import numpy as np
import sys
import SPHERLSanal_path
import disect_filename

def parseOptions():
  #setup command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END]"
    ,version="%prog 1.0",description="Generates a fourier transform and puts the data in an "\
    +"output file. The radial grid velocity is used, at the zone specified by \"--zone\" option.")
  parser.add_option("-o","--outputFile",dest="outputFile",default="out"
    ,help="Specifies the OUTPUTFILE. [default: %default]"
    ,metavar="OUTPUTFILE", type="string")
  parser.add_option("-z","--zone",type="int",dest="zone",default=0,help="Sets the radial zone from which the radial velocities will be extracted to compute the fourier transform of. Best to use the location of the node of the first overtone if looking for the period of the fundamental. This from starts the center of the star outward, so zero will be at the very bottom of the model. [default: %default]")
  parser.add_option("-k","--keep",action="store_true",dest="keep"
    ,help="Keeps distributed binary files [default].",default=True)
  parser.add_option("-r","--remove",action="store_false",dest="keep"
    ,help="Removes distributed binary files")
  parser.add_option("-m","--remake",action="store_true",dest="remake"
    ,help="Will remake profiles if they already exist. [not default].",default=False)
  #parse command line options
  return parser.parse_args()
def main():
  #parse command line options
  (options,args)=parseOptions()
  
  #get base file name
  [start,end,baseFileName]=disect_filename.disectFileName(args[0])
  
  #make sure that all the combined binary files have profiles made
  failedFiles=make_profiles.make_profiles(options.keep,args[0],options.remake,False)
  
  #get and sort files
  extension="_pro"+".txt"
  filesExistProfiles=glob.glob(baseFileName+"*"+extension)
  files=[]
  for file in filesExistProfiles:
    intOfFile=int(file[len(baseFileName):len(file)-len(extension)])
    if intOfFile>=start and intOfFile<end:
      files.append(file)
  
  files.sort()
  
  if len(files)<1:
    print __name__+":"+main.__name__+": no files found matching \""+baseFileName\
      +"\" in the range "+str(start)+"-"+str(end)
    return False
  
  #make list of radial grid velocity, and time
  radialGridVelocity=[]
  times=[]
  nCount=0
  for file in files:
    print "reading profile ",file," ..."
    fileData=datafile.DataFile()
    fileData.readFile(file)
    radialGridVelocity.append(fileData.fColumnValues[options.zone][20])
    fileHeader=fileData.sHeader.split()
    times.append(float(fileHeader[1]))
    nCount=nCount+1
    
  #print data to file
  print "writing to file ..."
  f=open("tmp_period.txt",'w')
  for i in range(0,len(radialGridVelocity)-1):
    line=str(times[i])+" "+str(radialGridVelocity[i])+"\n"
    f.write(line)
    
  f.close()
  success=os.system(SPHERLSanal_path.SPHERLSanalPath+' -tl tmp_period.txt')
  if success!=0:
    print "error computing fourier transform of tmp_period.txt"
  
  if __name__=="__main__":#keeps from redundently reporting errors
    #report failed files
    for file in failedFiles:
      print file
if __name__ == "__main__":
  main()