#!/usr/bin/env python
import os
import datafile
import optparse as op
import make_profiles
import glob
import math
import numpy as np
import sys
import paths
import disect_filename

def parseOptions():
  #setup command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END]"
    ,version="%prog 1.0",description="Generates a fourier transform from radial profile files "\
    +"and puts the data in an output file.")
  parser.add_option("-o","--outputFile",dest="outputFile",default="out"
    ,help="Specifies the OUTPUTFILE. [default: %default]"
    ,metavar="OUTPUTFILE", type="string")
  parser.add_option("-z","--zone",type="int",dest="zone",default=0,help="Sets the radial zone "\
    +"from which to extract the variable to use in computing the fourier transform"\
    +" [default: %default]")
  parser.add_option("-c","--column",type="int",dest="column",default=21,help="Sets the column "\
    +"in the radial profile file from which the fourier transform in computed. The default value "\
    +"is %default, the column for the radial grid velocity.")
  
  #add make_profile options to parser
  make_profiles.addParserOptions(parser)
  
  #parse command line options
  return parser.parse_args()
def main():
  #parse command line options
  (options,args)=parseOptions()
  
  #get base file name
  [start,end,baseFileName]=disect_filename.disectFileName(args[0])
  
  #make sure that all the combined binary files have profiles made
  failedFiles=make_profiles.make_fileSet(args[0],options)
  
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
    nNumZones=len(fileData.fColumnValues)
    radialGridVelocity.append(fileData.fColumnValues[nNumZones-options.zone-1][options.column-1])
    fileHeader=fileData.sHeader.split()
    times.append(float(fileHeader[1]))
    nCount=nCount+1
    
  #print data to file
  print "writing to file ..."
  outputFileName=options.outputFile+".txt"
  f=open(outputFileName,'w')
  for i in range(0,len(radialGridVelocity)-1):
    line=str(times[i])+" "+str(radialGridVelocity[i])+"\n"
    f.write(line)
    
  f.close()
  success=os.system(paths.SPHERLSanalPath+' -tl '+outputFileName)
  if success!=0:
    print "error computing fourier transform of \""+outputFileName+"\""
  
  if __name__=="__main__":#keeps from redundently reporting errors
    #report failed files
    for file in failedFiles:
      print file
if __name__ == "__main__":
  main()