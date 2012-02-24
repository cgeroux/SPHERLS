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
  parser=op.OptionParser(usage="Usage: %prog [options] INPUTFILE",version="%prog 1.0"
    ,description="Creates a list of times and periods. The expected input file is one created by "\
    +"average_PKE.py \"average_PKE.txt\" which it reads in and calculated periods based on the"\
    +"the time between every other peak, as there are two peaks in the kinetic energy every "\
    +"period.")
  
  # parser.add_option("-k","--keep",action="store_true",dest="keep"
    # ,help="Keeps distributed binary files [default].",default=True)
  # parser.add_option("-r","--remove",action="store_false",dest="keep"
    # ,help="Removes distributed binary files")
  # parser.add_option("-m","--remake",action="store_true",dest="remakeBins"
    # ,help="Will remake combined binary files from distributed binary files even if combined "
    # +"binary file already exists. [not default].",default=False)
    
  #parse command line options
  (options,args)=parser.parse_args()
    
  #should have at least one argument
  if len(args)!=1 :
    parser.error("Need an input file.")
  calculatePeriods(args[0],"periods.txt")
def calculatePeriods(inputFile,outputFile):
  f=open(inputFile)
  h=open(outputFile,"w")
  timeOld=None
  
  #through out first line of input file since it is a header
  line=f.readline()
  h.write("time[s] period[s]\n")
  numPeaks=0
  for line in f:
    columns=line.split()
    if columns[3]!='-':#there is a peak kinetic energy
      
      if timeOld==None:#grab first time
        timeOld=float(columns[1])
      
      if numPeaks==2:#already had one peak, numPeaks since last write is 2
        period=float(columns[1])-timeOld
        lineOut=str(columns[1])+" "+str(period)+"\n"
        h.write(lineOut)
        numPeaks=0
        timeOld=float(columns[1])
      numPeaks+=1
      
if __name__ == "__main__":
  main()
  