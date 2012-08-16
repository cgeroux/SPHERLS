#!/usr/bin/env python

import datafile
import optparse as op
import make_profiles
import glob
import math
import numpy as np
import sys
import disect_filename
import os

def parseOptions():
  #setup command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END]"
    ,version="%prog 1.0",description="Creates a file containing KE, Peak KE, and Average Peak KE.")
  parser.add_option("-k","--keep",action="store_true",dest="keep"
    ,help="Keeps distributed binary files [default].",default=True)
  parser.add_option("-r","--remove",action="store_false",dest="keep"
    ,help="Removes distributed binary files")
  parser.add_option("-m","--remake",action="store_true",dest="remake"
    ,help="Will remake profiles if they already exist. [not default].",default=False)
  parser.add_option("--remake-bins",action="store_true",dest="remakeBins"
    ,help="Will remake binaries if they already exist. [not default].",default=False)
  parser.add_option("--re-sum",action="store_true",dest="resum"
    ,help="Will re-sum all model profiles, usefull when files have problems being made from "
    +"corruption and have to be re-made. Other wise should not be used as it takes more time"
    +" [not default].",default=False)
  parser.add_option("-e",action="store",dest="eosFile"
    ,help="Overrides the equation of state file set in the model.",default=None)
  
  #parse command line options
  return parser.parse_args()
def main():
  
  #parse command line options
  (options,args)=parseOptions()
  
  #get base file name
  [start,end,baseFileName]=disect_filename.disectFileName(args[0])
  
  #make sure that all the combined binary files have profiles made
  failedFiles=make_profiles.make_profiles(options.keep,args[0],options.remake,options.remakeBins
    ,options.eosFile)
  
  #compute the average PKE
  averagePKE(start,end,baseFileName,options)
  
  #report failed files
  for file in failedFiles:
    print file
def averagePKE(start,end,baseFileName,options):

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
    print __name__+":"+averagePKE.__name__+": no files found matching \""+baseFileName+"\" in range given"
    return False
    
  if len(files)<2:
    print __name__+":"+averagePKE.__name__+": need at least 2 files!"
    return False
  
  
  print __name__+":"+averagePKE.__name__+":summing up kinetic energies in profiles ..."
  times=[]
  indexes=[]
  KE=[]
  #read in KE from average file, if file already exists
  lastIndex=-1
  directory=os.path.dirname(baseFileName)
  if directory=="":
    averagePKEFile="averagePKE.txt"
  else:
    averagePKEFile=directory+"/averagePKE.txt"
  if os.path.exists(averagePKEFile):
    if not options.resum:
      print "  \"",averagePKEFile,"\" already exists, not re-summing KE for entries already in file"
      f=open(averagePKEFile,'r')
      f.readline() #skip header
      for line in f:
        lineParts=line.split()
        indexes.append(int(lineParts[0]))
        times.append(float(lineParts[1]))
        KE.append(float(lineParts[2]))
      lastIndex=int(lineParts[0])
      f.close()
    else:
      print __name__+":"+averagePKE.__name__+":  \"",averagePKEFile,"\" already exists, but re-sum set so resumming KE entries in profiles"
  
  nColumKE=61
  for file in files:
    fileIndex=int(file[len(file)-16:len(file)-8])
    if fileIndex>lastIndex:#make sure that file hasn't already been done
      print __name__+":"+averagePKE.__name__+": summing KE of file \""+file+"\" ..."
      fileData=datafile.DataFile()
      fileData.readFile(file)
      fileHeader=fileData.sHeader.split()
      times.append(float(fileHeader[1]))
      
      #get file index
      indexes.append(int(fileIndex))
      
      #sum up all kinetic energies
      dKESum=0.0
      for i in range(len(fileData.fColumnValues)-1):
        dKESum=dKESum+fileData.fColumnValues[i][nColumKE]
      KE.append(dKESum)
  
  bIncreasing=True
  peakKESum=0.0
  numInKESum=0
  numPeaksInAve=6
  halfPeriodAve=0.0
  totalNumPeaks=0
  timeOfLastPeak=0.0
  fracHalfPeriodAcceptPeaks=0.35
  weightOfCurrentPeriod=0.25
  valueOfLastPeak=0
  f=open(averagePKEFile,'w')
  
  #write header to file
  f.write("index(1) t[s](2) KE(3) PKE(4) 3_Period_Ave_PKE(5)\n")
  print __name__+":"+averagePKE.__name__+":computing average Peak KE ..."
  
  #find first two peaks
  bContinue=True
  bContinue1=True
  i=0
  iOfLastPrint=-1
  while bContinue1 and bContinue:
    if KE[i+1]<KE[i] and bIncreasing:#stopped increasing
      bIncreasing=False
      if totalNumPeaks==0:#assume first peak is a peak to include in average
        
        #keep important info for peak finding
        timeOfLastPeak=times[i]
        totalNumPeaks=totalNumPeaks+1
        
        #add peak to summing
        peakKESum=peakKESum+KE[i]
        numInKESum=numInKESum+1
        
        #write first peak
        iOfLastPrint=i
        f.write(str(indexes[i])+" "+str(times[i])+" "+str(KE[i])+" "+str(KE[i])+" -\n")
      elif totalNumPeaks==1:#assume second peak is a peak to include in average
        
        #keep important info for peak finding
        halfPeriodAve=(times[i]-timeOfLastPeak)
        totalNumPeaks=totalNumPeaks+1
        timeOfLastPeak=times[i]
        
        #add peak to summing
        peakKESum=peakKESum+KE[i]
        numInKESum=numInKESum+1
        
        #write second peak
        iOfLastPrint=i
        f.write(str(indexes[i])+" "+str(times[i])+" "+str(KE[i])+" "+str(KE[i])+" -\n")
        bContinue1=False
    elif KE[i+1]>=KE[i]:
      bIncreasing=True
    if iOfLastPrint!=i:#if didn't print out a peak
      #write KE[i]
      iOfLastPrint=i
      f.write(str(indexes[i])+" "+str(times[i])+" "+str(KE[i])+" - -\n")
    
    i=i+1
    if i>=len(times)-1:#check for EOF
      bContinue=False
  
  #repeat search for peak until EOF reached
  while bContinue:
    
    #search expected range for maximum peak
    iOfMaxKEPeak=None
    KEMaxPeak=0.0
    lowerTime=timeOfLastPeak+halfPeriodAve-fracHalfPeriodAcceptPeaks*halfPeriodAve
    upperTime=timeOfLastPeak+halfPeriodAve+(fracHalfPeriodAcceptPeaks)*halfPeriodAve
    bIncreasing=True
    bContinue1=True
    while bContinue1 and bContinue:
      if KE[i+1]<KE[i] and bIncreasing:#stopped increasing
        #if in range of next expected peak
        if times[i]>=lowerTime and times[i]<=upperTime:
          if KE[i]>KEMaxPeak:
            KEMaxPeak=KE[i]
            iOfMaxKEPeak=i
        
        bIncreasing=False
      elif KE[i+1]>=KE[i]:
        bIncreasing=True
      
      if times[i]>upperTime:
        bContinue1=False
      
      #increament
      i=i+1
      if i>=len(times)-1:#if too large, stop everything
        bContinue=False
      
    #if no peaks found in expected range, use next peak found
    if iOfMaxKEPeak==None:
      bIncreasing=True
      bContinue1=True
      while bContinue1 and bContinue:
        
        if KE[i+1]<KE[i] and bIncreasing:#stopped increasing
          KEMaxPeak=KE[i]
          iOfMaxKEPeak=i
          bIncreasing=False
          bContinue1=False
        elif KE[i+1]>=KE[i]:
          bIncreasing=True
        
        i=i+1
        if i+1>=len(times)-1:#check for EOF
          bContinue=False
    
    #write out KE from but not including iOfLastPrint upto but not including iOfMaxKEPeak
    printStop=iOfMaxKEPeak
    if iOfMaxKEPeak==None and not bContinue:
      printStop=len(times)
    for j in range(iOfLastPrint+1,printStop):
      
      #write out KE[j]
      iOfLastPrint=j
      f.write(str(indexes[j])+" "+str(times[j])+" "+str(KE[j])+" - -\n")
    if iOfMaxKEPeak!=None:
      
      #we are garanteed to have a peak, either in the expected range or the next peak after it
      halfPeriod=(times[iOfMaxKEPeak]-timeOfLastPeak)
      halfPeriodAve=weightOfCurrentPeriod*halfPeriod+(1.0-weightOfCurrentPeriod)*halfPeriodAve  #compute a running average
      totalNumPeaks=totalNumPeaks+1
      timeOfLastPeak=times[iOfMaxKEPeak]
      if numInKESum==0:
        peakKESum=0.0
      if numInKESum==numPeaksInAve-1:#finish average
        numInKESum=numInKESum+1
        peakKESum=peakKESum+KEMaxPeak
        
        #write out average
        f.write(str(indexes[iOfMaxKEPeak])+" "+str(times[iOfMaxKEPeak])+" "+str(KE[iOfMaxKEPeak])+" "+str(KE[iOfMaxKEPeak])
          +" "+str(peakKESum/numInKESum)+"\n")
        iOfLastPrint=iOfLastPrint+1
        numInKESum=0
        peakKESum=0.0
      else:#keep averaging
        numInKESum=numInKESum+1
        peakKESum=peakKESum+KEMaxPeak
        #write out peak
        f.write(str(indexes[iOfMaxKEPeak])+" "+str(times[iOfMaxKEPeak])+" "+str(KE[iOfMaxKEPeak])+" "+str(KE[iOfMaxKEPeak])
          +" -\n")
        iOfLastPrint=iOfLastPrint+1
      i=i+1
      if i>=len(times)-2:#check for EOF
        bContinue=False
    
    #write out KE from but not including iOfLastPrint upto last value read in
    if not bContinue:
      printStop=len(times)
      for j in range(iOfLastPrint+1,printStop):
        
        #write out KE[j]
        iOfLastPrint=j
        f.write(str(indexes[j])+" "+str(times[j])+" "+str(KE[j])+" - -\n")
  f.close()
  
  print __name__+":"+averagePKE.__name__+":Peak Kinetic energy data saved in \"",averagePKEFile,"\""
if __name__ == "__main__":
  main()