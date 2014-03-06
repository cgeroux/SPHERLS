#!/usr/bin/env python

"""This script doesn't produce reliable results but might be used as a starting
point to develop something which could do this in a reasonable and reliable way.
"""

from calculate_residuals_of_light_curve_fit import Bin, BinnedData, DataFunction
import datafile
import optparse as op
import xml.etree.ElementTree as xml
import numpy as np
import scipy

def addParserOptions(parser):
  """Add parser options here
  """
  
  pass
def parseOptions():
  """Creates a command line parser, and parses command line arguments and 
  options
  """
  
  #create parser
  parser=op.OptionParser(usage="Usage: %prog [options] INPUTFILE"
    ,version="%prog 1.0",description=r"Fits input to a series of fits composed "
    +"of data points, the fit quality is evaluated using the algorithms from "
    +"calculate_residuals_of_light_curve_fit")
  
  #parse command line options
  return parser.parse_args()
def parseXML(fileName):
  """Parses the xml configuration file given by fileName
  """
  
  print "parsing \""+str(fileName)+"\" ..."
  
  #get root node
  tree=xml.parse(fileName)
  root=tree.getroot()
  
  #get all fit nodes
  fitElements=root.findall("fit")
  
  #read in fit settings
  fits=[]
  for fitElement in fitElements:
    fit={}
    
    #get input
    inputElements=fitElement.findall("input")
    fit["input"]=inputElements[0].text
    
    #get output
    inputElements=fitElement.findall("output")
    fit["output"]=inputElements[0].text
    
    #get DMInit
    inputElements=fitElement.findall("DMInit")
    fit["DMInit"]=float(inputElements[0].text)
    
    #get PSInit
    inputElements=fitElement.findall("PSInit")
    fit["PSInit"]=float(inputElements[0].text)
    
    #get data
    fitToElements=fitElement.findall("fit-to")
    count=0
    fitTo=[]
    for fitToElement in fitToElements:
      fitTo.append(fitToElement.text)
    fit["fit-to"]=fitTo
    fits.append(fit)
  return fits
def makeFits(fits):
  """Creates fits
  """
  
  for fitSettings in fits:
    print "fitting ",fitSettings['input'],
    
    #get input model to fit
    inputData=datafile.DataFile()
    inputData.readFile(fitSettings['input'])
    
    x1_0=inputData.fColumnValues[:,0]
    y1_0=inputData.fColumnValues[:,1]
    data=np.empty(inputData.fColumnValues.shape)
    DMInit=fitSettings['DMInit']
    PSInit=fitSettings['PSInit']
    
    DM=DMInit#distance modulus, this must be a good estimate, radius of convergence is small
    PS=PSInit
    
    f=open(fitSettings['output'],'w')
    f.write("Fit norm DM PS count DNDDMC DNDPSC\n")
    for fitTo in fitSettings['fit-to']:
      print "to "+fitTo
      #read in data to fit to
      fitData=datafile.DataFile()
      fitData.readFile(fitTo)
      
      PS=PSInit#reset phase shift initialization every fit (they aren't related)
      DM=DMInit#reset phase shift initialization every fit (they aren't related)
      DeltaDM=1.0e-7
      DeltaPS=1.0e-5
      epsDM=1e-4
      epsPS=1e-2
      fracDM=0.8#fraction of correction to add
      fracPS=0.2
      #setup the fit with the data
      bContinue=True
      fitCreated=True
      fitFailedMessage=""
      try:
        fit=Fit(inputData.fColumnValues,fitData.fColumnValues)
      except Exception as e:
        fitFailedMessage=e
        fitCreated=False
        bContinue=False
      
      #this will be a loop until updates to DM and PS get small
      PSCorrection=epsPS+1.0
      DMCorrection=epsDM+1.0
      count=0
      maxIter=20
      while(bContinue):
        
        #Compute first and second derivatives in DM direction
        
        #get three points to evaluate first and second derivatives with
        norm=fit.computeNormOfDifferencePointByPoint(DM,PS)#at centre point
        normmDM=fit.computeNormOfDifferencePointByPoint(DM-DeltaDM,PS)#at DM-Delta
        normpDM=fit.computeNormOfDifferencePointByPoint(DM+DeltaDM,PS)#at DM+Delta

        #compute first derivatives
        DNDDMC=(normpDM-normmDM)/(2.0*DeltaDM)
        DNDDML=(norm-normmDM)/(DeltaDM)
        DNDDMH=(normpDM-norm)/(DeltaDM)
        DNDDM2=(DNDDMH-DNDDML)/(DeltaDM)
        
        #Compute first and second derivatives in PS direction
        
        #get three points to evaluate first and second derivatives with
        normmPS=fit.computeNormOfDifferencePointByPoint(DM,PS-DeltaPS)#at PS-Delta
        normpPS=fit.computeNormOfDifferencePointByPoint(DM,PS+DeltaPS)#at PS+Delta
        
        #compute first derivatives
        DNDPSC=(normpPS-normmPS)/(2.0*DeltaPS)
        DNDPSL=(norm-normmPS)/(DeltaPS)
        DNDPSH=(normpPS-norm)/(DeltaPS)
        DNDPS2=(DNDPSH-DNDPSL)/(DeltaPS)
        
        if(abs(DNDDM2)>epsDM):
          if(DNDDM2!=0.0):
            DMCorrection=-1.0*DNDDMC/(DNDDM2)
          else:
            DMCorrection=1e-10
        else:
          DMCorrection=0.0
          
        #only correct if needed
        if(abs(DNDPSC)>epsPS):
          if(DNDPS2!=0.0):
            PSCorrection=-1.0*DNDPSC/(DNDPS2)
          else:
            PSCorrection=1e-10
        else:
          PSCorrection=0.0
        
        DM=DM+DMCorrection*fracDM
        PS=PS+PSCorrection*fracPS
        if(PS>1.0):
          PS=int(PS)-PS
        elif(PS<-1.0):
          PS=int(PS)+PS
        print "    norm=",norm," count=",count
        print "    DM=",DM," DMCorrection=",DMCorrection," DNDDMC=",DNDDMC," DNDDM2=",DNDDM2," DNDDML=",DNDDML," DNDDMH=",DNDDMH
        print "    PS=",PS," PSCorrection=",PSCorrection," DNDPSC=",DNDPSC," DNDPS2=",DNDPS2," DNDPSL=",DNDPSL," DNDPSH=",DNDPSH
        if(abs(DNDDMC)<epsDM):
          #print "DNDDMC<eps"
          if(abs(DNDPSC)<epsPS):
            #print "DNDPSC<eps"
            bContinue=False
        if(count>=maxIter):
          print "  **FAILED TO CONVERGE**",
          bContinue=False
        count+=1
      
      if fitCreated:
        norm=fit.computeNormOfDifferencePointByPoint(DM,PS)#at centre point
        print "  "+str(fitTo)+" norm=",norm," DM=",DM," PS=",PS," count=",count," DNDDMC=",DNDDMC," DNDPSC=",DNDPSC
        f.write(str(fitTo)+" "+str(norm)+" "+str(DM)+" "+str(PS)+" "+str(count)+" "+str(DNDDMC)+" "+str(DNDPSC)+"\n")
      else:
        print str(fitFailedMessage)
        f.write(str(fitFailedMessage))
      
      #reset initial guesses as, we just failed to converge
      if(count>=maxIter):
        DM=DMInit
        PS=PSInit
    f.close()
def printNorms(fits):
  """Basically a debugging print out
  """
  
  for fitSettings in fits:
    print "fitting ",fitSettings['input']
    
    #get input model to fit
    inputData=datafile.DataFile()
    inputData.readFile(fitSettings['input'])
    
    data=np.empty(inputData.fColumnValues.shape)
    for fitTo in fitSettings['fit-to']:
      
      print "to ",fitTo
      fileName="norms_DM_"+fitTo
      print "writting norms to \""+fileName+"\" ..."
      f=open(fileName,'w')
      #read in data to fit to
      fitData=datafile.DataFile()
      fitData.readFile(fitTo)
      
      #setup the fit with the data
      fit=Fit(inputData.fColumnValues,fitData.fColumnValues)
      
      numPoints=1000
      minDM=14.5
      maxDM=16.5
      minPS=-0.2
      maxPS=0.2
      deltaDM=(maxDM-minDM)/numPoints
      deltaPS=(maxPS-minPS)/numPoints
      PS=0.0
      norm2=fit.computeNormOfDifference(minDM-2*deltaDM,PS)
      norm1=fit.computeNormOfDifference(minDM-deltaDM,PS)
      f.write("DM norm0 norm1 norm2 DIVL DIVH DIV2")
      for i in range(0,numPoints):
        
        #Compute first and second derivatives in DM direction
        
        #get three points to evaluate first and second derivatives with
        DM=minDM+i*deltaDM
        norm0=fit.computeNormOfDifference(DM,PS)#at centre point
        DIVL=(norm1-norm2)/deltaDM
        DIVH=(norm0-norm1)/deltaDM
        DIV2=(DIVH-DIVL)/deltaDM
        f.write(str(DM)+" "+str(norm0)+" "+str(norm1)+" "+str(norm2)+" "+str(DIVL)+" "+str(DIVH)+" "+str(DIV2)+"\n")
        norm2=norm1
        norm1=norm0
      f.close()
class Fit:
  def __init__(self,data1,data2):
    """data 1 should have more points than data 2
    """
    
    self.data1=data1
    self.data1Deriv=np.empty(data1.shape[0])
    self.data1Mod=np.empty(data1.shape)
    self.data2=data2
    self.data2Deriv=np.empty(data2.shape[0])
    
    #Calculate derivatives for data1
    binData1=BinnedData()
    binData1.addEvenBins(0.0,1.0,30)

    binData1.binData(self.data1)
    centers1=binData1.getBinCenters()
    mean1=binData1.getMean()

    func1=DataFunction(centers1,mean1)
    dx=1e-7
    f=open("temp_bin1.txt",'w')
    for i in range(data1.shape[0]):
      
      #calculate derivative
      x=data1[i,0]
      y="-"
      if(x>func1.minRange and x<func1.maxRange):
        #print "[min,max]=["+str(func1.minRange)+","+str(func1.maxRange)+" x="+str(x)+" x+dx="+str(x+dx)+" x-dx="+str(x-dx)
        y=func1.getPointByCubicInt(x)

        if((x+dx)<func1.maxRange):
          dydx=(func1.getPointByCubicInt(x+dx)-func1.getPointByCubicInt(x))/dx
        elif((x-dx)>func1.minRange):
          dydx=(func1.getPointByCubicInt(x-dx)-func1.getPointByCubicInt(x))/dx
        else:
          dydx=0.0
        self.data1Deriv[i]=dydx
      else:
        self.data1Deriv[i]=0.0
      f.write(str(data1[i,0])+" "+str(data1[i,1])+" "+str(y)+" "+str(self.data1Deriv[i])+"\n")
    f.close()
    
    #calculate derivatives for data2
    binData2=BinnedData()
    binData2.addEvenBins(0.0,1.0,30)

    binData2.binData(self.data2)
    centers2=binData2.getBinCenters()
    mean2=binData2.getMean()

    func2=DataFunction(centers2,mean2)
    dx=1e-7
    f=open("temp_bin2.txt",'w')
    for i in range(data2.shape[0]):
      
      #calculate derivative
      x=data2[i,0]
      y="-"
      if(x>func2.minRange and x<func2.maxRange):
        #print "[min,max]=["+str(func2.minRange)+","+str(func2.maxRange)+" x="+str(x)+" x+dx="+str(x+dx)+" x-dx="+str(x-dx)
        y=func2.getPointByCubicInt(x)

        if((x+dx)<func2.maxRange):
          dydx=(func2.getPointByCubicInt(x+dx)-func2.getPointByCubicInt(x))/dx
        elif((x-dx)>func2.minRange):
          dydx=(func2.getPointByCubicInt(x-dx)-func2.getPointByCubicInt(x))/dx
        else:
          dydx=0.0
        self.data2Deriv[i]=dydx
      else:
        self.data2Deriv[i]=0.0
      f.write(str(data2[i,0])+" "+str(data2[i,1])+" "+str(y)+" "+str(self.data2Deriv[i])+"\n")
    f.close()
  def computeNormOfDifference(self,DM,PS,pow=2):
    """Computes the norm of the difference, the order of the norm is given by pow
    """
    
    #apply PS and DM adjustment
    self.data1Mod[:,0]=self.data1[:,0]+PS
    
    #wrap points greater than 1
    for i in range(self.data1Mod.shape[0]):
      if (self.data1Mod[i,0]>1.0):
        self.data1Mod[i,0]=self.data1Mod[i,0]-1.0
        
    self.data1Mod[:,1]=self.data1[:,1]+DM
    
    
    #set binning
    binData=BinnedData()
    binData.addEvenBins(0.0,1.0,50)
    #binData.addEvenBins(0.04,0.89,40)
    #binData.addEvenBins(0.89,1.0,10)

    binData.binData(self.data1Mod)
    centers=binData.getBinCenters()
    mean=binData.getMean()

    #calculate model points at observed points
    funcModel=DataFunction(centers,mean)
    
    #compute the norm
    i=0
    norm=0.0
    for x in self.data2[:,0]:
      if x<=funcModel.maxRange and x>funcModel.minRange:
        diff=self.data2[i,1]-funcModel.getPointByCubicInt(x)
        norm+=diff**pow
      i+=1
    return norm**(1.0/pow)
  def computeNormOfDifferencePointByPoint(self,DM,PS,pow=2):
    """Computes the norm of the difference, the order of the norm is given by pow
    """
    
    #apply PS and DM adjustment
    self.data1Mod[:,0]=self.data1[:,0]+PS
    
    #wrap points greater than 1
    for i in range(self.data1Mod.shape[0]):
      if (self.data1Mod[i,0]>1.0):
        self.data1Mod[i,0]=self.data1Mod[i,0]-1.0
        
    self.data1Mod[:,1]=self.data1[:,1]+DM
    
    #compute the norm
    norm=0.0
    count=0.0

    for i in range(self.data2.shape[0]):
      for j in range(self.data1Mod.shape[0]):
        if(self.data1Deriv[j]<-2.0 and self.data2Deriv[i]<-2.0):#weight steep gradients to match (improve matching rising light)
          weight=self.data2.shape[0]*self.data1Mod.shape[0]*1.0e1
        else:
          weight=1.0
        diff=(self.data2[i,0]-self.data1Mod[j,0])**pow+(self.data2[i,1]-self.data1Mod[j,1])**pow
        norm+=weight*diff
        count+=weight
    return (norm/count)**(1.0/pow)
def main():
  
  (options,args)=parseOptions()
  
  print "this script doesn't produce reliable results, but may be used a "\
    +"starting point "
  
  #parse xml file for fits to make
  fits=parseXML(args[0])
  
  makeFits(fits)
  
  #printNorms(fits)
if __name__ == "__main__":
  main()