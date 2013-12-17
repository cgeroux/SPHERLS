#!/usr/bin/env python

import datafile
import optparse as op
import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np

class Bin:
  def __init__(self,lowerBound,upperBound):
    """Set bounds of the bin and the center of the bin"""
    
    self.mean=None
    self.sigma=None
    self.lowerBound=0.0
    self.upperBound=0.0
    self.center=0.0
    self.dataAddedSinceMeanCal=False
    self.dataAddedSinceSigmaCal=False
    self.points=[]
    self.lowerBound=lowerBound
    self.upperBound=upperBound
    self.center=(lowerBound+upperBound)/2.0
  def addPoint(self,x,y):
    """Add a point to the bin"""
    
    self.points.append([x,y])
    self.dataAddedSinceMeanCal=True
    self.dataAddedSinceSigmaCal=True
  def getMean(self):
    """Returns the mean of the bin, calculating if needed"""
    
    if len(self.points)==0:
      raise Exception("no points in bin ["+str(self.lowerBound)+", "
        +str(self.upperBound)+"]")
    
    #if data added since mean calculated, recalculate
    if self.dataAddedSinceMeanCal:
      
      self.mean=0.0
      for point in self.points:
        self.mean+=point[1]
      self.mean=self.mean/float(len(self.points))
      dataAddedSinceMeanCal=False
    
    #return mean
    return self.mean
  def getSTDD(self):
    """Returns the standard deviation calculating if needed"""
    
    #if data added since standard deviation calculated, recalculate
    if self.dataAddedSinceSigmaCal:
      mean=self.getMean()
      
      sum=0.0
      for point in self.points:
        diff=point[1]-mean
        sum+=diff*diff
        
      self.sigma=math.sqrt(sum/float(len(self.points)))
      self.dataAddedSinceSigmaCal=False
    
    #return sigma
    return self.sigma
class BinnedData:
  def __init__(self):
    self.bins=[]
  def addEvenBins(self,domainMin,domainMax,numBins):
    """Sets the number and edges of the bins"""
    
    #set bin positions as evenly spaced, and initialize averages
    domainSize=float(domainMax)-float(domainMin)
    
    for i in range(numBins):
      lowerBound=float(domainMin)+float(i)*domainSize/float(numBins)
      upperBound=float(domainMin)+float(i+1)*domainSize/float(numBins)
      self.bins.append(Bin(lowerBound,upperBound))
  def binData(self,data):
    """Puts points into bins"""
    
    #bin Data
    print data.shape
    for i in range(len(data)):
      x=data[i][0]
      y=data[i][1]
      
      
      #put point into the correct bin
      binIndex=0
      for j in range(len(self.bins)):
        if x>=self.bins[j].lowerBound and x<self.bins[j].upperBound:
          #print "adding ["+str(x)+", "+str(y)+"] to bin with range ["+str(self.bins[j].lowerBound)+", "+str(self.bins[j].upperBound)+"]"
          self.bins[j].addPoint(x,y)
          break
  def getMean(self):
    """Returns a list of the mean values in each bin"""
    
    mean=[]
    for bin in self.bins:
      mean.append(bin.getMean())
    return mean
  def getBinCenters(self):
    """Returns a list of bin centers"""
    
    centers=[]
    for bin in self.bins:
      centers.append(bin.center)
    return centers
  def getSTDD(self):
    """Returns a list of standard deviations of each bin"""
    
    sigma=[]
    for bin in self.bins:
      sigma.append(bin.getSTDD())
    return sigma
class DataFunction:
  def __init__(self,x,y):
    self.x=x
    self.y=y
    self.maxRange=max(x)
    self.minRange=min(x)
  def getPointByLinearInt(self,x):
    
    #find bracketing values
    index=None
    for i in range(len(self.x)-1):
      if x>=self.x[i] and x<self.x[i+1]:
        index=i
        break
    
    if index!=None:
      #do linear interpolation
      slope=1.0
      try:
        slope=(self.y[index+1]-self.y[index])/(self.x[index+1]-self.x[index])
      except TypeError as e:
        print e
        print "index=",index,"self.y[index+1]=",self.y[index+1],"self.y[index]=",self.y[index],"self.x[index+1]=",self.x[index+1],"self.x[index]=",self.x[index], len(self.x)
        raise e
      return slope*(x-self.x[index])+self.y[index]
    else:
      raise Exception("\""+str(x)+"\" not found in data with range ["+str(max(self.x))+", "+str(min(self.x))+"]")
def parseOptions():
  #note: newlines are not respected in the optparse description string :(, maybe someday will use
  #argparse, which does allow for raw formating (repects indents, newlines etc.)
  
  #setup command line parser
  '''This is out of date, needs to be updated to reflect new additions'''
  parser=op.OptionParser(usage="Usage: %prog [options] INPUTFILE"
    ,version="%prog 1.0",description=r"Calculates the residual between two "
    "curves defined by a set of points in INPUTFILE.")
  parser.add_option("--x1",dest="x1",
    help=("Specify the column in INPUTFILE for the x-coordinate for curve 1 "
      +"starting at 0 for the first column in the INPUTFILE. "
      +"[default: %default]"),
    metavar="COLUMN",type="int",default=0)
  parser.add_option("--y1",dest="y1",
    help=("Specify the COLUMN number in INPUTFILE for the y-coordinate for "
      +"curve 1 starting at 0 for the first column in the INPUTFILE. "
      +"[default: %default]"),
    metavar="COLUMN",type="int",default=1)
  parser.add_option("--x2",dest="x2",
    help=("Specify the COLUMN number in INPUTFILE for the x-coordinate for "
      +"curve 2 starting at 0 for the first column in the INPUTFILE. "
      +"[default: %default]"),
    metavar="COLUMN",type="int",default=2)
  parser.add_option("--y2",dest="y2",
    help=("Specify the COLUMN number in INPUTFILE for the y-coordinate for "
      +"curve 2 starting at 0 for the first column in the INPUTFILE. "
      +"[default: %default]"),
    metavar="COLUMN",type="int",default=3)
    
  #parse command line options
  return parser.parse_args()
def main():
  
  #parse command line options
  (options,args)=parseOptions()
  
  if len(args)!=1:
    raise Exception("need one argument")
  
  #set x/y columns for first curve
  data1XColumn=options.x1
  data1YColumn=options.y1
  
  #x,y columns of second curve
  data2XColumn=options.x2
  data2YColumn=options.y2
  
  #read in data
  print "reading in file \""+args[0]+"\" ..."
  inputFile=datafile.DataFile()
  inputFile.readFile(args[0])
  
  #get first set of data, in our original case it was the model data
  x=inputFile.fColumnValues[:,data1XColumn]
  y=inputFile.fColumnValues[:,data1YColumn]
  data1=np.append(x.reshape(x.shape[0],1),y.reshape(y.shape[0],1),axis=1)
  
  #get second set of data, in our original case this was the observations
  #x=inputFile[:,data2XColumn]
  #y=inputFile[:,data2YColumn]
  #data2=np.append(x.reshape(x.shape[0],1),y.reshape(y.shape[0],1),axis=1)
  
  
  #bin data and get centers, means and standard deviations, of bins
  #this is the model data
  binData=BinnedData()
  binData.addEvenBins(0.0,0.04,5)
  binData.addEvenBins(0.04,0.89,40)
  binData.addEvenBins(0.89,1.0,10)
  binData.binData(data1)
  centers=binData.getBinCenters()
  mean=binData.getMean()
  sigma=binData.getSTDD()
  
  # for bin in binData.bins:
    # print bin.points, bin.lowerBound, bin.upperBound
    # print
  
  # print centers
  # print
  # print mean
  # print
  # print sigma
  
  #calculate model points at observed points
  funcModel=DataFunction(centers,mean)
  modelX=[]
  modelY=[]
  diffY=[]
  i=0
  j=0
  
  for x in inputFile.fColumnValues[:,data2XColumn]:
    
    if x<=funcModel.maxRange and x>funcModel.minRange:
      modelX.append(x)
      modelY.append(funcModel.getPointByLinearInt(x))
      diffY.append(inputFile.fColumnValues[i,data2YColumn]-modelY[j])
      j+=1
    i+=1
  
  #bin differences
  binDataRes=BinnedData()
  binDataRes.addEvenBins(0.0,1.0,20)
  data=np.transpose(np.array([modelX,diffY]))
  binDataRes.binData(data)
  centersRes=binDataRes.getBinCenters()
  meanRes=binDataRes.getMean()
  sigmaRes=binDataRes.getSTDD()
  
  #print len(binDataRes.bins)
  #print centersRes
  #print
  #print meanRes
  #print
  #print sigmaRes
  
  #write residuals of the fit to a file
  fileName=args[0][:len(args[0])-4]+"_residual.txt"
  print "writting residuals to \""+fileName+"\" ..."
  f=open(fileName,'w')
  header="phase residual standard_deviation"
  f.write(header)
  for i in range(len(centersRes)):
    line=str(centersRes[i])+" "+str(meanRes[i])+" "+str(sigmaRes[i])+"\n"
    f.write(line)
  for i in range(len(centersRes)):
    line=str(1+centersRes[i])+" "+str(meanRes[i])+" "+str(sigmaRes[i])+"\n"
    f.write(line)
  f.close()
  
  #plot stuff
  fig,axs=plt.subplots(nrows=2,ncols=1)
  #fig=plt.figure()
  ax=axs[0]
  #ax.plot(modelX,diffY,"o")
  ax.errorbar(centersRes,meanRes,yerr=sigmaRes)
  #ax.plot(centersRes,meanRes)
  ax.plot([0.0,1.0],[0.0,0.0])
  ax.plot([0.0,1.0],[0.05,0.05])
  ax.plot([0.0,1.0],[-0.05,-0.05])
  ax=axs[1]
  ax.plot(modelX,modelY,"o")
  #ax.plot(inputFile.fColumnValues[:,0],inputFile.fColumnValues[:,1],"o")
  ax.plot(inputFile.fColumnValues[:,data2XColumn],inputFile.fColumnValues[:,data2YColumn],"ro")
  #ax.errorbar(centers,mean,yerr=sigma)
  plt.ylim((16.5,14.5))
  plt.xlim((0,1.0))
  #ax.plot(modelX,modelY,"o")

  plt.show()
if __name__ == "__main__":
  main()