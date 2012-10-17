#!/usr/bin/env python
from scipy import interpolate
import xml.etree.ElementTree as xml
import optparse as op
import numpy as np
import math
import struct
import paths
import os
a=7.56591e-15 #radiation constant [erg cm^{-3} K^{-4}]
def parseOptions():
  #note: newlines are not respected in the optparse description string :(, maybe someday will use
  #argparse, which does allow for raw formating (repects indents, newlines etc.)
  
  #setup command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] XMLFILE\n"
    +"       %prog [options] EOSFILE1 EOSFILE2"
    ,version="%prog 1.0",description=r"Interpolates in a set of equation of state files. XMLFILE "
    +"gives all the settings needed to create multiple eos files used in SPHERLS. For an example of"
    +" a configuration file see eos_interp_reference.xml under docs/templateXML. In the second"
    +" usage case it can be used to compare two equation of state files.")
  parser.add_option("--rhoIndex1",dest="rho1",type="int",default=0
    ,help="Sets the index for density in table 1 [default: %default]")
  parser.add_option("--rhoIndex2",dest="rho2",type="int",default=0
    ,help="Sets the index for density in table 2 [default: %default]")
  parser.add_option("--numRho",dest="numRho",type="int",default=1
    ,help="Sets the number of densities to plot starting at the specified indices [default: %default]")
  parser.add_option("-o",dest="output",type="string",default=None
    ,help="Sets the output file name, so that the comparison plot is save to a file [default: %default]")
  #parse command line options
  (options,args)=parser.parse_args()
  if len(args)==0:
    raise Exception("need to specify configuration file")
  return (options,args)
def checkLineForMinusSplitError(line):
  """Checks for spaces being replaced with "-" as the seperator and properly splits the line if that
  is the case"""
  
  splitLine=line.split()
  
  #check that we can convert elements to a floats
  if len(splitLine)>0:
    try:
      for i in range(len(splitLine)):
        temp=float(splitLine[i])
    except ValueError:
      
      '''if not able to convert to a float it likely means that we have an issue with '-' take up
      ' ' so numbers aren't seperated by ' ', such a pain'''
      
      splitLineMinus=line.split("-")
      count=0
      splitLine=[]#empty list so we can refill it properly
      for splitTemp in splitLineMinus:
        splitLineSpace=splitTemp.split()
        if len(splitLineSpace)>0:
          if count>0:
            splitLine.append(str(float(splitLineSpace[0])*-1.0))
          else:
            splitLine.append(splitLineSpace[0])
          for i in range(1,len(splitLineSpace)):
            splitLine.append(splitLineSpace[i])
        count+=1
  return splitLine
class eosTable:
  """Holds equation of state data."""
  
  def load(self):
    """Reads in an OPAL equation of state file.
    
    It puts the resulting file info into:
    self.X: the hydrogen mass fraction
    self.Z: the metal mass fraction
    self.logD: numpy array of log density grid points [g/cm^3]
    self.logT: numpy array of log tempeature grid points [K]
    self.logE: numpy array of log energy [ergs/g]
    self.logP: numpy array of log pressure [dynes/cm^2]
    
    self.logD, self.logT, self.logE, and self.logP are all the same size numpy arrays, empty 
    emelents have logE and logP as nans.
    """
    
    self.status="loaded"
    print "loading \""+self.sFileName+"\" ..."
    
    f=open(self.sFileName,'r')
    
    #get composition
    line=f.readline()
    lineSplit=line.split()
    self.X=float(lineSplit[1])
    self.Z=float(lineSplit[3])
    Moles=float(lineSplit[5])
    MassPerMole=float(lineSplit[7])
    
    #skip 2 lines
    line=f.readline()
    line=f.readline()
    
    #make lists to hold data
    logD=[]
    logT=[]
    logE=[]
    logP=[]
    
    #loop over densities
    while(True):
      
      #read next line
      line=f.readline()
      lineSplit=line.split()
      
      #get number of temperatures to follow
      numTemps=int(lineSplit[1])
      
      #check if we are done
      if numTemps==0:#this will be zero if we are at the end of the file
        break
      
      #get density
      density=float(lineSplit[5])
      
      #get some stuff I don't understand, but allows us to convert to "normal" units?
      [tmassz,amoles,eground,fracz,fracz]=self.__gmass()
      
      #skip 2 lines
      line=f.readline()
      line=f.readline()
      
      #loop over all temperatures
      logDRow=[]
      logTRow=[]
      logERow=[]
      logPRow=[]
      for i in range(0,numTemps):
        
        #read in a temperature line
        line=f.readline()
        lineSplit=line.split()
        temperature=float(lineSplit[0])
        pressure=float(lineSplit[3])
        energy=float(lineSplit[4])
        #po=temperature*density
        #pressure=pressure*po#remove temperture*density scaling
        pressure=pressure*1.0e12#convert from MB to dynes/cm^2
        energy=(energy*MassPerMole-eground)/(MassPerMole)*1.0e12#set e=0 at T=0
        
        #save values
        logDRow.append(math.log10(density))
        logTRow.append(math.log10(temperature*1.0e6))#convert to K from MK
        logERow.append(math.log10(energy))
        logPRow.append(math.log10(pressure))
      
      #add to 2D lists
      logD.append(logDRow)
      logT.append(logTRow)
      logE.append(logERow)
      logP.append(logPRow)
      
      #skip three lines before the next density
      line=f.readline()
      line=f.readline()
      line=f.readline()
    
    #convert lists to numpy arrays and store them
    #find largest row dim
    #all rows will be the same for T,R, and K
    nMaxRow=0
    for i in range(len(logT)):
      if len(logT[i])>nMaxRow:
        nMaxRow=len(logT[i])
    
    #initialize to all nans
    self.logT=np.empty((len(logT),nMaxRow))
    self.logT[:]=np.NAN
    self.logD=np.empty((len(logT),nMaxRow))
    self.logD[:]=np.NAN
    self.logP=np.empty((len(logT),nMaxRow))
    self.logP[:]=np.NAN
    self.logE=np.empty((len(logT),nMaxRow))
    self.logE[:]=np.NAN
    
    #fill in elements that are not nans
    for i in range(len(logT)):
      for j in range(len(logT[i])):
        self.logD[i][j]=logD[i][j]
        self.logT[i][j]=logT[i][j]
        self.logE[i][j]=logE[i][j]
        self.logP[i][j]=logP[i][j]
    
    #fill in logD and logT nans to make a rectangular grid
    self.__fillInDepNans()
    
    '''if temperature decreases instead of increases flip arrays. This is the way the interpolation 
    routines like things.
    '''
    if self.logT[0][0]>self.logT[0][1]:
      self.logD=np.fliplr(self.logD)
      self.logT=np.fliplr(self.logT)
      self.logE=np.fliplr(self.logE)
      self.logP=np.fliplr(self.logP)
  def write(self,*args):
    """Generic write function that calls either writeToScreen, or writeToFiel depending on if a 
    file name is specified or not."""
    
    if len(args)==0:
      self.__writeToScreen()
    elif len(args)==1:
      self.__writeToFile(args[0])
  def plotLogE(self,otherTables=None,logDIndexList=None,wireFrame=True):
    """Plots LogE
    
    Keywords:
    otherTables: a list of other eosTables to include in the plot
    logDIndexList: a list of integers corresponding to which densities to plot the tables at
    wireFrame: if set to true (the default) and logDIndexList is set to None it will plot a 3D
      wireframe of logE.
    """
    
    if logDIndexList!=None:
      wireFrame=False
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib
    import matplotlib.pyplot as plt
    fig=plt.figure()
    if wireFrame:
      ax=fig.add_subplot(111,projection='3d')
      ax.set_xlabel("log10(rho) [g/cm^3]")
      ax.set_ylabel("log10(T) [K]")
      ax.set_zlabel("log10(E) [ergs/g]")
      ax.plot_wireframe(self.logD,self.logT,self.logE)
      if otherTables:
        for otherTable in otherTables:
          ax.plot_wireframe(otherTable.logD,otherTable.logT,otherTable.logE,color="green")
    else:
      ax=fig.add_subplot(111)
      l=logDIndexList[0]
      h=logDIndexList[0]+1
      print self.logD[l:h,0]
      ax.plot(self.logT[l:h,:][0],self.logE[l:h,:][0], "bo-")
      counter=1
      if otherTables:
        for otherTable in otherTables:
          l=logDIndexList[counter]
          h=logDIndexList[counter]+1
          print otherTable.logD[l:h,0]
          ax.plot(otherTable.logT[l:h,:][0],otherTable.logE[l:h,:][0], "go-")
          counter+=1
    plt.show()
  def plotLogP(self,otherTables=None,logDIndexList=None,wireFrame=True):
    """Plots LogP
    
    Keywords:
    otherTables: a list of other eosTables to include in the plot
    logDIndexList: a list of integers corresponding to which densities to plot the tables at
    wireFrame: if set to true (the default) and logDIndexList is set to None it will plot a 3D
      wireframe of logP.
    """
    
    if logDIndexList!=None:
      wireFrame=False
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib
    import matplotlib.pyplot as plt
    fig=plt.figure()
    if wireFrame:
      ax=fig.add_subplot(111,projection='3d')
      ax.set_xlabel("log10(rho) [g/cm^3]")
      ax.set_ylabel("log10(T) [K]")
      ax.set_zlabel("log10(P) [dynes/cm^2]")
      ax.plot_wireframe(self.logD,self.logT,self.logP)
      if otherTables:
        for otherTable in otherTables:
          ax.plot_wireframe(otherTable.logD,otherTable.logT,otherTable.logP,color="green")
    else:
      ax=fig.add_subplot(111)
      l=logDIndexList[0]
      h=logDIndexList[0]+1
      print self.logD[l:h,0]
      ax.plot(self.logT[l:h,:][0],self.logP[l:h,:][0], "bo-")
      counter=1
      if otherTables:
        for otherTable in otherTables:
          l=logDIndexList[counter]
          h=logDIndexList[counter]+1
          print otherTable.logD[l:h,0]
          ax.plot(otherTable.logT[l:h,:][0],otherTable.logP[l:h,:][0], "go-")
          counter+=1
    plt.show()
  def interpolate(self,gridConfig,setExtrapolatedToNan=True):
    """Interpolate from self's table to the griding specified by:
    
    logDMin: first (smallest) logD value of grid
    logDDel: spacing in logD
    numLogD: number of logD grid points
    logTMin: first (smallest) logT value of grid
    logTDel: spacing in logT
    numLogT: number of logT grid points
    """
    logDMin=gridConfig[0]
    logDDel=gridConfig[1]
    numLogD=gridConfig[2]
    logTMin=gridConfig[3]
    logTDel=gridConfig[4]
    numLogT=gridConfig[5]
    
    print "  interpolating eos to new grid ..."
    
    #create a new interpolated table
    interpolatedTable=self.__class__()
    interpolatedTable.X=self.X
    interpolatedTable.Z=self.Z
    
    #save some of the interpolation information
    interpolatedTable.logDMin=logDMin
    interpolatedTable.logDDel=logDDel
    interpolatedTable.logTMin=logTMin
    interpolatedTable.logTDel=logTDel
    
    #make space
    interpolatedTable.logD=np.empty((numLogD,numLogT))
    interpolatedTable.logT=np.empty((numLogD,numLogT))
    interpolatedTable.logE=np.empty((numLogD,numLogT))
    interpolatedTable.logP=np.empty((numLogD,numLogT))
    
    #set logD and logT grid points
    for i in range(numLogD):
      for j in range(numLogT):
        interpolatedTable.logD[i][j]=logDMin+float(i)*logDDel
        interpolatedTable.logT[i][j]=logTMin+float(j)*logTDel
    
    #remove nans from grid
    #this part is slow
    print "    filling in nans by extrapolating (this takes a while) ..."
    [logE,logP]=self.__fillDepNans()
    
    #interpolate to grid
    print "    creating interpolator ..."
    splineInterpLogP=interpolate.RectBivariateSpline(self.logD[:,0],self.logT[0,:],logP)
    splineInterpLogE=interpolate.RectBivariateSpline(self.logD[:,0],self.logT[0,:],logE)
    
    print "    interpolating to new grid ..."
    for i in range(numLogD):
      for j in range(numLogT):
        logPGas=splineInterpLogP(interpolatedTable.logD[i][j],interpolatedTable.logT[i][j])
        logEGas=splineInterpLogE(interpolatedTable.logD[i][j],interpolatedTable.logT[i][j])
        
        #print 10**interpolatedTable.logP[i][j]/10**interpolatedTable.logD[i][j],interpolatedTable.logE[i][j]
        
        #add radiation pressure
        PGas=10**logPGas
        EGas=10**logEGas
        T=10**interpolatedTable.logT[i][j]
        rho=10**interpolatedTable.logD[i][j]
        
        PRad=a*T**4/3.0
        ERad=3.0*PRad/rho
        interpolatedTable.logP[i][j]=math.log10(PGas+PRad)
        interpolatedTable.logE[i][j]=math.log10(EGas+ERad)
    
    #set values outside grid, and near nans on original grid equal to nan in interpolated grid
    if setExtrapolatedToNan:
      
      print "    setting nans for extrapolated values (this takes a while)..."
      
      #get range in logT, and logR
      logTRange=[min(self.logT[0,:]),max(self.logT[0,:])]
      logDRange=[min(self.logD[:,0]),max(self.logD[:,0])]
      np.set_printoptions(precision=8)
      for i in range(numLogD):
        for j in range(numLogT):
          
          #if logR and logT outside original table it is a nan
          bOutSide=False
          if interpolatedTable.logT[i][j]<logTRange[0] or logTRange[1]<interpolatedTable.logT[i][j]:
            bOutSide=True
          if interpolatedTable.logD[i][j]<logDRange[0] or logDRange[1]<interpolatedTable.logD[i][j]:
            bOutSide=True
          
          #if not outside see if it is near a nan in the oringal table
          if not bOutSide:
            
            #find location on original table in logR
            nIndexLogD=0
            for k in range(len(self.logD[:,0])-1):
              if self.logD[k][0]<=interpolatedTable.logD[i][j]\
                and interpolatedTable.logD[i][j]<self.logD[k+1][0]:
                nIndexLogD=k
                break
            
            #find location on original table in logT
            nIndexLogT=0
            for k in range(len(self.logT[0,:])-1):
              if self.logT[0][k]<=interpolatedTable.logT[i][j]\
                and interpolatedTable.logT[i][j]<self.logT[0][k+1]:
                nIndexLogT=k
                break
            
            #check if any of the nearest neighbors are nans if so current point should be a nan too
            bNearNan=False
            if np.isnan(self.logE[nIndexLogD][nIndexLogT]):
              bNearNan=True
            if nIndexLogD<self.logE.shape[1]-2:
              if np.isnan(self.logE[nIndexLogD+1][nIndexLogT]):
                bNearNan=True
            if nIndexLogT<self.logE.shape[0]-2:
              if np.isnan(self.logE[nIndexLogD][nIndexLogT+1]):
                bNearNan=True
            if nIndexLogT<self.logE.shape[0]-2 and nIndexLogD<self.logE.shape[1]-2:
              if np.isnan(self.logE[nIndexLogD+1][nIndexLogT+1]):
                bNearNan=True
            
            if bNearNan:
              interpolatedTable.logE[i][j]=np.nan
              interpolatedTable.logP[i][j]=np.nan
          else:
            interpolatedTable.logE[i][j]=np.nan
            interpolatedTable.logP[i][j]=np.nan
      
    return interpolatedTable
  def __init__(self,sFileName=None):
    """Returns a new instance of eosTable.
    
    If sFileName is set it will use that to set the filename to load the data from."""
    
    
    self.status="Haven't loaded"
    
    if sFileName!=None:
      self.sFileName=sFileName
  def __gmass(self):
    """This code has been reproduced in python from the fortran code from OPAL tables in the
    ZFS_interp_EOS5.f file."""
    
    #I appologize for the lack of explanation/comments below but it was just more or less a copy
    #of a funciton in ZFS_intper_EOS5.f
    
    x=self.X
    z=self.Z
    xc=0.247137766
    xn=0.0620782
    xo=0.52837118
    xne=0.1624188
    amas=[0.00054858,20.179,15.9994,14.0067,12.011,4.0026,1.0079]#7 numbers, atomic masses?
    anum=[10.0,8.0,7.0,6.0,2.0,1.0]#6 numbers, atomic numbers?
    eion=[-3394.873554,-1974.86545,-1433.92718,-993.326315,-76.1959403,-15.28698437]#ionization energies?
    fracz=z/(xc*amas[4]+xn*amas[3]+xo*amas[2]+xne*amas[1])
    xc2=fracz*xc
    xn2=fracz*xn
    xo2=fracz*xo
    xne2=fracz*xne
    xh=x/amas[6]
    xhe=(1.0-x-z)/amas[5]
    xtot=xh+xhe+xc2+xn2+xo2+xne2
    frac=[]
    frac.append(xne2/xtot)
    frac.append(xo2/xtot)
    frac.append(xn2/xtot)
    frac.append(xc2/xtot)
    frac.append(xhe/xtot)
    frac.append(xh/xtot)
    eground=0.0
    amoles=0.0
    for i in range(len(eion)):
      eground+=eion[i]*frac[i]
      #print eground,eion[i],frac[i]
      amoles+=(1.0+anum[i])*frac[i]
    anume=amoles-1.0
    gmass=0.0
    for i in range(1,6):
      gmass+=amas[i]*frac[i-1]
    return [gmass,amoles,eground,fracz,frac]
  def __writeToScreen(self):
    """writes eos table to screen. mostly used for debugging"""
    
    print "status: ",self.status
    print "X=",self.X
    print "Y=",self.Z
    print "density temperature energy pressure"
    for i in range(len(self.density)):
      line="{0:.15e} {1:.15e} {2:.15e} {3:.15e}".format(self.density[i],self.temperature[i]\
        ,self.energy[i],self.pressure[i])
      print line
  def __writeToFile(self,sFileName):
    """writes eos table to a file"""
    
    out=open(sFileName)
    out.write("X="+str(self.X))
    out.write("Y="+str(self.Z))
    out.write("density temperature energy pressure")
    for i in range(len(self.density)):
      line="{0:.15e} {1:.15e} {2:.15e} {3:.15e}".format(self.density[i],self.temeprature[i]\
        ,self.energy[i],self.pressure[i])
      out.write(line)
  def __fillInDepNans(self):
    """Fills in logD and logT values to make a rectangular grid"""
    
    
    '''fill in logT differently since it varies by columns instead of by rows as density does'''
    for k in range(self.logT.shape[1]):
      lastLogT=float('nan')
      for j in range(self.logT.shape[0]):
        
        #for logT
        if not np.isnan(self.logT[j][k]):
          lastLogT=self.logT[j][k]
        elif not np.isnan(lastLogT):
          self.logT[j][k]=lastLogT
          
    for j in range(len(self.logD)):
      lastLogD=float('nan')
      
      for k in range(len(self.logD[j])):
        
        #fill in logD nans
        if not np.isnan(self.logD[j][k]):
          lastLogD=self.logD[j][k]
        elif not np.isnan(lastLogD):
          self.logD[j][k]=lastLogD
  def __fillDepNans(self):
    """Extrapolates boundary grid points to Nan elements using an averaged weighted by a high power 
    of distance away from the boundary point. This function returns [logE,logP] free of 
    nans."""
    
    #fill in nans with extrapolation only for purposes of the interpolation
    logE=np.copy(self.logE)
    logP=np.copy(self.logP)
    
    '''find points along nan boundaries'''
    boundaryLogD=[]
    boundaryLogT=[]
    boundaryLogE=[]
    boundaryLogP=[]
    for j in range(len(self.logD)):
      for k in range(len(self.logD[j])):
        
        #only add non-nan points to boundary list
        if not np.isnan(self.logE[j][k]):
          
          '''if there is a nan either behind, or infront of the current point in either j or k 
          space it is on a boundary'''
          bOnBoundary=False
          if j>0:
            if np.isnan(self.logE[j-1][k]):
              bOnBoundary=True
          if j<len(self.logE)-1:
            if np.isnan(self.logE[j+1][k]):
              bOnBoundary=True
          if k>0:
            if np.isnan(self.logE[j][k-1]):
              bOnBoundary=True
          if k<len(self.logE[j])-1:
            if np.isnan(self.logE[j][k+1]):
              bOnBoundary=True
          if bOnBoundary:
            boundaryLogD.append(self.logD[j][k])
            boundaryLogT.append(self.logT[j][k])
            boundaryLogE.append(self.logE[j][k])
            boundaryLogP.append(self.logP[j][k])
    
    '''extrapolate E,P outside boundary using a weighted average of boundary points. The weight 
    is simply the inverse of the distnace to a high power from a given boundary point'''
    for j in range(len(self.logD)):
      lastLogE=float('nan')
      lastLogP=float('nan')
      
      for k in range(len(self.logE[j])):
        
        #for logE
        if np.isnan(self.logE[j][k]):
          sumE=0.0
          sumP=0.0
          sumd=0.0
          for n in range(len(boundaryLogD)):
            distanceLogD=(self.logD[j][k]-boundaryLogD[n])/2.0
            distanceLogT=(self.logT[j][k]-boundaryLogT[n])
            distance=distanceLogD*distanceLogD+distanceLogT*distanceLogT
            distance=distance*distance
            distance=distance*distance
            distance=distance*distance
            oneDdistance=1.0/distance
            sumE+=oneDdistance*boundaryLogE[n]
            sumP+=oneDdistance*boundaryLogP[n]
            sumd+=oneDdistance
          #print sumE,sumd
          logE[j][k]=sumE/sumd
          logP[j][k]=sumP/sumd
    
    ###########
    #DEBUG, to take a look at tables
    #from mpl_toolkits.mplot3d import axes3d
    #import matplotlib
    #import matplotlib.pyplot as plt
    #fig=plt.figure()
    #ax=fig.add_subplot(111,projection='3d')
    #ax.plot_wireframe(logD,logT,logE, color="green")
    #plt.show()
    ###########
    
    return [logE,logP]
class opacityTable:
  """Holds opacity table data.
  
  Initialize with a composition (X,Z), file name and weather the file name contains multiple.
  """
  
  def load(self):
    """Load from a file an opacity table for composition of the current opacity object. 
    It does this by advancing a file until the composition is matched and then calls 
    __loadTableFromFile to load the logR, logT, and logK values."""
    
    f=open(self.sFileName,'r')
    
    bContinue=True
    nCount=1
    while bContinue:
      
      #look for the right composition table
      line=f.readline()
      lineParts=line.split()
      nIndexX=-1
      nIndexZ=-1
      for i in range(len(lineParts)):
        if lineParts[i][0]=="X":# line has hydrogen mass fraction
          nIndexX=i
        if lineParts[i][0]=="Z":#line has metal mass fraction
          nIndexZ=i
      
      #if the line has a composition on it
      if nIndexX!=-1 and nIndexZ!=-1:
        
        #handle both the case where X=0.0 and X= 0.0, need to work with the space and without
        if len(lineParts[nIndexX])>2:
          xString=lineParts[nIndexX][2:]
        else:
          xString=lineParts[nIndexX+1]
        if len(lineParts[nIndexZ])>2:
          zString=lineParts[nIndexZ][2:]
        else:
          zString=lineParts[nIndexZ+1]
        
        #if there is a match between X and Z we have the right composition
        if float(xString)==self.X and float(zString)==self.Z:
          
          if nCount==2 and self.multitableFile:
            self.__loadTableFromFile(f)
            bContinue=False
          elif self.multitableFile==False:
            self.__loadTableFromFile(f)
            bContinue=False
          else:
            nCount=2
  def plotLogK(self,otherTables=None,logRIndex=None,wireFrame=True):
    """Plots opacity
    
    Keywords:
    otherTables: a list of opacity tables to also be ploted
    logRIndex: a list of integers used to indicate a specific logR index to plot 2D line plots at.
    """
    
    if logRIndex!=None:
      wireFrame=False
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib
    import matplotlib.pyplot as plt
    fig=plt.figure()
    if wireFrame:
      ax=fig.add_subplot(111,projection='3d')
      ax.plot_wireframe(self.logR,self.logT,self.logK)
      ax.set_xlabel("log10(R), R=rho/(T/1e6)^3 [g/cm^3]")
      ax.set_ylabel("log10(T) [K]")
      ax.set_zlabel("log10(K) [cm^2/g]")
      if otherTables:
        for otherTable in otherTables:
          ax.plot_wireframe(otherTable.logR,otherTable.logT,otherTable.logK,color="green")
    else:
      ax=fig.add_subplot(111)
      l=logRIndex[0]
      h=logRIndex[0]+1
      print self.logR[0][0]
      print self.logR[0,l:h]
      ax.plot(self.logT[:,l:h][0],self.logK[:,l:h][0], "bo-")
      counter=1
      ax.set_xlabel("log10(T)")
      if otherTables:
        for otherTable in otherTables:
          l=logRIndex[counter]
          h=logRIndex[counter]+1
          print otherTable.logR[0][0]
          print otherTable.logR[0,l:h]
          ax.plot(otherTable.logT[:,l:h][0],otherTable.logK[:,l:h][0], "go-")
          counter+=1
    plt.show()
  def interpolate(self,gridConfig,setExtrapolatedToNan=True):
    """Interpolate from self's table to the griding specified by:
    
    paramters:
    logDMin: first (smallest) logD value of grid
    logDDel: spacing in logD
    numLogD: number of logD grid points
    logTMin: first (smallest) logT value of grid
    logTDel: spacing in logT
    numLogT: number of logT grid points
    
    keyword:
    setExtrapolatedToNan: controls weather extrapolated points are set to nans (default is True)
    
    returns:
    an opacity table interpolated to the specified grid. In addition to the regular members of an
    opacity table logD is also included.
    """
    
    logDMin=gridConfig[0]
    logDDel=gridConfig[1]
    numLogD=gridConfig[2]
    logTMin=gridConfig[3]
    logTDel=gridConfig[4]
    numLogT=gridConfig[5]
    
    print "  interpolating opacity to new grid ..."
    
    #create a new interpolated table
    interpolatedTable=self.__class__()
    interpolatedTable.X=self.X
    interpolatedTable.Z=self.Z
    
    #save some of the interpolation information
    interpolatedTable.logDMin=logDMin
    interpolatedTable.logDDel=logDDel
    interpolatedTable.logTMin=logTMin
    interpolatedTable.logTDel=logTDel
    
    #make space
    interpolatedTable.logD=np.empty((numLogD,numLogT))
    interpolatedTable.logR=np.empty((numLogD,numLogT))
    interpolatedTable.logT=np.empty((numLogD,numLogT))
    interpolatedTable.logK=np.empty((numLogD,numLogT))
    
    #set logD and logT grid points
    for i in range(numLogD):
      for j in range(numLogT):
        interpolatedTable.logD[i][j]=logDMin+float(i)*logDDel
        interpolatedTable.logT[i][j]=logTMin+float(j)*logTDel
    
    #remove nans from grid
    #this part is slow
    print "    filling in nans by extrapolating ..."
    logK=self.__fillDepNans()
    
    #interpolate to grid
    print "    creating interpolator ..."
    splineInterpLogK=interpolate.RectBivariateSpline(self.logT[:,0],self.logR[0,:],logK)
    
    print "    interpolating to new grid ..."
    for i in range(numLogD):
      for j in range(numLogT):
        
        #calculate logR from logD and logT
        T6=10.0**interpolatedTable.logT[i][j]/1.0e6
        R=10.0**interpolatedTable.logD[i][j]/(T6**3)
        logR=math.log10(R)
        interpolatedTable.logR[i][j]=logR
        interpolatedTable.logK[i][j]=splineInterpLogK(interpolatedTable.logT[i][j],logR)
    
    #set values outside grid, and near nans on original grid equal to nan in interpolated grid
    if setExtrapolatedToNan:
      
      print "    setting nans for extrapolated values ..."
      
      #get range in logT, and logR
      logTRange=[min(self.logT[:,0]),max(self.logT[:,0])]
      logRRange=[min(self.logR[0,:]),max(self.logR[0,:])]
      np.set_printoptions(precision=8)
      for i in range(numLogD):
        for j in range(numLogT):
          
          #if logR and logT outside original table it is a nan
          bOutSide=False
          if interpolatedTable.logT[i][j]<logTRange[0] or logTRange[1]<interpolatedTable.logT[i][j]:
            bOutSide=True
          if interpolatedTable.logR[i][j]<logRRange[0] or logRRange[1]<interpolatedTable.logR[i][j]:
            bOutSide=True
          
          #if not outside see if it is near a nan in the oringal table
          if not bOutSide:
            
            #find location on original table in logR
            nIndexLogR=0
            for k in range(len(self.logR[0,:])-1):
              if self.logR[0][k]<=interpolatedTable.logR[i][j]\
                and interpolatedTable.logR[i][j]<self.logR[0][k+1]:
                nIndexLogR=k
                break
            
            #find location on original table in logT
            nIndexLogT=0
            for k in range(len(self.logT[:,0])-1):
              if self.logT[k][0]<=interpolatedTable.logT[i][j]\
                and interpolatedTable.logT[i][j]<self.logT[k+1][0]:
                nIndexLogT=k
                break
            
            #check if any of the nearest neighbors are nans if so current point should be a nan too
            bNearNan=False
            if np.isnan(self.logK[nIndexLogT][nIndexLogR]):
              bNearNan=True
            if nIndexLogR<self.logK.shape[1]-2:
              if np.isnan(self.logK[nIndexLogT][nIndexLogR+1]):
                bNearNan=True
            if nIndexLogT<self.logK.shape[0]-2:
              if np.isnan(self.logK[nIndexLogT+1][nIndexLogR]):
                bNearNan=True
            if nIndexLogT<self.logK.shape[0]-2 and nIndexLogR<self.logK.shape[1]-2:
              if np.isnan(self.logK[nIndexLogT+1][nIndexLogR+1]):
                bNearNan=True
            
            if bNearNan:
              interpolatedTable.logK[i][j]=np.nan
          else:
            interpolatedTable.logK[i][j]=np.nan
      
    #return interpolated table
    return interpolatedTable
  def __init__(self,X=None,Z=None,sFileName=None,multitableFile=None):
    """Initializes the opacity object.
    
    sets:
    self.X: the hydrogen mass fraction
    self.Z: the metal mass fraction
    self.sFileName: the file name to load the table from
    self.multitableFile: weather or not the file has more than one table in it
    """
    
    if X!=None:
      self.X=X
    if Z!=None:
      self.Z=Z
    if sFileName!=None:
      self.sFileName=sFileName
    if multitableFile!=None:
      self.multitableFile=multitableFile
  def __loadTableFromFile(self,f):
    """Loads an opacity table from a file object. This file object can contain multiple tables at 
    different compositions.
    
    This function should only be called from the load function.
    f : is the file object that is used to start reading the table.
    
    It loads the values into:
    self.logR: numpy array of log(densities[g/cm^3]/(T*1e-6)^3)
    self.logT: numpy array of log(T[K])
    self.logK: numpy array of log(opacities [cm^2/g])
    self.logR, self.logT, and self.logK are garanteed to be the same size. Array elements that are
    not defined have logK set to nans. This function assumes that there are only rectangular grids
    in logR and logT, with logK=9.999 for unknown opacities.
    """
    
    print "loading table with (X="+str(self.X)+",Z="+str(self.Z)+") from file, "\
      +self.sFileName+" ..."
    
    #skip the right number of lines, depends on the table
    line=f.readline()
    bContinue=True
    while bContinue:
      if line[0:5]=="log T":
        line=line[5:]
        bContinue=False
      elif  line[0:4]=="logT":
        line=line[4:]
        bContinue=False
      else:
        line=f.readline()
    
    #get Log(R); R=density/(T*1e-6)
    #line=f.readline()
    splitLine=line.split()
    logR1D=[]
    for entry in splitLine:
      logR1D.append(entry)
    
    #skip lines until they aren't empty
    bContinue=True
    line=""
    while bContinue:
      line=f.readline()
      
      #need to have more than a return characater, 2 makes it safe for windows return lines too.
      if len(line)>2:
        bContinue=False
    
    #read in opacities, temperatures, and densities
    bContinue=True
    logT1D=[]
    logK=[]
    splitLine=checkLineForMinusSplitError(line)
    
    while bContinue:
      
      #read line
      logT1D.append(float(splitLine[0]))
      KRow=[]
      for i in range(len(splitLine)-1):
        if 9.999!=float(splitLine[i+1]):#values of 9.99 are nans
          KRow.append(float(splitLine[i+1]))
        else:
          KRow.append(float('nan'))
      
      #add to 2D array
      logK.append(KRow)
      
      #get next line
      line=f.readline()
      splitLine=checkLineForMinusSplitError(line)
      
      #check if we are at the end
      if len(splitLine)<1:
        bContinue=False
    
    #find largest row dim
    #all rows sizes will be the same for T,R, and K
    #nMaxRow=0
    #for i in range(len(logK)):
    #  if len(logK[i])>nMaxRow:
    #    nMaxRow=len(logK[i])
    
    #allocate space for arrays
    self.logT=np.empty((len(logT1D),len(logR1D)))
    self.logR=np.empty((len(logT1D),len(logR1D)))
    self.logK=np.empty((len(logT1D),len(logR1D)))
    
    #fill in elements that are not nans
    for i in range(len(self.logR)):
      for j in range(len(self.logR[i])):
        self.logR[i][j]=logR1D[j]
        self.logT[i][j]=logT1D[i]
        if j<len(logK[i]):
          self.logK[i][j]=logK[i][j]
        else:
          self.logK[i][j]=np.nan
    
    #check order, AF tables are backwards
    if self.logT[0][0]>self.logT[1][0]:
      self.logR=np.flipud(self.logR)
      self.logT=np.flipud(self.logT)
      self.logK=np.flipud(self.logK)
  def __fillDepNans(self):
    """Extrapolates boundary grid points to Nan elements using an averaged weighted by a high power 
    of distance away from the boundary point. This function returns [logK] free of 
    nans."""
    
    #fill in nans with extrapolation only for purposes of the interpolation
    logK=np.copy(self.logK)
    
    '''find points along nan boundaries'''
    boundaryLogR=[]
    boundaryLogT=[]
    boundaryLogK=[]
    for j in range(len(self.logK)):
      for k in range(len(self.logK[j])):
        
        #only add non-nan points to boundary list
        if not np.isnan(self.logK[j][k]):
          
          '''if there is a nan either behind, or infront of the current point in either j or k 
          space it is on a boundary'''
          bOnBoundary=False
          if j>0:
            if np.isnan(self.logK[j-1][k]):
              bOnBoundary=True
          if j<len(self.logK)-1:
            if np.isnan(self.logK[j+1][k]):
              bOnBoundary=True
          if k>0:
            if np.isnan(self.logK[j][k-1]):
              bOnBoundary=True
          if k<len(self.logK[j])-1:
            if np.isnan(self.logK[j][k+1]):
              bOnBoundary=True
          if bOnBoundary:
            boundaryLogR.append(self.logR[j][k])
            boundaryLogT.append(self.logT[j][k])
            boundaryLogK.append(self.logK[j][k])
    
    '''extrapolate K outside boundary using a weighted average of boundary points. The weight 
    is simply the inverse of the distnace to a high power from a given boundary point'''
    for j in range(len(self.logK)):
      lastLogK=float('nan')
      
      for k in range(len(self.logK[j])):
        
        #for logK
        if np.isnan(self.logK[j][k]):
          sumK=0.0
          sumd=0.0
          for n in range(len(boundaryLogR)):
            distanceLogR=(self.logR[j][k]-boundaryLogR[n])/2.0
            distanceLogT=(self.logT[j][k]-boundaryLogT[n])
            distance=distanceLogR*distanceLogR+distanceLogT*distanceLogT
            distance=distance*distance
            distance=distance*distance
            distance=distance*distance
            oneDdistance=1.0/distance
            sumK+=oneDdistance*boundaryLogK[n]
            sumd+=oneDdistance
          #print sumE,sumd
          logK[j][k]=sumK/sumd
    
    ###########
    #DEBUG, to take a look at tables
    #from mpl_toolkits.mplot3d import axes3d
    #import matplotlib
    #import matplotlib.pyplot as plt
    #fig=plt.figure()
    #ax=fig.add_subplot(111,projection='3d')
    #ax.plot_wireframe(logD,logT,logE, color="green")
    #plt.show()
    ###########
    
    return logK
  def fillInDepNans(self):
    """Fills in logR and logT values to make a rectangular grid"""
    
    #find largest column in logT of non-nans
    columnLens=np.zeros(self.logT.shape[1])
    for i in range(self.logT.shape[0]):
      for j in range(self.logT.shape[1]):
        if not np.isnan(self.logT[i][j]):
          columnLens[j]+=1
    maxColumnIndex=columnLens.argmax()
    
    #fill in logT
    for i in range(self.logT.shape[0]):
      for j in range(self.logT.shape[1]):
        self.logT[i][j]=self.logT[i][maxColumnIndex]
    
    #find largest row in logR of non-nans
    rowLens=np.zeros(self.logR.shape[0])
    for i in range(self.logR.shape[0]):
      for j in range(self.logR.shape[1]):
        if not np.isnan(self.logR[i][j]):
          rowLens[i]+=1
    maxRowIndex=rowLens.argmax()
    
    #fill in logR
    for i in range(self.logR.shape[0]):
      for j in range(self.logR.shape[1]):
        self.logR[i][j]=self.logR[maxRowIndex][j]
class opacityTableManager:
  """Manages opacity files, including how they are interpolated between in composition."""
  
  def load(self):
    """Loads opacity files and merge files at duplicate compositions (i.e. merges low and high 
    temperature opacity tables).
    
    Sets the following:
    self.X: list of hydrogen mass fractions convered by opacity tables
    self.Z: list of metal mass fractions covered by opacity tables"""
    
    #load eos files
    nCount=0
    for opacityTable in self.opacityTables:
      #print nCount
      opacityTable.load()
      nCount+=1
    
    #sort list in X
    self.opacityTables.sort(key=lambda tempTable: tempTable.X)
    
    #sort list in Z
    self.opacityTables.sort(key=lambda tempTable: tempTable.Z)
    
    #combine opacity tables (low temperature AF tables, and OPAL tables)
    self.__merge()
    
    #make lists of composition
    self.__setCompLists()
  def interpComp(self,X,Z):
    """Interpolates a set of opacity files to the desired X and Z, and returns an the interpolated
    opacityTable.
    
    Parameters:
    X: hydrogen mass fraction
    Z: metal mass fraction
    """
    
    print "  interpolating opacity in composition to (X,Z)=("+str(X)+","+str(Z)+") ..."
    
    #option for printing all numpy arrays
    np.set_printoptions(threshold='nan')
    
    #if we have enough poitns in X and Z we can do 2D cubic spline interpolation
    if len(self.X)>=4 and len(self.Z)>=4:
      print "    interpolating in both X and Z ..."
      interpolatedOpacityTable=self.__bicubicSplineInXZ(X,Z)
      return interpolatedOpacityTable
    else:
      print "    interpolating with fewer than 4 X or Z points is not yet implemented"
      return None
  def plotGrids(self,opacityIndex):
    """Plot LogR and LogT points that form the opacity grid.
    
    Parameters:
    opacityIndex: a list of integers used to select which opacity tables will be plotted
    """
    
    #get eos
    import matplotlib.pyplot as plt
    fig=plt.figure()
    ax=fig.add_subplot(111)
    
    #loop over different eos files
    for index in opacityIndex:
      opacityTable=self.opacityTables[index]
    
      #make plot of rho and T
      ax.plot(opacityTable.logR,opacityTable.logT,'o')
    plt.show()
  def getTableFromComp(self,X,Z):
    """Returns a shallow copy of the opacity table with matching composition."""
    
    for opacityTableTemp in self.opacityTables:
      if opacityTableTemp.X==X and opacityTableTemp.Z==Z:
        return opacityTableTemp
    return None
  def __init__(self,opacityConfigFile=None):
    """Creates a new instance of opacityTableManager.
    
    If opacityConfigFile is set it will try to parse it for xml settings to get all the file names
    of the opacity files to include in the opacityTableManager."""
    
    if opacityConfigFile!=None:
      self.opacityConfigFileName=opacityConfigFile
      
      #get root element
      tree=xml.parse(self.opacityConfigFileName)
      root=tree.getroot()
      #switch to opacities node
      opacityElements=root.findall("opacity")
      if len(opacityElements)>1:
        print "WARNING: found "+str(len(opacityElements))+" \"opacity\" nodes, all but first node will "\
          +"be ignored."
      opacityElement=opacityElements[0]
      
      #get files element
      filesElements=opacityElement.findall("files")
      if len(filesElements)>1:#check that there is only one files element
        print "WARNING: found "+str(len(filesElements))+" \"files\" nodes under \"opacity\" node,"\
          +" all but first node will be ignored."
      
      #get opacity file names
      self.opacityFileNames={}
      fileElements=filesElements[0].findall("file")
      for fileElement in fileElements:
        id=fileElement.get("id")
        
        fileName=fileElement.text
        
        #if not absolute or relative add path assume it is in SPHERLS data path
        if fileName[0:2]!="./" and fileName[0:1]!="/":
          fileName=paths.SPHERLSDATA+fileName
        
        self.opacityFileNames[id]=fileName
      
      #get compositions in each file
      fileCompositions={}
      for id in self.opacityFileNames.keys():
        compositions=self.__getCompositions(self.opacityFileNames[id])
        fileCompositions[id]=compositions
      
      #get tables element
      tablesElement=opacityElement.findall("tables")
      if len(tablesElement)>1:#check that there is only one tables element
        print "WARNING: found "+str(len(tablesElement))+" \"tables\" nodes under \"opacity\" node,"\
          +" all but first node will be ignored."
      
      
      #create opacity tables
      self.opacityTables=[]
      tableElements=tablesElement[0].findall("table")
      for table in tableElements:
        fileID=table.get("fileID")
        
        #get X
        X=None
        try:
          X=float(table.get("X"))
        except ValueError:
          raise Exception("expecting a float in attribute \"X\" in \"table\" node under \"tables\""\
            +" node under \"opacity\" node . Could not convert \""+table.get("X")+"\" to a float.")
        
        #get Z
        Z=None
        try:
          Z=float(table.get("Z"))
        except ValueError:
          raise Exception("expecting a float in attribute \"Z\" in \"table\" node under \"tables\""\
            +" node under \"opacity\" node . Could not convert \""+table.get("Z")+"\" to a float.")
        
        #check that file with fileID has the correction composition X,Z in it
        if not ([X,Z] in fileCompositions[fileID]) :
          raise Exception("Unable to find the composition "+str([X,Z])+" in the file specified by"\
            +" fileID=\""+fileID+"\" ("+self.opacityFileNames[fileID]+").")
        
        self.opacityTables.append(opacityTable(X=X,Z=Z,sFileName=self.opacityFileNames[fileID]\
          ,multitableFile=(len(fileCompositions[fileID])>1)))
  def __merge(self):
    """Merges opacity tables of equal compositions together through out the whole set of tables."""
    
    #find two tables with composition X,Z
    table1Index=None
    index1=0
    for table in self.opacityTables:
      table1Index=index1
      
      #find matching composition, ordered so shouldn't have to go far
      index2=0
      table2Index=None
      for tableTest in self.opacityTables:
        
        if table.X==tableTest.X and table.Z==tableTest.Z and index2!=table1Index:
          table2Index=index2
          break
        index2+=1
      index1+=1
      self.opacityTables[table1Index]=self.__merge2files(table1Index,table2Index)
      del self.opacityTables[table2Index]
    return
  def __merge2files(self,table1Index,table2Index):
    """Merge two opacity tables together.
    
    returns an opacityTable containing the merged table
    """
    
    '''find start indices for logR of table 1 and table 2 and how many logR points the merged table
    will have in logR direction'''
    numLogR=0
    startIndexTable1LogR=0
    startIndexTable2LogR=0
    if self.opacityTables[table1Index].logR[0][0]<self.opacityTables[table2Index].logR[0][0]:
      for i in range(self.opacityTables[table1Index].logR.shape[1]):
        if self.opacityTables[table1Index].logR[0][i]==self.opacityTables[table2Index].logR[0][0]:
          numLogR=i+self.opacityTables[table2Index].logR.shape[1]
          startIndexTable2LogR=i
          break
    else:
      for i in range(self.opacityTables[table2Index].logR.shape[1]):
        if self.opacityTables[table2Index].logR[0][i]==self.opacityTables[table1Index].logR[0][0]:
          numLogR=i+self.opacityTables[table1Index].logR.shape[1]
          startIndexTable1LogR=i
          break
    
    '''find start indices for logT of table 1 and table 2 and how many logR points the merged table
    will have in logT direction'''
    numLogT=0
    startIndexTable1LogT=0
    startIndexTable2LogT=0
    if self.opacityTables[table1Index].logT[0][0]<self.opacityTables[table2Index].logT[0][0]:
      for i in range(self.opacityTables[table1Index].logT.shape[0]):
        if self.opacityTables[table1Index].logT[i][0]==self.opacityTables[table2Index].logT[0][0]:
          numLogT=i+self.opacityTables[table2Index].logT.shape[0]
          startIndexTable2LogT=i
          break
    else:
      for i in range(self.opacityTables[table2Index].logT.shape[0]):
        if self.opacityTables[table2Index].logT[i][0]==self.opacityTables[table1Index].logT[0][0]:
          numLogT=i+self.opacityTables[table1Index].logT.shape[0]
          startIndexTable1LogT=i
          break
    
    #set intependent variables of new grid
    mergedOpacityTable=opacityTable(X=self.opacityTables[table1Index].X\
      ,Z=self.opacityTables[table1Index].Z)
    
    #initialize logR, logT, and logK
    mergedOpacityTable.logR=np.empty((numLogT,numLogR))
    mergedOpacityTable.logR[:]=np.NAN
    mergedOpacityTable.logT=np.empty((numLogT,numLogR))
    mergedOpacityTable.logT[:]=np.NAN
    mergedOpacityTable.logK=np.empty((numLogT,numLogR))
    mergedOpacityTable.logK[:]=np.NAN
    
    #merge logR, logT, and logK
    for j in range(numLogR):
      for i in range(numLogT):
        
        bInTable1=False
        bInTable2=False
        
        #table set logR and logT from table 1 or table 2 depending which is in range
        if (i>=startIndexTable1LogT and i<self.opacityTables[table1Index].logR.shape[0]+startIndexTable1LogT)\
          and (j>=startIndexTable1LogR and j<self.opacityTables[table1Index].logR.shape[1]+startIndexTable1LogR):
          mergedOpacityTable.logR[i][j]=\
            self.opacityTables[table1Index].logR[i-startIndexTable1LogT][j-startIndexTable1LogR]
          mergedOpacityTable.logT[i][j]=\
            self.opacityTables[table1Index].logT[i-startIndexTable1LogT][j-startIndexTable1LogR]
          bInTable1=True
        if (i>=startIndexTable2LogT and i<self.opacityTables[table2Index].logR.shape[0]+startIndexTable2LogT)\
          and (j>=startIndexTable2LogR and j<self.opacityTables[table2Index].logR.shape[1]+startIndexTable2LogR):
          mergedOpacityTable.logR[i][j]=\
            self.opacityTables[table2Index].logR[i-startIndexTable2LogT][j-startIndexTable2LogR]
          mergedOpacityTable.logT[i][j]=\
            self.opacityTables[table2Index].logT[i-startIndexTable2LogT][j-startIndexTable2LogR]
          bInTable2=True
        
        #merge logK
        if bInTable1 and bInTable2:#merge logK in overlapping regin
          
          '''figure out values of logT at the edges of the overlap. Used in computing weight for averaging
          of two logK'''
          logT0=None#lower bound of overlap in T
          logT1=None#upper bound of overlap in T
          logTOverlapLowerEdgeTable1=False
          
          #set lower logT for table 1 at this density
          logT10=self.opacityTables[table1Index].logT[0][j-startIndexTable1LogR]
          logK10=self.opacityTables[table1Index].logK[0][j-startIndexTable1LogR]
          indexi=0
          while np.isnan(logK10):#remove any nans from logT10
            logK10=self.opacityTables[table1Index].logK[indexi][j-startIndexTable1LogR]
            logT10=self.opacityTables[table1Index].logT[indexi][j-startIndexTable1LogR]
            indexi+=1#move up
          
          #set upper logT for table 1 at this density
          logT11=self.opacityTables[table1Index].logT[self.opacityTables[table1Index].logT.shape[0]-1][j-startIndexTable1LogR]
          logK11=self.opacityTables[table1Index].logK[self.opacityTables[table1Index].logT.shape[0]-1][j-startIndexTable1LogR]
          indexi=self.opacityTables[table1Index].logT.shape[0]-1
          while np.isnan(logK11):#remove any nans from logT11
            logT11=self.opacityTables[table1Index].logT[indexi][j-startIndexTable1LogR]
            logK11=self.opacityTables[table1Index].logK[indexi][j-startIndexTable1LogR]
            indexi-=1#move down
          
          #set lower logT for table 2 at this density
          logT20=self.opacityTables[table2Index].logT[0][j-startIndexTable2LogR]
          logK20=self.opacityTables[table2Index].logK[0][j-startIndexTable2LogR]
          indexi=0
          while np.isnan(logK20):#remove any nans from logT10
            logT20=self.opacityTables[table2Index].logT[indexi][j-startIndexTable2LogR]
            logK20=self.opacityTables[table2Index].logK[indexi][j-startIndexTable2LogR]
            indexi+=1#move up
          
          #set upper logT for table 2 at this density
          logT21=self.opacityTables[table2Index].logT[self.opacityTables[table2Index].logT.shape[0]-1][j-startIndexTable2LogR]
          logK21=self.opacityTables[table2Index].logK[self.opacityTables[table2Index].logT.shape[0]-1][j-startIndexTable2LogR]
          indexi=self.opacityTables[table2Index].logT.shape[0]-1
          while np.isnan(logK21):#remove any nans from logT11
            logT21=self.opacityTables[table2Index].logT[indexi][j-startIndexTable2LogR]
            logK21=self.opacityTables[table2Index].logK[indexi][j-startIndexTable2LogR]
            indexi-=1#move down
          
          #pick upper and lower boundaries for overlap region
          if (logT10-logT11)<(logT20-logT21):
            logT0=logT10
            logT1=logT21
          else :
            logT0=logT20
            logT1=logT11
            logTOverlapLowerEdgeTable1=True
          
          #use simple weighted average
          w1=(logT0-mergedOpacityTable.logT[i][j])
          w1=w1*w1
          w2=(logT1-mergedOpacityTable.logT[i][j])
          w2=w2*w2
          
          #switch weights if lower bound of overlap is from table 2
          if logTOverlapLowerEdgeTable1:
            wTemp=w1
            w1=w2
            w2=wTemp
          
          #if one is nan use the other
          logK1=self.opacityTables[table1Index].logK[i-startIndexTable1LogT][j-startIndexTable1LogR]
          logK2=self.opacityTables[table2Index].logK[i-startIndexTable2LogT][j-startIndexTable2LogR]
          if not(np.isnan(logK1) and np.isnan(logK2)):#only set =0.0 if both aren't nans
            if np.isnan(logK1):
              logK1=0.0
              w1=0.0
            if np.isnan(logK2):
              logK2=0.0
              w2=0.0
          
          mergedOpacityTable.logK[i][j]=(w1*logK1+w2*logK2)/(w1+w2)
        elif bInTable1:#use table 1 logK
          
          #use table 1 logK
          mergedOpacityTable.logK[i][j]=\
            self.opacityTables[table1Index].logK[i-startIndexTable1LogT][j-startIndexTable1LogR]
        elif bInTable2:#use table 2 logK
          
          #use table 2 logK
          mergedOpacityTable.logK[i][j]=\
            self.opacityTables[table2Index].logK[i-startIndexTable2LogT][j-startIndexTable2LogR]
    
      '''#DEBUG, to take a look at tables
      bPlotWireFrame=False
      from mpl_toolkits.mplot3d import axes3d
      import matplotlib
      import matplotlib.pyplot as plt
      fig=plt.figure()
      if bPlotWireFrame:
        ax=fig.add_subplot(111,projection='3d')
        #ax.set_xlim3d([-3,1])
        ax.plot_wireframe(mergedOpacityTable.logR,mergedOpacityTable.logT,mergedOpacityTable.logK, color="green")
        ax.plot_wireframe(self.opacityTables[table1Index].logR,self.opacityTables[table1Index].logT,self.opacityTables[table1Index].logK, color="red")
        ax.plot_wireframe(self.opacityTables[table2Index].logR,self.opacityTables[table2Index].logT,self.opacityTables[table2Index].logK, color="blue")
        #ax.set_xlim3d([-1,1])
        #ax.set_ylim3d([3.5,4.5])
        #ax.set_zlim3d([0,4])
      else:
        ax=fig.add_subplot(111)
        l=3
        h=5
        #print "table 1 logR=",self.opacityTables[table1Index].logR[:,l-startIndexTable1LogR:h-startIndexTable1LogR]
        #print "table 2 logR=",self.opacityTables[table2Index].logR[:,l-startIndexTable2LogR:h-startIndexTable2LogR]
        #print "merged table logR=",mergedTable.logR[:,l:h]
        ax.plot(mergedOpacityTable.logT[:,l:h],mergedOpacityTable.logK[:,l:h], "g-o")
        ax.plot(self.opacityTables[table1Index].logT[:,l-startIndexTable1LogR:h-startIndexTable1LogR]\
          ,self.opacityTables[table1Index].logK[:,l-startIndexTable1LogR:h-startIndexTable1LogR],"r-o")
        ax.plot(self.opacityTables[table2Index].logT[:,l-startIndexTable2LogR:h-startIndexTable2LogR]\
          ,self.opacityTables[table2Index].logK[:,l-startIndexTable2LogR:h-startIndexTable2LogR], "b-o")
      plt.show()
      ###'''
    
    mergedOpacityTable.fillInDepNans()
    return mergedOpacityTable
  def __cubicSplineInX(self,eosManagerInterpZ,X):
    """Interpolates a set of equations of state which are at the same Z to a particular X using
    cubic spline interpolation."""
    
    #check if non-nan shapes are the same across all X's
    logDIsNan=np.isnan(eosManagerInterpZ.eosFiles[0].logD)
    for i in range(1,len(eosManagerInterpZ.eosFiles)):
      if not np.array_equal(logDIsNan,np.isnan(eosManagerInterpZ.eosFiles[i].logD)):
        print "WARNING: table 0, and table ",i," have different nan shapes, will introduce more nans"
      if eosManagerInterpZ.eosFiles[i].Z!=eosManagerInterpZ.eosFiles[0].Z:
        raise Exception("file 0 has Z of "+str(eosManagerInterpZ.eosFiles[0].Z)+" while file "\
          +str(i)+" has Z of "+str(eosManagerInterpZ.eosFiles[i].Z)+" should have same Z.")
    
    eosInterp=eosTables()
    eosInterp.Z=eosManagerInterpZ.eosFiles[0].Z
    eosInterp.X=X
    logE=[]
    logP=[]
    for i in range(eosManagerInterpZ.eosFiles[0].logD.shape[0]):
      logERow=[]
      logPRow=[]
      for j in range(eosManagerInterpZ.eosFiles[0].logD.shape[1]):
        
        #get values from all X tables
        xlist=[]
        logElist=[]
        logPlist=[]
        for eosFile in eosManagerInterpZ.eosFiles:
          xlist.append(eosFile.X)
          logElist.append(eosFile.logE[i][j])
          logPlist.append(eosFile.logP[i][j])
        
        #interpolate the energy in X
        x=np.array(xlist)
        y=np.array(logElist)
        tck=interpolate.splrep(x,y,s=0)
        logERow.append(interpolate.splev(X,tck,der=0))
        
        #interpolate the pressure in X
        y=np.array(logPlist)
        tck=interpolate.splrep(x,y,s=0)
        logPRow.append(interpolate.splev(X,tck,der=0))
      
      #add rows
      logE.append(logERow)
      logP.append(logPRow)
    
    #convert to numpy arrays and add
    eosInterp.logD=np.empty(eosManagerInterpZ.eosFiles[0].logD.shape)
    eosInterp.logD[:]=float('nan')
    eosInterp.logT=np.empty(eosManagerInterpZ.eosFiles[0].logD.shape)
    eosInterp.logT[:]=float('nan')
    eosInterp.logE=np.empty(eosManagerInterpZ.eosFiles[0].logD.shape)
    eosInterp.logE[:]=float('nan')
    eosInterp.logP=np.empty(eosManagerInterpZ.eosFiles[0].logD.shape)
    eosInterp.logP[:]=float('nan')
    for i in range(len(logE)):
      for j in range(len(logE[i])):
        if not np.isnan(logE[i][j]) and not np.isnan(logP[i][j]):
          eosInterp.logE[i][j]=logE[i][j]
          eosInterp.logP[i][j]=logP[i][j]
          eosInterp.logD[i][j]=eosManagerInterpZ.eosFiles[0].logD[i][j]
          eosInterp.logT[i][j]=eosManagerInterpZ.eosFiles[0].logT[i][j]
    
    return eosInterp
  def __cubicSplineInZ(self,eosManagerInterpZ,Z):
    """Interpolates a set of equations of state which are at the same Z to a particular X using
    cubic spline interpolation."""
    
    #check if non-nan shapes are the same across all X's
    logDIsNan=np.isnan(eosManagerInterpZ.eosFiles[0].logD)
    for i in range(1,len(eosManagerInterpZ.eosFiles)):
      if not np.array_equal(logDIsNan,np.isnan(eosManagerInterpZ.eosFiles[i].logD)):
        print "WARNING: table 0, and table ",i," have different nan shapes, will introduce more nans"
      if eosManagerInterpZ.eosFiles[i].Z!=eosManagerInterpZ.eosFiles[0].Z:
        raise Exception("file 0 has Z of "+str(eosManagerInterpZ.eosFiles[0].Z)+" while file "\
          +str(i)+" has Z of "+str(eosManagerInterpZ.eosFiles[i].Z)+" should have same Z.")
    
    eosInterp=eosTables()
    eosInterp.Z=eosManagerInterpZ.eosFiles[0].Z
    eosInterp.X=X
    logE=[]
    logP=[]
    for i in range(eosManagerInterpZ.eosFiles[0].logD.shape[0]):
      logERow=[]
      logPRow=[]
      for j in range(eosManagerInterpZ.eosFiles[0].logD.shape[1]):
        
        #get values from all X tables
        xlist=[]
        logElist=[]
        logPlist=[]
        for eosFile in eosManagerInterpZ.eosFiles:
          xlist.append(eosFile.X)
          logElist.append(eosFile.logE[i][j])
          logPlist.append(eosFile.logP[i][j])
        
        #interpolate the energy in X
        x=np.array(xlist)
        y=np.array(logElist)
        tck=interpolate.splrep(x,y,s=0)
        logERow.append(interpolate.splev(X,tck,der=0))
        
        #interpolate the pressure in X
        y=np.array(logPlist)
        tck=interpolate.splrep(x,y,s=0)
        logPRow.append(interpolate.splev(X,tck,der=0))
      
      #add rows
      logE.append(logERow)
      logP.append(logPRow)
    
    #convert to numpy arrays and add
    eosInterp.logD=np.empty(eosManagerInterpZ.eosFiles[0].logD.shape)
    eosInterp.logD[:]=float('nan')
    eosInterp.logT=np.empty(eosManagerInterpZ.eosFiles[0].logD.shape)
    eosInterp.logT[:]=float('nan')
    eosInterp.logE=np.empty(eosManagerInterpZ.eosFiles[0].logD.shape)
    eosInterp.logE[:]=float('nan')
    eosInterp.logP=np.empty(eosManagerInterpZ.eosFiles[0].logD.shape)
    eosInterp.logP[:]=float('nan')
    for i in range(len(logE)):
      for j in range(len(logE[i])):
        if not np.isnan(logE[i][j]) and not np.isnan(logP[i][j]):
          eosInterp.logE[i][j]=logE[i][j]
          eosInterp.logP[i][j]=logP[i][j]
          eosInterp.logD[i][j]=eosManagerInterpZ.eosFiles[0].logD[i][j]
          eosInterp.logT[i][j]=eosManagerInterpZ.eosFiles[0].logT[i][j]
    
    return eosInterp
  def __bicubicSplineInXZ(self,X,Z):
    """Performs 2D cubic spline interpolation in X and Z to produce an opaicty table at the 
    specified composition.
    
    Parameters:
    X: hydrogen mass fraction to interpolate to
    Z: metal mass fraction to interpolate to
    """
    
    interpOpacityTable=opacityTable(X,Z)
    
    #for each point on logT, logR grid

    logK=[]
    #for opacityTableTemp in self.opacityTables:
    #  print "X=",opacityTableTemp.X," Z=",opacityTableTemp.Z
    #print len(self.opacityTables)
    for i in range(self.opacityTables[0].logR.shape[0]):
      logKRow=[]
      for j in range(self.opacityTables[0].logR.shape[1]):
        #x=np.empty((len(self.Z),len(self.X)))
        #z=np.empty(x.shape)
        k=np.empty((len(self.Z),len(self.X)))
        nTable=0
        #print "i,j=",i,j
        bGotNan=False
        for l in range(len(self.Z)):
          for m in range(len(self.X)):
            #x[l][m]=self.X[m]
            #z[l][m]=self.Z[l]
            #print "l=",l," m=",m," i=",i," j=",j," nTable=",nTable," X=",self.X[m]," Z=",self.Z[l]
            k[l][m]=self.opacityTables[nTable].logK[i][j]
            if np.isnan(k[l][m]):
              bGotNan=True
              break
            nTable+=1
          if bGotNan:
            break
        if not bGotNan:
          
          ############
          #DEBUG
          #from mpl_toolkits.mplot3d import axes3d
          #import matplotlib
          #import matplotlib.pyplot as plt
          #fig=plt.figure()
          #ax=fig.add_subplot(111,projection='3d')
          #ax.plot_wireframe(x,z,k)
          #plt.show()
          ##########
          
          #print "  creating interpolator ..."
          splineInterpLogK=interpolate.RectBivariateSpline(self.Z,self.X,k)
          logKNew=splineInterpLogK(Z,X)
          logKRow.append(logKNew)
          
        else:
          logKRow.append(np.nan)
      logK.append(logKRow)
    
    #set values of new opacity table
    interpOpacityTable.logR=self.opacityTables[0].logR
    interpOpacityTable.logT=self.opacityTables[0].logT
    interpOpacityTable.logK=np.empty(self.opacityTables[0].logR.shape)
    #print logK
    for i in range(interpOpacityTable.logR.shape[0]):
      for j in range(interpOpacityTable.logR.shape[1]):
        interpOpacityTable.logK[i][j]=logK[i][j]
    
    return interpOpacityTable
  def __getCompositions(self,sFileName):
    """Returns a list of compositions found in the given opacity file.
    
    Parameters:
    sFileName: the name of the file to search for (X,Z) compositions.
    """
    
    f=open(sFileName,'r')
    
    nCurrentTable=0
    compositions=[]
    #get composition
    for line in f:
      
      lineParts=line.split()
      nIndexX=-1
      nIndexZ=-1
      for i in range(len(lineParts)):
        if lineParts[i][0]=="X":# line has hydrogen mass fraction
          nIndexX=i
        if lineParts[i][0]=="Z":#line has metal mass fraction
          nIndexZ=i
      
      #if the line has a composition on it
      if nIndexX!=-1 and nIndexZ!=-1:
        
        #if it is a table entry
        if lineParts[0][0:5]=="TABLE":
          
          #check to make sure this table index isn't less than current table index
          nTempTable=0
          
          #combine second and third parts of the line if needed
          if lineParts[1]=="#":
            nTempTable=int(lineParts[2])
          else:
            nTempTable=int(lineParts[1][1:])
          
          if nTempTable<nCurrentTable:#getting second set of compositions from the actual tables
            break;
          else:
            nCurrentTable=nTempTable
        
        #handle both the case where X=0.0 and X= 0.0, need to work with the space and without
        if len(lineParts[nIndexX])>2:
          xString=lineParts[nIndexX][2:]
        else:
          xString=lineParts[nIndexX+1]
        if len(lineParts[nIndexZ])>2:
          zString=lineParts[nIndexZ][2:]
        else:
          zString=lineParts[nIndexZ+1]
          
        #add composition to list
        compositions.append([float(xString),float(zString)])
        
    return compositions
  def __setCompLists(self):
    """
    Sets composition lists based on the compositions of the tables.
    
    Sets:
    self.X: hydrogen mass fraction
    self.Z: mettal mass fraction"""
    
    #require a rectangular grid
    self.Z=[]
    self.X=[]
    zCurrent=self.opacityTables[0].Z
    self.Z.append(zCurrent)
    nXIndex=0
    for table in self.opacityTables:
      
      #at next Z
      if table.Z>zCurrent:
        zCurrent=table.Z;
        self.Z.append(zCurrent)
        nXIndex=0
      
      #if still on first Z
      if len(self.Z)==1:
        self.X.append(table.X)
      
      #if not on first Z check X's against those at first Z. Expecting rectangular grid so X's 
      #should be the same
      else:
        if table.X!=self.X[nXIndex]:
          raise Exception("current opacity table \""+table.sFileName+"\" has X of "+str(table.X)\
            +" expecting "+str(self.X[nXIndex])+" for a rectangular grid.")
        nXIndex+=1
class eosTableManager:
  """Manages equation of state files, including how they are interpolated between."""
  
  def load(self):
    """Loads eos files.
    
    Sets the following:
    self.Z: a list of Z (metal mass fraction) values of the equation of state files
    self.X: a list of X (hydrogen mass fraction) values of the equation of state files
    """
    
    #load eos files
    for eosTable in self.eosTables:
      eosTable.load()
    
    #sort list in X
    self.eosTables.sort(key=lambda tempEOS: tempEOS.X)
    
    #sort list in Z
    self.eosTables.sort(key=lambda tempEOS: tempEOS.Z)
    
    #require a rectangular grid
    self.Z=[]
    self.X=[]
    zCurrent=self.eosTables[0].Z
    self.Z.append(zCurrent)
    nXIndex=0
    for eosFile in self.eosTables:
      
      #at next Z
      if eosFile.Z>zCurrent:
        zCurrent=eosFile.Z;
        self.Z.append(zCurrent)
        nXIndex=0
      
      #if still on first Z
      if len(self.Z)==1:
        self.X.append(eosFile.X)
      
      #if not on first Z check X's against those at first Z. Expecting rectangular grid so X's 
      #should be the same
      else:
        if eosFile.X!=self.X[nXIndex]:
          raise Exception("current eosFile \""+eosFile.sFileName+"\" has X of "+str(eosFile.X)+" expecting "\
            +str(self.X[nXIndex])+" for a rectangular grid.")
        nXIndex+=1
  def interpComp(self,X,Z):
    """Interpolates a set of eos files and opacities to the desired X and Z, and returns an 
    eosManager with this new set of files which can then be interpolated to the desired rho and 
    T's."""
    
    print "  interpolating eos in composition to (X,Z)=("+str(X)+","+str(Z)+") ..."
    
    #option for printing all numpy arrays
    np.set_printoptions(threshold='nan')
    
    #interpolate to requested Z
    print "    interpolating in Z ..."
    if len(self.Z)==3:
      
      '''if self.Z has 3 entries use a quadratic interpolation to interpolate a set of tables to the 
      new Z.
      '''
      eosManagerInterpZ=self.__class__()
      eosManagerInterpZ.Z=[Z]
      eosManagerInterpZ.X=self.X
      eosManagerInterpZ.eosTables=[]
      
      #create a new eos file for each X
      for i in range(len(self.X)):
        index1=i
        index2=i+len(self.X)
        index3=i+2*len(self.X)
        eosManagerInterpZ.eosTables.append(self.__quadInterpInZ(Z,index1,index2,index3))
    else:
      raise Exception("interpolation between more than or less than 3 different Z values is not "\
        +"yet support")
    
    #now to interpolate to new X
    '''to do cubic spline interp can't have nans or will propegate throughout the calculation.
    So, keep track of the rho and T they are at, and use constant extrapolation in (rho or T?) 
    to fill in nans until after interpolation, will just effect the value of the derivatives a 
    little near the nans. non-nan elements of logD, logT, logE, and logP are garanteed to be the
    same shape for a given X because of the way Z interpolation was done in __quadInterpInZ.
    only need to keep shape of one. This might not be the case across X's.'''
    print "    interpolating in X ..."
    return self.__cubicSplineInX(eosManagerInterpZ,X)
  def plotGrid(self,eosIndex):
    """Plot rho and T points that form the grid"""
    
    #get eos
    import matplotlib.pyplot as plt
    fig=plt.figure()
    ax=fig.add_subplot(111)
    
    #loop over different eos files
    for index in eosIndex:
      eosData=self.eosTables[index]
    
      #make plot of rho and T
      ax.plot(eosData.density,eosData.temperature,'o')
    plt.show()
  def getTableFromComp(self,X,Z):
    """Returns a shallow copy of the eos table with matching composition. If none found it returns
    None."""
    
    for eosTableTemp in self.eosTables:
      if eosTableTemp.X==X and eosTableTemp.Z==Z:
        return eosTableTemp
    return None
  def __init__(self,eosFileName=None):
    """Returns a new instance of eosTableManager.
    
    if eosFileName is set it will call __initFromFile to load settings from a file to initialize
    the new eosTableManager."""
    if eosFileName!=None:
      self.__initFromFile(eosFileName)
  def __initFromFile(self,eosFileName):
    """parse xml settings file to get all the file names of the eos files
    
    sets the following:
    self.eosFileName: the name of the xml file name specifying all the eos files
    self.eosFiles: a list of eos objects for each of the of the files given in self.eosFileName
    
    """
    
    self.eosFileName=eosFileName
    
    #get root element
    tree=xml.parse(self.eosFileName)
    root=tree.getroot()
    
    #switch to first eos element
    eosElements=root.findall("eos")
    if len(eosElements)>1:
      print "WARNING: found "+str(len(eosElements))+" \"eos\" nodes, all but first node will be"\
        +" ignored."
    eosElement=eosElements[0]
    
    #get files element
    filesElements=eosElement.findall("files")
    if len(filesElements)>1:#check that there is only one files element
      print "WARNING: found "+str(len(filesElements))+" \"files\" nodes under \"eos\" node, all"\
        +" but first node will be ignored."
    
    #get file names, and create eos objects for each
    self.eosTables=[]
    fileElements=filesElements[0].findall("file")
    for fileElement in fileElements:
      fileName=fileElement.text
      
      #if not absolute or relative add path assume it is in SPHERLS data path
      if fileName[0:2]!="./" and fileName[0:1]!="/":
        fileName=paths.SPHERLSDATA+fileName
      self.eosTables.append(eosTable(fileName))
  def __cubicSplineInX(self,eosManagerInterpZ,X):
    """Interpolates a set of equations of state which are at the same Z to a particular X using
    cubic spline interpolation."""
    
    #check if non-nan shapes are the same across all X's
    logDIsNan=np.isnan(eosManagerInterpZ.eosTables[0].logD)
    for i in range(1,len(eosManagerInterpZ.eosTables)):
      if not np.array_equal(logDIsNan,np.isnan(eosManagerInterpZ.eosTables[i].logD)):
        print "WARNING: table 0, and table ",i," have different nan shapes, will introduce more nans"
      if eosManagerInterpZ.eosTables[i].Z!=eosManagerInterpZ.eosTables[0].Z:
        raise Exception("file 0 has Z of "+str(eosManagerInterpZ.eosTables[0].Z)+" while file "\
          +str(i)+" has Z of "+str(eosManagerInterpZ.eosTables[i].Z)+" should have same Z.")
    
    eosInterp=eosTable()
    eosInterp.Z=eosManagerInterpZ.eosTables[0].Z
    eosInterp.X=X
    logE=[]
    logP=[]
    for i in range(eosManagerInterpZ.eosTables[0].logD.shape[0]):
      logERow=[]
      logPRow=[]
      for j in range(eosManagerInterpZ.eosTables[0].logD.shape[1]):
        
        #get values from all X tables
        xlist=[]
        logElist=[]
        logPlist=[]
        for eosFile in eosManagerInterpZ.eosTables:
          xlist.append(eosFile.X)
          logElist.append(eosFile.logE[i][j])
          logPlist.append(eosFile.logP[i][j])
        
        #interpolate the energy in X
        x=np.array(xlist)
        y=np.array(logElist)
        tck=interpolate.splrep(x,y,s=0)
        logERow.append(interpolate.splev(X,tck,der=0))
        
        #interpolate the pressure in X
        y=np.array(logPlist)
        tck=interpolate.splrep(x,y,s=0)
        logPRow.append(interpolate.splev(X,tck,der=0))
      
      #add rows
      logE.append(logERow)
      logP.append(logPRow)
    
    #convert to numpy arrays and add
    eosInterp.logD=np.empty(eosManagerInterpZ.eosTables[0].logD.shape)
    eosInterp.logD[:]=float('nan')
    eosInterp.logT=np.empty(eosManagerInterpZ.eosTables[0].logD.shape)
    eosInterp.logT[:]=float('nan')
    eosInterp.logE=np.empty(eosManagerInterpZ.eosTables[0].logD.shape)
    eosInterp.logE[:]=float('nan')
    eosInterp.logP=np.empty(eosManagerInterpZ.eosTables[0].logD.shape)
    eosInterp.logP[:]=float('nan')
    for i in range(len(logE)):
      for j in range(len(logE[i])):
        
        eosInterp.logD[i][j]=eosManagerInterpZ.eosTables[0].logD[i][j]
        eosInterp.logT[i][j]=eosManagerInterpZ.eosTables[0].logT[i][j]
        eosInterp.logE[i][j]=logE[i][j]
        eosInterp.logP[i][j]=logP[i][j]
    
    return eosInterp
  def __interpInZ(self,z):
    """Figures out which equation of state files bracket z  and determines which 3 equation of state
    files to use for quadtratic interpolation. Only needed if there are more than 3 Z tables"""
    
    #make sure z is bigger than zero
    if z<self.Z[0]:
      raise Exception("the value of z to interpolate to "+str(z)
        +" is less than the lowest Z in the table "+str(self.Z[0]))
    
    #find bracketing Z indices, and check for out of bounds Z
    zIndexHigh=0
    for i in range(1,len(self.Z)):
      if self.Z[i]>z:#if grid Z bigger than z interpolating to then we found the upper bound
        zIndexHigh=i
        break
    
    #check that it isn't outside upper bound of table
    if zIndexHigh==0:#out side upper bound of Z in table
      raise Exception("the value of z to interpolate to "+str(z)
        +" larger than the highest Z in the table "+str(self.Z[len(self.Z)-1]))
    
    #low index is simply one before upper index
    zIndexLow=zIndexHigh-1
    
    if zIndexLow==0:#at inner edge of table in Z
      zIndex1=0
      zIndex2=1
      zIndex3=2
    elif zIndexHigh==len(self.Z)-1:#at outer edge of table in Z
      zIndex1=len(self.Z)-4
      zIndex2=len(self.Z)-3
      zIndex3=len(self.Z)-2
    else:#in between, just use the next point
      zIndex1=zIndexLow
      zIndex2=zIndexHigh
      zIndex3=zIndex+1
    
    #create an interpolated eos object for each X
    interpInZEOSs=[]
    for i in range (len(self.X)):
      fileIndex1=i+zIndex1*len(self.X)
      fileIndex2=i+zIndex2*len(self.X)
      fileIndex3=i+zIndex3*len(self.X)
      interpInZEOSs.append(self.__quadInterpInZ(z,fileIndex1,fileIndex2,fileIndex3))
    return interpInZEOSs
  def __quadInterpInZ(self,z,fileIndex1,fileIndex2,fileIndex3):
    """Interpolates between three equations of state quadratically in z to produce a new equation of
    state"""
    
    #check that they all have the same X
    if self.eosTables[fileIndex1].X!=self.eosTables[fileIndex2].X:
      raise Exception("file \""+self.eosTables[fileIndex1].sFileName+"\" has X= "
        +str(self.eosTables[fileIndex1].X)+" while X of \""\
        +self.eosTables[fileIndex2].sFileName+"\" is "+str(self.eosTables[fileIndex2].X)\
        +" must have the same X!")
    if self.eosTables[fileIndex1].X!=self.eosTables[fileIndex3].X:
      raise Exception("file \""+self.eosTables[fileIndex1].sFileName+"\" has X= "
        +str(self.eosTables[fileIndex1].X)+" while X of \""\
        +self.eosTables[fileIndex3].sFileName+"\" is "+str(self.eosTables[fileIndex3].X)\
        +" must have the same X!")
    
    #check that the three equation of states have the same number of points
    if self.eosTables[fileIndex1].logD.shape!=self.eosTables[fileIndex2].logD.shape:
      raise Exception("shape of \""+self.eosTables[fileIndex1].sFileName+"\" is "
        +str(self.eosTables[fileIndex1].logD.shape)+" while shape of \""\
        +self.eosTables[fileIndex2].sFileName+"\" is "+str(self.eosTables[fileIndex2].logD.shape)\
        +" shape of the two tables must be the same!")
    if self.eosTables[fileIndex1].logD.shape!=self.eosTables[fileIndex3].logD.shape:
      raise Exception("shape of \""+self.eosTables[fileIndex1].sFileName+"\" is "
        +str(self.eosTables[fileIndex1].logD.shape)+" while shape of \""\
        +self.eosTables[fileIndex3].sFileName+"\" is "+str(self.eosTables[fileIndex3].logD.shape)\
        +" shape of the two tables must be the same!")
    
    #check that the three eos's have the same density/temperature points
    for i in range(self.eosTables[fileIndex1].logD.shape[0]):
      for j in range(self.eosTables[fileIndex1].logD.shape[1]):
        
        #check densities
        if not np.isnan(self.eosTables[fileIndex1].logD[i][j])\
          and not np.isnan(self.eosTables[fileIndex2].logD[i][j])\
          and not np.isnan(self.eosTables[fileIndex3].logD[i][j]):
          if self.eosTables[fileIndex1].logD[i][j]!=self.eosTables[fileIndex2].logD[i][j]:
            #print self.eosTables[fileIndex1].logD[i][j],self.eosTables[fileIndex2].logD[i][j]
            raise Exception("\""+str(self.eosTables[fileIndex1].sFileName)
              +" has a different density point than \""+str(self.eosTables[fileIndex2].sFileName)+"\"")
          if self.eosTables[fileIndex1].logD[i][j]!=self.eosTables[fileIndex3].logD[i][j]:
            #print self.eosTables[fileIndex1].logD[i][j],self.eosTables[fileIndex3].logD[i][j]
            raise Exception("\""+str(self.eosTables[fileIndex1].sFileName)
              +" has a different density point than \""+str(self.eosTables[fileIndex3].sFileName)+"\"")
            
        #check tempeartures
        if not np.isnan(self.eosTables[fileIndex1].logT[i][j])\
          and not np.isnan(self.eosTables[fileIndex2].logT[i][j])\
          and not np.isnan(self.eosTables[fileIndex3].logT[i][j]):
          if self.eosTables[fileIndex1].logT[i][j]!=self.eosTables[fileIndex2].logT[i][j]:
            #print self.eosTables[fileIndex1].logT[i][j],self.eosTables[fileIndex2].logT[i][j]
            raise Exception("\""+str(self.eosTables[fileIndex1].sFileName)
              +" has a different temperature point than \""
              +str(self.eosTables[fileIndex2].sFileName)+"\"")
          if self.eosTables[fileIndex1].logT[i][j]!=self.eosTables[fileIndex3].logT[i][j]:
            #print self.eosTables[fileIndex1].logT[i][j],self.eosTables[fileIndex3].logT[i][j]
            raise Exception("\""+str(self.eosTables[fileIndex1].sFileName)
              +" has a different temperature point than \""
              +str(self.eosTables[fileIndex3].sFileName)+"\"")
    
    #next create a new equation of state object, and interpolate quadtatically between these 
    #3 equations of state in energy and pressure
    eosInterp=eosTable()
    eosInterp.Z=z
    eosInterp.X=self.eosTables[fileIndex1].X
    logE=[]
    logP=[]
    for i in range(self.eosTables[fileIndex1].logD.shape[0]):
      logERow=[]
      logPRow=[]
      for j in range(self.eosTables[fileIndex1].logD.shape[1]):
      
        #interpolate energy
        logETemp=self.__quad(z,self.eosTables[fileIndex1].Z,self.eosTables[fileIndex2].Z
          ,self.eosTables[fileIndex3].Z,self.eosTables[fileIndex1].logE[i][j]
          ,self.eosTables[fileIndex2].logE[i][j],self.eosTables[fileIndex3].logE[i][j])
        logERow.append(logETemp)
        
        #interpolate pressure
        logPTemp=self.__quad(z,self.eosTables[fileIndex1].Z,self.eosTables[fileIndex2].Z
          ,self.eosTables[fileIndex3].Z,self.eosTables[fileIndex1].logP[i][j]
          ,self.eosTables[fileIndex2].logP[i][j],self.eosTables[fileIndex3].logP[i][j])
        logPRow.append(logPTemp)
        
      #add rows to 2D arrays
      logE.append(logERow)
      logP.append(logPRow)
    
    #convert to numpy arrays and add
    eosInterp.logD=np.empty(self.eosTables[fileIndex1].logD.shape)
    eosInterp.logD[:]=float('nan')
    eosInterp.logT=np.empty(self.eosTables[fileIndex1].logD.shape)
    eosInterp.logT[:]=float('nan')
    eosInterp.logE=np.empty(self.eosTables[fileIndex1].logD.shape)
    eosInterp.logE[:]=float('nan')
    eosInterp.logP=np.empty(self.eosTables[fileIndex1].logD.shape)
    eosInterp.logP[:]=float('nan')
    for i in range(len(logE)):
      for j in range(len(logE[i])):
        eosInterp.logE[i][j]=logE[i][j]
        eosInterp.logP[i][j]=logP[i][j]
        eosInterp.logD[i][j]=self.eosTables[fileIndex1].logD[i][j]
        eosInterp.logT[i][j]=self.eosTables[fileIndex1].logT[i][j]
    
    #return equation of state interpolated to new Z
    return eosInterp
  def __quad(self,x,x1,x2,x3,y1,y2,y3):
    """Interpolates quadratically between 3 points, to the location of x, and returns the result"""
    
    x1Sq=x1*x1
    x2Sq=x2*x2
    x3Sq=x3*x3
    x1Sqmx2Sq=(x1Sq-x2Sq)
    a=((y1-y2)*(x2-x3)-(y2-y3)*(x1-x2))/(x1Sqmx2Sq*(x2-x3)-(x2Sq-x3Sq)*(x1-x2))
    b=(y1-y2-a*x1Sqmx2Sq)/(x1-x2)
    c=y1-a*x1Sq-b*x1
    return a*x*x+b*x+c
class interpTable:
  """This class reads in and holds data for an equations of state and opacities from a file formated
  in the same was as read to and written by the class defined in eos.h, and implemented in eos.cpp."""
  
  def interpolate(self,eosSet,opacitySet,withoutNans=None):
      """creates the interpolated table and writes it out"""
      
      #if not overriden, use value set in xml, or the default which is False
      if withoutNans==None:
        withoutNans=(not self.setNans)
      
      #interpolate in composition
      self.eosAtNewComp=eosSet.interpComp(self.X,self.Z)
      self.opacityAtNewComp=opacitySet.interpComp(self.X,self.Z)
      
      #interpolate in logD and logT
      self.eosTable=self.eosAtNewComp.interpolate(self.gridConfig,setExtrapolatedToNan=(not withoutNans))
      self.opacityTable=self.opacityAtNewComp.interpolate(self.gridConfig,setExtrapolatedToNan=(not withoutNans))
      
      #write out table
      print "  saving table to file \""+self.outputFile+"\" ..."
      self.__writeCompleteEOS()
      
      #plot final table for inspection
      if self.plot:
      
        print "plotting final tables, logP, logE, and logK in that order ..."
        self.eosTable.plotLogE()
        self.eosTable.plotLogP()
        self.opacityTable.plotLogK()
  def read(self,sFilename):
    """Reads in an interpolated table"""
    
    self.sFileName=sFilename
    
    f=open(sFilename,'rb')
    data=f.read()
    
    self.numLogR=-1
    sizeInt=struct.calcsize('i')
    sizeDouble=struct.calcsize('d')
    pos=0
    
    #get numLogD
    numLogD=struct.unpack('i',data[pos:pos+sizeInt])[0]
    pos+=sizeInt
    
    #get numLogT
    numLogT=struct.unpack('i',data[pos:pos+sizeInt])[0]
    pos+=sizeInt
    
    #get X
    self.X=struct.unpack('d',data[pos:pos+sizeDouble])[0]
    pos+=sizeDouble
    
    #get Y
    Y=struct.unpack('d',data[pos:pos+sizeDouble])[0]
    pos+=sizeDouble
    
    #calculate Z
    self.Z=1.0-self.X-Y
    
    #get minLogD
    minLogD=struct.unpack('d',data[pos:pos+sizeDouble])[0]
    pos+=sizeDouble
    
    #get delLogD
    delLogD=struct.unpack('d',data[pos:pos+sizeDouble])[0]
    pos+=sizeDouble
    
    #get minLogT
    minLogT=struct.unpack('d',data[pos:pos+sizeDouble])[0]
    pos+=sizeDouble
    
    #get delLogT
    delLogT=struct.unpack('d',data[pos:pos+sizeDouble])[0]
    pos+=sizeDouble
    
    #put everything into grid config
    self.gridConfig=[minLogD,delLogD,numLogD,minLogT,delLogT,numLogT]
    
    #create arrays to store data
    self.logD=np.empty((numLogD,numLogT))
    self.logT=np.empty((numLogD,numLogT))
    self.logP=np.empty((numLogD,numLogT))
    self.logE=np.empty((numLogD,numLogT))
    self.logK=np.empty((numLogD,numLogT))
    for i in range(numLogD):
      
      #get logP at all temps at this density
      fmt=str(numLogT)+'d'
      logPSet=struct.unpack(fmt,data[pos:pos+numLogT*sizeDouble])
      pos+=numLogT*sizeDouble
      
      #get logE at all temps at this density
      logESet=struct.unpack(fmt,data[pos:pos+numLogT*sizeDouble])
      pos+=numLogT*sizeDouble
      
      #get logK at all temps at this density
      logKSet=struct.unpack(fmt,data[pos:pos+numLogT*sizeDouble])
      pos+=numLogT*sizeDouble
      
      for j in range(numLogT):
        
        #set logD and logT
        self.logD[i][j]=minLogD+float(i)*delLogD
        self.logT[i][j]=minLogT+float(j)*delLogT
        
        #read in pressure
        self.logP[i][j]=logPSet[j]
        
        #read in energy
        self.logE[i][j]=logESet[j]
        
        #read in opacity
        self.logK[i][j]=logKSet[j]
        
    f.close()
  def plotLogE(self,otherTables=None,logDIndexList=None,logDRangeList=None,wireFrame=True,rstride=1
    ,cstride=1,outputfile=None):
    """Plots LogE
    
    Keywords:
    otherTables: a list of other eosTables to include in the plot
    logDIndexList: a list of integers corresponding to which densities to plot the tables at
    wireFrame: if set to true (the default) and logDIndexList is set to None it will plot a 3D
      wireframe of logE.
    """
    
    if logDIndexList!=None:
      wireFrame=False
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib
    if outputfile!=None:
      matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig=plt.figure()
    lines=[]
    lables=[]
    if wireFrame:
      ax=fig.add_subplot(111,projection='3d')
      ax.plot_wireframe(self.logD,self.logT,self.logE,rstride=rstride,cstride=cstride)
      if otherTables:
        for otherTable in otherTables:
          ax.plot_wireframe(otherTable.logD,otherTable.logT,otherTable.logE,color="green")
    else:
      ax=fig.add_subplot(111)
      l=logDIndexList[0]
      h=logDIndexList[0]+logDRangeList[0]
      print self.sFileName," logD=",self.logD[l:h,0]
      if logDRangeList[0]>1:
        temp=ax.plot(np.transpose(self.logT[l:h,:]),np.transpose(self.logE[l:h,:]), "b-")
      else:
        temp=ax.plot(self.logT[l:h,:][0],self.logE[l:h,:][0], "b-")
      lines.append(temp)
      lables.append(self.sFileName)
      counter=1
      if otherTables:
        for otherTable in otherTables:
          l=logDIndexList[counter]
          h=logDIndexList[counter]+logDRangeList[counter]
          print otherTable.sFileName," logD=",otherTable.logD[l:h,0]
          if logDRangeList[counter]>1:
            temp=ax.plot(np.transpose(otherTable.logT[l:h,:]),np.transpose(otherTable.logE[l:h,:]), "g-")
          else:
            temp=ax.plot(otherTable.logT[l:h,:][0],otherTable.logE[l:h,:][0], "g-")
          counter+=1
          lines.append(temp)
          lables.append(otherTable.sFileName)
    
    #make ledgend
    #if len(lines)>0:
      #ax.legend(lines,lables)
    
    #show results
    if outputfile==None:
      plt.show()
    else:
      [path,ext]=os.path.splitext(outputfile)
      supportedFileTypes=["png", "pdf", "ps", "eps", "svg"]
      if ext[1:] not in supportedFileTypes:
        print "File type \""+ext[1:]+"\" not suported. Supported types are ",supportedFileTypes," please choose one of those"
        quit()
      print __name__+":"+main.__name__+": saving figure to file \""+outputfile
      fig.savefig(path+"_E"+ext,format=ext[1:],transparent=False)#save to file
  def plotLogP(self,otherTables=None,logDIndexList=None,logDRangeList=None,wireFrame=True
    ,outputfile=None):
    """Plots LogP
    
    Keywords:
    otherTables: a list of other eosTables to include in the plot
    logDIndexList: a list of integers corresponding to which densities to plot the tables at
    wireFrame: if set to true (the default) and logDIndexList is set to None it will plot a 3D
      wireframe of logP.
    """
    
    if logDIndexList!=None:
      wireFrame=False
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib
    if outputfile!=None:
      matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig=plt.figure()
    lines=[]
    lables=[]
    if wireFrame:
      ax=fig.add_subplot(111,projection='3d')
      ax.plot_wireframe(self.logD,self.logT,self.logP)
      if otherTables:
        for otherTable in otherTables:
          ax.plot_wireframe(otherTable.logD,otherTable.logT,otherTable.logP,color="green")
    else:
      ax=fig.add_subplot(111)
      l=logDIndexList[0]
      h=logDIndexList[0]+logDRangeList[0]
      print self.sFileName," logD=",self.logD[l:h,0]
      if logDRangeList[0]>1:
        temp=ax.plot(np.transpose(self.logT[l:h,:]),np.transpose(self.logP[l:h,:]), "b-")
      else:
        temp=ax.plot(self.logT[l:h,:][0],self.logP[l:h,:][0], "b-")
      lines.append(temp)
      lables.append(self.sFileName)
      counter=1
      if otherTables:
        for otherTable in otherTables:
          l=logDIndexList[counter]
          h=logDIndexList[counter]+logDRangeList[counter]
          print otherTable.sFileName," logD=",otherTable.logD[l:h,0]
          if logDRangeList[counter]>1:
            temp=ax.plot(np.transpose(otherTable.logT[l:h,:]),np.transpose(otherTable.logP[l:h,:]), "g-")
          else:
            temp=ax.plot(otherTable.logT[l:h,:][0],otherTable.logP[l:h,:][0], "g-")
          counter+=1
          lines.append(temp)
          lables.append(otherTable.sFileName)
    #if len(lines)>0:
      #ax.legend(lines,lables)
    
    #show results
    if outputfile==None:
      plt.show()
    else:
      [path,ext]=os.path.splitext(outputfile)
      supportedFileTypes=["png", "pdf", "ps", "eps", "svg"]
      if ext[1:] not in supportedFileTypes:
        print "File type \""+ext[1:]+"\" not suported. Supported types are ",supportedFileTypes," please choose one of those"
        quit()
      print __name__+":"+main.__name__+": saving figure to file \""+outputfile
      fig.savefig(path+"_P"+ext,format=ext[1:],transparent=False)#save to file
  def plotLogK(self,otherTables=None,logDIndexList=None,logDRangeList=None,wireFrame=True
    ,outputfile=None):
    """Plots opacity
    
    Keywords:
    otherTables: a list of opacity tables to also be ploted
    logDIndex: a list of integers used to indicate a specific logR index to plot 2D line plots at.
    """
    
    if logDIndexList!=None:
      wireFrame=False
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib
    if outputfile!=None:
      matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig=plt.figure()
    lines=[]
    lables=[]
    if wireFrame:
      ax=fig.add_subplot(111,projection='3d')
      ax.plot_wireframe(self.logD,self.logT,self.logK)
      if otherTables:
        for otherTable in otherTables:
          ax.plot_wireframe(otherTable.logD,otherTable.logT,otherTable.logK,color="green")
    else:
      ax=fig.add_subplot(111)
      l=logDIndexList[0]
      h=logDIndexList[0]+logDRangeList[0]
      print self.sFileName," logD=",self.logD[l:h,0],logDRangeList[0],l,h
      if logDRangeList[0]>1:
        #fig.suptitle("logD=",str(self.logD[l:h,0]))
        temp=ax.plot(np.transpose(self.logT[l:h,:]),np.transpose(self.logK[l:h,:]), "b-")
      else:
        temp=ax.plot(self.logT[l:h,:][0],self.logK[l:h,:][0], "b-")
      lines.append(temp)
      lables.append(self.sFileName)
      counter=1
      if otherTables:
        for otherTable in otherTables:
          l=logDIndexList[counter]
          h=logDIndexList[counter]+logDRangeList[counter]
          print otherTable.sFileName," logD=",otherTable.logD[l:h,0],logDRangeList[counter],l,h
          if logDRangeList[counter]>1:
            temp=ax.plot(np.transpose(otherTable.logT[l:h,:]),np.transpose(otherTable.logK[l:h,:]), "g-",)
          else:
            temp=ax.plot(otherTable.logT[l:h,:][0],otherTable.logK[l:h,:][0], "g-")
          counter+=1
          lines.append(temp)
          lables.append(otherTable.sFileName)
    #if len(lines)>0:
      #ax.legend(lines,lables)
    
    if outputfile==None:
      plt.show()
    else:
      [path,ext]=os.path.splitext(outputfile)
      supportedFileTypes=["png", "pdf", "ps", "eps", "svg"]
      if ext[1:] not in supportedFileTypes:
        print "File type \""+ext[1:]+"\" not suported. Supported types are ",supportedFileTypes," please choose one of those"
        quit()
      print __name__+":"+main.__name__+": saving figure to file \""+outputfile
      fig.savefig(path+"_kappa"+ext,format=ext[1:],transparent=False)#save to file
  def __init__(self,tableElement=None):
    """Reads in an interpolation table info from from the xml element tableElement."""
    
    if tableElement!=None:
      
      #get hydrogen mass fraction
      element=tableElement.findall("X")
      if len(element)>1:
        print "WARNING: more than one \"X\" found using only first node found"
      elif len(element)==0:
        raise Exception("no \"X\" node found")
      try:
        self.X=float(element[0].text)
      except (TypeError, ValueError):
        raise Exception("\"X\" node expected to be a float got \""+str(element[0].text)+"\"")
        
      #get metal mass fraction
      element=tableElement.findall("Z")
      if len(element)>1:
        print "WARNING: more than one \"Z\" found using only first node found"
      elif len(element)==0:
        raise Exception("no \"Z\" node found")
      try:
        self.Z=float(element[0].text)
      except (TypeError, ValueError):
        raise Exception("\"Z\" node expected to be a float got \""+str(element[0].text)+"\"")
      
      #get grid configuration of new table
      self.gridConfig=[]
      
      #get minLogD
      element=tableElement.findall("minLogD")
      if len(element)>1:
        print "WARNING: more than one \"minLogD\" found using only first node found"
      elif len(element)==0:
        raise Exception("no \"minLogD\" node found")
      try:
        self.gridConfig.append(float(element[0].text))
      except (TypeError, ValueError):
        raise Exception("\"minLogD\" node expected to be a float got \""+str(element[0].text)+"\"")
      
      #get delLogD
      element=tableElement.findall("delLogD")
      if len(element)>1:
        print "WARNING: more than one \"delLogD\" found using only first node found"
      elif len(element)==0:
        raise Exception("no \"delLogD\" node found")
      try:
        self.gridConfig.append(float(element[0].text))
      except (TypeError, ValueError):
        raise Exception("\"delLogD\" node expected to be a float got \""+str(element[0].text)+"\"")
      
      #get numLogD
      element=tableElement.findall("numLogD")
      if len(element)>1:
        print "WARNING: more than one \"numLogD\" found using only first node found"
      elif len(element)==0:
        raise Exception("no \"numLogD\" node found")
      try:
        self.gridConfig.append(int(element[0].text))
      except (TypeError, ValueError):
        raise Exception("\"numLogD\" node expected to be an int got \""+str(element[0].text)+"\"")
      
      #get minLogT
      element=tableElement.findall("minLogT")
      if len(element)>1:
        print "WARNING: more than one \"minLogT\" found using only first node found"
      elif len(element)==0:
        raise Exception("no \"minLogT\" node found")
      try:
        self.gridConfig.append(float(element[0].text))
      except (TypeError, ValueError):
        raise Exception("\"minLogT\" node expected to be a float got \""+str(element[0].text)+"\"")
      
      #get delLogT
      element=tableElement.findall("delLogT")
      if len(element)>1:
        print "WARNING: more than one \"delLogT\" found using only first node found"
      elif len(element)==0:
        raise Exception("no \"delLogT\" node found")
      try:
        self.gridConfig.append(float(element[0].text))
      except (TypeError, ValueError):
        raise Exception("\"delLogT\" node expected to be a float got \""+str(element[0].text)+"\"")
      
      #get numLogT
      element=tableElement.findall("numLogT")
      if len(element)>1:
        print "WARNING: more than one \"numLogT\" found using only first node found"
      elif len(element)==0:
        raise Exception("no \"numLogT\" node found")
      try:
        self.gridConfig.append(int(element[0].text))
      except (TypeError, ValueError):
        raise Exception("\"numLogT\" node expected to be an int got \""+str(element[0].text)+"\"")
      
      #get outputFile
      element=tableElement.findall("outputFile")
      if len(element)>1:
        print "WARNING: more than one \"outputFile\" found using only first node found"
      elif len(element)==0:
        raise Exception("no \"outputFile\" node found")
      self.outputFile=element[0].text
      
      #get wether to plot, if not set, set to False
      element=tableElement.findall("plot")
      if len(element)>1:
        print "WARNING: more than one \"plot\" found using only first node found"
      elif len(element)==0:
        self.plot=False
      if element[0].text.lower() in ['true','1','t','y','yes']:
        self.plot=True
      else:
        self.plot=False
      
      #get wether to set nans, if not set, set to False
      element=tableElement.findall("setNans")
      if len(element)>1:
        print "WARNING: more than one \"setNans\" found using only first node found"
      if len(element)==0:
        self.setNans=False
      elif element[0].text.lower() in ['true','1','t','y','yes']:
        self.setNans=True
      else:
        self.setNans=False
  def __writeCompleteEOS(self):
    """Writes out an eosTable and and an opacity Table in a form that SPHERLS can read in
    """
    
    #check that eosTable and opacityTable are compatible
    if self.eosTable.X!=self.opacityTable.X:
      raise Exception("eosTable has X="+str(self.eosTable.X)+" opacityTable has X="+str(self.opacityTable.X)\
        +" must be the same!")
    elif self.eosTable.Z!=self.opacityTable.Z:
      raise Exception("eosTable has Z="+str(self.eosTable.Z)+" opacityTable has Z="+str(self.opacityTable.Z)\
        +" must be the same!")
    elif self.eosTable.logD.shape!=self.opacityTable.logD.shape:
      raise Exception("eosTable has shape="+str(self.eosTable.shape)+" opacityTable has shape="\
        +str(self.opacityTable.shape)+" must be the same!")
    elif self.eosTable.logDMin!=self.opacityTable.logDMin:
      raise Exception("eosTable has logDMin="+str(self.eosTable.logDMin)+" opacityTable has logDMin="\
        +str(self.opacityTable.logDMin)+" must be the same!")
    elif self.eosTable.logDDel!=self.opacityTable.logDDel:
      raise Exception("eosTable has logDDel="+str(self.eosTable.logDDel)+" opacityTable has logDDel="\
        +str(self.opacityTable.logDDel)+" must be the same!")
    elif self.eosTable.logTMin!=self.opacityTable.logTMin:
      raise Exception("eosTable has logTMin="+str(self.eosTable.logTMin)+" opacityTable has logTMin="\
        +str(self.opacityTable.logTMin)+" must be the same!")
    elif self.eosTable.logTDel!=self.opacityTable.logTDel:
      raise Exception("eosTable has logTDel="+str(self.eosTable.logTDel)+" opacityTable has logTDel="\
        +str(self.opacityTable.logTDel)+" must be the same!")
    
    f=open(self.outputFile,'wb')
    
    #write out header info
    Y=1.0-self.eosTable.X-self.eosTable.Z
    data=struct.pack("2i6d",self.eosTable.logD.shape[0],self.eosTable.logD.shape[1],self.eosTable.X,Y,\
      self.eosTable.logDMin,self.eosTable.logDDel,self.eosTable.logTMin,self.eosTable.logTDel)
    f.write(data)
    
    for i in range(self.eosTable.logD.shape[0]):
      for j in range(self.eosTable.logD.shape[1]):
        data=struct.pack("d",self.eosTable.logP[i][j])
        f.write(data)
      for j in range(self.eosTable.logD.shape[1]):
        data=struct.pack("d",self.eosTable.logE[i][j])
        f.write(data)
      for j in range(self.eosTable.logD.shape[1]):
        data=struct.pack("d",self.opacityTable.logK[i][j])
        f.write(data)
      
    f.close()
class interpTableManager:
  def createTables(self,withoutNans=None):
    """Creates interpolated tables and write them out."""
    
    for interpTableTemp in self.tables:
      print "creating table \""+interpTableTemp.outputFile+"\" ..."
      interpTableTemp.interpolate(self.eosSet,self.opacitySet,withoutNans)
  def __init__(self,configFile=None):
    """Initializes interpTableManager from the given configuration file."""
    
    self.configFile=configFile
    if configFile!=None:
      
      #process equation of state
      self.eosSet=eosTableManager(configFile)
      self.eosSet.load()
      
      #process opacity
      self.opacitySet=opacityTableManager(configFile)
      self.opacitySet.load()
      
      #get infos about table to create
      self.__readInterpTableConfigs()
  def __readInterpTableConfigs(self):
    """Finds table elements from the self.confFile xml file and passes the xml nodes to """
    
    tree=xml.parse(self.configFile)
    root=tree.getroot()
    
    #get interpolatedTables node
    interpolatedTables=root.findall("interpolatedTables")
    if len(interpolatedTables)>1:
      print "WARNING: more than one \"interpolatedTables\" found, using only first node found"
    elif len(interpolatedTables)==0:
      raise Exception("no \"interpolatedTables\" node found")
    
    #get all tables nodes
    tables=interpolatedTables[0].findall("table")
    if len(tables)==0:
      raise Exception("no \"tables\" node found under \"interpolatedTables\" node.")
    
    self.tables=[]
    for tableElement in tables:
      self.tables.append(interpTable(tableElement))
def createTable(args):
  if len(args)==1:#could add a test to make sure it is an xml file
    interp=interpTableManager(args[0])
    #interp.createTables(withoutNans=True)
    interp.createTables()
def compareTables(args,options):
  if len(args)==2:
    table1=interpTable()
    table1.read(args[0])
    
    table2=interpTable()
    table2.read(args[1])
    rhoIndex1=options.rho1
    rhoIndex2=options.rho2
    numRho=options.numRho
    
    table1.plotLogP([table2],[rhoIndex1,rhoIndex2],[numRho,numRho],outputfile=options.output)
    table1.plotLogE([table2],[rhoIndex1,rhoIndex2],[numRho,numRho],outputfile=options.output)
    table1.plotLogK([table2],[rhoIndex1,rhoIndex2],[numRho,numRho],outputfile=options.output)
def printTablePart(args):
  #could add a test to see if it isn't an xml file, then assume it is an equation of state
  table=interpTable()
  table.read(args[0])
  for i in range(table.logD.shape[0]):#loop over density
    print table.logD[i][0]
    if table.logD[i][0]<=-7 and table.logD[i][0]>=-9:
      f=open(table.sFileName+"_opacity_LogD"+str(table.logD[i][0])+".txt",'w')
      for j in range(table.logD.shape[1]):#loopover temperature
        line=str(table.logT[i][j])+" "+str(table.logK[i][j])+"\n"
        f.write(line)
      f.close()
def main():
  
  #parse command line options
  (options,args)=parseOptions()
  
  #assume reading a configuration file to make a new table
  createTable(args)
    
  #assume comparing two pre-existing tables
  compareTables(args,options)
  
  #print out part of a table
  #printTablePart(args)
if __name__ == "__main__":
  main()