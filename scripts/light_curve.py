#!/usr/bin/env python

import datafile
import optparse as op
import make_profiles
import glob
import numpy as np
import sys
import os
import disect_filename
import parser
from math import *
import xml.etree.ElementTree as xml
import parse_formula
import mywarnings
import warnings
import paths

luminosityIndex=59
tempIndex=49
massIndex=1
gridUIndex=20
radiusIndex=4
GravConst=6.67259e-08
M_bol_sun=4.75
withAccelInGrav=False

def parseOptions():
  
  #setup command line parser
  '''This is out of date, needs to be updated to reflect new additions'''
  parser=op.OptionParser(usage="Usage: %prog [options] XMLFILE"
    ,version="%prog 1.0",description=r"Reads in the XMLFILE which specifies "
    +"how the light curve is to be made. see light_curve_reference.xml for a "
    +"template xml file to be used with this script. Found in the SPHERLS "
    +"directory docs/templateXML.")
  parser.add_option("--write-diagnostics",dest="writeDiagnostics"
    ,action="store_true"
    ,help="If set it will write out Logg,g,accel,T_eff,BC to the output file in"
    +" that order in addition to the time and mag in the first and second "
    +"columns",default=False)
  parser.add_option("-d",dest="double"
    ,action="store_true"
    ,help="If set it will double the data in the light curve for easier "
    +"viewing of the light curve shape.",default=False)
  
  make_profiles.addParserOptions(parser)
  
  #parse command line options
  return parser.parse_args()
class LightCurve:
  def __init__(self,element):
    
    #get file of bolo metric corrections
    boloCorrElements=element.findall("boloCorr")
    if len(boloCorrElements)>1:
      warnings.warn("more than one \"boloCorr\" node in \"lightCurve\" node"
        +", ignoring all but first node")
    if len(boloCorrElements)>0:
      self.boloCorrFile=boloCorrElements[0].text
    else:
      raise Exception("Must have one \"boloCorr\" input file per \"lightCurve\" node")
    
    #get column for bolometric correction
    try:
      self.columnBC=int(boloCorrElements[0].get("columnBC"))
    except Exception as e:
      raise Exception(e.message+"\n Must have a \"columnBC\" attribute set to the column number for"
        +" the bolometric correction. The first column in the file is column 0")
    
    #get column for bolometric correction
    if boloCorrElements[0].get("withAcceleration").lower() in ["true","yes","y","1","t"]:
      self.withAcceleration=True
    else:
      self.withAcceleration=False
    
    #get input set
    inputFilesElements=element.findall("inputFiles")
    if len(inputFilesElements)>1:
      warnings.warn("more than one \"inputFiles\" node in \"lightCurve\" node"
        +", ignoring all but first node")
    if len(inputFilesElements)>0:
      self.inputFileRange=inputFilesElements[0].text
    else:
      raise Exception("Must have one \"inputFiles\" node per \"lightCurve\" node")
    
    #get frequency
    if inputFilesElements[0].get("frequency")!=None:
      try:
        self.frequency=int(inputFilesElements[0].get("frequency"))
      except ValueError as e:
        raise Exception(e.message+"\n \"frequency\" attribute in \"inputFiles\" node must be an "
          +"integer")
    else:
      self.frequency=1
    
    #split up file set string
    [self.start,self.end,self.baseFileName]=disect_filename.disectFileName(self.inputFileRange)
    
    #get equation of state file if one specified
    eosFileElement=element.findall("eosFile")
    if len(eosFileElement)>1:
      warnings.warn("more than one \"eosFile\" node, ignoring all but first node")
    if len(eosFileElement)>0:
      self.eosFile=eosFileElement[0].text
    else:
      self.eosFile=None
    
    #get get number of zones in from the surface, 0=surface
    zonesFromSurfElement=element.findall("zonesFromSurf")
    if len(zonesFromSurfElement)>1:
      warnings.warn("more than one \"zonesFromSurf\" node, ignoring all but first node")
    if len(zonesFromSurfElement)>0:
      try:
        self.zonesFromSurf=int(zonesFromSurfElement[0].text)
      except ValueError as e:
        raise Exception(e.message+"\nExpecting an integer for \"zonesFromSurf\" node")
    else:
      self.zonesFromSurf=0
    
    #get outputFile 
    outputFileElements=element.findall("outputFile")
    if len(outputFileElements)>1:
      warnings.warn("more than one \"outputFile\" node in \"lightCurve\" node"
        +", ignoring all but first node")
    if len(outputFileElements)>0:
      self.outputFile=outputFileElements[0].text
    else:
      raise Exception("Must have one \"outputFile\" node per \"lightCurve\" node")
      
    #get period
    periodElements=element.findall("period")
    if len(periodElements)>1:
      warnings.warn("more than one \"period\" node in \"lightCurve\" node"
        +", ignoring all but first node")
    if len(periodElements)>0:
      self.period=float(periodElements[0].text)
    else:
      self.period=None
  def create(self,options):
    
    #read in bolometric corrections
    self.readBoloCorr()
    
    #read in needed data from radial profile set
    self.readProfiles(options)
    
    #calculate light curve
    curve=self.calculateCurve()
    
    #if a period is given phase the light curve
    if self.period!=None:
      curve=self.phase(curve,double=options.double)
    
    #write out curve
    self.write(curve,writeDiagnostics=options.writeDiagnostics)
  def readProfiles(self,options):
    '''Reads the needed data to create the light curve from the radial profile files'''
    
    #make sure needed radial profiles are made
    fileName=self.baseFileName+"["+str(self.start)+"-"+str(self.end)+"]"
    failedFiles=make_profiles.make_fileSet(fileName,options)#default is to make 
                                                            #a profile set
    
    #let user know if there were any failed files
    if len(failedFiles)>0:
      for failedFile in failedFiles:
        print failedFile
        
    #get and sort profiles within range of dataset
    extension="_pro"+".txt"
    filesExistProfiles=glob.glob(self.baseFileName+"*"+extension)
    filesExistProfiles.sort()
    files=[]
    for i in range(0,len(filesExistProfiles),self.frequency):
      file=filesExistProfiles[i]
      intOfFile=int(file[len(self.baseFileName):len(file)-len(extension)])
      if intOfFile>=self.start and intOfFile<self.end:
        files.append(file)
    if len(files)==0:
      raise Exception("no files found in range")
    
    self.nNumFiles=len(files)
    
    #get time, Luminosity, Teff, M_r and grid velocity
    self.luminosity=[]
    self.temperature=[]
    self.interiorMass=[]
    self.time=[]
    self.gridVelocity=[]
    self.radius=[]
    count=1
    
    for file in files:
      
      #read data file
      print "reading in \""+file+"\" "+str(count)+"/"+str(self.nNumFiles)+" ..."
      fileData=datafile.DataFile()
      fileData.readFile(file)
      
      #get time
      #these are required in case the time isn't space seperated
      fileHeader=fileData.sHeader.split("=")
      indexBracket=fileHeader[1].find("[s]")
      self.time.append(float(fileHeader[1][0:indexBracket]))
      
      #get index for surface zone, minus zoneFromSurf
      rowIndex=len(fileData.fColumnValues)-1-self.zonesFromSurf
      
      #get luminosity
      self.luminosity.append(fileData.fColumnValues[rowIndex][luminosityIndex])
      
      #get Teff, of outer interface
      temp=(fileData.fColumnValues[rowIndex-1][tempIndex])*(2**(0.25))
      self.temperature.append(temp)
      
      #get M_r
      self.interiorMass.append(fileData.fColumnValues[rowIndex][massIndex])
      
      #get radius
      self.radius.append(fileData.fColumnValues[rowIndex][radiusIndex])
      
      #get grid velocity
      self.gridVelocity.append(fileData.fColumnValues[rowIndex][gridUIndex])
      
      #increment file counter
      count=count+1
  def readBoloCorr(self):
    '''Reads in the bolometric correction table'''
    
    #assume path is relative to atmospheres directory
    if (self.boloCorrFile[0:2] in ["./","/"]) or (self.boloCorrFile[0] in ["./","/"]):
      boloCorrFile=self.boloCorrFile
    else:
      boloCorrFile=os.path.join(paths.SPHERLSDataPath,"atmospheres/"+self.boloCorrFile)
    
    #figure out zoning, assuming first column is teff, and second collumn is logg
    f=open(boloCorrFile)
    self.TMin=1e300
    TMax=-1e300
    self.loggMin=1e300
    loggMax=-1e300
    self.TDel=None
    self.loggDel=None
    lineNum=1
    for line in f:
      columns=line.split()
      tempT=float(columns[0])
      tempLogg=float(columns[1])
      
      #set old values equal to current values for first itteration only
      if lineNum==1:
        tempTOld=tempT
        tempLoggOld=tempLogg
      
      #save smallest T
      if tempT<self.TMin:
        self.TMin=tempT
      
      #save largest T
      if tempT>TMax:
        TMax=tempT
      
      #save smallest logg
      if tempLogg<self.loggMin:
        self.loggMin=tempLogg
      
      #save largest logg
      if tempLogg>loggMax:
        loggMax=tempLogg
      
      #check to see if T changed at the same logg
      if tempTOld!=tempT and tempLoggOld==tempLogg:
        if self.TDel==None:#set delta first time temperatures change
          self.TDel=abs(tempT-tempTOld)
        else:#check that the change is the same
          TDelTest=abs(tempT-tempTOld)
          if TDelTest!=self.TDel:
            raise Exception("not a rectangular table, Del T="+str(self.TDel)+", but T changed by "
              +str(TDelTest)+" in one step "+str(lineNum))
      
      #check to see if logg changed at the same T
      if tempLoggOld!=tempLogg:
        if self.loggDel==None:#set delta first time temperatures change 
          self.loggDel=abs(tempLogg-tempLoggOld)
        else:#check that the change is the same
          loggDelTest=abs(tempLogg-tempLoggOld)
          if loggDelTest!=self.loggDel:
            raise Exception("not a rectangular table, Del logg="+str(self.loggDel)
              +", but logg changed by "+str(loggDelTest)+" in one step on line "+str(lineNum))
      
      #save old values
      tempTOld=tempT
      tempLoggOld=tempLogg
      lineNum=lineNum+1
    self.numLogg=int((loggMax-self.loggMin)/self.loggDel)+1
    self.numT=int((TMax-self.TMin)/self.TDel)+1
    #print self.loggMin,loggMax,self.loggDel,self.numLogg, self.TMin,TMax,self.TDel,self.numT
    
    #read in bolometric corrections
    f.close()
    f=open(boloCorrFile)
    self.BC=np.empty([self.numLogg,self.numT],float)
    self.BC.fill(np.nan)
    for line in f:
      
      #split up line and get temperature and logg
      columns=line.split()
      tempT=float(columns[0])
      tempLogg=float(columns[1])
      
      #figure out array indices based on T and logg
      j=int((tempT-self.TMin)/self.TDel)
      i=int((tempLogg-self.loggMin)/self.loggDel)
      
      #put BC in the right array location
      self.BC[i][j]=float(columns[self.columnBC])
    
    #print self.BC
    f.close()
  def calculateCurve(self):
    '''Creates the light curve by converting luminosity to bolometric magnitude and then appling a
    bolometric correction and returns a 2D list of times and light curve magnitudes.'''
    
    #loop over time
    curve=[]
    for n in range(len(self.time)-1):
      
      #calculate logg
      radSq=self.radius[n]*self.radius[n]
      accel=(self.gridVelocity[n+1]-self.gridVelocity[n])/(self.time[n+1]-self.time[n])
      g=GravConst*self.interiorMass[n]/(radSq)
      T=self.temperature[n]
      if self.withAcceleration:
        logg=log10(g+accel)
      else:
        logg=log10(g)
      
      #calculate indices at lower logg and lower T
      j=int((T-self.TMin)/self.TDel)
      i=int((logg-self.loggMin)/self.loggDel)
      
      #check that we are within table
      loggMax=self.loggMin+float(self.numLogg-1)*self.loggDel
      TMax=self.TMin+float(self.numT-1)*self.TDel
      if i>self.numLogg-2:
        raise Exception("logg="+str(logg)+" is out side the table with largest logg="+str(loggMax))
      if j>self.numT-2:
        raise Exception("T="+str(T)+" is out side the table with largest T="+str(TMax))
      
      #calculate bracketing T's and fraction of the distance between them
      Tj=self.TMin+float(j)*self.TDel
      Tjp1=self.TMin+float(j+1)*self.TDel
      Tfrac=(T-Tj)/(Tjp1-Tj)
      
      #calculate bracketing logg's and the fraction of the distance between them
      loggi=self.loggMin+float(i)*self.loggDel
      loggip1=self.loggMin+float(i+1)*self.loggDel
      loggFrac=(logg-loggi)/(loggip1-loggi)
      
      #linearly interpolate BC in T at lower logg
      BCI=(self.BC[i][j+1]-self.BC[i][j])*Tfrac+self.BC[i][j]
      
      #linearly interpolate BC in T at upper logg
      BCIP1=(self.BC[i+1][j+1]-self.BC[i+1][j])*Tfrac+self.BC[i+1][j]
      
      #linearly interpolate BC in logg
      BC=(BCIP1-BCI)*loggFrac+BCI
      if np.isnan(BC):
        raise Exception("BC became a nan when interpolated to logg="+str(logg)+" and T="+str(T)
          +" this indicates that it was inside the rectangular table, but used a (logg,T)"
          +" combination that doesn't have a BC defined")
      
      #compute the magnitude
      mag=-2.5*log10(self.luminosity[n])+M_bol_sun-BC
      
      #print self.time[n], self.luminosity[n],self.temperature[n]\
      #,self.interiorMass[n],self.gridVelocity[n],g,accel,logg,BC,self.BC[i][j]
      
      temp=[self.time[n],mag,logg,g,accel,T,BC]
      curve.append(temp)
    return curve
  def phase(self,curve,double=False):
    """If user specified a period convert the time to a phase, and fold it
    
    Parameters:
      curve: the light curve to phase
    Keywords:
      double: by default it is False, but if set to true it will double the data
        to provide two phases of the light curve to make it easier to see the
        shape
    """
    
    #insure that we have a period, otherwise there is nothing to do
    if self.period!=None:
      
      for i in range(len(curve)):
        
        #force phase to be between 0-1
        phase=(curve[i][0]-int(curve[i][0]/self.period)*self.period)/self.period
        curve[i][0]=phase
      
      #shift max light to phase 0
      curveTmp=np.array(curve)
      indexOfMax=curveTmp[:,1].argmin()
      phaseOfMax=curveTmp[indexOfMax][0]
      for i in range(len(curve)):
        phase=curve[i][0]-phaseOfMax
        if phase<0.0:
          curve[i][0]=phase+1
        else:
          curve[i][0]=phase
      
      #double the phase
      for i in range(len(curve)):
        curve.append( [curve[i][0]+1.0,curve[i][1]] )
    else:
      raise Exception("No period given unable to phase light curve")
    return curve
  def write(self,curve,writeDiagnostics=False):
    """Writes out the light curve to the specified output file.
    
    Parameters:
      curve: a list containing the time, magnitude, logg, g,accel,T for a 
      particular light curve. The logg,g,accel, T are only needed if these
      quantities are to be written out. But the order must remain, e.g. if T is 
      to be printed out, but the other quantities listed before T aren't present
      it will cause problems.
    Keywords:
      writeDiagnostics: if true will write out the additional information of:
        logg: log of the effective gravtity with or without pulsational
          acceleration included
        g: the gravintational acceleration not including the pulsational 
          acceleration
        accel: the acceleration due to the pulsation just d(U_r0)/dt
        T_eff: the effective temaperature
        BC: the calculated bolometric correction from T_eff, and logg
    """
    
    print "writting light curve to file \""+self.outputFile+"\" ..."
    
    #open file
    f=open(self.outputFile,'w')
    
    #write out light curve
    for i in range(len(curve)):
      line=str(curve[i][0])+" "+str(curve[i][1])
      if writeDiagnostics:
        line+=" "+str(curve[i][2])
        line+=" "+str(curve[i][3])
        line+=" "+str(curve[i][4])
        line+=" "+str(curve[i][5])
        line+=" "+str(curve[i][6])
      line+="\n"
      f.write(line)
    f.close()
def main():
  
  #parse command line options
  (options,args)=parseOptions()
  
  if len(args)==0:
    print "must supply an input file"
    quit()
  
  #get root element
  tree=xml.parse(args[0])
  root=tree.getroot()
  
  #get all light curve initialization data
  lightCurveElements=root.findall("lightCurve")
  lightCurves=[]
  for lightCurveElement in lightCurveElements:
    lightCurves.append(LightCurve(lightCurveElement))
  
  #create the light curves
  for lightCurve in lightCurves:
    lightCurve.create(options)
  
  return
if __name__ == "__main__":
  main()