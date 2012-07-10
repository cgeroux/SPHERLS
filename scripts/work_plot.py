#!/usr/bin/env python

import datafile
import optparse as op
import make_profiles
import glob
import math
import numpy as np
import sys
import os
import disect_filename
import average_PKE
import xml.etree.ElementTree as xml
import warnings
import mywarnings#imports some custom warning output formating
import xmlParseFunctions#some commonly used xmlParsing functions
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec
supportedFormats=["png","pdf"]

class WorkPlotSettings:
  def __init__(self):
    self.minTemp=0.0
    self.ylim=[None,None]
    self.grid=True#will put grid on work plots
    self.points=False#will point points on work plots
    self.lines=True#if true will plot lines on work plots
    self.plotPdVCurves=False#if true will plot PdV curves
    self.startZone=True#zone to start create PdV curves
    self.temperatureProfileFile=None
    self.format="pdf"
    self.outputFile="workPlot"
  def parseXML(self,element):
    
    #get outputFile
    temp=xmlParseFunctions.getOneChildElementText(element,"output-file-name")
    if temp!=None:
      self.outputFile=temp
    
    #check format based on outputFile
    outputFileExtension=self.outputFile[len(self.outputFile)-3:len(self.outputFile)]
    if outputFileExtension not in supportedFormats:
      raise Exception("output file extension of \""+outputFileExtension
        +"\" not a supported format "+str(supportedFormats))
    self.format=outputFileExtension
    self.outputFile=self.outputFile[0:len(self.outputFile)-4]
    
    #get temperature profile file
    temp=xmlParseFunctions.getOneChildElementText(element,"temp-profile-file")
    if temp!=None:
      self.temperatureProfileFile=temp
    
    #get grid
    temp=xmlParseFunctions.getOneChildElementTextToBool(element,"grid")
    if temp!=None:#if None don't override default
      self.grid=temp
    
    #get points
    temp=xmlParseFunctions.getOneChildElementTextToBool(element,"with-points")
    if temp!=None:#if None don't override default
      self.points=temp
    
    #get lines
    temp=xmlParseFunctions.getOneChildElementTextToBool(element,"with-lines")
    if temp!=None:#if None don't override default
      self.lines=temp
    
    #get minTemp
    temp=xmlParseFunctions.getOneChildElementTextToFloat(element,"min-temp")
    if temp!=None:#if None don't override default
      self.minTemp=temp
    
    #get ymin
    temp=xmlParseFunctions.getOneChildElementTextToFloat(element,"y-min")
    if temp!=None:#if None don't override default
      self.ylim[0]=temp
    
    #get ymax
    temp=xmlParseFunctions.getOneChildElementTextToFloat(element,"y-max")
    if temp!=None:#if None don't override default
      self.ylim[1]=temp
class PdVPlotSettings:
  def __init__(self):
    self.startZone=2
    self.points=False
    self.grid=True
    self.format="pdf"
    self.outputFile="workPlot"
    self.lines=True
    self.show=False
  def parseXML(self,element):
    
    #get outputFile
    temp=xmlParseFunctions.getOneChildElementText(element,"output-file-name")
    if temp!=None:
      self.outputFile=temp
    
    #check format based on outputFile
    outputFileExtension=self.outputFile[len(self.outputFile)-3:len(self.outputFile)]
    if outputFileExtension not in supportedFormats:
      raise Exception("output file extension of \""+outputFileExtension
        +"\" not a supported format "+str(supportedFormats))
    self.format=outputFileExtension
    self.outputFile=self.outputFile[0:len(self.outputFile)-4]
    
    #get start zone
    temp=xmlParseFunctions.getOneChildElementTextToInt(element,"start-zone")
    if temp!=None:
      self.startZone=temp
    
    #get points
    temp=xmlParseFunctions.getOneChildElementTextToBool(element,"with-points")
    if temp!=None:#if None don't override default
      self.points=temp
    
    #get lines
    temp=xmlParseFunctions.getOneChildElementTextToBool(element,"with-lines")
    if temp!=None:#if None don't override default
      self.lines=temp
    
    #get grid
    temp=xmlParseFunctions.getOneChildElementTextToBool(element,"with-grid")
    if temp!=None:#if None don't override default
      self.grid=temp
    
    #get show
    temp=xmlParseFunctions.getOneChildElementTextToBool(element,"show-plots")
    if temp!=None:#if None don't override default
      self.show=temp
class Settings:
  def __init__(self,oldColumns=False):
    """Initialize settings"""
    
    if not oldColumns:
      #newer format
      self.pColumn=62#the column number for P_ave in the profile file, starting at 0 for first column
      self.pColumnHeader="P_ave[dynes/cm^2](63)"
      self.tColumn=49#the column number for T_ave in the profile file
      self.tColumnHeader="T_ave[K](50)"
    else:
      #older profile format
      self.tColumn=45#the column number for T_ave in the profile file
      self.tColumnHeader="T_ave[K](46)"
      self.pColumn=58#the column number for P_ave in the profile file, starting at 0 for first column
      self.pColumnHeader="P_ave[dynes/cm^2](59)"
    
    #unchanged column references
    self.rhoColumn=6#the column number for Rho_ave in the profile file
    self.rhoColumnHeader="D_ave[g/cm^3](7)"
    self.QColumn=35
    self.QColumnHeader="Q[dynes/cm^2](36)"
    self.deltaMColumn=2#the column number for delta M in the profile file, start at 0 for the first column
    self.deltaMColumnHeader="DM_r[g](3)"
    self.AV=True#if true will include artificial viscosity in pressure
    self.outputFile="out.pdf"
    self.plotPdVCurves=False
  def parseXML(self,fileName):
    """Get user settings from XML file"""
    
    tree=xml.parse(fileName)
    root=tree.getroot()
    
    #get file set
    self.files=xmlParseFunctions.getOneChildElementText(root,"files")
    
    #get if using AV
    temp=xmlParseFunctions.getOneChildElementTextToBool(root,"include-av")
    if temp!=None:
      self.AV=temp
    
    #get workplot settings
    workPlotsElement=xmlParseFunctions.getOneChildElement(root,"work-plots")
    
    #if there is a workplot element
    if workPlotsElement!=None:
      
      self.workPlotSettings=WorkPlotSettings()
      self.workPlotSettings.parseXML(workPlotsElement)
        
    #get PdV plot settings
    pdVPlotsElement=xmlParseFunctions.getOneChildElement(root,"PdV-plots")
    
    #if there is a pdVPlotsElement element
    if pdVPlotsElement!=None:
    
      #there is this element we want to make PdV curves
      self.plotPdVCurves=True
      
      self.PdVPlotSettings=PdVPlotSettings()
      self.PdVPlotSettings.parseXML(pdVPlotsElement)
def getPeriodRanges(baseFileName,start,end,options):
  #check that the averagePKE.txt file is there, and make sure it has all entries
  averagePKEFile=os.path.dirname(baseFileName)+"/averagePKE.txt"
  average_PKE.averagePKE(start,sys.maxint,baseFileName,options)
  
  #get periods in range indicated
  averagePKEFileData=datafile.DataFile()
  print __name__+": reading file \"",averagePKEFile,"\" ..."
  averagePKEFileData.readFile(averagePKEFile)
  periodCount=0
  periodRange=[]
  nSet=0
  time=[]
  for i in range(len(averagePKEFileData.fColumnValues)-1):
    if averagePKEFileData.fColumnValues[i][0]>=start and averagePKEFileData.fColumnValues[i][0]<end:
      if averagePKEFileData.fColumnValues[i][3]!=None and nSet==0:
        periodRange.append([averagePKEFileData.fColumnValues[i][0],None])
        nSet=1
      elif averagePKEFileData.fColumnValues[i][3]!=None and nSet==1:
        nSet=2
      elif averagePKEFileData.fColumnValues[i][3]!=None and nSet==2:
        if periodCount>0:
          periodRange.append([periodRange[periodCount-1][1],None])
        time.append(averagePKEFileData.fColumnValues[i][1])
        periodRange[periodCount][1]=averagePKEFileData.fColumnValues[i][0]
        nSet=1
        periodCount+=1
  return [periodRange,time]
def parseOptions():
  #setup command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END]"
    ,version="%prog 1.0",description="Computes the work done in the model over each period between "
    +"START and END")
  parser.add_option("-k","--keep",action="store_true",dest="keep"
    ,help="Keeps distributed binary files [default].",default=True)
  parser.add_option("-r","--remove",action="store_false",dest="keep"
    ,help="Removes distributed binary files")
  parser.add_option("-m","--remake",action="store_true",dest="remake"
    ,help="Will remake profiles if they already exist. [not default].",default=False)
  parser.add_option("--re-sum",action="store_true",dest="resum"
    ,help="Will re-sum all model profiles, usefull when files have problems being made from "
    +"corruption and have to be re-made. Other wise should not be used as it takes more time"
    +" [not default].",default=False)
  parser.add_option("-s","--show",dest="show",action="store_true",default=False
    ,help="Display plot to x11 display rather than saving to a file.")
  
  
  #parse command line options
  (options,args)=parser.parse_args()
  
  if len(args)==0:
    raise Exception("need to specify configuration file")
  
  #parse command line options
  return (options,args)
def main():
  
  #parse command line options
  (options,args)=parseOptions()
  
  #get XML settings
  currentSettings=Settings()
  currentSettings.parseXML(args[0])
  
  #get base file name
  [start,end,baseFileName]=disect_filename.disectFileName(currentSettings.files)
  
  #make sure that all the combined binary files have profiles made
  failedFiles=make_profiles.make_profiles(options.keep,currentSettings.files,options.remake,False)
  
  #get period ranges and times for each period from PKE file 
  [periodRange,time]=getPeriodRanges(baseFileName,start,end,options)
  
  dWSum=[]
  dWSumError=[]
  dClosestToTurn=[]
  
  #get all possible file names
  extension="_pro"+".txt"
  filesExistProfiles=glob.glob(baseFileName+"*"+extension)
  filesExistProfiles.sort()
  
  #check for the t=0 model
  if currentSettings.workPlotSettings.temperatureProfileFile!=None:#if set use user set value
    firstFile=currentSettings.workPlotSettings.temperatureProfileFile
  else:
    firstFile=baseFileName+"00000000"
  failedFiles2=make_profiles.make_profiles(options.keep,firstFile,options.remake,False)
  if not os.path.isfile(firstFile+extension):
    warnings.warn("didn't find profile or dump at \""+firstFile+extension+
      "\" using, \""+filesExistProfiles[0]+"\" instead.")
    firstFile=filesExistProfiles[0]
  else:
    firstFile+=extension
  
  fileData=datafile.DataFile()
  fileData.readFile(firstFile)
  Log10T=fileData.fColumnValues[0:len(fileData.fColumnValues)-1,currentSettings.tColumn]
  Log10T=np.array(Log10T,dtype=np.float64)#force double precision
  Log10T=np.log10(Log10T)
  
  #find out how close to the surface to go
  nEndZone=len(fileData.fColumnValues)-1
  for j in range(len(fileData.fColumnValues)-1):#position
    if Log10T[j]<currentSettings.workPlotSettings.minTemp:
      nEndZone=j
      break
  
  for n in range(len(periodRange)):
    
    #get and sort files
    files=[]
    for file in filesExistProfiles:
      intOfFile=int(file[len(baseFileName):len(file)-len(extension)])
      if intOfFile>=periodRange[n][0] and intOfFile<periodRange[n][1]:
        files.append(file)
    
    #check to make sure we have a start and end for the period
    if periodRange[n][1]==None:
      raise Exception("file range index range "+str(start)+"-"+str(end)+" should contain at least"
        +" one period as indicated by PKE peaks, but does not.")
    
    if len(files)<3:
      raise Exception("need more than 3 files to compute the work, likely alot more with file"
        +" indices in the range  ("+str(periodRange[n][0])+","+str(periodRange[n][1])+")")
    
    #for first model dump
    fileData=datafile.DataFile()
    print "reading file ",files[0]," ..."
    fileData.readFile(files[0])
    
    #read in p, and 1/rho for first file
    p=np.zeros( (len(fileData.fColumnValues)-1,len(files)+1) )
    rhoInvert=np.zeros( (len(fileData.fColumnValues)-1,len(files)+1) )
    deltaM=np.zeros( (len(fileData.fColumnValues)-1,len(files)+1) )
    temperature=np.zeros( (len(fileData.fColumnValues)-1,len(files)+1) )
    maxP=np.empty(len(fileData.fColumnValues)-1)
    maxP.fill(-1.0*sys.float_info.max)
    minP=np.empty(len(fileData.fColumnValues)-1)
    minP.fill(sys.float_info.max)
    for i in range(len(fileData.fColumnValues)-1):
      p[i][0]=fileData.fColumnValues[i][currentSettings.pColumn]
      p[i][len(files)]=fileData.fColumnValues[i][currentSettings.pColumn]
      if currentSettings.AV:
        p[i][0]+=fileData.fColumnValues[i][currentSettings.QColumn]
        p[i][len(files)]+=fileData.fColumnValues[i][currentSettings.QColumn]
        
      #get max P
      if p[i][0]>=maxP[i]:
        maxP[i]=p[i][0]
      
      #get min P
      if p[i][0]<=minP[i]:
        maxP[i]=p[i][0]
      
      rhoInvert[i][0]=1.0/fileData.fColumnValues[i][currentSettings.rhoColumn]
      rhoInvert[i][len(files)]=1.0/fileData.fColumnValues[i][currentSettings.rhoColumn]
      deltaM[i][0]=fileData.fColumnValues[i][currentSettings.deltaMColumn]
      deltaM[i][len(files)]=fileData.fColumnValues[i][currentSettings.deltaMColumn]
      temperature[i][0]=fileData.fColumnValues[i][currentSettings.tColumn]
      temperature[i][len(files)]=fileData.fColumnValues[i][currentSettings.tColumn]
    
    #read all files in for current period
    for i in range(1,len(files)):#for each dump
      
      print "reading file ",files[i]," ..."
      
      fileData.readFile(files[i])
      
      #check headers for used columns
      if fileData.sColumnNames[currentSettings.pColumn]!=currentSettings.pColumnHeader:
        warings.warn("file \""+files[i]+" has pressure column header as \""\
          +fileData.sColumnNames[currentSettings.pColumn]+" expected something like \""\
          +currentSettings.pColumnHeader+"\".")
      if currentSettings.AV:
        if fileData.sColumnNames[currentSettings.QColumn]!=currentSettings.QColumnHeader:
          warnings.warn("file \""+files[i]+" has A.V. column header as \""\
            +fileData.sColumnNames[currentSettings.QColumn]+" expected something like \""\
            +currentSettings.QColumnHeader+"\".")
      if fileData.sColumnNames[currentSettings.rhoColumn]!=currentSettings.rhoColumnHeader:
        warnings.warn("file \""+files[i]+" has density column header as \""\
          +fileData.sColumnNames[currentSettings.rhoColumn]+" expected something like \""\
          +currentSettings.rhoColumnHeader+"\".")
      if fileData.sColumnNames[currentSettings.deltaMColumn]!=currentSettings.deltaMColumnHeader:
        warnings.warn("file \""+files[i]+" has delta M_r column header as \""\
          +fileData.sColumnNames[currentSettings.deltaMColumn]+" expected something like \""\
          +currentSettings.deltaMColumnHeader+"\".")
      
      for j in range(len(fileData.fColumnValues)-1):#for each zone
        p[j][i]=fileData.fColumnValues[j][currentSettings.pColumn]
        if currentSettings.AV:
          p[j][i]+=fileData.fColumnValues[j][currentSettings.QColumn]
          
        #get max P
        if p[j][i]>=maxP[j]:
          maxP[j]=p[j][i]
        
        #get min P
        if p[j][i]<=minP[j]:
          maxP[j]=p[j][i]
        
        deltaM[j][i]=fileData.fColumnValues[j][currentSettings.deltaMColumn]
        rhoInvert[j][i]=1.0/fileData.fColumnValues[j][currentSettings.rhoColumn]
    
    #compute work
    dW=np.zeros(len(fileData.fColumnValues)-1)
    for i in range(1,len(files)+1):#time
      for j in range(0,nEndZone):#position
        
        dW[j]+=(0.5*(rhoInvert[j][i]-rhoInvert[j][i-1])*(p[j][i]+p[j][i-1]))*deltaM[j][i]
    
    dWMax=0.0
    dWMin=0.0
    dWError=np.zeros(len(fileData.fColumnValues)-1)
    dClosenessToEdge=np.zeros(len(fileData.fColumnValues)-1)
    dWSum.append(0.0)
    dClosestToTurnThisPeriod=sys.float_info.max
    for j in range(0,nEndZone):#position
      dWSum[n]+=dW[j]
      
      #
      y1=p[j][len(files)]
      y0=p[j][len(files)-1]
      ym1=p[j][len(files)-2]
      x1=rhoInvert[j][len(files)]
      x0=rhoInvert[j][len(files)-1]
      xm1=rhoInvert[j][len(files)-2]
      deltaX=(x1-x0)
      y1py0=y1+y0
      tiny=1e-300#fixes the case where things don't move, e.g. near the core
      m=(y0-ym1)/((x0-xm1)+tiny)
      y1p=y0+m*deltaX
      y1ppy0=y1p+y0
      #how far off is our guess from assuming a linear extrapolation
      dW1=0.5*deltaX*y1py0#work from setting first point equal to last point
      dW1p=0.5*deltaX*y1ppy0#work from linear extrapolation
      dWError[j]=abs(dW1-dW1p)*deltaM[j][len(files)]
      dWMax+=dW[j]+dWError[j]
      dWMin+=dW[j]-dWError[j]
      
      #how close to a turn is the start/end position
      dtotal=maxP[j]-minP[i]
      dmax=max(abs(p[j][0]-maxP[j]),abs(p[j][0]-minP[j]))
      distanceFromEdge=dmax/dtotal#will be 1/2 if near middle, will be near 0 if near edge
      
      #keep value of closest start/end position to turn
      if distanceFromEdge<dClosestToTurnThisPeriod:
        dClosestToTurnThisPeriod=distanceFromEdge
    dClosestToTurn.append(dClosestToTurnThisPeriod)
    dWSumError.append((dWMax-dWMin)*0.5)
    
    #make plots of P and 1/rho for each zone
    fig=plt.figure(figsize=(13,8))
    ax1 = plt.subplot2grid((3,3), (0,0),colspan=3,rowspan=3)
    if currentSettings.plotPdVCurves:
      for j in range(currentSettings.PdVPlotSettings.startZone,len(fileData.fColumnValues)-1):
        #print rhoInvert[j]
        make_PdV_plot(rhoInvert[j],p[j],j,n,fig,ax1,currentSettings.PdVPlotSettings)
    
    #make work plot
    make_Work_plot(Log10T,dW,dWError,fig,ax1
      ,currentSettings.workPlotSettings,n,time[n],"Log10(T)")
    make_Work_plot(range(len(fileData.fColumnValues)-1),dW,dWError,fig,ax1
      ,currentSettings.workPlotSettings,n,time[n],"zone #")
    
  #print out total work done by model
  f=open("work_per_period.txt",'w')
  f.write("period time[s]      work[ergs] work_uncertianty[ergs] relative_distance_from_turn\n")
  for i in range(len(dWSum)):
    f.write(str(i)+" "+str(time[i])+" "+str(dWSum[i])+" "+str(dWSumError[i])+" "
      +str(dClosestToTurn[i])+"\n")
  f.close()
def make_PdV_plot(x,y,zone,period,fig,ax1,settings):
  if not settings.grid:
    ax1.grid()
  plotString='-'
  if settings.lines:
    plotString='o'
  if settings.points:
    plotString=plotString+'o'
  line1=ax1.plot(x,y,plotString,markersize=3)
  line1=ax1.plot(x[0],y[0],'x',markersize=15)
  line1=ax1.plot(x[len(x)-1],y[len(x)-1],'x',markersize=15)
  ax1.set_xlabel("<1/rho>")
  ax1.set_ylabel("<P>")
  
  sOutFileName=settings.outputFile+"_period"+str(period)+"_zone"+str(zone)+"."+settings.format
  print __name__+":"+main.__name__+": creating plot \""+sOutFileName+"\" ..."
  fig.savefig(sOutFileName,format=settings.format,transparent=False)#save to file
  
  ax1.cla()#clear plot
  return True
def make_Work_plot(x,y,yerr,fig,ax1,settings,period,time,xlabel):
  if settings.grid:
    ax1.grid(True,which="both")
  plotString='-'
  if settings.lines and settings.points:
    plotString='-o'
  elif settings.points:
    plotString='o'
  
  line1=ax1.errorbar(x,y,yerr=yerr,fmt=plotString,markersize=5,ecolor='r')
  
  ax1.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
  ax1.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
  
  fig.suptitle("Work plot for period "+str(period)+" at time "+str(time)+"[s]")
  
  #set yrange
  if settings.ylim!=None:
    ax1.set_ylim(settings.ylim)
  
  ax1.set_ylabel("Work [ergs]")
  ax1.set_ylabel(xlabel)
  
  sOutFileName=""
  if xlabel=="zone #":
    sOutFileName=settings.outputFile+"_zones_period"+str(period)+"."+settings.format
  elif xlabel=="Log10(T)":
    sOutFileName=settings.outputFile+"_temp_period"+str(period)+"."+settings.format
  
  print __name__+":"+main.__name__+": creating plot \""+sOutFileName+" ..."
  fig.savefig(sOutFileName,format=settings.format,transparent=False)#save to file
  
  ax1.cla()#clear plot
  return True
if __name__ == "__main__":
  main()