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

def parseOptions():
  #note: newlines are not respected in the optparse description string :(, maybe someday will use
  #argparse, which does allow for raw formating (repects indents, newlines etc.)
  
  #setup command line parser
  '''This is out of date, needs to be updated to reflect new additions'''
  parser=op.OptionParser(usage="Usage: %prog [options] XMLFILE"
    ,version="%prog 1.0",description=r"Reads in the XMLFILE which specifies how the plot is to be "
    +"made and generates the plot of plots accordingly.")
    
  #these options apply globaly
  parser.add_option("-o","--outputFile",dest="outputFile",default="out"
    ,help="Specifies the OUTPUTFILE. [default: %default]"
    ,metavar="OUTPUTFILE", type="string")
  parser.add_option("-f","--format",dest="format",default="png", type="choice",
    help="Sets the FMT of the OUTPUTFILE, availble formats are 'png', 'pdf', 'ps', 'eps', and 'svg'."
    +"[default: %default]", metavar="FMT",choices=('png','pdf','ps','eps','svg'))
  parser.add_option("-s","--show",dest="show",action="store_true",default=False
    ,help="Display plot to x11 display rather than saving to a file.")
  parser.add_option("--dpi",dest="dpi",type="float",default=100
    ,help="Sets the dots per inch of the figure [default: %default]")
  parser.add_option("--space-axis-evenly",action="store_true",dest="spaceAxisEvenly"
    ,help="Will give the same amount of space per x-axis, instead of per plot. [not default]."
    ,default=False)
  
  #should likely be per plot, at the momement they are across all plots
  parser.add_option("--fig-width",dest="figWidth",type="float",default=15
    ,help="Set the width of the figure [default: %default]")
  parser.add_option("--fig-height",dest="figHeight",type="float",default=9
    ,help="Set the height of the figure [default: %default]")
  
  #parse command line options
  return parser.parse_args()
class Curve:
  '''This class holds all the information for a curve on a plot.'''
  def __init__(self,element):
    
    '''This method initilizes a curve object, the type parameter allows checking curve syntax 
    against axis syntax to see if they match.'''
    
    self.nColumnX=None
    self.nColumnY=None
    self.y=[]
    self.x=[]
    self.index=[]
    self.formulaOrigY=None
    self.formulaOrigX=None
    self.formulaX=None
    self.formulaY=None
    self.codeY=None
    self.codeX=None
    self.style=""
    self.color=""
    self.markersize=3.0
    self.linewidth=1.0
    self.label=None
    self.fileReference=None
    
    #set curve style
    if element.get("style")!=None:
      self.style=element.get("style")
    
    #set curve color
    if element.get("color")!=None:
      self.color=element.get("color")
      
    #set curve label
    if element.get("label")!=None:
      self.label=element.get("label")
    
    #get marker size
    if element.get("markersize")!=None:
      if isFloat(element.get("markersize")):
        self.markersize=float(element.get("markersize"))
      else:
        print "\"markersize\" must be a float, got \"",element.get("markersize"),"\""
    
    #get linewidth
    if element.get("linewidth")!=None:
      if isFloat(element.get("linewidth")):
        self.linewidth=float(element.get("linewidth"))
      else:
        print "\"linewidth\" must be a float, got \"",element.get("linewidth"),"\""
    
    #get file references
    if element.get("file")!=None and element.get("file"):
      self.fileReference=element.get("file")
    else:
      print "Curve must have a \"file\" set to the name of a file form which the data is to be taken."
      quit()
    
    #get xcolumns
    if element.get("xcolumn").isdigit():
      self.nColumnX=int(element.get("xcolumn"))-1
    else:
      
      #if there is a foluma get list of columns
      self.formulaOrigX=element.get("xcolumn")
      columnsTemp=self.formulaOrigX.split('$')
      self.formulaX=columnsTemp[0]
      if len(columnsTemp)>1:#turn into a list
        self.nColumnX=[]
        for i in range(1,len(columnsTemp)):
          columnAndFormula=splitFirstInt(columnsTemp[i])
          self.formulaX+="a["+str(i-1)+"]"+columnAndFormula[1]
          self.nColumnX.append(columnAndFormula[0]-1)
        self.codeX=parser.expr(self.formulaX).compile()
      else:
        self.formulaX=None
        if element.get("xcolumn") == None or element.get("xcolumn")=="":
          print "\"xcolumn\" number given for curve, must be either an integer, or"\
            +" a mathatmical formula containing column references prefixed with a \"$\"."
          quit()
    
    #get ycolumn
    if element.get("ycolumn").isdigit():
      self.nColumnY=int(element.get("ycolumn"))-1
    else:
      
      #if there is a foluma get list of columns
      self.formulaOrigY=element.get("ycolumn")
      columnsTemp=self.formulaOrigY.split('$')
      self.formulaY=columnsTemp[0]
      if len(columnsTemp)>1:#turn into a list
        self.nColumnY=[]
        for i in range(1,len(columnsTemp)):
          columnAndFormula=splitFirstInt(columnsTemp[i])
          self.formulaY+="a["+str(i-1)+"]"+columnAndFormula[1]
          self.nColumnY.append(columnAndFormula[0]-1)
        self.codeY=parser.expr(self.formulaY).compile()
      else:
        self.formulaY=None
        if element.get("ycolumn") == None or element.get("ycolumn")=="":
          print "\"xcolumn\" number given for curve, must be either an integer, or"\
            +" a mathatmical formula containing column references prefixed with a \"$\"."
          quit()
  def load(self,files,options):
    '''This method adds a y value and index to the curve for the current fileData.'''
    
    #get y's from file
    for i in range(len(files[self.fileReference].fColumnValues)):
      y=getY(self.nColumnY,files[self.fileReference],self.codeY,i)
      x=getY(self.nColumnX,files[self.fileReference],self.codeX,i)
      if x!=None and y!=None:
        self.y.append(y)
        self.x.append(x)
class Plot:
  '''This class holds all the information for a single plot, namely the list of curves for that plot.'''
  
  def __init__(self,element):
    '''This method initlizes the plot object'''
    
    self.ylabel=None
    self.curves=[]
    self.limits=None
    self.grid=None
    self.bMinorTics=False
    self.legendloc=1
    
    #check for grid setting
    if element.get("grid")!=None:
      if element.get("grid").lower() in ["major","minor","both"]:
        self.grid=element.get("grid").lower()
      else:
        print "Grid type for plot unrecognized \"",element.get("grid"),"\" should be \"major\","\
          +" \"minor\", or \"both\"."
        quit()
    
    #check if using minor tics
    if element.get("yminortics")!=None:
      if element.get("yminortics").lower() in ["true","1","t","yes","y"]:
        self.bMinorTics=True
    
    #get ylabel
    if element.get("ylabel")!=None and element.get("ylabel")!="":
      self.ylabel=element.get("ylabel")
    
    #get range
    yMin=None
    if element.get("ymin")!=None and element.get("ymin")!="":
      yMin=float(element.get("ymin"))
    yMax=None
    if element.get("ymax")!=None and element.get("ymax")!="":
      yMax=float(element.get("ymax"))
    self.limits=[yMin,yMax]
    
    #get ledgend location
    if element.get("legendloc")!=None and element.get("legendloc")!="":
      self.legendloc=int(element.get("legendloc"))
    
    #add curves to plot
    curveElements=element.findall("curve")
    for curveElement in curveElements:
      self.curves.append(Curve(curveElement))
  def load(self,files,options):
    '''loads the data for a plot, y-data is stored in the curves, and sets the ylabel from the first
    file read in'''
    
    for curve in self.curves:
      curve.load(files,options)
class Axis:
  '''This class holds all the information needed for a particular x-axis. An axis can either be 
  either of time, or of some column in the data files.'''
  
  def __init__(self,element,options):
    '''This function initizalizes the axis object.'''
    
    self.plots=[]
    self.xlabel=None
    self.limits=None
    self.bMinorTics=False
    
    #check for grid setting
    if element.get("grid")!=None:
      if element.get("grid").lower() in ["major","minor","both"]:
        self.grid=element.get("grid").lower()
      else:
        print "Grid type for plot unrecognized \"",element.get("grid"),"\" should be \"major\","\
          +" \"minor\", or \"both\"."
        quit()
    
    #check if using minor tics
    if element.get("xminortics")!=None:
      if element.get("xminortics").lower() in ["true","1","t","yes","y"]:
        self.bMinorTics=True
    
    #get xlabel
    self.xlabel=element.get("xlabel")
    if self.xlabel==None and self.bTime:
      if self.period!=None:
        self.xlabel="phase"
      else:
        self.xlabel="t [s]"
    
    #get xrange
    xMin=None
    if element.get("xmin")!=None and element.get("xmin")!="":
      xMin=float(element.get("xmin"))
    xMax=None
    if element.get("xmax")!=None and element.get("xmax")!="":
      xMax=float(element.get("xmax"))
    self.limits=[xMin,xMax]
    
    #add plots to axis
    plotElements=element.findall("plot")
    for plotElement in plotElements:
      self.plots.append(Plot(plotElement))
  def load(self,files,options):
    '''This function loads the values needed for the x-axis data from the fileData argument'''
    
    #load plots
    for plot in self.plots:
      plot.load(files,options)
class DataSet:
  '''This class holds all the information for a single dataSet, which includes the baseFileName of 
  the dataset, the range of the dataSet (start-end), the times and phases of the files within the 
  range of the dataSet, and the plots made from the dataSet.'''
  
  def __init__(self,element,options):
    '''Initilizes the dataSet by setting baseFileName, start, end, and intilizing plots from an xml 
    element'''
    
    #set some initial values
    self.axes=[]
    self.files={}
    
    #load files
    fileElements=element.findall("file")
    for fileElement in fileElements:
      fileData=datafile.DataFile()
      fileData.readFile(fileElement.text)
      self.files[fileElement.get("name")]=fileData
      
    #add axes to dataset
    axisElements=element.findall("axis")
    for axisElement in axisElements:
      axis=Axis(axisElement,options)
      self.axes.append(axis)
  def load(self,options):
    '''Loads the dataSet, this means that it sets, time, phases, and plots data'''
    for axis in self.axes:
      axis.load(self.files,options)
  def getCurve(self,ID):
    '''Returns a curve object that has ID, ID'''
    nCount=0
    for axis in self.axes:
      for plot in axis.plots:
        for curve in plot.curves:
          if nCount==ID:
            return curve
          nCount=nCount+1
def plot(dataSets,options,title):
  
  #set import options based on weather it will saved to a file, or sent to x11
  #print "interactive backends=",matplotlib.rcsetup.interactive_bk
  #print "non-interactive backends=",matplotlib.rcsetup.non_interactive_bk
  #print "all backends=",matplotlib.rcsetup.all_backends
  
  import matplotlib
  if not options.show:
    matplotlib.use("Agg")
  import matplotlib.pyplot as plt
  
  from matplotlib.gridspec import GridSpec
  
  #create a figure object, to be reused for each figure, could be 1 or more
  fig=plt.figure(figsize=(options.figWidth,options.figHeight))
  
  #add title
  if title!="":
    fig.suptitle(title)
  
  #basic plot spacing options
  #axisSpacing=0.05
  #figBottom=0.05
  #figTop=0.95
  axisSpacing=0.05
  figBottom=0.1
  figTop=0.93
  
  #count number of axes in all plots, this number will be the same for all plots.
  nNumAxes=0
  for dataSet in dataSets:
    nNumAxes=nNumAxes+len(dataSet.axes)
  heightAxis=((figTop-figBottom)-axisSpacing*(nNumAxes-1.0))/nNumAxes
  
  #count number of plots
  nNumPlots=0
  for dataSet in dataSets:
    for axis in dataSet.axes:
      nNumPlots=nNumPlots+len(axis.plots)
  heightPlot=((figTop-figBottom)-axisSpacing*(nNumAxes-1.0))/nNumPlots
  
  #create plots
  fileCount=0
  fileStart=0
  axisCount=0
  nTotalPlotCount=1
  dataSetCount=0
  ax=[]
  gs=[]
  top=figTop
  for dataSet in dataSets:
    for axisMine in dataSet.axes:
      gs.append(GridSpec(len(axisMine.plots),1))
      bottom=top-heightPlot*(len(axisMine.plots))
      if options.spaceAxisEvenly:
        top=figTop-axisCount*(heightAxis+axisSpacing)
        bottom=figBottom+(nNumAxes-1.0-axisCount)*(heightAxis+axisSpacing)
      gs[axisCount].update(top=top,bottom=bottom,hspace=0.0)
      
      #for each plot
      nPlotCount=0
      for plot in axisMine.plots:
        
        #create a subplot
        ax.append(plt.subplot(gs[axisCount][nPlotCount,0]))
        
        #set minor axis tics
        if axisMine.bMinorTics:
          ax[nTotalPlotCount-1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        if plot.bMinorTics:
          ax[nTotalPlotCount-1].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        
        #set grid settings
        if plot.grid!=None:
          ax[nTotalPlotCount-1].grid(True, which=plot.grid)
        else:
          ax[nTotalPlotCount-1].grid(False)
        
        #for each curve in the plot
        lines=[]
        labels=[]
        curveCount=0
        for curve in plot.curves:
          
          #plot the curve
          temp=ax[nTotalPlotCount-1].plot(curve.x,curve.y,str(curve.color)+str(curve.style)
            ,markersize=curve.markersize,linewidth=curve.linewidth)
          if curve.label!=None and curve.label!="":
            lines.append(temp)
            labels.append(curve.label)
          curveCount=curveCount+1
      ax[nTotalPlotCount-1].set_xlim(axisMine.limits)
      ax[nTotalPlotCount-1].set_ylim(plot.limits)
      
      #set legend
      if len(lines)>0:
        ax[nTotalPlotCount-1].legend(lines,labels,loc=plot.legendloc)
      
      #remove x-axis labels
      if nPlotCount!=len(axisMine.plots)-1:
        plt.setp(plt.gca(), 'xticklabels', [])
      else:
        ax[nTotalPlotCount-1].set_xlabel(axisMine.xlabel)
      
      #set y-axis labels
      ax[nTotalPlotCount-1].set_ylabel(plot.ylabel)
      
      #remove top and bottom y-axis tic labels
      ax[nTotalPlotCount-1].yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(prune='both'))
      
      #adjust hspace to be 0 so plots for the same x-axis will be tight
      nPlotCount=nPlotCount+1#resets every axis
      nTotalPlotCount=nTotalPlotCount+1#continues to increase for all plots
      
      top=bottom-axisSpacing
      axisCount=axisCount+1
  
  #done making plot
  if options.show:
    print "ploting file ",i," to screen, close for next plot"
    plt.show()
  else:
    [path,ext]=os.path.splitext(options.outputFile)
    supportedFileTypes=["png", "pdf", "ps", "eps", "svg"]
    if ext[1:] not in supportedFileTypes:
      print "File type \""+ext[1:]+"\" not suported. Supported types are ",supportedFileTypes," please choose one of those"
      quit()
    print __name__+":"+main.__name__+": saving figure to file \""+options.outputFile
    fig.savefig(options.outputFile,format=ext[1:],transparent=False,dpi=options.dpi)#save to file
def splitFirstInt(str):
  '''Returns the integer which starts at the very first character of str, and drops everything else 
  from str'''
  
  strInt=''
  i=0
  while str[i].isdigit() or (str[i]=='-' and i==0):
    strInt=strInt+str[i]
    i=i+1
    if i>=len(str):#check that we haven't passed the end of the string
      break
  return [int(strInt),str[i:]]
def splitFirstFloat(str):
  '''Returns the integer which starts at the very first character of str, and drops everything else 
  from str'''
  
  strInt=''
  i=0
  nExponentIndex=-1
  nDecimalCount=0
  while str[i].isdigit()\
    or ((str[i]=='-' or str[i]=='+') and i==0)\
    or ((str[i]=="e" or str[i]=="E") and nExponentIndex==-1)\
    or ((str[i]=="-" or str[i]=="+") and nExponentIndex+1==i)\
    or (str[i]=="." and nDecimalCount==0) :
    if str[i]=="e" or str[i]=="E":
      nExponentIndex=i
    if str[i]==".":
      nDecimalCount=nDecimalCount+1
    strInt=strInt+str[i]
    i=i+1
    if i>=len(str):#check that we haven't passed the end of the string
      break
  return [float(strInt),str[i:]]
def getY(nColumn,fileData,code,i):
  if isinstance(nColumn,int):
    #normal simple 1 column data
    return fileData.fColumnValues[i][nColumn]
  else:
    #new column operation
    a=[]
    for j in range(len(nColumn)):
      if fileData.fColumnValues[i][nColumn[j]]==None:
        return None
      else:
        a.append(fileData.fColumnValues[i][nColumn[j]])
    return eval(code)
def isFloat(str):
  try:
    float(str)
    return True
  except:
    return False
def main():
  
  #parse command line options
  (options,args)=parseOptions()
  
  #get root element
  tree=xml.parse(args[0])
  root=tree.getroot()
  
  #set plot title
  title=""
  if root.get("title")!=None and root.get("title")!="":
    title=root.get("title")
  
  #set plot output
  if root.get("outputfile")!=None and root.get("outputfile")!="":
    options.outputFile=root.get("outputfile")
    [path,ext]=os.path.splitext(options.outputFile)
    supportedFileTypes=["png", "pdf", "ps", "eps", "svg"]
    if ext[1:] not in supportedFileTypes:
      print "File type \""+ext[1:]+"\" not suported. Supported types are ",supportedFileTypes," please choose one of those"
      quit()
  
  #get list of files elements
  dataSetElements=root.findall("dataSet")
  
  #initialize dataSets
  dataSets=[]
  for dataSetElement in dataSetElements:
    dataSets.append(DataSet(dataSetElement,options))
  
  #load datasets
  i=0
  for dataSet in dataSets:
    print "reading in data set ",i
    dataSet.load(options)
    i=i+1
  
  #plot datasets
  plot(dataSets,options,title)
  
  return
if __name__ == "__main__":
  main()