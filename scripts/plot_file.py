#!/usr/bin/env python

import datafile
import optparse as op
import make_profiles
import glob
import numpy as np
import sys
import os
import parser
from math import *
import xml.etree.ElementTree as xml
import parse_formula
figureNumber=0
def checkForRecogAttrib(element,recogAttrib):
  """Checks for attributes in element that arn't in recogAttrib
  
  Elements not in recogAttrib are added to a list of attributes not recognized
  and are returned
  """
  
  notRecogAttrib=[]
  for key in element.attrib.keys():
    if key not in recogAttrib:
      notRecogAttrib.append(key)
  return notRecogAttrib
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
class Text:
  """This class holds informatin for a text object on a plot.
  
  """
  
  def __init__(self,element):
    """This method initializest a text object from an xml element
    
    """
    
    self.x=0
    self.y=0
    self.text=None
    
    #check for unrecognized attributes
    recogAttribs=[
      "x",
      "y",
      ]
    notRecogAttribs=checkForRecogAttrib(element,recogAttribs)
    if len(notRecogAttribs)>0:
      raise Exception("attributes in text "+str(notRecogAttribs)
      +" are not in list of recognized text attributes "+str(recogAttribs))
    
    
    #get x position
    if element.get("x")!=None:
      if isFloat(element.get("x")):
        self.x=float(element.get("x"))
      else:
        print "\"x\" must be a float, got \"",element.get("x"),"\""
    
    #get y position
    if element.get("y")!=None:
      if isFloat(element.get("y")):
        self.y=float(element.get("y"))
      else:
        print "\"y\" must be a float, got \"",element.get("y"),"\""
    
    #set text
    self.text=element.text
class Curve:
  """Holds all the information for a curve on a plot.
  
  """
  
  def __init__(self,element):
    """Initilizes a curve object
    
    """
    
    self.nColumnX=None
    self.nColumnY=None
    self.nColumnErr=None
    self.y=[]
    self.x=[]
    self.error=[]
    self.index=[]
    self.formulaOrigY=None
    self.formulaOrigX=None
    self.formulaOrigErr=None
    self.formulaX=None
    self.formulaY=None
    self.formulaErr=None
    self.codeY=None
    self.codeX=None
    self.codeErr=None
    self.style=" "
    self.color="b"
    self.markerfacecolor="b"
    self.markeredgecolor="b"
    self.markersize=2.0
    self.linewidth=1.0
    self.label=None
    self.fileReference=None
    self.nRowShiftErr=None
    self.nRowShiftX=None
    self.nRowShiftY=None
    self.marker=None
    self.ecolor="red"
    self.elinewidth=1.0
    self.capsize=1.0
    
    #check for unrecognized attributes
    recogAttribs=[
      "line",
      "marker",
      "color",
      "markerfacecolor",
      "markeredgecolor",
      "label",
      "markersize",
      "capsize",
      "elinewidth",
      "ecolor",
      "linewidth",
      "file",
      "errcolumn",
      "xcolumn",
      "ycolumn",
      ]
    notRecogAttribs=checkForRecogAttrib(element,recogAttribs)
    if len(notRecogAttribs)>0:
      raise Exception("attributes in cuve "+str(notRecogAttribs)
      +" are not in list of recognized curve attributes "+str(recogAttribs))
    
    #set curve style
    if element.get("line")!=None:
      self.style=element.get("line")
      
    #set curve marker
    if element.get("marker")!=None:
      self.marker=element.get("marker")
      
    #set curve color
    if element.get("color")!=None:
      self.color=element.get("color")
    
    #set curve marker face color
    if element.get("markerfacecolor")!=None:
      self.markerfacecolor=element.get("markerfacecolor")
    elif element.get("color")!=None:#if not set use the color element
      self.markerfacecolor=element.get("color")
      
    #set curve marker edge color
    if element.get("markeredgecolor")!=None:
      self.markeredgecolor=element.get("markeredgecolor")
    elif element.get("color")!=None:#if not set use the color element
      self.markeredgecolor=element.get("color")
      
    #set curve label
    if element.get("label")!=None:
      self.label=element.get("label")
    
    #get marker size
    if element.get("markersize")!=None:
      if isFloat(element.get("markersize")):
        self.markersize=float(element.get("markersize"))
      else:
        print "\"markersize\" must be a float, got \"",element.get("markersize"),"\""
    
    #get capsize
    if element.get("capsize")!=None:
      if isFloat(element.get("capsize")):
        self.capsize=float(element.get("capsize"))
      else:
        print "\"capsize\" must be a float, got \"",element.get("capsize"),"\""
    
    #get elinewidth
    if element.get("elinewidth")!=None:
      if isFloat(element.get("elinewidth")):
        self.elinewidth=float(element.get("elinewidth"))
      else:
        print "\"elinewidth\" must be a float, got \"",element.get("elinewidth"),"\""
        
    #get ecolor
    if element.get("ecolor")!=None:
      if isFloat(element.get("ecolor")):
        self.ecolor=float(element.get("ecolor"))
      else:
        print "\"ecolor\" must be a float, got \"",element.get("ecolor"),"\""
    
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
    
    #get xcolumn
    [self.formulaOrigX,self.formulaX,self.nColumnX,self.nRowShiftX,self.codeX]=\
      parse_formula.getFormula(element.get("xcolumn"))
    
    #get ycolumn
    [self.formulaOrigY,self.formulaY,self.nColumnY,self.nRowShiftY,self.codeY]=\
      parse_formula.getFormula(element.get("ycolumn"))
    
    #get errcolumns
    if element.get("errcolumn")!=None:
      [self.formulaOrigErr,self.formulaErr,self.nColumnErr,self.nRowShiftErr,self.codeErr]=\
        parse_formula.getFormula(element.get("errcolumn"))
  def load(self,files,options):
    """Adds a y-value and index to the curve for the current fileData."""
    
    #get y's from file
    for i in range(len(files[self.fileReference].fColumnValues)):
      y=parse_formula.getY(self.nRowShiftY,self.nColumnY,files[self.fileReference],self.codeY,i)
      x=parse_formula.getY(self.nRowShiftX,self.nColumnX,files[self.fileReference],self.codeX,i)
      if self.nColumnErr!=None:
        error=parse_formula.getY(self.nRowShiftErr,self.nColumnErr,files[self.fileReference]
          ,self.codeErr,i)
        if error!=None:
          self.error.append(error)
      if x!=None and y!=None:
        self.y.append(y)
        self.x.append(x)
class Plot:
  """Holds all the information for a single plot
  
  namely the list of curves for that plot.
  
  """
  
  def __init__(self,element):
    """Initlizes the plot object
    
    """
    
    self.ylabel=None
    self.curves=[]
    self.texts=[]
    self.limits=None
    self.grid=None
    self.bMinorTics=False
    self.legendloc=1
    self.numpoints=None
    self.weightHeight=1.0
    self.autoLimits=[True,True]
    
    recogAttribs=[
      "grid",
      "yminortics",
      "ticks",
      "ylabel",
      "ymin",
      "ymax",
      "weightHeight",
      "legendloc",
      "numlegendpoints"
      ]
    notRecogAttribs=checkForRecogAttrib(element,recogAttribs)
    if len(notRecogAttribs)>0:
      raise Exception("attributes in plot "+str(notRecogAttribs)
      +" are not in list of recognized plot attributes "+str(recogAttribs))
    
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
    
    self.ticks=None
    
    #check for manual ticks
    if element.get("ticks")!=None:
      tickStrings=element.get("ticks").split()
      if len(tickStrings)>0:
        self.ticks=[]
      for tickString in tickStrings:
        self.ticks.append(float(tickString))
    
    #get ylabel
    #if element.get("ylabel")!=None and element.get("ylabel")!="":
    self.ylabel=element.get("ylabel")
    
    #get range
    yMin=None
    if element.get("ymin")!=None and element.get("ymin")!="":
      yMin=float(element.get("ymin"))
      self.autoLimits[0]=False
    yMax=None
    if element.get("ymax")!=None and element.get("ymax")!="":
      yMax=float(element.get("ymax"))
      self.autoLimits[1]=False
    self.limits=[yMin,yMax]
    
    #get weight to use when deciding plot height
    if element.get("weightHeight")!=None and element.get("weightHeight")!="":
      self.weightHeight=float(element.get("weightHeight"))
    
    #get legend location
    if element.get("legendloc")!=None and element.get("legendloc")!="":
      self.legendloc=int(element.get("legendloc"))
    
    #get numlegendpoints
    if element.get("numlegendpoints")!=None and element.get("numlegendpoints")!="":
      self.numpoints=int(element.get("numlegendpoints"))
    
    #add curves to plot
    curveElements=element.findall("curve")
    for curveElement in curveElements:
      self.curves.append(Curve(curveElement))
    
    #add text objects to plot
    textElements=element.findall("text")
    for textElement in textElements:
      self.texts.append(Text(textElement))
  def load(self,files,options):
    """Loads the data for a plot
    
    y-data is stored in the curves, and sets the ylabel from the first
    file read in
    
    """
    
    #load the plots
    for curve in self.curves:
      curve.load(files,options)
  def setLimits(self,xlimits):
    """Sets the y limits
    
    from the maximum, and minimum y values in the curves on
    the plot in the xrange specified, over all files.
    
    """
    
    #if both y limits were already set in configuraiton file, don't set anything
    if self.autoLimits[0]==False and self.autoLimits[1]==False:
      return
    
    #Only upper bound set
    if xlimits[1]==None:
      xLimitsTemp=[xlimits[0],sys.float_info.max]
      
    #Only lower bound set
    elif xlimits[0]==None:
      xLimitsTemp=[-1.0*sys.float_info.max,xlimits[1]]
    
    #Both Bounds set
    else:
      xLimitsTemp=xlimits
    
    minyInCurves=sys.float_info.max
    maxyInCurves=-1.0*sys.float_info.max
    
    for curve in self.curves:
      
      #check that num x, and num y match for each curve, for each file
      if len(curve.x)!=len(curve.y):
        raise Exception("number of x points for axis, and number of y points"\
          +" for curve don't match")
      
      for i in range(len(curve.y)):
        
        #if in xrange
        if xLimitsTemp[0]<=curve.x[i] and curve.x[i]<=xLimitsTemp[1]:
          
          #keep smallest y
          if minyInCurves> curve.y[i] and curve.y[i]!=None:
            minyInCurves=curve.y[i]
          
          #keep largest y
          if maxyInCurves< curve.y[i] and curve.y[i]!=None:
            maxyInCurves=curve.y[i]
      
    #set only limits which are not set in configuration file
    if self.autoLimits[0]==True and self.autoLimits[1]==True:
      self.limits=[minyInCurves,maxyInCurves]
    elif self.autoLimits[0]==True:
      self.limits[0]=minyInCurves
    elif self.autoLimits[1]==True:
      self.limits[1]=maxyInCurves
class Axis:
  """Holds all the information needed for a particular x-axis.
  
  """
  
  def __init__(self,element,options):
    """Initizalizes the axis object.
    
    """
    
    self.plots=[]
    self.xlabel=None
    self.limits=None
    self.bMinorTics=False
    self.ticks=None
    
    #check for unrecognized attributes
    recogAttribs=[
      "ticks",
      "grid",
      "xminortics",
      "ticks",
      "xlabel",
      "xmin",
      "xmax",
      ]
    notRecogAttribs=checkForRecogAttrib(element,recogAttribs)
    if len(notRecogAttribs)>0:
      raise Exception("attributes in axis "+str(notRecogAttribs)
      +" are not in list of recognized axis attributes "+str(recogAttribs))
    
    #check for manual ticks
    if element.get("ticks")!=None:
      tickStrings=element.get("ticks").split()
      if len(tickStrings)>0:
        self.ticks=[]
      for tickString in tickStrings:
        self.ticks.append(float(tickString))
    
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
    
    #get plot weightHeight
    self.plotHeightWeights=[]
    for plot in self.plots:
      self.plotHeightWeights.append(plot.weightHeight)
  def load(self,files,options):
    """Loads the values needed for the x-axis points from the fileData argument
    
    """
    
    #load plots
    for plot in self.plots:
      plot.load(files,options)
    
    #set plot limits
    for plot in self.plots:
      plot.setLimits(self.limits)
class DataSet:
  """Holds all the information for a single dataSet
  
  which includes the baseFileName of the dataset, the range of the dataSet
  (start-end), the times and phases of the files within the range of the 
  dataSet, and the plots made from the dataSet.
  
  """
  
  def __init__(self,element,options):
    """Initilizes the dataSet

    by setting baseFileName, start, end, and intilizing plots from an xml 
    element
    
    """
    
    #set some initial values
    self.axes=[]
    self.files={}
    
    #check for unrecognized attributes
    recogAttribs=[]
    notRecogAttribs=checkForRecogAttrib(element,recogAttribs)
    if len(notRecogAttribs)>0:
      raise Exception("attributes in dataset "+str(notRecogAttribs)
      +" are not in list of recognized dataset attributes "+str(recogAttribs))
    
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
    """Loads the dataSet
    
    this means that it sets, time, phases, and plots data
    
    """
    
    for axis in self.axes:
      axis.load(self.files,options)
  def getCurve(self,ID):
    """Returns a curve object with the given ID
    
    """
    
    nCount=0
    for axis in self.axes:
      for plot in axis.plots:
        for curve in plot.curves:
          if nCount==ID:
            return curve
          nCount=nCount+1
def plot(dataSets,options,title):
  """Creates the plot
  
  """
  
  import matplotlib
  matplotlib.rc('text', usetex=True)  #use latex for fonts
  matplotlib.rc('font',size=options.fontSize) #set font size
  
  if not options.show and figureNumber==0:
    matplotlib.use("Agg")
  import matplotlib.pyplot as plt
  
  from matplotlib.gridspec import GridSpec
  
  #create a figure object, to be reused for each figure, could be 1 or more
  fig=plt.figure(figsize=(options.figWidth,options.figHeight))
  basicMatPlotColors=['b','g','r','c','m','y','k','w']
  #add title
  if title!="":
    fig.suptitle(title)
  
  #basic plot spacing options
  axisSpacing=options.axisSpacing
  figBottom=options.figBottom
  figTop=options.figTop
  figLeft=options.figLeft
  figRight=options.figRight
  
  #count number of axes in all figures, this number will be the same for all
  #figures.
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
  nDataSetCount=0
  for dataSet in dataSets:
    for axisMine in dataSet.axes:
      gs.append(GridSpec(len(axisMine.plots),1
        ,height_ratios=axisMine.plotHeightWeights))
      bottom=top-heightPlot*(len(axisMine.plots))
      if options.spaceAxisEvenly:
        top=figTop-axisCount*(heightAxis+axisSpacing)
        bottom=figBottom+(nNumAxes-1.0-axisCount)*(heightAxis+axisSpacing)
      gs[axisCount].update(left=figLeft,right=figRight,top=top,bottom=bottom,hspace=0.0)
      
      #for each plot
      nPlotCount=0
      for plot in axisMine.plots:
        
        #create a subplot
        ax.append(plt.subplot(gs[axisCount][nPlotCount,0]))
        
        #set ticks if manual
        if axisMine.ticks!=None:
          ax[nTotalPlotCount-1].xaxis.set_ticks(axisMine.ticks)
          
        #set ticks if manual
        if plot.ticks!=None:
          ax[nTotalPlotCount-1].yaxis.set_ticks(plot.ticks)
        
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
          if curve.color in basicMatPlotColors or isHexColor(curve.color):
            if len(curve.error)==0:
              temp=ax[nTotalPlotCount-1].plot(curve.x,curve.y
                ,linestyle=curve.style
                ,color=curve.color
                ,marker=curve.marker
                ,markerfacecolor=curve.markerfacecolor
                ,markeredgecolor=curve.markeredgecolor
                ,markersize=curve.markersize
                ,linewidth=curve.linewidth)
            else:
              temp=ax[nTotalPlotCount-1].errorbar(curve.x,curve.y,yerr=curve.error
                ,linestyle=curve.style
                ,color=curve.color
                ,marker=curve.marker
                ,markerfacecolor=curve.markerfacecolor
                ,markeredgecolor=curve.markeredgecolor
                ,markersize=curve.markersize
                ,linewidth=curve.linewidth
                ,ecolor=curve.ecolor
                ,elinewidth=curve.elinewidth
                ,capsize=curve.capsize
                )
          else:
            raise Exception("only hex strings are accepted for colors other than "
              +str(basicMatPlotColors)+" string given was \""+curve.color+"\"")
          if curve.label!=None and curve.label!="":
            lines.append(temp[0])
            labels.append(curve.label)
          curveCount=curveCount+1
        ax[nTotalPlotCount-1].set_xlim(axisMine.limits)
        ax[nTotalPlotCount-1].set_ylim(plot.limits)
        
        #for each text in the plot
        for text in plot.texts:
          ax[nTotalPlotCount-1].text(text.x,text.y,text.text)
        
        #set legend
        if len(lines)>0:
          ax[nTotalPlotCount-1].legend(lines,labels,loc=plot.legendloc,numpoints=plot.numpoints)
        
        #remove x-axis labels
        if nPlotCount!=len(axisMine.plots)-1:
          plt.setp(plt.gca(), 'xticklabels', [])
        else:
          ax[nTotalPlotCount-1].set_xlabel(axisMine.xlabel)
        
        #set y-axis labels
        ax[nTotalPlotCount-1].set_ylabel(plot.ylabel)
        
        #remove top and bottom y-axis tic labels
        if plot.ticks==None:
          ax[nTotalPlotCount-1].yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(prune='both'))
        
        #adjust hspace to be 0 so plots for the same x-axis will be tight
        nPlotCount=nPlotCount+1#resets every axis
        nTotalPlotCount=nTotalPlotCount+1#continues to increase for all plots
      
      top=bottom-axisSpacing
      axisCount=axisCount+1
      
    nDataSetCount+=1
  #done making plot
  if options.show:
    print "ploting file to screen, close for next plot"
    plt.show()
  else:
    [path,ext]=os.path.splitext(options.outputFile)
    supportedFileTypes=["png", "pdf", "ps", "eps", "svg"]
    if ext[1:] not in supportedFileTypes:
      print "File type \""+ext[1:]+"\" not suported. Supported types are ",supportedFileTypes," please choose one of those"
      quit()
    print __name__+":"+main.__name__+": saving figure to file \""+options.outputFile+"\" ..."
    fig.savefig(options.outputFile,format=ext[1:],transparent=False,dpi=options.dpi)#save to file
def isFloat(str):
  """Returns true if str can be converted to a float
  
  """
  
  try:
    float(str)
    return True
  except:
    return False
def isHexColor(str):
  """Returns True if str is a hex color else returns False
  
  """
  
  #split up string
  numberSign=str[0]
  redHex=str[1:3]
  greenHex=str[3:5]
  blueHex=str[5:7]
  
  #test for the right format
  if numberSign!="#":
    return False
  try:
    redDec=int("0x"+redHex,0)
    greenDec=int("0x"+greenHex,0)
    blueDec=int("0x"+blueHex,0)
  except ValueError as e:
    return False
  
  if not (0<=redDec and redDec<=255 ):
    return False
  if not (0<=greenDec and redDec<=255 ):
    return False
  if not (0<=blueDec and redDec<=255 ):
    return False
  
  #if we got here it passed all the tests
  return True
def main():
  
  #parse command line options
  (options,args)=parseOptions()
  
  if len(args)==0:
    print "must supply an input file"
    quit()
    
  #get root element
  tree=xml.parse(args[0])
  root=tree.getroot()
  
  #get all figures
  figureElements=root.findall("figure")
  if len(figureElements)==0:
    raise Exception("No \"figure\" nodes found under \"figures\" node! Nothing to be done.")
  global figureNumber
  figureNumber=0
  for figureElement in figureElements:
    
    #set plot title
    title=""
    if figureElement.get("title")!=None and figureElement.get("title")!="":
      title=figureElement.get("title")
    
    #set spacing between axies
    options.axisSpacing=0.05
    if figureElement.get("axisSpacing")!=None and figureElement.get("axisSpacing")!="":
      options.axisSpacing=float(figureElement.get("axisSpacing"))
    
    #set location of top of the plot area
    options.figTop=0.95
    if figureElement.get("figTop")!=None and figureElement.get("figTop")!="":
      options.figTop=float(figureElement.get("figTop"))
    
    #set location of bottom of the plot area
    options.figBottom=0.05
    if figureElement.get("figBottom")!=None and figureElement.get("figBottom")!="":
      options.figBottom=float(figureElement.get("figBottom"))
      
    #set location of left of the plot area
    options.figLeft=0.05
    if figureElement.get("figLeft")!=None and figureElement.get("figLeft")!="":
      options.figLeft=float(figureElement.get("figLeft"))
    
    #set location of left of the plot area
    options.figRight=0.95
    if figureElement.get("figRight")!=None and figureElement.get("figRight")!="":
      options.figRight=float(figureElement.get("figRight"))
    
    #set figure height
    options.figHeight=11.0
    if figureElement.get("figHeight")!=None and figureElement.get("figHeight")!="":
      options.figHeight=float(figureElement.get("figHeight"))
    
    #set figure width
    options.figWidth=8.5
    if figureElement.get("figWidth")!=None and figureElement.get("figWidth")!="":
      options.figWidth=float(figureElement.get("figWidth"))
      
    #set font size
    options.fontSize=12.0
    if figureElement.get("fontSize")!=None and figureElement.get("fontSize")!="":
      options.fontSize=float(figureElement.get("fontSize"))
      
    #set figure dpi
    options.dpi=300
    if figureElement.get("dpi")!=None and figureElement.get("dpi")!="":
      options.dpi=int(figureElement.get("dpi"))
    
    #set plot output
    if figureElement.get("outputfile")!=None and figureElement.get("outputfile")!="":
      options.outputFile=figureElement.get("outputfile")
      [path,ext]=os.path.splitext(options.outputFile)
      supportedFileTypes=["png", "pdf", "ps", "eps", "svg"]
      if ext[1:] not in supportedFileTypes:
        print "File type \""+ext[1:]+"\" not suported. Supported types are ",supportedFileTypes," please choose one of those"
        quit()
    
    #get list of files elements
    dataSetElements=figureElement.findall("dataSet")
    
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
    figureNumber+=1
  return
if __name__ == "__main__":
  main()