#!/usr/bin/env python

#import getopt
import argparse
import sys
import os
from math import *
import glob
import numpy as np
import make_2DSlices
import time
import disect_filename
import xml.etree.ElementTree as xml
import parser
import re
import matplotlib
matplotlib.use('Agg')#prevents a window from poping up
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d

nNumCoords=4

class File2DSlice:
  def load(self,fileName):
    '''sets:
    fileName, file name of the 2D slice
    planeType, type of the 2D slice ("rt","rp", "tp")
    eosFile, file name of the equition of state file, if using a gamma-law gas it is None
    gamma, value of gamma for a gamma-law gass, if using an equation of state table it is None
    coordinateNames, Names of the coordinates
    coordinates, values of the coordinates
    dataNames, names of the data columns
    data, the data columns
    '''
    self.fileName=fileName
    if not os.access(fileName,os.F_OK|os.R_OK):
      print "error opening ",fileName," for reading"
      return
    f=open(fileName,'r')
    
    #set time step index
    temp=fileName[fileName.find("_t")+2:]
    self.index=temp[:temp.find("_2D")]
    
    #set type of plane
    self.planeType=f.readline().strip()
    
    #read time
    line=f.readline()
    words=line.split()
    self.time=float(words[1])
    
    #read gamma/eos
    if words[2]!="0":
      self.eosFile=words[3]
      self.gamma=None
    else:
      self.eosFile=None
      self.gamma=float(words[3])
    
    #read coordinates
    self.coordinateNames=[]
    self.coordinates=[]
    for i in range(nNumCoords):
      line=f.readline()
      words=line.split()
      self.coordinates.append([])
      self.coordinateNames.append(words[0])
      for j in range(1,len(words)):
        self.coordinates[i].append(float(words[j]))
    
    #read in data names
    line=f.readline()
    words=line.split()
    self.dataNames=words
    
    #read in data
    self.data=[]
    line=f.readline()
    nCount=0
    for word in line.split():
      self.data.append([])
      if word!="-":
        self.data[nCount].append(float(word))
      nCount+=1
    for line in f:
      nCount=0
      for word in line.split():
        if word!="-":
          self.data[nCount].append(float(word))
        nCount+=1
def isFloat(testFloat):
  '''try to convert to a float, if an exception is thrown it isn't a float otherwise it is'''
  try:
    float(testFloat)
  except:
    return False
  else:
    return True
def isInt(testInt):
  '''try to convert to a integer, if an exception is thrown it isn't an integer otherwise it is'''
  try:
    int(testInt)
  except:
    return False
  else:
    return True
def main():
  
  #set parser options
  parser=argparse.ArgumentParser(description="Plots 2D slices from xml configuration file "\
    +"specified by \"fileName\". A description of the xml configuration file can be found in the"\
    +" reference file \"plot_2DSlices_reference.xml\" found with this script")
  parser.add_argument('fileName',action="store",type=str,help="Name of the xml configuration file")
  parser.add_argument('-m',action="store_true",default=False,help="Remake 2D slices")
  parser.add_argument('-b',action="store_true",default=False,help="Re-combine binary files, at "\
    +"this option momement doesn't do anything")
  parser.add_argument('-r',action="store_true",default=False,help="Remove distributed binary files")
  parser.add_argument('-c',action="store_true",default=False,help="Show codes as the will be "\
    +"executed, mostly a debugging option")
  
  #parse arguments
  parsed=parser.parse_args()
  
  #get xml settings
  settings=parseXMLFile(parsed.fileName)
  
  #set codes
  setCodes(settings,parsed)
  
  #create plots
  createPlots(settings,parsed)
def parseXMLFile(fileName):
  
  #get root element
  tree=xml.parse(fileName)
  root=tree.getroot()
  settings={}
  
  if root==None:
    print "Require a \"figure\" root node"
  
  #get output file and extension
  if root.get("outputFile")!=None and root.get("outputFile")!="":
    outputFile=root.get("outputFile")
    [path,ext]=os.path.splitext(outputFile)
    supportedFileTypes=["png", "pdf", "ps", "eps", "svg"]
    if ext[1:] not in supportedFileTypes:
      print "File type \""+ext[1:]+"\" not suported. Supported types are ",supportedFileTypes,\
        " please choose one of those"
      quit()
    settings['outputFile']=path
    settings['fileFormat']=ext[1:]
  else:#ues defaults
    settings['outputFile']="2Dslice"
    settings['fileFormat']="png"
  
  #get figure width
  if root.get("width")!=None and root.get("width")!="":
    if isFloat(root.get("width")):
      settings['figWidth']=float(root.get("width"))
    else:
      print "Expecting a float for \"width\" attribute of node \"figure\""
      quit()
  else:
    settings['figWidth']=20.0
    
  #get figure height
  if root.get("height")!=None and root.get("height")!="":
    if isFloat(root.get("height")):
      settings['figHeight']=float(root.get("height"))
    else:
      print "Expecting a float for \"height\" attribute of node \"figure\""
      quit()
  else:
    settings['figHeight']=10.0
    
  #get figure dpi
  if root.get("dpi")!=None and root.get("dpi")!="":
    if isInt(root.get("dpi")):
      settings['figDpi']=int(root.get("dpi"))
    else:
      print "Expecting an integer for \"dpi\" attribute of node \"figure\""
      quit()
  else:
    settings['figDpi']=90
    
  #get starting index for output file names
  if root.get("startIndex")!=None and root.get("startIndex")!="":
    if isInt(root.get("startIndex")):
      settings['startIndex']=int(root.get("startIndex"))
    else:
      print "Expecting an integer for \"startIndex\" attribute of node \"figure\""
      quit()
  else:
    settings['startIndex']=0
  
  #get figure title
  if root.get("title")!=None:
    settings['title']=root.get("title")
  else:
    settings['title']=""
  
  #get baseFileName
  if root.findtext("inputFileName")==None or root.findtext("inputFileName")=="":
    print "Requires an \"inputFileName\" node"
    quit()
  settings['inputFileName']=root.findtext("inputFileName")
  
  # get plane element
  planeElement=root.find("plane")
  if planeElement==None:
    print "Requires a plane node"
    quit()
  
  #get plane type
  allowedPlaneTypes=["rt","rp","tp"]
  if planeElement.get("type") not in allowedPlaneTypes:
    print "Requires a plane \"type\" attribute from one of the following",allowedPlaneTypes
    quit()
  settings['planeType']=planeElement.get("type")
  
  #get plane index
  if not isInt(planeElement.get("index")):
    print "Requires a plane \"index\" attribute that indicates the location of the plane in the "\
      +"model, should be an integer indicating the zone in the direction perpendicular to the"\
      +" plane."
    quit()
  settings['planeIndex']=int(planeElement.get("index"))
  
  #get grid
  if planeElement.get("grid")!=None:
    if planeElement.get("grid").lower() in ["both","minor","major","none"]:
      settings['grid']=planeElement.get("grid").lower()
    else:#use default
      settings['grid']=None
  
  #get x-axis
  xaxisElement=root.find("xaxis")
  if xaxisElement==None:
    print "Requires an \"xaxis\" node"
    quit()
 
 #get xmin
  if xaxisElement.get("min")!=None:
    if not isFloat(xaxisElement.get("min")):
      print "Expecting a float for xaxis \"min\" attribute"
      quit()
    settings['xMin']=float(xaxisElement.get("min"))
  else:#use default
    settings['xMin']=None
  
  #get xmax
  if xaxisElement.get("max")!=None:
    if not isFloat(xaxisElement.get("max")):
      print "Expecting a float for xaxis \"max\" attribute"
      quit()
    settings['xMax']=float(xaxisElement.get("max"))
  else:#use default
    settings['xMax']=None
    
  #get xaxis formula
  if xaxisElement.text==None or xaxisElement.text=="":
    print "Must specify an xaxis varible"
    quit()
  settings['xFormula']=xaxisElement.text
  
  #get xlabel
  if xaxisElement.get("label")!=None:
    settings['xLabel']=xaxisElement.get("label")
  else:#use default
    settings['xLabel']="x axis label"
  
  #get x minorTics
  if xaxisElement.get("minortics")!=None:
    if xaxisElement.get("xminortics").lower() in ["true","1","t","yes","y"]:
      settings['xminortics']=True
  else:#use default
    settings['xminortics']=False
  
  #get y-axis
  yaxisElement=root.find("yaxis")
  if yaxisElement==None:
    print "Requires an \"yaxis\" node"
    quit()
 
 #get ymin
  if yaxisElement.get("min")!=None:
    if not isFloat(yaxisElement.get("min")):
      print "Expecting a float for yaxis \"min\" attribute"
      quit()
    settings['yMin']=float(yaxisElement.get("min"))
  else:#use default
    settings['yMin']=None
  
  #get xmax
  if yaxisElement.get("max")!=None:
    if not isFloat(yaxisElement.get("max")):
      print "Expecting a float for yaxis \"max\" attribute"
      quit()
    settings['yMax']=float(yaxisElement.get("max"))
  else:#use default
    settings['yMax']=None
    
  #get xaxis formula
  if yaxisElement.text==None or yaxisElement.text=="":
    print "Must specify an yaxis varible"
    quit()
  settings['yFormula']=yaxisElement.text
  
  #get ylabel
  if yaxisElement.get("label")!=None:
    settings['yLabel']=yaxisElement.get("label")
  else:#use default
    settings['yLabel']="y axis label"
  
  #get y minorTics
  if yaxisElement.get("minortics")!=None:
    if yaxisElement.get("yminortics").lower() in ["true","1","t","yes","y"]:
      settings['yminortics']=True
  else:#use default
    settings['yminortics']=False
  
  #get scalor
  scalorElement=root.find("scalor")
  if scalorElement==None:
    settings['scalorFormula']=None
    settings['scalorMin']=sys.float_info.min
    settings['scalorMax']=sys.float_info.max
    settings['scalorPallet']="jet"
    settings['scalorLabel']=None
  else:
    
    #get scalorMin
    if scalorElement.get("min")!=None:
      if not isFloat(scalorElement.get("min")):
        print "Expecting a float for scalor \"min\" attribute"
        quit()
      settings['scalorMin']=float(scalorElement.get("min"))
    else:#use default
      settings['scalorMin']=None
    
    #get scalorMax
    if scalorElement.get("max")!=None:
      if not isFloat(scalorElement.get("max")):
        print "Expecting a float for scalor \"max\" attribute"
        quit()
      settings['scalorMax']=float(scalorElement.get("max"))
    else:#use default
      settings['scalorMax']=None
    
    #get scalor formula
    if scalorElement.text==None or scalorElement.text=="":
      print "Must specify an scalor varible"
      quit()
    settings['scalorFormula']=scalorElement.text
    
    #get scalor label
    if scalorElement.get("label")!=None:
      settings['scalorLabel']=scalorElement.get("label")
    else:#use default
      settings['scalorLabel']="scalor label"
    
    #get scalor pallet
    if scalorElement.get("pallet")!=None:
      settings['scalorPallet']=scalorElement.get("pallet")
      if settings['scalorPallet']=="stellar":
        if scalorElement.get("palletFocus")!=None:
          if not isFloat(scalorElement.get("palletFocus")):
            print "Expecting a float for scalor \"palletFocus\" attribute"
            quit()
          settings['palletFocus']=float(scalorElement.get("palletFocus"))
        else:#use default
          print "Scalor pallet \"stellar\" also needs the attribute \"palletFocus\" to be set."\
            " This indicates which value of the scale the pallet should set to white."
    else:#use default
      settings['scalorPallet']="jet"
  
  #read in vectors
  vectorElements=root.findall("vector")
  nCount=0
  settings['numVectors']=len(vectorElements)
  for vectorElement in vectorElements:
    
    vector={}
    
    #get xfrequency
    vector['xfrequency']=1
    if isInt(vectorElement.get("xfrequency")):
      vector['xfrequency']=int(vectorElement.get("xfrequency"))
    
    #get yfrequency
    vector['yfrequency']=1
    if isInt(vectorElement.get("yfrequency")):
      vector['yfrequency']=int(vectorElement.get("yfrequency"))
    
    #get scale
    vector['scale']=1.0
    if isFloat(vectorElement.get("scale")):
      vector['scale']=float(vectorElement.get("scale"))
    
    #get color
    vector['color']='k'
    if vectorElement.get("color")!=None and vectorElement.get("color")!="":
      vector['color']=vectorElement.get("color")
    
    #get label
    labelElement=vectorElement.find("label")
    if labelElement!=None:
      
      #get label xposition
      if labelElement.get("xpos")!=None:
        if isFloat(labelElement.get("xpos")):
          vector['labelXPos']=labelElement.get("xpos")
        else:
          print "Expecting a float for \"xPos\" in vector label of vector",nCount
      else:#default
        vector['labelXPos']=0.1
        
      #get label yposition
      if labelElement.get("ypos")!=None:
        if isFloat(labelElement.get("ypos")):
          vector['labelYPos']=labelElement.get("ypos")
        else:
          print "Expecting a float for \"ypos\" in vector label of vector",nCount
      else:#default
        vector['labelYPos']=0.92
    
      #get text
      if labelElement.text!=None:
        vector['label']=labelElement.text
      else:#use default
        vector['label']=""
    else:
      vector['label']=""
      vector['labelYPos']=0.0
      vector['labelXPos']=0.0
    
    #get xposition formula
    if vectorElement.findtext("xposition")==None or vectorElement.findtext("xposition")=="":
      print "Must specify an xposition varible for vector"
      quit()
    vector['xposition']=vectorElement.findtext("xposition")
      
    #get yposition formula
    if vectorElement.findtext("yposition")==None or vectorElement.findtext("yposition")=="":
      print "Must specify an yposition varible for vector"
      quit()
    vector['yposition']=vectorElement.findtext("yposition")
    
    #get xcomponent formula
    if vectorElement.findtext("xcomponent")==None or vectorElement.findtext("xcomponent")=="":
      print "Must specify an xcomponent varible for vector"
      quit()
    vector['xcomponent']=vectorElement.findtext("xcomponent")
    
    #get ycomponent formula
    if vectorElement.findtext("ycomponent")==None or vectorElement.findtext("ycomponent")=="":
      print "Must specify an ycomponent varible for vector"
      quit()
    vector['ycomponent']=vectorElement.findtext("ycomponent")
    key='vector'+str(nCount)
    settings[key]=vector
    nCount+=1
  return settings
def setCodes(settings,parsed):
  '''
  Checks formulas for calculating coordinates, scalors, and vectors. It makes sure the references to
  intreface varialbes and coordinates are properly averaged so that all data will have the correct
  spacing to match either grid centeres (in the case of scalors) or grid interfaces or centers (in
  the case of vectors). It does this in a way that allows the greatest amount of freedom while still
  allowing as much freedom as possilbe in plotting combinations.
  '''
  
  #set allowed coordinates,vectors, and scalors. This makes sure that the zoning always works out 
  #correctly
  if settings['planeType']=="rt":
    allowedCoordinate1=[0,1]
    allowedCoordinate2=[2]
    allowedVectors1=[4,5]
    allowedVectors2=[6]
    allowedScalors=[7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
  if settings['planeType']=="rp":
    allowedCoordinate1=[0,1]
    allowedCoordinate2=[3]
    allowedVectors1=[4,5]
    allowedVectors2=[7]
    allowedScalors=[6,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
  if settings['planeType']=="tp":
    allowedCoordinate1=[2]
    allowedCoordinate2=[3]
    allowedVectors1=[6]
    allowedVectors2=[7]
    allowedScalors=[4,5,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
  
  #get all coordinate references
  coordRefs=[]
  xCoordRefs=[]
  yCoordRefs=[]
  [temp,tokens]=splitStrAtList(settings['xFormula'],['$','#'])
  for i in range(1,len(temp)):
    xCoordRefs.append(splitFirstInt(temp[i])[0])
  [temp,tokens]=splitStrAtList(settings['yFormula'],['$','#'])
  for i in range(1,len(temp)):
    yCoordRefs.append(splitFirstInt(temp[i])[0])
  coordRefs=xCoordRefs+yCoordRefs
  
  #check to make sure we have one coordinate for each direction at least
  coordinate1=None
  coordinate2=None
  for coordinate in allowedCoordinate1:
    if coordinate in coordRefs:
      coordinate1=coordinate
  for coordinate in allowedCoordinate2:
    if coordinate in coordRefs:
      coordinate2=coordinate
  if coordinate1==None:
    print "Need to have a refernece to one of the following coordinates", allowedCoordinate1\
      ," in one or both of the axis nodes"
    quit()
  if coordinate2==None:
    print "Need to have a refernece to one of the following coordinates", allowedCoordinate2\
      ," in one or both of the axis nodes"
    quit()
  settings['coordinate1']=coordinate1
  settings['coordinate2']=coordinate2
  
  #set xCode
  [temp,tokens]=splitStrAtList(settings['xFormula'],['$','#'])
  string=temp[0]
  gotCoordinateInX=False
  for i in range(1,len(temp)):
    splits=splitFirstInt(temp[i])
    if splits[0]==coordinate1:
      if tokens[i]=="#":#add as average
        print "In caculation of x-coordinate can not average the independent coordinate \""\
          +str(splits[0])+"\" to a zone center, must be at interfaces."
        quit()
      gotCoordinateInX=True
      string+="slice2D.coordinates["+str(coordinate1)+"][i]"
    elif splits[0]==coordinate2:
      if tokens[i]=="#":#add as average
        print "In caculation of x-coordinate can not average the independent coordinate \""\
          +str(splits[0])+"\" to a zone center, must be at interfaces."
        quit()
      gotCoordinateInX=True
      string+="slice2D.coordinates["+str(coordinate2)+"][j]"
    elif not splits[0]<nNumCoords:
      print "Need to pick from the coordinate rows (0-"+str(nNumCoords)+")"
      quit()
    elif splits[0] not in allowedCoordinate1 or splits[0] not in allowedCoordinate2:#assume it is a constant
      if tokens[i]=="$":#add as is
        print "The coordinate reference",splits[0]," in yaxis is not in the plane, and thus it"\
          +" must be averaged by replacing \"$\" with a \"#\"."
        quit()
      if tokens[i]=="#":#add as average
        string+="(slice2D.coordinates["+str(splits[0])+"][0]+slice2D.coordinates["\
        +str(splits[0])+"][1])*0.5"
    string+=splits[1]
  if not gotCoordinateInX:
    print "Need a reference to at least one independent coordinate for x-axis formula, should be "\
      +"one or more of ",allowedCoordinate1+allowedCoordinate2
    quit()
  if parsed.c:
    print "xCode=",string
  settings['xCode']=parser.expr(string).compile()
  
  #set yCode
  [temp,tokens]=splitStrAtList(settings['yFormula'],['$','#'])
  string=temp[0]
  gotCoordinateInY=False
  for i in range(1,len(temp)):
    splits=splitFirstInt(temp[i])
    if splits[0]==coordinate1:#if it is coordinate 1
      if tokens[i]=="#":#add as average
        print "In caculation of x-coordinate can not average the independent coordinate \""\
          +str(splits[0])+"\" to a zone center, must be at interfaces."
        quit()
      gotCoordinateInY=True
      string+="slice2D.coordinates["+str(coordinate1)+"][i]"
    elif splits[0]==coordinate2:#if it is coordinate 2
      if tokens[i]=="#":#add as average
        print "In caculation of x-coordinate can not average the independent coordinate \""\
          +str(splits[0])+"\" to a zone center, must be at interfaces."
        quit()
      gotCoordinateInY=True 
      string+="slice2D.coordinates["+str(coordinate2)+"][j]"
    elif splits[0]>=nNumCoords:#if it is not a coordinate
      print "Need to pick from the coordinate rows (0-"+str(nNumCoords)+")"
      quit()
    elif splits[0] not in allowedCoordinate1 or splits[0] not in allowedCoordinate2:#assume it is a constant
      if tokens[i]=="$":#must be averaged
        print "The coordinate reference",splits[0]," in yaxis is not in the plane, and thus it"\
          +" must be averaged by replacing \"$\" with a \"#\"."
        quit()
      if tokens[i]=="#":#add as average
        string+="(slice2D.coordinates["+str(splits[0])+"][0]+slice2D.coordinates["\
        +str(splits[0])+"][1])*0.5"
    string+=splits[1]
  if not gotCoordinateInY:
    print "Need a reference to at least one independent coordinate for y-axis formula, should be "\
      +"one or more of ",allowedCoordinate1+allowedCoordinate2
    quit()
  if parsed.c:
    print "yCode=",string
  settings['yCode']=parser.expr(string).compile()
  
  #set scalorCode
  if settings['scalorFormula']!=None:
    [temp,tokens]=splitStrAtList(settings['scalorFormula'],['$','#'])
    string=temp[0]
    for i in range(1,len(temp)):
      splits=splitFirstInt(temp[i])
      if splits[0]>=nNumCoords:
        ref=splits[0]-nNumCoords
        if splits[0] in allowedVectors1:#it is a vector at coordinate1 interfaces
          if tokens[i]=="#":
            string+="(slice2D.data["+str(ref)\
              +"][j+i*shape[1]]+slice2D.data["+str(ref)+"][j+(i+1)*shape[1]])*0.5"
          else:
            print "The reference to vector",splits[0],"should be averaged by using \"#\" in place "\
              +"of \"$\" so that it is at zone centers where other scalors are defined."
            quit()
        elif splits[0] in allowedVectors2:#it is a vector at coordinate2 interfaces
          if tokens[i]=="#":
            string+="(slice2D.data["+str(ref)\
              +"][j+i*shape[1]]+slice2D.data["+str(ref)+"][(j+1)+i*shape[1]])*0.5"
          else:
            print "The reference to vector",splits[0],"should be averaged by using \"#\" in place "\
              +"of \"$\" so that it is at zone centers where other scalors are defined."
            quit()
        else:#it is a scalor
          string+="slice2D.data["+str(ref)+"][j+i*shape[1]]"
      else:#it is a coordinate
        if splits[0] in allowedCoordinate1:#it is a coordinate at coordinate1 interfaces
          if tokens[i]=="#":
            string+="(slice2D.coordinates["+str(splits[0])+"][i]"\
              +"+slice2D.coordinates["+str(splits[0])+"][i+1])*0.5"
          else:
            print "The reference to coordinate",splits[0],"in scalor formula should be averaged"\
              +" by using \"#\" in place of \"$\" so that it is at zone centers where other"\
              +" scalors are defined."
            quit()
        if splits[0] in allowedCoordinate2:#it is a coordinate at coordinate2 interfaces
          if tokens[i]=="#":
            string+="(slice2D.coordinates["+str(splits[0])+"][j]"\
              +"+slice2D.coordinates["+str(splits[0])+"][j+1])*0.5"
          else:
            print "The reference to coordinate",splits[0],"in scalor formula should be averaged"\
              +" by using \"#\" in place of \"$\" so that it is at zone centers where other"\
              +" scalors are defined."
            quit()
      string+=splits[1]
    if parsed.c:
      print "scalorCode=",string
    settings['scalorCode']=parser.expr(string).compile()
  
  #for each vector, make code
  for j in range(settings['numVectors']):
    
    #make xposition code
    key='vector'+str(j)
    [temp,tokens]=splitStrAtList(settings[key]['xposition'],['$','#'])
    string=temp[0]
    coordinate1Averaged=None
    coordinate2Averaged=None
    for i in range(1,len(temp)):
      splits=splitFirstInt(temp[i])
      if splits[0]==coordinate1:
        if tokens[i]=="#":#add as average
          coordinate1Averaged=True
          string+="(slice2D.coordinates["+str(coordinate1)+"][i]"\
            +"+slice2D.coordinates["+str(coordinate1)+"][i+1])*0.5"
        else:
          coordinate1Averaged=False
          string+="slice2D.coordinates["+str(coordinate1)+"][i]"
      elif splits[0]==coordinate2:
        if tokens[i]=="#":#add as average
          coordinate2Averaged=True
          string+="(slice2D.coordinates["+str(coordinate2)+"][j]"\
            +"+slice2D.coordinates["+str(coordinate2)+"][j+1])*0.5"
        else:
          coordinate2Averaged=False
          string+="slice2D.coordinates["+str(coordinate2)+"][j]"
      elif splits[0]>=nNumCoords:
        print "In vector",j,"need to pick from the coordinate rows (0-"+str(nNumCoords)+")"
        quit()
      elif splits[0] not in allowedCoordinate1 or splits[0] not in allowedCoordinate2:#assume it is a constant
        if tokens[i]=="$":#add as is
          print "The coordinate reference",splits[0]," in xcomponent of vector",j," is not in the "\
            +"plane, and thus it must be averaged by replacing \"$\" with a \"#\"."
          quit()
        if tokens[i]=="#":#add as average
          string+="(slice2D.coordinates["+str(splits[0])+"][0]+slice2D.coordinates["\
          +str(splits[0])+"][1])*0.5"
      string+=splits[1]
    if coordinate1Averaged==None and coordinate2Averaged==None:
      print "In vector",j,"need a reference to at least one independent coordinate "\
        +"for xposition formula, should be one or more of ",allowedCoordinate1+allowedCoordinate2
      quit()
    if parsed.c:
      print "vector ",j,"'s xCode=",string
    settings[key]['xCode']=parser.expr(string).compile()
    
    #make yposition code
    [temp,tokens]=splitStrAtList(settings[key]['yposition'],['$','#'])
    string=temp[0]
    for i in range(1,len(temp)):
      splits=splitFirstInt(temp[i])
      if splits[0]==coordinate1:
        if tokens[i]=="#":#add as average
          if coordinate1Averaged!=None:#means we have gotten coordinate 1 already
            if not coordinate1Averaged:# check for consistancy of averaging
              print "In vector",j,"coordinate \""+str(coordinate1)\
                +"\" not averaged in xposition, must also not be averaged in yposition"
              quit()
          else:
            coordinate1Averaged=True
          string+="(slice2D.coordinates["+str(settings['coordinate1'])+"][i]"\
            +"+slice2D.coordinates["+str(settings['coordinate1'])+"][i+1])*0.5"
        else:#add as interface
          if coordinate1Averaged!=None:#means we have gotten coordinate 1 already
            if coordinate1Averaged:#check for consistancy of averaging
              print "In vector ",j," coordinate \""+str(coordinate1)\
                +"\" averaged in xposition, must also be averaged in yposition"
              quit()
          else:
            coordinate1Averaged=False
          string+="slice2D.coordinates["+str(settings['coordinate1'])+"][i]"
      elif splits[0]==coordinate2:
        if tokens[i]=="#":#add as average
          if coordinate2Averaged!=None:#means we have gotten coordinate 2 already
            if not coordinate2Averaged:#check for consistancy of averaging
              print "In vector",j,"coordinate \""+str(coordinate2)\
                +"\" not averaged in xposition, must also not be averaged in yposition"
              quit()
          else:
            coordinate2Averaged=True
          string+="(slice2D.coordinates["+str(settings['coordinate2'])+"][j]"\
            +"+slice2D.coordinates["+str(settings['coordinate2'])+"][j+1])*0.5"
        else:
          if coordinate2Averaged!=None:#means we have gotten coordinate 2 already
            if coordinate2Averaged:#check for consistancy of averaging
              print "In vector",j,"Coordinate \""+str(coordinate2)\
                +"\" averaged in xposition, must also be averaged in yposition"
              quit()
          else:
            coordinate2Averaged=False
          string+="slice2D.coordinates["+str(settings['coordinate2'])+"][j]"
      elif splits[0]>=nNumCoords:
        print "In vector",j,"need to pick from the coordinate rows (0-"+str(nNumCoords)+")"
        quit()
      elif splits[0] not in allowedCoordinate1 or splits[0] not in allowedCoordinate2:#assume it is a constant
        if tokens[i]=="$":#add as is
          print "The coordinate reference",splits[0]," in ycomponent of vector",j," is not in the "\
            +"plane, and thus it must be averaged by replacing \"$\" with a \"#\"."
          quit()
        if tokens[i]=="#":#add as average
          string+="(slice2D.coordinates["+str(splits[0])+"][0]+slice2D.coordinates["\
          +str(splits[0])+"][1])*0.5"
      string+=splits[1]
    if coordinate1Averaged==None or coordinate2Averaged==None:
      print "In vector",j,"need a reference to at least one independent coordinate for "\
        +"yposition formula, should be one or more of ",allowedCoordinate1+allowedCoordinate2
      quit()
    if parsed.c:
      print "vector ",j,"'s yCode=",string
    settings[key]['yCode']=parser.expr(string).compile()
    settings[key]['coordinate1Averaged']=coordinate1Averaged
    settings[key]['coordinate2Averaged']=coordinate2Averaged
    
    #make xcomponent code
    [temp,tokens]=splitStrAtList(settings[key]['xcomponent'],['$','#'])
    string=temp[0]
    for i in range(1,len(temp)):
      splits=splitFirstInt(temp[i])
      if splits[0]>=nNumCoords:
        ref=splits[0]-nNumCoords
        if splits[0] in allowedVectors1:#it is a vector at coordinate1 interfaces
          if tokens[i]=="#":
            if not coordinate1Averaged:#coordinate not averaged so vector cannot be averaged
              print "Vector",splits[0],"in xcomponent can not be averaged unless the xposition"\
                " cooridnate",coordinate1,"is also averaged."
              quit()
            elif not coordinate2Averaged:
              print "Vector",splits[0],"in xcomponent must have yposition"\
                " cooridnate",coordinate2," averaged."
              quit()
            else:#coordinate is averaged so the vector can be also
              string+="(slice2D.data["+str(ref)\
                +"][j+i*shape1[1]]+slice2D.data["+str(ref)+"][j+(i+1)*shape1[1]])*0.5"
          else:
            if coordinate1Averaged:
              print "Vector",splits[0],"in xcomponent must be averaged since the coorisponding"\
                " cooridnate",coordinate1,"is also averaged."
              quit()
            elif not coordinate2Averaged:
              print "Vector",splits[0],"in xcomponent must have yposition"\
                " cooridnate",coordinate2," averaged."
              quit()
            if not coordinate2Averaged:#coordinate not averaged so vector cannot be averaged
              print "Vector",splits[0],"in xcomponent must have yposition"\
                " cooridnate",coordinate2," averaged."
              quit()
            else:
              string+="slice2D.data["+str(ref)+"][j+i*shape1[1]]"
        elif splits[0] in allowedVectors2:#it is a vector at coordinate2 interfaces
          if tokens[i]=="#":
            if not coordinate2Averaged:#coordinate not averaged so vector cannot be averaged
              print "Vector",splits[0],"in xcomponent can not be averaged unless the coorisponding"\
                " cooridnate",coordinate2,"is also averaged."
              quit()
            elif not coordinate1Averaged:
              print "Vector",splits[0],"in xcomponent must have xposition"\
                " cooridnate",coordinate1," averaged."
              quit()
            else:#coordinate is averaged so the vector can be also
              string+="(slice2D.data["+str(ref)\
                +"][j+i*shape2[1]]+slice2D.data["+str(ref)+"][(j+1)+i*shape2[1]])*0.5"
          else:
            if coordinate2Averaged:
              print "Vector",splits[0],"in xcomponent must be averaged since the coorisponding"\
                " cooridnate",coordinate2,"is also averaged."
              quit()
            elif not coordinate1Averaged:
              print "Vector",splits[0],"in xcomponent must have xposition"\
                " cooridnate",coordinate1," averaged."
              quit()
            else:
              string+="slice2D.data["+str(ref)+"][j+i*shape2[1]]"
      else:#it is a coordinate
        if splits[0] in allowedCoordinate1:#it is a coordinate at coordinate1 interfaces
          if tokens[i]=="#":
            if not coordinate1Averaged:
              print "Coordinate ",splits[0],"in xcomponent can not be averaged unless the"\
                +" coorisponding cooridnate",coordinate1,"is also averaged."
              quit()
            else:
              string+="(slice2D.coordinates["+str(splits[0])+"][i]"\
                +"+slice2D.coordinates["+str(splits[0])+"][i+1])*0.5"
          else:
            if coordinate1Averaged:
              print "Coordinate ",splits[0],"in xcomponent must be averaged if the"\
                +" coorisponding cooridnate",coordinate1,"is also averaged."
              quit()
            else:
              string+="slice2D.coordinates["+str(splits[0])+"][i]"
        elif splits[0] in allowedCoordinate2:#it is a coordinate at coordinate2 interfaces
          if tokens[i]=="#":
            if not coordinate2Averaged:
              print "Coordinate ",splits[0],"in xcomponent can not be averaged unless the"\
                +" coorisponding cooridnate",coordinate2,"is also averaged."
              quit()
            else:
              string+="(slice2D.coordinates["+str(splits[0])+"][j]"\
                +"+slice2D.coordinates["+str(splits[0])+"][j+1])*0.5"
          else:
            if coordinate2Averaged:
              print "Coordinate ",splits[0],"in xcomponent must be averaged if the"\
                +" coorisponding cooridnate",coordinate2,"is also averaged."
              quit()
            else:
              string+="slice2D.coordinates["+str(splits[0])+"][j]"
        else:
          if tokens[i]!="#":
            print "Coordinate ",splits[0],"in xcomponent must be averaged since it isn't an "\
              +"allowed coordinate (",allowedCoordinate1+allowedCoordinate2,") for this plane type."
            quit()
          else:
            string+="(slice2D.coordinates["+str(splits[0])+"][0]"\
              +"+slice2D.coordinates["+str(splits[0])+"][1])*0.5"
      string+=splits[1]
    if parsed.c:
      print "vector ",j,"'s uCode=",string
    settings[key]['uCode']=parser.expr(string).compile()
    
    #make ycomponent code
    [temp,tokens]=splitStrAtList(settings[key]['ycomponent'],['$','#'])
    string=temp[0]
    for i in range(1,len(temp)):
      splits=splitFirstInt(temp[i])
      if splits[0]>=nNumCoords:
        ref=splits[0]-nNumCoords
        if splits[0] in allowedVectors1:#it is a vector at coordinate1 interfaces
          if tokens[i]=="#":
            if not coordinate1Averaged:#coordinate not averaged so vector cannot be averaged
              print "Vector",splits[0],"in ycomponent can not be averaged unless the coorisponding"\
                " cooridnate",coordinate1,"is also averaged."
              quit()
            elif not coordinate2Averaged:
              print "Vector",splits[0],"in ycomponent must have yposition"\
                " cooridnate",coordinate2," averaged."
              quit()
            else:#coordinate is averaged so the vector can be also
              string+="(slice2D.data["+str(ref)\
                +"][j+i*shape1[1]]+slice2D.data["+str(ref)+"][j+(i+1)*shape1[1]])*0.5"
          else:
            if coordinate1Averaged:
              print "Vector",splits[0],"in ycomponent must be averaged since the coorisponding"\
                " cooridnate",coordinate1,"is also averaged."
              quit()
            elif not coordinate2Averaged:
              print "Vector",splits[0],"in ycomponent must have yposition"\
                " cooridnate",coordinate2," averaged."
              quit()
            else:
              string+="slice2D.data["+str(ref)+"][j+i*shape1[1]]"
        elif splits[0] in allowedVectors2:#it is a vector at coordinate2 interfaces
          if tokens[i]=="#":
            if not coordinate2Averaged:#coordinate not averaged so vector cannot be averaged
              print "Vector",splits[0],"in ycomponent can not be averaged unless the coorisponding"\
                " cooridnate",coordinate2,"is also averaged."
              quit()
            elif not coordinate1Averaged:
              print "Vector",splits[0],"in ycomponent must have xposition"\
                " cooridnate",coordinate1," averaged."
              quit()
            else:#coordinate is averaged so the vector can be also
              string+="(slice2D.data["+str(ref)\
                +"][j+i*shape2[1]]+slice2D.data["+str(ref)+"][(j+1)+i*shape2[1]])*0.5"
          else:
            if coordinate2Averaged:
              print "Vector",splits[0],"in ycomponent must be averaged since the coorisponding"\
                " cooridnate",coordinate2,"is also averaged."
              quit()
            elif not coordinate1Averaged:
              print "Vector",splits[0],"in ycomponent must have xposition"\
                " cooridnate",coordinate1," averaged."
              quit()
            else:
              string+="slice2D.data["+str(ref)+"][j+i*shape2[1]]"
      else:#it is a coordinate
        if splits[0] in allowedCoordinate1:#it is a coordinate at coordinate1 interfaces
          if tokens[i]=="#":
            if not coordinate1Averaged:
              print "Coordinate ",splits[0],"in ycomponent can not be averaged unless the"\
                +" coorisponding cooridnate",coordinate1,"is also averaged."
              quit()
            else:
              string+="(slice2D.coordinates["+str(splits[0])+"][i]"\
                +"+slice2D.coordinates["+str(splits[0])+"][i+1])*0.5"
          else:
            if coordinate1Averaged:
              print "Coordinate ",splits[0],"in ycomponent must be averaged if the"\
                +" coorisponding cooridnate",coordinate1,"is also averaged."
              quit()
            else:
              string+="slice2D.coordinates["+str(splits[0])+"][i]"
        elif splits[0] in allowedCoordinate2:#it is a coordinate at coordinate2 interfaces
          if tokens[i]=="#":
            if not coordinate2Averaged:
              print "Coordinate ",splits[0],"in ycomponent can not be averaged unless the"\
                +" coorisponding cooridnate",coordinate2,"is also averaged."
              quit()
            else:
              string+="(slice2D.coordinates["+str(splits[0])+"][j]"\
                +"+slice2D.coordinates["+str(splits[0])+"][j+1])*0.5"
          else:
            if coordinate2Averaged:
              print "Coordinate ",splits[0],"in ycomponent must be averaged if the"\
                +" coorisponding cooridnate",coordinate2,"is also averaged."
              quit()
            else:
              string+="slice2D.coordinates["+str(splits[0])+"][j]"
        else:
          if tokens[i]!="#":
            print "Coordinate ",splits[0],"in xcomponent must be averaged since it isn't an "\
              +"allowed coordinate (",allowedCoordinate1+allowedCoordinate2,") for this plane type."
            quit()
          else:
            string+="(slice2D.coordinates["+str(splits[0])+"][0]"\
              +"+slice2D.coordinates["+str(splits[0])+"][1])*0.5"
      string+=splits[1]
    if parsed.c:
      print "vector ",j,"'s vCode=",string
    settings[key]['vCode']=parser.expr(string).compile()
def createPlots(settings,parsed):
  
  #get base file name
  [start,end,baseFileName]=disect_filename.disectFileName(settings['inputFileName'])
  
  #make sure that all combined binary files have 2D slices made
  if settings['planeType']=="rt":
    nPlaneID=0
    planeID="k"
  if settings['planeType']=="tp":
    nPlaneID=1
    planeID="i"
  if settings['planeType']=="rp":
    nPlaneID=2
    planeID="j"
  make_2DSlices.make_2DSlices(not(parsed.r),settings['inputFileName'],nPlaneID,settings['planeIndex']
    ,parsed.m)
    
  #get and sort files
  extension="_2D"+planeID+"="+str(settings['planeIndex'])+".txt"
  filesExistSlices=glob.glob(baseFileName+"*"+extension)
  files=[]
  for file in filesExistSlices:
    intOfFile=int(file[len(baseFileName):len(file)-len(extension)])
    if intOfFile>=start and intOfFile<end:
      files.append(file)
  files.sort()
  
  if len(files)==0:
    print "no files found in range"
    quit()
  
  #set stellar pallet only once if both max/min are set, otherwise will have to reset for each plot
  if settings['scalorPallet']=="stellar" and (settings['scalorMin']!=None and settings['scalorMax']!=None):
    setStellarPallet(settings['scalorMax'],settings['scalorMin'],settings['palletFocus'])
  
  nCount=settings['startIndex']
  for file in files:
    fig=plt.figure(figsize=(settings['figWidth'],settings['figHeight']),frameon=False)
    ax=fig.add_subplot(111)
    plot_plane(file,nCount,fig,ax,settings)
    nCount+=1
    plt.close(fig)
def splitFirstInt(str):
  '''Returns the integer which starts at the very first character of str, and drops everything else 
  from str'''
  
  strInt=''
  i=0
  while str[i].isdigit() or (str[i]=='-' and i==0):
    strInt+=str[i]
    i=i+1
    if i>=len(str):#check that we haven't passed the end of the string
      break
  return [int(strInt),str[i:]]
def plot_plane(fileName,nCount,fig,ax,settings):
  
  #load data
  slice2D=File2DSlice()
  slice2D.load(fileName)
  
  #plot scalor data if it is set
  if settings['scalorFormula']!=None:
    [X,Y,C,max,min]=makeScalorPlotData(slice2D,settings)
    cmap=settings['scalorPallet']
    if cmap==None:
      cmap="jet"
    scalerMap=ax.pcolor(X,Y,C,cmap=cmap,edgecolors='none',antialiased=False
      ,vmin=min,vmax=max)
    cb=fig.colorbar(scalerMap,pad=0.01,fraction=0.07)
    cb.set_label(settings['scalorLabel'])
    ax.set_xlabel(settings['xLabel'])
    ax.set_ylabel(settings['yLabel'])
    
  #plot vectors
  for i in range(settings['numVectors']):
    key='vector'+str(i)
    [X,Y,U,V,scale,maxMag,diag]=makeVectorPlotData(slice2D,settings,i)
    scaleString=format(maxMag,"0.1e")
    scaleNumber=float(scaleString)
    label=settings[key]['label']
    labelScale=format(scaleNumber,"0.1e")
    
    label=label.replace("\scale",labelScale)
    uvec=ax.quiver(X,Y,U,V,pivot='middle',scale=scale,units='xy',angles='xy'
      ,width=diag*0.03*settings[key]['scale'],color=settings[key]['color'])
    qk=ax.quiverkey(uvec,settings[key]['labelXPos'],settings[key]['labelYPos']
      ,scaleNumber,label)
  
  #set limits
  ax.set_xlim([settings['xMin'],settings['xMax']])
  ax.set_ylim([settings['yMin'],settings['yMax']])
  
  #set title
  title=settings['title']
  timeStr=format(slice2D.time,"0.4e")
  title=title.replace("\\time",timeStr)
  title=title.replace("\index",str(slice2D.index))
  fig.suptitle(title)
  
  #save figure
  fname=settings['outputFile']+"_"+str(nCount)+"."+settings['fileFormat']
  print "saving figure \""+fname+"\" made from 2D slice \""+fileName+"\" ..."
  fig.savefig(fname,format=settings['fileFormat'],transparent=False,dpi=settings['figDpi'])
def makeScalorPlotData(slice2D,settings):
  
  #calculate grid coordinates
  shape=(len(slice2D.coordinates[settings['coordinate1']])
    ,len(slice2D.coordinates[settings['coordinate2']]))
  X=np.zeros(shape)
  Y=np.zeros(shape)
  for i in range(shape[0]):
    for j in range(shape[1]):
      X[i][j]=eval(settings['xCode'])
      Y[i][j]=eval(settings['yCode'])
  
  #calculate scalor data
  shape=((shape[0]-1),(shape[1]-1))
  C=np.zeros(shape)
  max=-1.0*sys.float_info.max
  min=sys.float_info.max
  
  #if max/min not set, set them to something reasonable
  xMin=settings['xMin']
  xMax=settings['xMax']
  yMin=settings['yMin']
  yMax=settings['yMax']
  if xMin==None:
    xMin=-1.0*sys.float_info.max
  if xMax==None:
    xMax=sys.float_info.max
  if yMin==None:
    yMin=-1.0*sys.float_info.max
  if yMax==None:
    yMax=sys.float_info.max
  scalorsInRange=False
  for i in range(shape[0]):
    for j in range(shape[1]):
      C[i][j]=eval(settings['scalorCode'])
      if xMin<=X[i][j]<=xMax and yMin<=Y[i][j]<=yMax:
        scalorsInRange=True
        if C[i][j]>max:
          max=C[i][j]
        if C[i][j]<min:
          min=C[i][j]
  if not scalorsInRange:
    print "No scalor data found in spedified xrange [",xMin,":",xMax,"] and yrange [",yMin,":",yMax,"]"
    quit()
  
  #allow user setting values to overide these values
  if settings['scalorMin']!=None:
    min=settings['scalorMin']
  if settings['scalorMax']!=None:
    max=settings['scalorMax']
  
  #if both scalorMin and scalorMax are set, then can just define the colormap once
  if settings['scalorPallet']=="stellar":
    setStellarPallet(max,min,settings['palletFocus'])
  return [X,Y,C,max,min]
def setStellarPallet(max,min,focus):
  
  range=max-min
  cfocus=(focus-min)/range
  cdict = {
    'red':  (
    (0.0,  1.0, 1.0),      #1
    (cfocus*0.5,  1.0, 1.0), #2
    (cfocus,  1.0, 1.0),     #3
    (cfocus*1.5,  0.0, 0.0), #4
    (1,  0.0, 0.0),        #5
    ),
    'green':(
    (0.0,  0.0, 0.0),      #1
    (cfocus*0.5,  1.0, 1.0), #2
    (cfocus,  1.0, 1.0),     #3
    (cfocus*1.5,  1.0, 1.0), #4
    (1.0,  0.0, 0.0)       #5
    ),
    'blue': (
    (0.0,  0.0, 0.0),      #1
    (cfocus*0.5,  0.0, 0.0), #2
    (cfocus,  1.0, 1.0),     #3
    (cfocus*1.5,  1.0, 1.0), #4
    (1.0,  1.0, 1.0)       #5
    ),}
  stellar = matplotlib.colors.LinearSegmentedColormap('stellar', cdict)
  plt.register_cmap(cmap=stellar)
def makeVectorPlotData(slice2D,settings,vectorNumber):
  
  #calculate grid coordinates
  key='vector'+str(vectorNumber)
  xSize=len(slice2D.coordinates[settings['coordinate1']])
  ySize=len(slice2D.coordinates[settings['coordinate2']])
  shape1=(xSize,ySize-1)
  shape2=(xSize-1,ySize)
  if settings[key]['coordinate1Averaged']:
    xSize-=1
  if settings[key]['coordinate2Averaged']:
    ySize-=1
  shapeSmaller=(xSize,ySize)
  X=np.zeros(shapeSmaller)
  Y=np.zeros(shapeSmaller)
  U=np.zeros(shapeSmaller)
  V=np.zeros(shapeSmaller)
  scale=1.0
  diag=sys.float_info.max
  mag=-1.0*sys.float_info.max
  
  #calcuate position and components
  for i in range(0,shapeSmaller[0],settings[key]['xfrequency']):
    for j in range(0,shapeSmaller[1],settings[key]['yfrequency']):
      X[i][j]=eval(settings[key]['xCode'])
      Y[i][j]=eval(settings[key]['yCode'])
      U[i][j]=eval(settings[key]['uCode'])
      V[i][j]=eval(settings[key]['vCode'])
  #if max/min not set, set them to something reasonable
  xMin=settings['xMin']
  xMax=settings['xMax']
  yMin=settings['yMin']
  yMax=settings['yMax']
  if xMin==None:
    xMin=-1.0*sys.float_info.max
  if xMax==None:
    xMax=sys.float_info.max
  if yMin==None:
    yMin=-1.0*sys.float_info.max
  if yMax==None:
    yMax=sys.float_info.max
  
  #find smallest diagonal
  bVectorsInRange=False
  for i in range(settings[key]['xfrequency'],shapeSmaller[0],settings[key]['xfrequency']):
    for j in range(settings[key]['yfrequency'],shapeSmaller[1],settings[key]['yfrequency']):
      if xMin<=X[i][j]<=xMax and yMin<=Y[i][j]<=yMax:
        bVectorsInRange=True
        delX0=X[i][j]-X[i-settings[key]['xfrequency']][j]
        delX1=X[i][j]-X[i][j-settings[key]['yfrequency']]
        delY0=Y[i][j]-Y[i][j-settings[key]['yfrequency']]
        delY1=Y[i][j]-Y[i-settings[key]['xfrequency']][j]
        delX0Sq=delX0**2
        delX1Sq=delX1**2
        delY0Sq=delY0**2
        delY1Sq=delY1**2
        delXSq=delX0Sq
        delYSq=delY0Sq
        if delX1Sq>delXSq:
          delXSq=delX1Sq
        if delY1Sq>delYSq:
          delYSq=delY1Sq
        temp=sqrt(delXSq+delYSq)
        if temp<diag:
          diag=temp
  if not bVectorsInRange:
    print "No vector data found in spedified xrange [",xMin,":",xMax,"] and yrange [",yMin,":",yMax,"]"
    quit()
  
  #find largest magnitude
  for i in range(0,shapeSmaller[0],settings[key]['xfrequency']):
    for j in range(0,shapeSmaller[1],settings[key]['yfrequency']):
      if xMin<=X[i][j]<=xMax and yMin<=Y[i][j]<=yMax:
        temp=sqrt((U[i][j]**2+V[i][j]**2))
        if temp>mag:
          mag=temp
  if mag>0.0 and diag!=0.0 and settings[key]['scale']!=0.0:
    scale=mag/diag/settings[key]['scale']
  else:
    scale=1.0
    mag=0.0
  return [X,Y,U,V,scale,mag,diag]
def splitStrAtList(strToSplit,splitList):
  splits=[]
  splitTokens=[None]
  indexOfLastSplit=0
  for i in range(len(strToSplit)):
    for j in range(len(splitList)):
      if strToSplit[i]==splitList[j]:
        splits.append(strToSplit[indexOfLastSplit:i])
        splitTokens.append(splitList[j])
        indexOfLastSplit=i+1
  splits.append(strToSplit[indexOfLastSplit:])
  return [splits,splitTokens]
if __name__ == "__main__":
  main()
  