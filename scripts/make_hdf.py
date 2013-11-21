#!/usr/bin/env python
import sys
import paths
sys.path.append(paths.srcPath+"/pythonextensions/lib/python/")
import hdf
import dump
import optparse as op
import xml.etree.ElementTree as xml
import parse_formula
import warnings
import mywarnings#sets a custom warning print out
import combine_bins
import disect_filename
import glob
import numpy
import parser
from math import *

class hdfFile:
  def __init__(self,variables,interpVars,dump):
    self.variables=variables
    self.interpVars=interpVars
    print "getting data form dump file \""+dump.fileName+"\" ..."
    self.__getDataFromDump(dump)
    print "  calculating interpolation variables on grid ..."
    #print self.dataMax
    #print self.dataMin
    self.__makeInterpVars()
    #print self.dataMax
    #print self.dataMin
    print "  creating variables ..."
    self.__makeVars()
    quit()
  def __makeInterpVars(self):
    """This function calculates the formulas for the variables used for interpolation on the
    original model dump grid"""
    
    varID=len(self.varNames)-1
    for key in self.interpVars.keys():
      
      #create formula
      formula=self.interpVars[key].formula
      #print key,formula
      for var in self.varNames:
        if var!=None:
          formula=formula.replace(var,"self.data["+str(self.varIDs[var])+"][i][j][k]")
      code=parser.expr(formula).compile()
      
      #set data
      max=-1.0*sys.float_info.min
      min=sys.float_info.max
      for i in range(self.data.shape[1]):
        for j in range(self.data.shape[2]):
          for k in range(self.data.shape[3]):
            self.data[varID][i][j][k]=eval(code)
            if self.data[varID][i][j][k]>max:
              max=self.data[varID][i][j][k]
            if self.data[varID][i][j][k]<min:
              min=self.data[varID][i][j][k]
      
      self.dataMin[varID]=min
      self.dataMax[varID]=max
      self.interpVars[key].max=max
      self.interpVars[key].min=min
      self.interpVars[key].delta=(max-min)/float(self.interpVars[key].numPoints)
      self.varNames.append(key)
      self.varIDs[key]=varID
      varID+=1
  def __makeVars(self):
    
    for n in range(len(self.variables)):
      
      #determine dimension of interpolation
      nNumDims=0
      shape=[]
      nCurrentDim=0
      bIsInterpVar=False
      formula=[]
      for indepVar in self.variables[n].indep:
        if indepVar!=None:
          nNumDims+=1
          if indepVar in self.interpVars.keys():
            shape.append(int(self.interpVars[indepVar].numPoints))
            bIsInterpVar=True
          elif indepVar in self.varNames:
            shape.append(numPoints=int(self.data.shape[nCurrentDim+1]))
        else:
          shape.append(1)
      
      formula=self.variables[n].formula
      for var in self.varNames:
        if var!=None:
          formula=formula.replace(var,"self.data["+str(self.varIDs[var])+"][i][j][k]")
          
      #print "new formula=",formula
      
      print nNumDims,shape
      code=parser.expr(formula).compile()
      location=[0.0,0.0,0.0]
      indepNames=self.variables[n].indep
      indepIndices=[]
      indepIndices.append(self.varIDs[indepNames[0]])
      indepIndices.append(self.varIDs[indepNames[1]])
      indepIndices.append(self.varIDs[indepNames[2]])
      for i in range(shape[0]):
        for j in range(shape[1]):
          for k in range(shape[2]):
            
            location[0]=self.__getXFromInterp(indepNames[0],i)
            location[1]=self.__getXFromInterp(indepNames[1],j)
            location[2]=self.__getXFromInterp(indepNames[2],k)
            print location
            self.__interp3DIJK(location,indepIndices,code)
      
      '''
      if self.variables[n].indep[0]!=None:
        print "    interpolating to independent variable \""\
          +self.variables[n].indep[0]+"\" ..."
        tmp=self.__interpolateToIndep(0,n,self.data[n],code)
      if tmp!=None and self.variables[n].indep[1]!=None:
        print "    interpolating to independent variable \""\
          +self.variables[n].indep[1]+"\" ..."
        code=parser.expr("tmpIn[i][j][k]").compile()
        tmp=self.__interpolateToIndep(1,n,tmp,code)
      elif self.variables[n].indep[1]!=None:
        print "    interpolating to independent variable \""\
          +self.variables[n].indep[1]+"\" ..."
        tmp=self.__interpolateToIndep(1,n,self.data[n],code)
      if tmp!=None and self.variables[n].indep[2]!=None:
        print "    interpolating to independent variable \""\
          +self.variables[n].indep[2]+"\" ..."
        code=parser.expr("tmpIn[i][j][k]").compile()
        tmp=self.__interpolateToIndep(2,n,tmp,code)
      elif self.variables[n].indep[2]!=None:
        print "    interpolating to independent variable \""\
          +self.variables[n].indep[2]+"\" ..."
        tmp=self.__interpolateToIndep(2,n,self.data[n],code)
      
      #check that we had at least one independent variable set
      if tmp==None:
        raise Exception("no indpendent variable specified for variable with "\
          +"formula \""+self.variables[n].formula+"\"")
      
      self.variables[n].data=tmp'''
  def __interp3DIJK(self,location,indep,code):
    #find i,j,k indices
    for i in range(self.data.shape[1]-1):
      for j in range(self.data.shape[2]-1):
        for k in range(self.data.shape[3]-1):
          
          #use AABB for extents of volume element, should be pretty good for
          #usual geometries
          minPoint=[]
          maxPoint=[]
          lower=[None,None,None]
          for l in range(3):
            points=[
              self.data[indep[l]][i  ][j  ][k  ],
              self.data[indep[l]][i  ][j  ][k+1],
              self.data[indep[l]][i  ][j+1][k  ],
              self.data[indep[l]][i  ][j+1][k+1],
              self.data[indep[l]][i+1][j  ][k  ],
              self.data[indep[l]][i+1][j  ][k+1],
              self.data[indep[l]][i+1][j+1][k  ],
              self.data[indep[l]][i+1][j+1][k+1]]
            minPoint.append(min(points))
            maxPoint.append(max(points))
          
          print minPoint
          print maxPoint,location
          
          if minPoint[0]<=location[0]<=maxPoint[0] and\
            minPoint[1]<=location[1]<=maxPoint[1] and\
            minPoint[2]<=location[2]<=maxPoint[2]:
            lower=[i,j,k]
            quit()
    quit()
          #bracketed 
          #if self.data[indep[0]][i][j][k]
  def __getXFromData(self,indepVar,index):
    """indepVar is just a placeholder"""
  def __getXFromInterp(self,indepVar,index):
    """indepVar should be a string identifying the interpolation variable"""
    return self.interpVars[indepVar].min+float(index)*self.interpVars[indepVar].delta
  def __interpolateToIndep(self,l,n,tmpIn,code):
    
    #interpolate to new indep0
    #print self.variables[n].indep[l]
    
    #get min of indep0
    min=self.dataMin[self.varIDs[self.variables[n].indep[l]]]
    
    #get max of indep0
    max=self.dataMax[self.varIDs[self.variables[n].indep[l]]]
    
    #get num points of indep0
    if self.variables[n].indep[l]==None:
      numPoints=1
    elif self.variables[n].indep[l] in self.interpVars.keys():
      numPoints=int(self.interpVars[self.variables[n].indep[l]].numPoints)
    elif self.variables[n].indep[l] in self.varNames:
      numPoints=int(tmpIn.shape[l])
    
    #get independent variable ID
    indepID=self.varIDs[self.variables[n].indep[l]]
    
    #get delta of indep0
    delta=(max-min)/float(numPoints)
          
    #find shape of new array
    if l==0:
      shape=[numPoints,tmpIn.shape[1],tmpIn.shape[2]]
    elif l==1:
      shape=[tmpIn.shape[0],numPoints,tmpIn.shape[2]]
    elif l==2:
      shape=[tmpIn.shape[0],tmpIn.shape[1],numPoints]
    
    print min,max,numPoints,delta,shape
    tmp=numpy.empty(shape)
    return tmp
    lowerLast=None
    
    #loop over new array and interpolate
    #print 
    if l==0:
      for i in range(shape[0]):
        #print ".",
        for j in range(shape[1]):
          for k in range(shape[2]):
            x=min+float(i)*delta
            [value,lowerLast]=self.__interpolateLinearInI(tmpIn,code,x,indepID
              ,lowerLast,j,k)
            tmp[i][j][k]=value
    elif l==1:
      for i in range(shape[0]):
        #print ".",
        for j in range(shape[1]):
          for k in range(shape[2]):
            x=min+float(j)*delta
            [value,lowerLast]=self.__interpolateLinearInJ(tmpIn,code,x,indepID
              ,lowerLast,i,k)
            tmp[i][j][k]=value
    elif l==2:
      for i in range(shape[0]):
        #print ".",
        for j in range(shape[1]):
          for k in range(shape[2]):
            x=min+float(k)*delta
            [value,lowerLast]=self.__interpolateLinearInK(tmpIn,code,x,indepID
              ,lowerLast,i,j)
            tmp[i][j][k]=value
          #if value==None:#if outside grid use fill value
          #  tmp[i][j][k]=float(self.variables[n].fillValue)
          #else:
    return tmp
  def __interpolateLinearIn1DI(self,tmpIn,code,x,indepID,iLowerLast,j,k):
    
    #if value at iLowerLast too large, move in
    if iLowerLast==None:
      iLowerLast=0
    else:
      #check that iLowerLast is lower than x if not move in
      #this is a little tricky as it depends on which way things are increasing
      while(self.data[indepID][iLowerLast][j][k]>x and iLowerLast>0):
        #print self.data[indepID][iLowerLast][j][k],x,iLowerLast
        iLowerLast-=1
    #print iLowerLast
    
    #find bracketing values of x
    iLower=None
    for i in range(iLowerLast,self.data.shape[1]-1):
      #print i,self.data[indepID][i+1][j][k],x,self.data[indepID][i][j][k]
      if (self.data[indepID][i+1][j][k]+sys.float_info.epsilon>=x\
        and self.data[indepID][i][j][k]-sys.float_info.epsilon<=x)\
        or (self.data[indepID][i+1][j][k]-sys.float_info.epsilon<=x \
        and self.data[indepID][i][j][k]+sys.float_info.epsilon>=x) :
        iLower=i
        #print iLowerLast,iLower
        break
      if (self.data[indepID][i][j][k]>x):#this is a speed up for monotonicly increasing values
        break
    
    if iLower==None:
      #print "no bracketing value found for ",x,indepID,j,k,self.dataMin[indepID],self.dataMax[indepID]
      return [None,None]
    
    iLowerLast=iLower
    #do 1D linear interpolation
    i=iLower
    y0=eval(code)
    i=iLower+1
    y1=eval(code)
    x0=self.data[indepID][iLower][j][k]
    x1=self.data[indepID][iLower+1][j][k]
    y=(y1-y0)/(x1-x0)*(x-x0)+y0
    #print x0,x,x1,y0,y,y1,self.data[10][iLower][j][k]\
    #  ,self.data[10][iLower+1][j][k],iLower,iLowerLast
    return [y,iLowerLast]
  def __interpolateLinearIn1DJ(self,tmpIn,code,x,indepID,jLowerLast,i,k):
    
    #if value at jLowerLast too large, move in
    if jLowerLast==None:
      jLowerLast=0
    else:
      #check that jLowerLast is lower than x if not move in
      #this is a little tricky as it depends on which way things are increasing
      while(self.data[indepID][i][jLowerLast][k]>x and jLowerLast>0):
        #print self.data[indepID][i][jLowerLast][k],x,jLowerLast
        jLowerLast-=1
    #print jLowerLast
    
    #find bracketing values of x
    jLower=None
    for j in range(jLowerLast,self.data.shape[2]-1):
      #print j,self.data[indepID][i+1][j][k],x,self.data[indepID][i][j][k]
      if (self.data[indepID][i][j+1][k]+sys.float_info.epsilon>=x\
        and self.data[indepID][i][j][k]-sys.float_info.epsilon<=x)\
        or (self.data[indepID][i][j+1][k]-sys.float_info.epsilon<=x \
        and self.data[indepID][i][j][k]+sys.float_info.epsilon>=x) :
        jLower=j
        #print jLowerLast,jLower
        break
      if (self.data[indepID][i][j][k]>x):#this is a speed up for monotonicly increasing values
        break
    
    if jLower==None:
      #print "no bracketing value found for ",x,indepID,j,k,self.dataMin[indepID],self.dataMax[indep0ID]
      return [None,None]
    
    jLowerLast=jLower
    #do 1D linear interpolation
    j=jLower
    y0=eval(code)
    j=jLower+1
    y1=eval(code)
    x0=self.data[indepID][i][jLower][k]
    x1=self.data[indepID][i][jLower+1][k]
    y=(y1-y0)/(x1-x0)*(x-x0)+y0
    #print x0,x,x1,y0,y,y1,self.data[10][i][jLower][k]\
    #  ,self.data[10][i][jLower+1][k],jLower,jLowerLast
    return [y,jLowerLast]
  def __interpolateLinearIn1DK(self,tmpIn,code,x,indepID,kLowerLast,i,j):
    
    #if value at kLowerLast too large, move in
    if kLowerLast==None:
      kLowerLast=0
    else:
      #check that kLowerLast is lower than x if not move in
      #this is a little tricky as it depends on which way things are increasing
      while(self.data[indepID][i][j][kLowerLast]>x and kLowerLast>0):
        #print self.data[indepID][i][j][kLowerLast],x,kLowerLast
        kLowerLast-=1
    #print kLowerLast
    
    #find bracketing values of x
    kLower=None
    print kLowerLast,self.data[indepID].shape,self.data.shape[3]-1,tmpIn.shape,kLowerLast,self.data.shape[3]-1
    for k in range(kLowerLast,self.data.shape[3]-1):
      print i,j,k,self.data[indepID][i][j][k+1],x,self.data[indepID][i][j][k]
      if (self.data[indepID][i][j][k+1]+sys.float_info.epsilon>=x\
        and self.data[indepID][i][j][k]-sys.float_info.epsilon<=x)\
        or (self.data[indepID][i][j][k+1]-sys.float_info.epsilon<=x \
        and self.data[indepID][i][j][k]+sys.float_info.epsilon>=x) :
        kLower=k
        #print kLowerLast,kLower
        break
      if (self.data[indepID][i][j][k]>x):#this is a speed up for monotonicly increasing values
        break
    
    if kLower==None:
      #print "no bracketing value found for ",x,indepID,j,k,self.dataMin[indepID],self.dataMax[indepID]
      return [None,None]
    
    jLowerLast=kLower
    #do 1D linear interpolation
    k=kLower
    y0=eval(code)
    k=kLower+1
    y1=eval(code)
    x0=self.data[indepID][i][j][kLower]
    x1=self.data[indepID][i][j][kLower+1]
    y=(y1-y0)/(x1-x0)*(x-x0)+y0
    #print x0,x,x1,y0,y,y1,self.data[10][i][kLower][k]\
    #  ,self.data[10][i][kLower+1][k],kLower,jLowerLast
    return [y,jLowerLast]
  def __getDataFromDump(self,dump):
    
    self.varNames=dump.varNames
    self.varIDs=dump.varIDs
    shape=[]
    
    #number of data variables
    tmp=len(dump.vars)+len(self.interpVars)
    shape.append(tmp)
    
    shape.append(len(dump.vars[dump.varIDs["M_r"]])-1)
    if dump.numDims>1:
      shape.append(len(dump.vars[dump.varIDs["theta"]][0]))
    else:
      shape.append(1)
    if dump.numDims>2:
      shape.append(len(dump.vars[dump.varIDs["phi"]][0][0]))
    else:
      shape.append(1)
    
    #make data all 3D, copying/interpolating where needed
    
    self.data=numpy.empty(shape)
    self.dataMax=numpy.empty(shape[0])
    self.dataMax.fill(-1.0*sys.float_info.min)
    self.dataMin=numpy.empty(shape[0])
    self.dataMin.fill(sys.float_info.max)
    
    for n in range(len(dump.vars)):
      dim0=len(dump.vars[n])
      dim1=len(dump.vars[n][dim0-1])
      dim2=len(dump.vars[n][dim0-1][dim1-1])
      
      for i in range(self.data.shape[1]):
        #the below has assumed that the model is periodic in theta & phi, but not in radial direction
        
        if dump.varInfo[n][0:3]==[1,-1,-1]:#radial only interface variable
          formula="(dump.vars[n][i+1][0][0]+dump.vars[n][i][0][0])*0.5"
        if dump.varInfo[n][0:3]==[-1,1,-1]:#theta only interface variable
          formula="(dump.vars[n][0][j-1][0]+dump.vars[n][0][j][0])*0.5"
        if dump.varInfo[n][0:3]==[-1,-1,1]:#phi only interface variable
          formula="(dump.vars[n][0][0][k-1]+dump.vars[n][0][0][k])*0.5"
        if dump.varInfo[n][0:3]==[0,-1,-1]:#radial only centered variable
          formula="dump.vars[n][i][0][0]"
        if dump.varInfo[n][0:3]==[0,-1,-1]:#radial only centered variable
          formula="dump.vars[n][i][0][0]"
        
        #centered variable in 3D region
        if dump.varInfo[n][0:3]==[0,0,0] and i>=dump.num1DZones+dump.numGhostCells:
          formula="dump.vars[n][i][j][k]"
        
        #centered variable in 1D region
        if dump.varInfo[n][0:3]==[0,0,0] and i<dump.num1DZones+dump.numGhostCells:
          formula="dump.vars[n][i][0][0]"
        
        #radial interface, theta/phi centered variable 3D region
        if dump.varInfo[n][0:3]==[1,0,0] and i>dump.num1DZones+dump.numGhostCells:
          formula="(dump.vars[n][i+1][j][k]+dump.vars[n][i][j][k])*0.5"
        
        #radial interface, theta/phi centered variable 1D region
        if dump.varInfo[n][0:3]==[1,0,0] and i<=dump.num1DZones+dump.numGhostCells:
          formula="(dump.vars[n][i+1][0][0]+dump.vars[n][i][0][0])*0.5"
        
        #theta interface radial/phi centered variable 3D region
        if dump.varInfo[n][0:3]==[0,1,0] and i>=dump.num1DZones+dump.numGhostCells:
          formula="(dump.vars[n][i][j-1][k]+dump.vars[n][i][j][k])*0.5"
        
        #theta interface radial/phi centered variable 1D region
        if dump.varInfo[n][0:3]==[0,1,0] and i<dump.num1DZones+dump.numGhostCells:
          formula="dump.vars[n][i][0][0]"
        
        #theta interface radial/phi centered variable 3D region
        if dump.varInfo[n][0:3]==[0,0,1] and i>=dump.num1DZones+dump.numGhostCells:
          formula="(dump.vars[n][i][j][k-1]+dump.vars[n][i][j][k])*0.5"
        
        #theta interface radial/phi centered variable 1D region
        if dump.varInfo[n][0:3]==[0,0,1] and i<dump.num1DZones+dump.numGhostCells:
          formula="dump.vars[n][i][0][0]"
        
        if n==1:#it is theta, need to do something special for first zone 
          formula1="(dump.vars[n][0][j][0]-(dump.vars[n][0][j+1][0]-dump.vars[n][0][j][0])*0.5)"
        else:
          formula1=formula
        if n==2:#this is phi, need to do something special for first zone 
          formula2="(dump.vars[n][0][0][k]-(dump.vars[n][0][0][k+1]-dump.vars[n][0][0][k])*0.5)"
        else:
          formula2=formula
        
        code=parser.expr(formula).compile()
        code1=parser.expr(formula1).compile()
        code2=parser.expr(formula2).compile()
        for j in range(self.data.shape[2]):
          for k in range(self.data.shape[3]):
            
            #set values
            if j==0 and n==1:
              self.data[n][i][j][k]=eval(code1)
            if k==0 and n==2:
              self.data[n][i][j][k]=eval(code2)
            else:
              self.data[n][i][j][k]=eval(code)
            
            #save min/max
            if self.data[n][i][j][k]>self.dataMax[n]:
              self.dataMax[n]=self.data[n][i][j][k]
            if self.data[n][i][j][k]<self.dataMin[n]:
              self.dataMin[n]=self.data[n][i][j][k]
  def printVarToScreen(self, n):
    for i in range(self.data.shape[1]):
      for j in range(self.data.shape[2]):
        for k in range(self.data.shape[3]):
          print self.varNames[n]+" ("+str(i)+","+str(j)+","+str(k)+")="+str(self.data[n][i][j][k])
        print ""
      print ""
  def write(self):
    """this function writes the data specified in the configuration file to a new hdf file. It does
    this by interpolating where nessacary to get data at the right location"""
    
    #for i in range(len(self.variables)):
    #  print self.variables[i].name
    warnings.warn("hdfFile.write is not yet implemented")
class variable:
  def __init__(self,element):
    
    #get indpendent variables
    self.indep=[]
    self.indep.append(element.get("indep0"))
    self.indep.append(element.get("indep1"))
    self.indep.append(element.get("indep2"))
    self.fillValue=[]
    self.fillValue.append(element.get("fillValue"))
    #get variable name
    if element.text=="" or element.text==None:
      raise Exception("variable reference empty")
    
    self.formula=element.text
class interpVar:
  def __init__(self,element):
    
    self.numPoints=int(element.get("numPoints"))
    self.name=element.get("name")
    
    self.formula=element.text
    if element.text=="" or element.text==None:
      raise Exception("interpolation variable formula is empty")
    
    self.formula=element.text
class fileSet:
  def __init__(self,element):
    """Initialize an fileSet from an xml node"""
    
    #define allowed attributes for nodes
    self.__setSupportedNodeAttributes()
    
    #warn about unsupported attributes in node
    self.__checkSuppotedNodeAttributes(element)
    
    #get file range
    self.fileRange=element.get("fileRange")
    if self.fileRange==None or self.fileRange=="":
      raise Exception("Need a fileRange attribute in fileSet Node")
    [self.start,self.end,self.baseFileName]=disect_filename.disectFileName(self.fileRange)
    
    #get timeFile Name
    self.timeFile=element.get("timeFile")
    
    #get dataPerFile node
    dataPerFile=element.findall("dataPerFile")
    
    #check to see if there are more than one dataPerFile node
    if len(dataPerFile)>1:
      message="more than one \"dataPerFile\" node found only using first"\
        +" node and ignoring the rest"
      warnings.warn(message)
    
    #warn about unsupported attributes in node
    self.__checkSuppotedNodeAttributes(dataPerFile[0])
    
    #get list of variables to include in file
    self.variables=[]
    variableElements=dataPerFile[0].findall("variable")
    for variableElement in variableElements:
      
      #check that attributes are supported
      self.__checkSuppotedNodeAttributes(variableElement)
      self.variables.append(variable(variableElement))
    
    #get list of interpolation variables 
    self.interpVars={}
    interpVarElements=element.findall("interpVar")
    for interpVarElement in interpVarElements:
      
      #check that attributes of node are supported
      self.__checkSuppotedNodeAttributes(interpVarElement)
      tmp=interpVarElement.get("name")
      if tmp==None or tmp=="":
        raise Exception("\"interpVar\" must have a \"name\" attribute set")
      
      self.interpVars[tmp]=interpVar(interpVarElement)
  def __checkSuppotedNodeAttributes(self,element):
    """check for attributes not supported in dataPerFile element"""
    
    if element.tag not in self.supportedNodeAttributes.keys():
      message="node type \""+element.tag+"\" does not have supporting attributes listed"
      warnings.warn(message)
    
    for attribute in element.attrib.keys():
      if attribute not in self.supportedNodeAttributes[element.tag]:
        message="attribute \""+attribute+"\" in a \""+element.tag\
          +"\" node has no effect"
        warnings.warn(message)
  def __setSupportedNodeAttributes(self):
    """Sets supported attributes for the different nodes"""
    
    self.supportedNodeAttributes={}
    
    #for a fileSet Node
    self.supportedNodeAttributes["fileSet"]=["fileRange","timeFile"]
    
    #for a dataPerFile Node
    self.supportedNodeAttributes["dataPerFile"]=[]
    
    #for variable node
    self.supportedNodeAttributes["variable"]=["indep0","indep1","indep2"\
      ,"fillValue"]
    
    #for interpVar
    self.supportedNodeAttributes["interpVar"]=["numPoints","name"]
  def makeHDFFiles(self,options):
    """Makes HDF files specified by settings"""
    
    #check to make sure that all combined binary files are made
    fileName=self.baseFileName+"["+str(self.start)+"-"+str(self.end)+"]"
    failedFiles=combine_bins.combine_bin_files(options.keep,fileName,options.remakeBins)
    
    #get a list of files in our range
    filesExistCombBin=glob.glob(self.baseFileName+"[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]")
    files=[]
    for file in filesExistCombBin:
      intOfFile=int(file[len(self.baseFileName):len(file)])
      if intOfFile>=self.start and intOfFile<self.end:
        files.append(file)
    
    files.sort()
    
    #loop over files to make HDF files
    if len(files)==0:
      warnings.warn("no files found in \""+fileName+"\"")
    
    for file in files:
      tmp=dump.Dump(file)
      hdfTmp=self.convertDumpToHDF(tmp)
      hdfTmp.write()
  def convertDumpToHDF(self,dump):
    """Converts a dump file to an hdf file
    
    formated in the way sepcified in the xml configuration file
    
    """
    
    #make interpolated variables
    tmpHDF=hdfFile(self.variables,self.interpVars,dump)
    
    #create variables from dump
    
    return tmpHDF
def parseOptions():
  
  #setup command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] XMLFILE"
    ,version="%prog 1.0",description=r"Creates one or more HDF files. XMLFILE "
    +"gives all the settings needed. For an example of a configuration file see"
    +" make_hdf_reference.xml under docs/templateXML. The make_hdf2.py is "
    +"working, while this version make_hdf.py is not. This version is far more "
    +"ambitious in terms of flexability and capability, but alas does not "
    +"work. So use make_hdf2.py instead.")
  
  parser.add_option("-k","--keep",action="store_true",dest="keep"
    ,help="Keeps distributed binary files [default].",default=True)
  parser.add_option("-r","--remove",action="store_false",dest="keep"
    ,help="Removes distributed binary files")
  parser.add_option("--remake-bins",action="store_true",dest="remakeBins"
    ,help="Will remake binaries even if they already exist. [not default].",default=False)
  
  #parse command line options
  (options,args)=parser.parse_args()
  if len(args)==0:
    raise Exception("must supply an xml configure file")
  return (options,args)
def main():
  
  #parse command line options
  (options,args)=parseOptions()
  
  #get root element
  tree=xml.parse(args[0])
  root=tree.getroot()
  
  #for each fileSet get settings, and make HDF files
  fileSetElements=root.findall("fileSet")
  for fileSetElement in fileSetElements:
    tmp=fileSet(fileSetElement)
    tmp.makeHDFFiles(options)
if __name__ == "__main__":
  main()
  
#dumpFile=dump.Dump("./T5700_5x5_t00000000")
#dumpFile.printHeader()
#dumpFile.printVar(5)

'''
dims=[2,2,2]
data=[
[[1.0,1.1], [2.0,2.1]],
[[1.0,1.1], [2.0,2.1]]
]
hdf.open("test.hdf")
hdf.openData("test",hdf.DFNT_FLOAT64,dims)
hdf.writeData(data)
hdf.closeData()
hdf.openData("test2",hdf.DFNT_FLOAT64,dims)
hdf.closeData()
hdf.close()
'''