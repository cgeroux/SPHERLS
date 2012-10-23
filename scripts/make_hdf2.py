#!/usr/bin/env python
import sys
import paths
import optparse as op
import xml.etree.ElementTree as xml
import parser
import numpy
import disect_filename
import glob
import warnings
import mywarnings#sets a custom warning print out
import dump
import combine_bins
from math import *
sys.path.append(paths.srcPath+"/pythonextensions/lib/python/")
import hdf
import os.path

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
      raise Exception("Need a \"fileRange\" attribute in \"fileSet\" Node")
    [self.start,self.end,self.baseFileName]=disect_filename.disectFileName(self.fileRange)
    
    #get timeFile Name
    self.timeFile=element.get("timeFile")
    if self.timeFile==None:
      self.timeFile="timeFile.hdf"
    
    #get frequency
    self.frequency=1
    self.frequency=int(element.get("frequency"))
    
    #get output path
    self.outputPath=element.get("outputPath")
    if self.outputPath==None or self.outputPath=="":
      raise Exception("Need a \"outputPath\" attribute in \"fileSet\" Node")
    
    #get radialCutZone element text
    elementRadialCutZone=element.findall("radialCutZone")
    if len(elementRadialCutZone)>1:
      warnings.warn("more than one \"radialCutZone\" node ignoring all but first node")
    elif len(elementRadialCutZone)==1:
      self.__checkSuppotedNodeAttributes(elementRadialCutZone[0])
      if elementRadialCutZone[0].text!=None or elementRadialCutZone[0].text!="":
        self.radialCutZone=int(elementRadialCutZone[0].text)
    else:
      raise Exception("must have a \"radialCutZone\" node in a \"fileSet\" node with an integer value")
    
    #get includeBoundaries element text
    elementIncludeBoundaries=element.findall("includeBoundaries")
    if len(elementIncludeBoundaries)>1:
      warnings.warn("more than one \"includeBoundaries\" node ignoring all but first node")
    elif len(elementIncludeBoundaries)==1:
      self.__checkSuppotedNodeAttributes(elementIncludeBoundaries[0])
      if elementIncludeBoundaries[0]!=None:
        if elementIncludeBoundaries[0].text in ["true","yes","y","t","1"]:
          self.includeBoundaries=True
        elif elementIncludeBoundaries[0].text in ["false","no","n","n","0"]:
          self.includeBoundaries=False
        else:
          raise Exception("\"includeBoundaries\" node expects \"true\" or \"false\"")
    else:
      raise Exception("must have a \"includeBoundaries\" node in a \"fileSet\" node")
    
    #get numRInterp element text
    elementnumRInterp=element.findall("numRInterp")
    if len(elementnumRInterp)>1:
      warnings.warn("more than one \"numRInterp\" node ignoring all but first node")
    elif len(elementnumRInterp)==1:
      self.__checkSuppotedNodeAttributes(elementnumRInterp[0])
      if elementnumRInterp[0].text!=None or elementnumRInterp[0]!="":
        self.numRInterp=int(elementnumRInterp[0].text)
    else:
      raise Exception("must have a \"numRInterp\" node in a \"fileSet\" node with an integer value")
      
    '''
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
    '''
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
    self.supportedNodeAttributes["fileSet"]=["fileRange","timeFile","outputPath","frequency"]
    
    #for a dataPerFile Node
    self.supportedNodeAttributes["radialCutZone"]=[]
    
    #for variable node
    self.supportedNodeAttributes["includeBoundaries"]=[]
    
    #for interpVar
    self.supportedNodeAttributes["numRInterp"]=[]
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
    
    times=[]
    numFiles=len(files)-1
    for i in range(0,len(files),self.frequency):
      file=files[i]
      print "reading file \""+file+"\" "+str(i)+"/"+str(numFiles)+" ..."
      tmp=dump.dump(file)
      times.append(tmp.time)
      self.convertDumpToHDF(tmp)
    
    #write out time
    hdf.open(self.timeFile)
    hdf.openData("time",hdf.DFNT_FLOAT64,[len(times)])
    hdf.writeData(times)
    hdf.closeData()
    hdf.close()
  def convertDumpToHDF(self,dump):
    """Converts a dump ifle to an hdf file formated in the way sepcified in the 
    xml configuration file"""
    
    nRadialCutZone=self.radialCutZone
    bIncludeBoundaries=self.includeBoundaries
    numRInterp=self.numRInterp
    
    self.getDataFromDump(dump,nRadialCutZone,bIncludeBoundaries)
    self.setAdditionalVariables(numRInterp)
    if self.outputPath!=None:
      file=os.path.join(self.outputPath,os.path.basename(dump.fileName)+".hdf")
    else:
      file=dump.fileName+".hdf"
    
    hdf.open(file)
    n=0
    print "writting to hdf file \""+file+"\" ..."
    for data in self.data:
      #print "  writting \""+self.dataNames[n]+"\" to file ..."
      hdf.openData(self.dataNames[n],hdf.DFNT_FLOAT64,self.dataShape[n])
      hdf.writeData(data)
      hdf.closeData()
      n+=1
    hdf.close()
  def setAdditionalVariables(self,numRInterp):
    
    #make scaled theta
    nTheta=self.dataIDs["theta"]
    thetaScale=0.1/(self.dataMax[nTheta]-self.dataMin[nTheta])
    self.dataNames.append("theta scaled between 0 and 0.1")
    dataI=[]
    for i in range(self.dataShape[nTheta][0]):
      dataJ=[]
      for j in range(self.dataShape[nTheta][1]):
        dataK=[]
        for k in range(self.dataShape[nTheta][2]):
          dataK.append((self.data[nTheta][i][j][k]-self.dataMin[nTheta])*thetaScale)
        dataJ.append(dataK)
      dataI.append(dataJ)
    self.data.append(dataI)
    self.dataShape.append(self.dataShape[nTheta])
    
    #make scaled phi
    nPhi=self.dataIDs["phi"]
    phiScale=0.1/(self.dataMax[nPhi]-self.dataMin[nPhi])
    self.dataNames.append("phi scaled between 0 and 0.1")
    dataI=[]
    for i in range(self.dataShape[nPhi][0]):
      dataJ=[]
      for j in range(self.dataShape[nPhi][1]):
        dataK=[]
        for k in range(self.dataShape[nPhi][2]):
          dataK.append((self.data[nPhi][i][j][k]-self.dataMin[nPhi])*phiScale)
        dataJ.append(dataK)
      dataI.append(dataJ)
    self.data.append(dataI)
    self.dataShape.append(self.dataShape[nPhi])
    
    #make scaled v
    #needs some work, not sure the units are correct, should think about scaling factor some more
    nV=self.dataIDs["v"]
    nR=self.dataIDs["r"]
    nU=self.dataIDs["u"]
    nU0=self.dataIDs["u_0"]
    '''self.dataNames.append("v scaled to units of scaled theta/s")
    dataI=[]
    for i in range(self.dataShape[nV][0]):
      dataJ=[]
      for j in range(self.dataShape[nV][1]):
        dataK=[]
        for k in range(self.dataShape[nV][2]):
          dataK.append((self.data[nV][i][j][k])*thetaScale/(self.data[nR][i][0][0]))
        dataJ.append(dataK)
      dataI.append(dataJ)
    self.data.append(dataI)
    self.dataShape.append(self.dataShape[nV])'''
    
    #make scaled w
    #needs some work, not sure the units are correct, should think about scaling factor some more
    nW=self.dataIDs["w"]
    '''self.dataNames.append("w scaled to units of scaled phi/s")
    dataI=[]
    for i in range(self.dataShape[nW][0]):
      dataJ=[]
      for j in range(self.dataShape[nW][1]):
        dataK=[]
        for k in range(self.dataShape[nW][2]):
          dataK.append(self.data[nW][i][j][k]
            /(self.data[nR][i][0][0]*sin(self.data[nTheta][0][j][0]))*phiScale)
        dataJ.append(dataK)
      dataI.append(dataJ)
    self.data.append(dataI)
    self.dataShape.append(self.dataShape[nW])'''
    
    #make interpolation independent variable r/R_surf, T and rho,u-u0,v,w
    rScale=1.0/self.dataMax[nR]
    delta=(self.dataMax[nR]-self.dataMin[nR])/float(numRInterp)*rScale
    start=self.dataMin[nR]*rScale
    dataI=[]
    dataIT=[]
    dataIRho=[]
    dataIUmU0Scaled=[]
    dataIUScaled=[]
    dataIVScaled=[]
    dataIWScaled=[]
    dataIUmU0=[]
    dataIU=[]
    dataIV=[]
    dataIW=[]
    nT=self.dataIDs["T"]
    nRho=self.dataIDs["rho"]
    shapeR=[numRInterp,1,1]
    self.dataNames.append("r/R to evenly spaced intervals")
    self.dataNames.append("T interpolated to r/R")
    self.dataNames.append("Rho interpolated to r/R")
    self.dataNames.append("u-u_0 interpolated to r/R")
    self.dataNames.append("u interpolated to r/R")
    self.dataNames.append("v interpolated to r/R")
    self.dataNames.append("w interpolated to r/R")
    self.dataNames.append("scaled u-u_0 interpolated to r/R")
    self.dataNames.append("scaled u interpolated to r/R")
    self.dataNames.append("scaled v interpolated to r/R")
    self.dataNames.append("scaled w interpolated to r/R")
    
    if self.dataShape[nT]!=self.dataShape[nRho]:
      raise Exception("something strange happened, self.dataShape[nT]="+str(self.dataShape[nT])
        +" while self.dataShape[nRho]="+str(self.dataShape[nRho])+" they should be the same")
    shapeT=[numRInterp,self.dataShape[nT][1],self.dataShape[nT][2]]
    iLowerLast=None
    for i in range(shapeT[0]):
      r=delta*float(i)+start
      dataJ=[]
      for j in range(shapeR[1]):
        dataK=[]
        for k in range(shapeR[2]):
          dataK.append(r)
        dataJ.append(dataK)
      dataI.append(dataJ)
      dataJT=[]
      dataJRho=[]
      dataJUmU0=[]
      dataJU=[]
      dataJV=[]
      dataJW=[]
      dataJUmU0Scaled=[]
      dataJUScaled=[]
      dataJVScaled=[]
      dataJWScaled=[]
      
      for j in range(shapeT[1]):
        dataKT=[]
        dataKRho=[]
        dataKUmU0=[]
        dataKU=[]
        dataKV=[]
        dataKW=[]
        dataKUmU0Scaled=[]
        dataKUScaled=[]
        dataKVScaled=[]
        dataKWScaled=[]
        for k in range(shapeT[2]):
          
          #T
          [value,iLowerLast]=self.__interpolateLinearIn1DI(r,nR,nT,rScale
            ,iLowerLast,j,k)
          dataKT.append(value)
          
          #rho
          [value,iLowerLast]=self.__interpolateLinearIn1DI(r,nR,nRho,rScale
            ,iLowerLast,j,k)
          dataKRho.append(value)
          
          #U-U0 Scaled
          [valueU,iLowerLast]=self.__interpolateLinearIn1DI(r,nR,nU,rScale
            ,iLowerLast,j,k)
          [valueU0,iLowerLast]=self.__interpolateLinearIn1DI(r,nR,nU0,rScale
            ,iLowerLast,0,0)
          dataKUmU0Scaled.append((valueU-valueU0)*rScale)
          
          #U Scaled
          dataKUScaled.append(valueU*rScale)
          
          #V Scaled
          [valueV,iLowerLast]=self.__interpolateLinearIn1DI(r,nR,nV,rScale
            ,iLowerLast,j,k)
          dataKVScaled.append(valueV*thetaScale/r*rScale)
          
          #W Scaled
          [valueW,iLowerLast]=self.__interpolateLinearIn1DI(r,nR,nW,rScale
            ,iLowerLast,j,k)
          dataKWScaled.append(valueW/(r*sin(self.data[nTheta][0][j][0]))*phiScale*rScale)
          
          #U-U0 
          [valueU,iLowerLast]=self.__interpolateLinearIn1DI(r,nR,nU,rScale
            ,iLowerLast,j,k)
          [valueU0,iLowerLast]=self.__interpolateLinearIn1DI(r,nR,nU0,rScale
            ,iLowerLast,0,0)
          dataKUmU0.append((valueU-valueU0))
          
          #U
          dataKU.append(valueU)
          
          #V 
          dataKV.append(valueV)
          
          #W 
          dataKW.append(valueW)
          
        dataJT.append(dataKT)
        dataJRho.append(dataKRho)
        dataJUmU0.append(dataKUmU0)
        dataJU.append(dataKU)
        dataJV.append(dataKV)
        dataJW.append(dataKW)
        dataJUmU0Scaled.append(dataKUmU0Scaled)
        dataJUScaled.append(dataKUScaled)
        dataJVScaled.append(dataKVScaled)
        dataJWScaled.append(dataKWScaled)
      dataIT.append(dataJT)
      dataIRho.append(dataJRho)
      dataIUmU0.append(dataJUmU0)
      dataIU.append(dataJU)
      dataIV.append(dataJV)
      dataIW.append(dataJW)
      dataIUmU0Scaled.append(dataJUmU0Scaled)
      dataIUScaled.append(dataJUScaled)
      dataIVScaled.append(dataJVScaled)
      dataIWScaled.append(dataJWScaled)
    
    self.data.append(dataI)
    self.dataShape.append(shapeR)
    self.data.append(dataIT)
    self.dataShape.append(shapeT)
    self.data.append(dataIRho)
    self.dataShape.append(shapeT)
    self.data.append(dataIUmU0)
    self.dataShape.append(shapeT)
    self.data.append(dataIU)
    self.dataShape.append(shapeT)
    self.data.append(dataIV)
    self.dataShape.append(shapeT)
    self.data.append(dataIW)
    self.dataShape.append(shapeT)
    self.data.append(dataIUmU0Scaled)
    self.dataShape.append(shapeT)
    self.data.append(dataIUScaled)
    self.dataShape.append(shapeT)
    self.data.append(dataIVScaled)
    self.dataShape.append(shapeT)
    self.data.append(dataIWScaled)
    self.dataShape.append(shapeT)
  def __interpolateLinearIn1DI(self,x,indepVarID,varID,rScale,iLowerLast,j,k):
    
    #if value at iLowerLast too large, move in
    if iLowerLast==None:
      iLowerLast=0
    else:
      #check that iLowerLast is lower than x if not move in
      #this is a little tricky as it depends on which way things are increasing
      
      while(self.data[indepVarID][iLowerLast][0][0]*rScale>x and iLowerLast>0):
        print self.data[indepVarID][iLowerLast][0][0]*rScale,x,iLowerLast
        iLowerLast-=1
    #print iLowerLast
    
    #find bracketing values of x
    iLower=None
    for i in range(iLowerLast,self.dataShape[indepVarID][0]-1):
      #print i,self.data[indepID][i+1][j][k],x,self.data[indepID][i][j][k]
      if (self.data[indepVarID][i+1][0][0]*rScale+sys.float_info.epsilon>=x\
        and self.data[indepVarID][i][0][0]*rScale-sys.float_info.epsilon<=x) :
        iLower=i
        #print iLowerLast,iLower
        break
      if (self.data[indepVarID][i][0][0]*rScale>x):#this is a speed up for monotonicly increasing values
        break
    
    if iLower==None:
      tmp="no bracketing value found for "+str(x)+" "\
        +self.dataNames[indepVarID]+" "+str(j)+" "+str(k)+" "\
        +str(self.dataMin[indepVarID]*rScale)+" "+str(self.dataMax[indepVarID]\
        *rScale)+" "+str(rScale)
      warnings.warn(tmp)
      return [None,None]
    
    iLowerLast=iLower
    #do 1D linear interpolation
    y0=self.data[varID][iLower][j][k]
    y1=self.data[varID][iLower+1][j][k]
    x0=self.data[indepVarID][iLower][0][0]*rScale
    x1=self.data[indepVarID][iLower+1][0][0]*rScale
    y=(y1-y0)/(x1-x0)*(x-x0)+y0
    #print x0,x,x1,y0,y,y1,self.data[10][iLower][j][k]\
    #  ,self.data[10][iLower+1][j][k],iLower,iLowerLast
    return [y,iLowerLast]
  def getDataFromDump(self,dump,nRadialCutZone,bIncludeBoundaries):
    
    #create variables from dump
    self.data=[]
    self.dataMax=[]
    self.dataMin=[]
    self.dataShape=[]
    self.dataNames=[]
    self.dataIDs={}
    for n in range(len(dump.vars)):
      
      #get size of variable
      shape=[]
      self.dataNames.append(dump.varNames[n])
      self.dataIDs[dump.varNames[n]]=n
      
      #print "  converting \""+dump.varNames[n]+"\" to desired hdf format ..."
      
      #radial direction
      if dump.varInfo[n][0]==1:#radial interface
        shape.append(len(dump.vars[n])-1-nRadialCutZone)
      elif dump.varInfo[n][0]==0:#radial centered
        shape.append(len(dump.vars[n])-nRadialCutZone)
      elif dump.varInfo[n][0]==-1:#not defined
        shape.append(1)
      
      #theta direction
      if dump.varInfo[n][1]==1:#radial interface
        index=0
        if shape[0]!=1:
          index=shape[0]-1+nRadialCutZone
        shape.append(len(dump.vars[n][index]))
      elif dump.varInfo[n][1]==0:#radial centered
        index=0
        if shape[0]!=1:
          index=shape[0]-1+nRadialCutZone
        shape.append(len(dump.vars[n][index]))
      elif dump.varInfo[n][1]==-1:#not defined
        shape.append(1)
      
      #phi direction
      if dump.varInfo[n][2]==1:#radial interface
        index=0
        if shape[0]!=1:
          index=shape[0]-1+nRadialCutZone
        shape.append(len(dump.vars[n][index][shape[1]-1]))
      elif dump.varInfo[n][2]==0:#radial centered
        index=0
        if shape[0]!=1:
          index=shape[0]-1+nRadialCutZone
        shape.append(len(dump.vars[n][index][shape[1]-1]))
      elif dump.varInfo[n][2]==-1:#not defined
        shape.append(1)
      
      #adjust shape to include boundaries or not
      startY=0
      startZ=0
      if not bIncludeBoundaries:
        startY=dump.numGhostCells
        startZ=dump.numGhostCells
        if dump.varInfo[n][1]==-1:
          startY=0
        else:
          shape[1]-=2*startY
        if dump.varInfo[n][2]==-1:
          startZ=0
        else:
          shape[2]-=2*startZ
      
      self.dataShape.append(shape)
      
      minVar=sys.float_info.max
      maxVar=-1.0*sys.float_info.max
      dataTemp=[]
      
      for i in range(nRadialCutZone,shape[0]+nRadialCutZone):
        
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
        
        if n==1 and bIncludeBoundaries:#it is theta, need to do something special for first zone 
          formula1="(dump.vars[n][0][j][0]-(dump.vars[n][0][j+1][0]-dump.vars[n][0][j][0])*0.5)"
        else:
          formula1=formula
        if n==2 and bIncludeBoundaries:#this is phi, need to do something special for first zone 
          formula2="(dump.vars[n][0][0][k]-(dump.vars[n][0][0][k+1]-dump.vars[n][0][0][k])*0.5)"
        else:
          formula2=formula
        
        code=parser.expr(formula).compile()
        code1=parser.expr(formula1).compile()
        code2=parser.expr(formula2).compile()
        tmpDataJ=[]
        for j in range(startY,shape[1]+startY):
          tmpDataK=[]
          for k in range(startZ,shape[2]+startZ):
            
            #set values
            if j==0 and n==1:
              value=eval(code1)
            elif k==0 and n==2:
              value=eval(code2)
            else:
              value=eval(code)
            
            #save min/max
            minVar=min(minVar,value)
            maxVar=max(maxVar,value)
            tmpDataK.append(value)
          tmpDataJ.append(tmpDataK)
        dataTemp.append(tmpDataJ)
      self.data.append(dataTemp)
      self.dataMax.append(maxVar)
      self.dataMin.append(minVar)
def parseOptions():
  
  #setup command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] XMLFILE"
    ,version="%prog 1.0",description=r"Creates one or more HDF files. XMLFILE "
    +"gives all the settings needed. For an example of a configuration file see"
    +" make_hdf_reference.xml under docs/templateXML.")
  
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
