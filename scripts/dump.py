#!/usr/bin/env python
import struct
import sys
import numpy as np
import math
import paths
import os
from timeitDec import timeitDec

#The blow two lines should not needed when a proper install is done
sys.path.append(paths.srcPath+"/pythonextensions/lib/python/")
sys.path.append(paths.srcPath+"/pythonextensions/lib/python/evtk")

import evtk.hl
import eos

eosFileName=None
eosTable=None

class Dump:
  """Allows manipulation of SPHERLS binary and ascii dump files.
  
  \todo should probably add methods to write out binary/ascii dump files in the
  version that SPHERLS knows how to read.
  """
  
  def __init__(self,fileName=None):
    """Initilizes the dump by reading in a binary file.
    
    """
    
    self.rectVars=[]
    
    #if a file name was given read in a file
    if fileName!=None:
      self.read(fileName)
  def _addVarID(self,name,type):
    """Set ID for a new variable called name
    
    This should be called every time a new variable is appended to vars. And 
    should be called with the names in the same order the new variables are
    appended
    
    name: name of the new variable
    type: hybrid=uses a 1D-muti-D grid, rect=uses a rectangular grid, both
    """
    
    self._varIDs[name]=len(self.varNames)
    self.varNames.append(name)
    self.varType[name]=type
  def _setVarIDs(self):
    """Sets names for the interger values of the grid varibles
    
    The array indices are purely set by the order in which the dumps have
    organized the variables, the names below refelect this order.
    """
    
    self._varIDs={}
    self.varType={}#hybrid, rect, both
    self.varNames=[]
    if self.gamma!=None: #using gamma law gas
      if self.numDims==1:
        self.varNames=["M_r","Delta_M_r","r","rho","u","u_0","e"]
        i=0
        for varName in self.varNames:
          self.varIDs[varName]=i
          i+=1
        
        #set data type
        self.varType=dict.fromkeys(self.varNames,"hybrid")
        
      elif self.numDims==2:
        self._varIDs["M_r"]=0
        self._varIDs["theta"]=1
        self._varIDs["Delta_M_r"]=2
        self._varIDs["r"]=3
        self._varIDs["rho"]=4
        self._varIDs["u"]=5
        self._varIDs["u_0"]=6
        self._varIDs["v"]=7
        self._varIDs["e"]=8
        
        #set data type
        self._varIDs["M_r"]="hybrid"
        self._varIDs["theta"]="hybrid"
        self._varIDs["Delta_M_r"]="hybrid"
        self._varIDs["r"]="hybrid"
        self._varIDs["rho"]="hybrid"
        self._varIDs["u"]="hybrid"
        self._varIDs["u_0"]="hybrid"
        self._varIDs["v"]="hybrid"
        self._varIDs["e"]="hybrid"
        
        self.varNames.append("M_r")
        self.varNames.append("theta")
        self.varNames.append("Delta_M_r")
        self.varNames.append("r")
        self.varNames.append("rho")
        self.varNames.append("u")
        self.varNames.append("u_0")
        self.varNames.append("v")
        self.varNames.append("e")

      elif self.numDims==3:
        self._varIDs["M_r"]=0
        self._varIDs["theta"]=1
        self._varIDs["phi"]=2
        self._varIDs["Delta_M_r"]=3
        self._varIDs["r"]=4
        self._varIDs["rho"]=5
        self._varIDs["u"]=6
        self._varIDs["u_0"]=7
        self._varIDs["v"]=8
        self._varIDs["w"]=9
        self._varIDs["e"]=10
        
        self.varNames.append("M_r")
        self.varNames.append("theta")
        self.varNames.append("phi")
        self.varNames.append("Delta_M_r")
        self.varNames.append("r")
        self.varNames.append("rho")
        self.varNames.append("u")
        self.varNames.append("u_0")
        self.varNames.append("v")
        self.varNames.append("w")
        self.varNames.append("e")
    else: #using a tabulated equation of state
      if self.numDims==1:
        self._varIDs["M_r"]=0
        self._varIDs["Delta_M_r"]=1
        self._varIDs["r"]=2
        self._varIDs["rho"]=3
        self._varIDs["u"]=4
        self._varIDs["u_0"]=5
        self._varIDs["T"]=6
        
        self.varNames.append("M_r")
        self.varNames.append("Delta_M_r")
        self.varNames.append("r")
        self.varNames.append("rho")
        self.varNames.append("u")
        self.varNames.append("u_0")
        self.varNames.append("T")

      elif self.numDims==2:
        self.varNames=["M_r","theta","Delta_M_r","r","rho","u","u_0","v","T"]
        i=0
        for varName in self.varNames:
          self._varIDs[varName]=i
          i+=1
        
        #set data type
        self.varType=dict.fromkeys(self.varNames,"hybrid")
        
        # self._varIDs["M_r"]=0
        # self._varIDs["theta"]=1
        # self._varIDs["Delta_M_r"]=2
        # self._varIDs["r"]=3
        # self._varIDs["rho"]=4
        # self._varIDs["u"]=5
        # self._varIDs["u_0"]=6
        # self._varIDs["v"]=7
        # self._varIDs["T"]=8
        
        # self.varNames.append("M_r")
        # self.varNames.append("theta")
        # self.varNames.append("Delta_M_r")
        # self.varNames.append("r")
        # self.varNames.append("rho")
        # self.varNames.append("u")
        # self.varNames.append("u_0")
        # self.varNames.append("v")
        # self.varNames.append("T")

      elif self.numDims==3:
        self._varIDs["M_r"]=0
        self._varIDs["theta"]=1
        self._varIDs["phi"]=2
        self._varIDs["Delta_M_r"]=3
        self._varIDs["r"]=4
        self._varIDs["rho"]=5
        self._varIDs["u"]=6
        self._varIDs["u_0"]=7
        self._varIDs["v"]=8
        self._varIDs["w"]=9
        self._varIDs["T"]=10
        
        self.varNames.append("M_r")
        self.varNames.append("theta")
        self.varNames.append("phi")
        self.varNames.append("Delta_M_r")
        self.varNames.append("r")
        self.varNames.append("rho")
        self.varNames.append("u")
        self.varNames.append("u_0")
        self.varNames.append("v")
        self.varNames.append("w")
        self.varNames.append("T")
  def _readHeaderBinary(self,eosFile=None):
    """Reads a header from a binary file, after the type has been read in.
    
    eosFile: allows one to override the equation of state file name given in the
      model.
    """
    
    #Unpack binary data
    self.version=struct.unpack('i',self.f.read(4))[0]#file version integer
    self.time=struct.unpack('d',self.f.read(8))[0]#time
    self.timeStepIndex=struct.unpack('i',self.f.read(4))[0]#time step index
    self.delta_t_nm1half=struct.unpack('d',self.f.read(8))[0]#delta t nm1half
    self.delta_t_np1half=struct.unpack('d',self.f.read(8))[0]#delta t np1half
    self.alpha=struct.unpack('d',self.f.read(8))[0]#alpha
    self.eosStringLen=struct.unpack('i',self.f.read(4))[0]#size of equation of
                                                          #state string
    if self.eosStringLen>0:
      self.eosString=""
      self.gamma=None
      
      #read in the equation of state string character by character
      for i in range(self.eosStringLen):
        self.eosString+=struct.unpack('s',self.f.read(1))[0]
      
      if eosFile!=None:#allows overriding eos file name from model
        self.eosString=eosFile
    else:
      self.eosString=None
      self.gamma=struct.unpack('d',self.f.read(8))[0]#get value of gamma
    self.av=struct.unpack('d',self.f.read(8))[0]#artificial viscosity
    self.avthreshold=struct.unpack('d',self.f.read(8))[0]#artificial viscosity threshold
    self.globalDims=[]
    self.globalDims.append(struct.unpack('i',self.f.read(4))[0])#x0 dim
    self.globalDims.append(struct.unpack('i',self.f.read(4))[0])#x1 dim
    self.globalDims.append(struct.unpack('i',self.f.read(4))[0])#x2 dim
    self.boundaryConditions=[]
    self.boundaryConditions.append(struct.unpack('i',self.f.read(4))[0])#x0 BC
    self.boundaryConditions.append(struct.unpack('i',self.f.read(4))[0])#x1 BC
    self.boundaryConditions.append(struct.unpack('i',self.f.read(4))[0])#x2 BC
    self.num1DZones=struct.unpack('i',self.f.read(4))[0]# number of 1D zones
    self.numGhostCells=struct.unpack('i',self.f.read(4))[0]#number of ghost cells at boundaries
    self.numVars=struct.unpack('i',self.f.read(4))[0]# number of variables
    
    #get variable infos, and set variable sizes
    self.varInfo=[]
    self.varSize=[]
    for i in range(self.numVars):
      tmp=[0,0,0,0]
      tmp[0]=struct.unpack('i',self.f.read(4))[0]
      tmp[1]=struct.unpack('i',self.f.read(4))[0]
      tmp[2]=struct.unpack('i',self.f.read(4))[0]
      tmp[3]=struct.unpack('i',self.f.read(4))[0]
      self.varInfo.append(tmp)
      tmp2=[0,0,0]
      for l in range(3):
        if self.globalDims[l]==1:
          #make sure if global dimensions are 1 in direction l, that varible 
          #isn't defined in that direction
          self.varInfo[i][l]==-1
        if self.varInfo[i][l]==-1:#not defined in direction l
          tmp2[l]=1
        elif self.varInfo[i][l]==1 and self.boundaryConditions[l]==0:
          #interface variable, and boundary is not periodic
          tmp2[l]=self.globalDims[l]+1
        else:
          tmp2[l]=self.globalDims[l]
      self.varSize.append(tmp2)
    
    #set number of dimensions
    self.numDims=0
    if self.globalDims[0]>1:
      self.numDims+=1
    if self.globalDims[1]>1:
      self.numDims+=1
    if self.globalDims[2]>1:
      self.numDims+=1
    
    #make a list for future variables
    self.vars=[]
  def _readHeaderAscii(self,eosFile=None):
    """Reads a header from a ascii file, after the type has been read in.
    
    eosFile: allows one to override the equation of state file name given in the
      model.
    """
    
    self.version        =int(self.f.readline())#file version integer
    lineSplit=self.f.readline().split()
    self.time           =float(lineSplit[0])#time
    self.timeStepIndex  =int(lineSplit[1])#time step index
    self.delta_t_nm1half=float(self.f.readline())#delta t nm1half
    self.delta_t_np1half=float(self.f.readline())#delta t np1half
    self.alpha          =float(self.f.readline())#alpha
    lineSplit=self.f.readline().split()
    self.eosStringLen   =int(lineSplit[0])#size of equation of state string
    
    if self.eosStringLen>0:
      self.gamma=None
      self.eosString=lineSplit[1]#get a character of the eos string
      
      if eosFile!=None:#allows overriding eos file name from model
        self.eosString=eosFile
    else:
      self.eosString=None
      self.gamma=float(lineSplit[1])#get value of gamma
    self.av             =float(self.f.readline())#artificial viscosity
    self.avthreshold    =float(self.f.readline())#artificial viscosity threshold
    lineSplit=self.f.readline().split()
    self.globalDims=[int(lineSplit[0]),int(lineSplit[1]),int(lineSplit[2])]
    lineSplit=self.f.readline().split()
    self.boundaryConditions=[int(lineSplit[0]),int(lineSplit[1]),int(lineSplit[2])]
    self.num1DZones=int(self.f.readline())# number of 1D zones
    self.numGhostCells=int(self.f.readline())# number of ghost cells at boundaries
    self.numVars=int(self.f.readline())# number of variables
    
    #get variable infos, and set variable sizes
    self.varInfo=[]
    self.varSize=[]
    lineSplit=self.f.readline().split()
    for i in range(self.numVars):
      tmp=[0,0,0,0]
      tmp[0]=int(lineSplit[i*4])
      tmp[1]=int(lineSplit[i*4+1])
      tmp[2]=int(lineSplit[i*4+2])
      tmp[3]=int(lineSplit[i*4+3])
      self.varInfo.append(tmp)
      tmp2=[0,0,0]
      for l in range(3):
        if self.globalDims[l]==1:
          #make sure if global dimensions are 1 in direction l, that varible isn't defined in that
          #direction
          self.varInfo[i][l]==-1
        if self.varInfo[i][l]==-1:#not defined in direction l
          tmp2[l]=1
        elif self.varInfo[i][l]==1 and self.boundaryConditions[l]==0:
          #interface variable, and boundary is not periodic
          tmp2[l]=self.globalDims[l]+1
        else:
          tmp2[l]=self.globalDims[l]
      self.varSize.append(tmp2)
    
    #set number of dimensions
    self.numDims=0
    if self.globalDims[0]>1:
      self.numDims+=1
    if self.globalDims[1]>1:
      self.numDims+=1
    if self.globalDims[2]>1:
      self.numDims+=1
    
    #make a list for future variables
    self.vars=[]
  def _readBinaryVar(self,var):
    """Read in a variable from a binary dump file.
    
    Must be called in order with var increasing from 0 to self.numVars.
    
    Arguments:
    var: variable to read, it is an integer index ranging from 0 to self.numVars
    """
    
    #set ghost cells based on which directions the variable is defined in
    ghostCellsInX0=1
    if self.varInfo[var][0]==-1:
      ghostCellsInX0=0
    ghostCellsInX1=1
    if self.varInfo[var][1]==-1:
      ghostCellsInX1=0
    ghostCellsInX2=1
    if self.varInfo[var][2]==-1:
      ghostCellsInX2=0
    
    #read 1D part
    sizeX01=ghostCellsInX0*(self.num1DZones+self.numGhostCells)
    if self.varInfo[var][0]==1 and self.boundaryConditions[0]==0:
      sizeX01=ghostCellsInX0*(self.num1DZones+self.numGhostCells+1)
    sizeX1=1
    sizeX2=1
    varTmp=[]
    for i in range(sizeX01):
      tmpj=[]
      for j in range(sizeX1):
        tmpk=[]
        for k in range(sizeX2):
          tmpk.append(struct.unpack('d',self.f.read(8))[0])
        tmpj.append(tmpk)
      varTmp.append(tmpj)
    
    #read multi-D part
    sizeX02=self.varSize[var][0]+ghostCellsInX0*2*self.numGhostCells
    sizeX1=self.varSize[var][1]+ghostCellsInX1*2*self.numGhostCells
    sizeX2=self.varSize[var][2]+ghostCellsInX2*2*self.numGhostCells
    for i in range(sizeX01,sizeX02):
      tmpj=[]
      for j in range(sizeX1):
        tmpk=[]
        for k in range(sizeX2):
          tmpk.append(struct.unpack('d',self.f.read(8))[0])
        tmpj.append(tmpk)
      varTmp.append(tmpj)
    self.vars.append(varTmp)
  def _readAsciiVar(self,var):
    """Read in a variable from an ascii dump file.
    
    Must be called in order with var increasing from 0 to self.numVars.
    
    Arguments:
    var: variable to read, it is an integer index ranging from 0 to self.numVars
    """
    
    #set ghost cells based on which directions the variable is defined in
    ghostCellsInX0=1
    if self.varInfo[var][0]==-1:
      ghostCellsInX0=0
    ghostCellsInX1=1
    if self.varInfo[var][1]==-1:
      ghostCellsInX1=0
    ghostCellsInX2=1
    if self.varInfo[var][2]==-1:
      ghostCellsInX2=0
    
    #read 1D part
    sizeX01=ghostCellsInX0*(self.num1DZones+self.numGhostCells)
    if self.varInfo[var][0]==1 and self.boundaryConditions[0]==0:
      sizeX01=ghostCellsInX0*(self.num1DZones+self.numGhostCells+1)
    sizeX1=1
    sizeX2=1
    varTmp=[]
    for i in range(sizeX01):
      tmpj=[]
      for j in range(sizeX1):
        tmpk=[]
        for k in range(sizeX2):
          line=self.f.readline()
          tmpk.append(float(line))
        tmpj.append(tmpk)
      self.f.readline()#toss a junk line
      varTmp.append(tmpj)
    
    #read multi-D part
    sizeX02=self.varSize[var][0]+ghostCellsInX0*2*self.numGhostCells
    sizeX1=self.varSize[var][1]+ghostCellsInX1*2*self.numGhostCells
    sizeX2=self.varSize[var][2]+ghostCellsInX2*2*self.numGhostCells
    for i in range(sizeX01,sizeX02):
      tmpj=[]
      for j in range(sizeX1):
        tmpk=[]
        line=self.f.readline()
        lineSplit=line.split()
        for k in range(sizeX2):
          tmpk.append(float(lineSplit[k]))
          tmpj.append(tmpk)
      varTmp.append(tmpj)
      self.f.readline()#toss a junk line
    self.vars.append(varTmp)
    self.f.readline()#toss a junk line
  def _printVarByID(self,var,out):
    """Writes a non-rectangular variable to out."""
    
    #set ghost cells based on which directions the variable is defined in
    ghostCellsInX0=1
    if self.varInfo[var][0]==-1:
      ghostCellsInX0=0
    ghostCellsInX1=1
    if self.varInfo[var][1]==-1:
      ghostCellsInX1=0
    ghostCellsInX2=1
    if self.varInfo[var][2]==-1:
      ghostCellsInX2=0
    
    #print out 1D part
    sizeX01=ghostCellsInX0*(self.num1DZones+self.numGhostCells)
    if self.varInfo[var][0]==1 and self.boundaryConditions[0]==0:
      sizeX01=ghostCellsInX0*(self.num1DZones+self.numGhostCells+1)
    sizeX1=1
    sizeX2=1
    for i in range(sizeX01):
      for j in range(sizeX1):
        for k in range(sizeX2):
          out.write("var["+str(i)+"]["+str(j)+"]["+str(k)+"]="
            +str(self.vars[var][i][j][k])+"\n")
    
    #print out MD part
    sizeX02=self.varSize[var][0]+ghostCellsInX0*2*self.numGhostCells
    sizeX1=self.varSize[var][1]+ghostCellsInX1*2*self.numGhostCells
    sizeX2=self.varSize[var][2]+ghostCellsInX2*2*self.numGhostCells
    for i in range(sizeX01,sizeX02):
      for j in range(sizeX1):
        for k in range(sizeX2):
          out.write("var["+str(i)+"]["+str(j)+"]["+str(k)+"]="
            +str(self.vars[var][i][j][k])+"\n")
  def _adjustForPeriodBC(self):
    """Adjust for periodic boundary conditions in theta and phi directions
    
    In the theta and phi directions the theta and phi velocity componenets are
    copies of those at the other boundary; thus even though they are interface
    centered quantities they have the same number of elements as a zone centered
    quantity in that direction. Thus function copies the boundaries values of
    the velocities and cooridnates so this doesn't need to be worried about
    """
    
    for n in range(len(self.varInfo)):
      if self.varInfo[n][1]==1:
        
        if self.varNames[n]!="theta":
          
          #copy j=0 to j=N
          shape=self.rectVars[n].shape
          self.rectVars[n].resize((shape[0],shape[1]+1,shape[2]))
          shape=self.rectVars[n].shape
          for i in range(shape[0]):
            for k in range(shape[2]):
              self.rectVars[n][i][shape[1]-1][k]=self.rectVars[n][i][0][k]
        else:
          
          #copy j=0 to j=N
          shape=self.rectVars[n].shape
          self.rectVars[n].resize((shape[0],shape[1]+1,shape[2]))
          shape=self.rectVars[n].shape
          dTheta=(self.rectVars[n][0][1][0]-self.rectVars[n][0][0][0])
          for i in range(shape[0]):
            for k in range(shape[2]):
              self.rectVars[n][i][shape[1]-1][k]=(
                self.rectVars[n][i][shape[1]-2][k]+dTheta)
      if self.varInfo[n][2]==1:

        if self.varNames[n]!="phi":
          
          #copy k=0 to k=N
          shape=self.rectVars[n].shape
          self.rectVars[n].resize((shape[0],shape[1],shape[2]+1))
          shape=self.rectVars[n].shape
          for i in range(shape[0]):
            for j in range(shape[1]):
              self.rectVars[n][i][j][shape[2]-1]=self.rectVars[n][i][j][0]
        else:
          
          #copy j=0 to j=N
          shape=self.rectVars[n].shape
          self.rectVars[n].resize((shape[0],shape[1],shape[2]+1))
          shape=self.rectVars[n].shape
          dPhi=(self.rectVars[n][0][0][1]-self.rectVars[n][0][0][0])
          for i in range(shape[0]):
            for j in range(shape[1]):
              self.rectVars[n][i][j][shape[2]-1]=(
                self.rectVars[n][i][j][shape[2]-2]+dPhi)
  def readHeader(self,eosFile=None):
    """Reads header information from binary dump file.
    
    This version calls either the _readHeaderAscii or the _readHeaderBinary
    eosFile: allows one to override the equation of state file name given in the
      model.
    """
    
    #read header
    self.type=struct.unpack('c',self.f.read(1))[0]#file type, either a or b
    if self.type=='a':
      self._readHeaderAscii(eosFile=eosFile)
    else:
      self._readHeaderBinary(eosFile=eosFile)
    
    #set shapes
    includeGhostR=1
    includeGhostTheta=0
    includeGhostPhi=0
    if self.numDims==2:#2D
      includeGhostTheta=1
      includeGhostPhi=0
    elif self.numDims==3:#3D
      includeGhostTheta=1
      includeGhostPhi=1
    self.cellCenteredShape=(
      self.globalDims[0]+includeGhostR*(2*self.numGhostCells)
      ,self.globalDims[1]+includeGhostTheta*(2*self.numGhostCells)
      ,self.globalDims[2]+includeGhostPhi*(2*self.numGhostCells))
    self.meshShape=(
      self.globalDims[0]+includeGhostR*(2*self.numGhostCells+1)
      ,self.globalDims[1]+includeGhostTheta*(2*self.numGhostCells+1)
      ,self.globalDims[2]+includeGhostPhi*(2*self.numGhostCells)+1)
  def read(self,fileName,eosFile=None):
    """Reads in a combined dump file and puts the variables into the vars list
    
    It also converts the hybrid 1D-Multi-grid into a rectangular grid
    eosFile: allows one to override the equation of state file name given in the
      model.
    """
    
    self.fileName=fileName
    self.f=open(fileName,'rb')
    self.readHeader(eosFile=eosFile)
    if self.type=='b':
      for i in range(self.numVars):
        self._readBinaryVar(i)
    elif self.type=='a':
      for i in range(self.numVars):
        self._readAsciiVar(i)
    self._setVarIDs()#set variable names/id
    self.setRectVars()#create rectangular variables
    self._adjustForPeriodBC()#adjust interface variables effected by periodic BC
  def printHeader(self,out):
    """Writes the header of a binary dump file to out.
    
    Arguments:
    out: an object supporting the write() function
    """
    
    #print header
    out.write("type="+str(self.type)+"\n")
    out.write("version="+str(self.version)+"\n")
    out.write("time="+str(self.time)+"\n")
    out.write("timeStepIndex="+str(self.timeStepIndex)+"\n")
    out.write("delta_t_nm1half="+str(self.delta_t_nm1half)+"\n")
    out.write("delta_t_np1half="+str(self.delta_t_np1half)+"\n")
    out.write("alpha="+str(self.alpha)+"\n")
    out.write("eosStringLen="+str(self.eosStringLen)+"\n")
    if self.eosStringLen>0:
      out.write("eosString="+str(self.eosString)+"\n")
    else:
      out.write("gamma="+str(self.gamma)+"\n")
    out.write("av="+str(self.av)+"\n")
    out.write("avthreshold="+str(self.avthreshold)+"\n")
    out.write("globalDims="+str(self.globalDims)+"\n")
    out.write("boundaryConditions="+str(self.boundaryConditions)+"\n")
    out.write("num1DZones="+str(self.num1DZones)+"\n")
    out.write("numGhostCells="+str(self.numGhostCells)+"\n")
    out.write("numVars="+str(self.numVars)+"\n")
    out.write("varInfo="+str(self.varInfo)+"\n")
    out.write("varSize="+str(self.varSize)+"\n")
    out.write("numDims="+str(self.numDims)+"\n")
  def printDumpToSTDOut(self):
    """Prints dump to standard output"""
    
    self.printHeader(sys.stdout)
    for i in range(self.numVars):
      self._printVarByID(i,sys.stdout)
  def getAllVarNames(self):
    """Returns a list of variable names that are availble."""
    
    return self._varIDs.keys()
  def getRectVarNames(self):
    """Returns a list of variable names that are rectangular.
    
    That have a rectangular representation set
    """
    
    #make a list of all variable keys that aren't None
    varNames=[]
    for key in self._varIDs.keys():
      if self.varType[key]=="rect" or self.varType[key]=="both":
        varNames.append(key)
    return varNames
  def getHybridVarNames(self):
    """Returns a list of variable names that are hybrid 
    
    i.e. they have a 1D and a multi-D grid.
    
    """
    
    #make a list of all variable keys that aren't None
    varNames=[]
    for key in self._varIDs.keys():
      if self.varType[key]=="hybrid" or self.varType[key]=="both":
        varNames.append(key)
    return varNames
  def getVarID(self,var):
    """Returns the array index (ID) of a variable given by name var
    
    var: a string name for variable
    """
    
    return self._varIDs[var]
  def printVarToOut(self,var,out):
    """Print var to out"""
    
    varNames=self.getVarNames()
    if var not in varNames:
      raise Exception(str(var)+" not in available variable names "+str(varNames))
    
    self._printVarByID(self._varIDs[var],out)
  def printVarToSTDOut(self,var):
    """Prints varible to standard output"""
    
    self.printVarToOut(var,sys.stdout)
  def setRectVars(self):
    """Creates rectangular representations of variables for later use
    
    """
    
    for i in range(len(self.vars)):
      self.rectVars.append(self.getRectVar(i))
      self.varType[self.varNames[i]]="both"
  def getRectVar(self,var):
    """Returns a rectangular numpy array version of a varible.
    
    Variables are stored as a 1D part plus a 2D or 3D part. This function
    returns a variable that has the 1D part copied to match the 2D or 3D part.
    
    var: an index or name of a variable
    """
    
    strType=type("")
    if type(var)==strType:#if it is a string assume it is a variable name and
                          #get the ID
      var=self.getVarID(var)
    
    #if the length of the recVars array is larger than array index var don't
    #need to re-rectangularize the variable
    if len(self.rectVars)>var:
      return self.rectVars[var]
      
    #set ghost cells based on which directions the variable is defined in
    ghostCellsInX0=1
    if self.varInfo[var][0]==-1:
      ghostCellsInX0=0
    ghostCellsInX1=1
    if self.varInfo[var][1]==-1:
      ghostCellsInX1=0
    ghostCellsInX2=1
    if self.varInfo[var][2]==-1:
      ghostCellsInX2=0
    
    #set sizes of 1D grid (sizeXD1, D=0,1,2)
    sizeX01=ghostCellsInX0*(self.num1DZones+self.numGhostCells)
    if self.varInfo[var][0]==1 and self.boundaryConditions[0]==0:
      sizeX01=ghostCellsInX0*(self.num1DZones+self.numGhostCells+1)
    sizeX11=1
    sizeX21=1
    
    #set sizes of multi-D grid (sizeXD2, D=0,1,2)
    sizeX02=self.varSize[var][0]+ghostCellsInX0*2*self.numGhostCells
    sizeX12=self.varSize[var][1]+ghostCellsInX1*2*self.numGhostCells
    sizeX22=self.varSize[var][2]+ghostCellsInX2*2*self.numGhostCells
    
    #calculate rectangular size
    shapeRec=[sizeX02,sizeX12,sizeX22]
    
    #allocate numpy array
    rectVar=np.empty(shapeRec)
    
    #convert 1D part to 2D/3D
    for i in range(sizeX01):
      for j in range(sizeX12):
        for k in range(sizeX22):
          rectVar[i][j][k]=self.vars[var][i][0][0]
    
    #print out MD part
    for i in range(sizeX01,sizeX02):
      for j in range(sizeX12):
        for k in range(sizeX22):
          rectVar[i][j][k]=self.vars[var][i][j][k]
    return rectVar
  def getVarSlice(self,var,rIndexMin=0,rIndexMax=None,tIndexMin=0
                  ,tIndexMax=None,pIndexMin=0,pIndexMax=None):
    """Returns a 3D numpy array for variable slice.
    
    Returns a 3D numpy array for the named variable within the given range
    
    Arguments:
    var: variable name to get
    
    Keyword Arguements:
    rIndexMin: minimum radial index to include
    rIndexMax: maximum radial index-1 to include
    tIndexMin: minimum theta index to include
    tIndexMax: maximum theta index-1 to include
    pIndexMin: minimum theta index to include
    pIndexMax: maximum theta index-1 to include
    """
    
    #get a rectangular version
    variableToSlice=self.getRecVar(self.getVarID(var))
    
    #set default shape options, if max's are Nones use full extent
    shape=variableToSlice.shape
    if rIndexMax==None:
      rIndexMax=shape[0]
    if tIndexMax==None:
      tIndexMax=shape[1]
    if pIndexMax==None:
      pIndexMax=shape[2]
    
    return variableToSlice[rIndexMin:rIndexMax,tIndexMin:tIndexMax
      ,pIndexMin:pIndexMax]
  def getHorAveVarSlice(self,var,rIndexMin=0,rIndexMax=None
                                    ,tIndexMin=0 ,tIndexMax=None,pIndexMin=0
                                    ,pIndexMax=None):
    """Returns a 1D numpy array for the horizontal average of the variable slice.
    
    Returns a 1D numpy array for the named variable within the given range
    averaged in the horizontal direction so that the 1D array is a function of
    radius only.
    
    Arguments:
    var: variable name to get
    
    Keyword Arguements:
    rIndexMin: minimum radial index to include
    rIndexMax: maximum radial index-1 to include
    tIndexMin: minimum theta index to include
    tIndexMax: maximum theta index-1 to include
    pIndexMin: minimum theta index to include
    pIndexMax: maximum theta index-1 to include
    """
    
    #get a variable slice
    varSliced=self.getVarSlice(var,rIndexMin=rIndexMin,rIndexMax=rIndexMax
                  ,tIndexMin=tIndexMin,tIndexMax=tIndexMax,pIndexMin=pIndexMin
                  ,pIndexMax=pIndexMax)
    
    #Perform a horizontal average
    varHorAved=varSliced.sum(axis=1)/varSliced.shape[1]
    
    #make 3D numpy array to be more universally usable
    temp=varHorAved[:,:,np.newaxis]
    
    return temp
  def printVarToRadCol(self,temp,var,out,tIndexMin=0,pIndexMin=0):
    """Writes the numpy array temp to out
    
    It formats the write to out in a column wise manner with the first index
    written as the first column, and the second two as a column header seperated
    by a comma.
    
    Arguments:
    temp: 3D numpy array to print out
    var: string varaible name
    out: object supporting the write function
    
    Keyword Arguments:
    tIndexMin: integer, mimimum theta index of slice
    pIndexMin: integer, mimimum phi index of slice
    """
    
    #write temp to out
    out.write("varible=\""+var+"\"\n")
    columnWidthFloat=28
    columnWidthInt=8
    precision=15
    
    #write out header
    line=""
    columnFormatStringFloatWidth="{0: >"+str(columnWidthFloat)+"}"
    columnFormatStringIntWidth="{0: >"+str(columnWidthInt)+"}"
    columnFormatFloat="{0: >"+str(columnWidthFloat)+"."+str(precision)+"e}"
    columnFormatInt="{0: >"+str(columnWidthInt)+"d}"
    line+=columnFormatStringIntWidth.format("zone#")
    for j in range(temp.shape[1]):
      for k in range(temp.shape[2]):
        line+=columnFormatStringFloatWidth.format( (str(j+tIndexMin)+","
          +str(k+pIndexMin)) )
    out.write(line+"\n")
    
    #write out variable slice
    for i in range(temp.shape[0]):
      line=""
      line+=columnFormatInt.format(i)
      for j in range(temp.shape[1]):
        for k in range(temp.shape[2]):
          line+=columnFormatFloat.format(temp[i][j][k])
      out.write(line+"\n")
  def printHorAveVarSliceToOutInRadCol(self,var,out,rIndexMin=0,rIndexMax=None
                                        ,tIndexMin=0,tIndexMax=None,pIndexMin=0
                                        ,pIndexMax=None):
    """
    Prints a horizontally averaged variable slice to out in a radial column
    format.
    
    Arguments:
    var: variable name to write in column format
    out: object with a write method, such as a file object or sys.stdout
    
    Keyword Arguements:
    rIndexMin: minimum radial index to include
    rIndexMax: maximum radial index-1 to include
    tIndexMin: minimum theta index to include
    tIndexMax: maximum theta index-1 to include
    pIndexMin: minimum theta index to include
    pIndexMax: maximum theta index-1 to include
    """
    
    temp=self.getHorAveVarSlice(var,rIndexMin=rIndexMin,rIndexMax=rIndexMax
                  ,tIndexMin=tIndexMin,tIndexMax=tIndexMax,pIndexMin=pIndexMin
                  ,pIndexMax=pIndexMax)
    
    self.printVarToRadCol(temp,var,out)
  def printVarSliceToOutInRadCol(self,var,out,rIndexMin=0,rIndexMax=None
                                        ,tIndexMin=0,tIndexMax=None,pIndexMin=0
                                        ,pIndexMax=None):
    """
    Prints a variable slice to out in a radial column format.
    
    The columns will be for different j (theta) and k (phi) values while the
    rows cover the i (radial) indices.
    
    Arguments:
    var: variable name to write in column format
    out: object with a write method, such as a file object or sys.stdout
    
    Keyword Arguements:
    rIndexMin: minimum radial index to include
    rIndexMax: maximum radial index-1 to include
    tIndexMin: minimum theta index to include
    tIndexMax: maximum theta index-1 to include
    pIndexMin: minimum theta index to include
    pIndexMax: maximum theta index-1 to include
    """
    
    temp=self.getVarSlice(var,rIndexMin=rIndexMin,rIndexMax=rIndexMax
                  ,tIndexMin=tIndexMin,tIndexMax=tIndexMax,pIndexMin=pIndexMin
                  ,pIndexMax=pIndexMax)
    
    self.printVarToRadCol(temp,var,out,tIndexMin=tIndexMin,pIndexMin=pIndexMin)
  def writeVTKFile(self,fileName,withGhost=True,eosFile=None,curvature=True
    ,excludedScalars=[]
    ,includeScalars=None
    ,vectors=None
    ,coords=None
    ):
    """Writes out a vtk file from the model dump
    
    parameters:
      fileName: model dump to create vtk file from
    keywords:
      withGhost: if true will output ghost cells, if false it won't. Currently
        not implemented
      excludedScalars: a list of scalars to exclude from vtk file
      includeScalars: a list of scalar varaibles to include, if not set will
        output all available scalars. Possible scalars are: "dlnPdlnT"
        ,"dlnPdlnD","dEdT","p","e","kappa","gamma","cp","c","gamma1","DelAd"
        ,"cv"
      vectors: A dictonary with a vector name as a key, and a 3 component list
        of strings indicating the names of the components to be used for the
        vector. Default value depends on curvature and is either 
        {"vel":["vx_cen","vy_cen","vz_cen"]
        ,"vel_con":["vx_con_cen","vy_con_cen","vz_con_cen"]}
        or
        {"vel":["vr_cen","vt_cen","vp_cen"]
        ,"vel_con":["vr_con_cen","vt_con_cen","vp_con_cen"]}
      coords: A dictonary with a x,y,z keys, each indicating the name of the
        variable to be used for that component of the coordinate. Default value
        depends on the curvature and is either
        {"x":"x_mesh","y":"y_mesh","z":"z_mesh"}
        or
        {"x":"r_mesh","y":"theta_mesh","z":"phi_mesh"}
    """
    
    #set vtk file vectors and coordinates to reasonable defaults if not
    #explicitly set
    if curvature:
      if vectors==None:
        vectors={"vel":["vx_cen","vy_cen","vz_cen"]
          ,"vel_con":["vx_con_cen","vy_con_cen","vz_con_cen"]}
      if coords==None:
        coords={"x":"x_mesh","y":"y_mesh","z":"z_mesh"}
    else:
      if vectors==None:
        vectors={"vel":["vr_cen","vt_cen","vp_cen"]
          ,"vel_con":["vr_con_cen","vt_con_cen","vp_con_cen"]}
      if coords==None:
        coords={"x":"r_mesh","y":"theta_mesh","z":"phi_mesh"}
    
    #creates centered, r,theta, phi coordinates
    self.calCenCoords()
    
    #creates mesh coordinates either x,y,z or r,theta,phi depending on curvature
    self.calMesh(curvature=curvature)
    
    #create cell centered velocity
    self.calZoneCenteredVelocity(curvature=curvature)
    if curvature==True and ("vr_cen" in includeScalars
      or "vt_cen" in includeScalars
      or "vp_cen" in includeScalars):
      self.calZoneCenteredVelocity(curvature=False)
    if curvature==False and ("vx_cen" in includeScalars
      or "vy_cen" in includeScalars
      or "vz_cen" in includeScalars):
      self.calZoneCenteredVelocity(curvature=True)
      
    #create cell centered convective velocity
    self.calZoneCenteredConvVelocity(curvature=curvature)
    
    #if variables of opposite curvature are requested
    if curvature==True and ("vr_con_cen" in includeScalars
      or "vt_con_cen" in includeScalars
      or "vp_con_cen" in includeScalars):
      self.calZoneCenteredConvVelocity(curvature=False)
    if curvature==False and ("vx_con_cen" in includeScalars
      or "vy_con_cen" in includeScalars
      or "vz_con_cen" in includeScalars):
      self.calZoneCenteredConvVelocity(curvature=True)
      
    
    #set thermodynamic variables from equation of state/opacity table
    self.setAdditionalVars(varsToSet=includeScalars)
    
    #set coordinates
    x=self.rectVars[self.getVarID(coords["x"])]
    y=self.rectVars[self.getVarID(coords["y"])]
    z=self.rectVars[self.getVarID(coords["z"])]
    
    #get cell centered vars
    cellCenteredVars=self.getCellCenteredRect()
    if includeScalars==None:
      includeScalars=cellCenteredVars
    
    #add cell centered data
    cellData={}
    
    #Add Scalars
    for var in cellCenteredVars:
      includeScalar=True
      for vector in vectors:
        if var in vectors[vector]:
          includeScalar=False
          break
      if var in excludedScalars or var not in includeScalars:
        includeScalar=False
      if includeScalar:#add it
        cellData[var]=self.rectVars[self.getVarID(var)]
    
    #add vectors
    for vector in vectors:
      cellData[vector]=(self.rectVars[self.getVarID(vectors[vector][0])]
        ,self.rectVars[self.getVarID(vectors[vector][1])]
        ,self.rectVars[self.getVarID(vectors[vector][2])])
    
    evtk.hl.gridToVTK(fileName, x, y, z, cellData=cellData, pointData = {})
  def getMeshCoords(self):
    """Returns a list of rectangular variables defined at all cell corners
    
    """
    
    cellCornerRectVars=[]
    for i in range(len(self.rectVars)):
      if self.rectVars[i].shape==self.meshShape:
        cellCornerRectVars.append(self.varNames[i])
    return cellCornerRectVars
  def getCellCenteredRect(self):
    """Returns a list of rectangular variables defined at cell centres
    
    """
    
    cellCenteredRectVars=[]
    for i in range(len(self.rectVars)):
      if self.rectVars[i].shape==self.cellCenteredShape:
        cellCenteredRectVars.append(self.varNames[i])
    return cellCenteredRectVars
  def calCenCoords(self):
    """Calculates cell centered coordinates
    
    """
    
    if self.numDims>0:
      
      #add r_cen
      self._addVarID("r_cen","rect")
      varTmp=np.empty((self.cellCenteredShape[0],1,1))
      for i in range(self.cellCenteredShape[0]):
        varTmp[i][0][0]=(self.rectVars[self.getVarID("r")][i][0][0]
          +self.rectVars[self.getVarID("r")][i+1][0][0])*0.5
      self.rectVars.append(varTmp)
    if self.numDims>1:
      
      #add theta_cen
      self._addVarID("theta_cen","rect")
      varTmp=np.empty((1,self.cellCenteredShape[1],1))
      for j in range(self.cellCenteredShape[1]):
        varTmp[0][j][0]=(self.rectVars[self.getVarID("theta")][0][j][0]
          +self.rectVars[self.getVarID("theta")][0][j+1][0])*0.5
      self.rectVars.append(varTmp)
    if self.numDims>2:
      
      #add phi_cen
      self._addVarID("phi_cen","rect")
      varTmp=np.empty((1,1,self.cellCenteredShape[2]))
      for k in range(self.cellCenteredShape[2]):
        varTmp[0][0][k]=(self.rectVars[self.getVarID("phi")][0][0][k]
          +self.rectVars[self.getVarID("phi")][0][0][k+1])*0.5
      self.rectVars.append(varTmp)
  def calMesh(self,curvature=True):
    
    if curvature:
      
      #register the variables
      self._addVarID("x_mesh","rect")
      self._addVarID("y_mesh","rect")
      self._addVarID("z_mesh","rect")
      if self.numDims==1:
        
        #add x_cen
        varTmpX=np.empty(self.meshShape)
        for i in range(self.meshShape[0]):
          for j in range(self.meshShape[1]):
            for k in range(self.meshShape[2]):
              varTmpX[i][j][k]=self.rectVars[self.getVarID("r")][i][0][0]
        self.rectVars.append(varTmpX)
        
        #add y_cen
        self.rectVars.append(np.zeros(self.meshShape))
        
        #add z_cen
        self.rectVars.append(np.zeros(self.meshShape))
        
      elif self.numDims==2:
        
        #add x_cen and y_cen
        varTmpX=np.empty(self.meshShape)
        varTmpY=np.empty(self.meshShape)
        for i in range(self.meshShape[0]):
          for j in range(self.meshShape[1]):
            for k in range(self.meshShape[2]):
              r=self.rectVars[self.getVarID("r")][i][0][0]
              theta=self.rectVars[self.getVarID("theta")][0][j][0]
              varTmpX[i][j][k]=r*math.sin(theta)
              varTmpY[i][j][k]=r*math.cos(theta)
        self.rectVars.append(varTmpX)
        self.rectVars.append(varTmpY)
        
        #add z_cen
        self.rectVars.append(np.zeros(self.meshShape))
      elif self.numDims==3:
        
        #add x_cen, y_cen and z_cen
        varTmpX=np.empty(self.meshShape)
        varTmpY=np.empty(self.meshShape)
        varTmpZ=np.empty(self.meshShape)
        for i in range(self.meshShape[0]):
          for j in range(self.meshShape[1]):
            for k in range(self.meshShape[2]):
              r=self.rectVars[self.getVarID("r")][i][0][0]
              theta=self.rectVars[self.getVarID("theta")][0][j][0]
              phi=self.rectVars[self.getVarID("phi")][0][0][k]
              varTmpX[i][j][k]=r*math.sin(theta)*math.cos(phi)
              varTmpY[i][j][k]=r*math.sin(theta)*math.sin(phi)
              varTmpZ[i][j][k]=r*math.cos(theta)
        self.rectVars.append(varTmpX)
        self.rectVars.append(varTmpY)
        self.rectVars.append(varTmpZ)
    else:
      
      #register the variables
      self._addVarID("r_mesh","rect")
      self._addVarID("theta_mesh","rect")
      self._addVarID("phi_mesh","rect")
      if self.numDims==1:
        
        #add x_cen
        varTmpX=np.empty(self.meshShape)
        for i in range(self.meshShape[0]):
          for j in range(self.meshShape[1]):
            for k in range(self.meshShape[2]):
              varTmpX[i][j][k]=self.rectVars[self.getVarID("r")][i][0][0]
        self.rectVars.append(varTmpX)
        
        #add y_cen
        self.rectVars.append(np.zeros(self.meshShape))
        
        #add z_cen
        self.rectVars.append(np.zeros(self.meshShape))
        
      elif self.numDims==2:
        
        #add x_cen and y_cen
        varTmpX=np.empty(self.meshShape)
        varTmpY=np.empty(self.meshShape)
        for i in range(self.meshShape[0]):
          for j in range(self.meshShape[1]):
            for k in range(self.meshShape[2]):
              r=self.rectVars[self.getVarID("r")][i][0][0]
              theta=self.rectVars[self.getVarID("theta")][0][j][0]
              varTmpX[i][j][k]=r
              varTmpY[i][j][k]=theta
        self.rectVars.append(varTmpX)
        self.rectVars.append(varTmpY)
        
        #add z_cen
        self.rectVars.append(np.zeros(self.meshShape))
      elif self.numDims==3:
        
        #add x_cen, y_cen and z_cen
        varTmpX=np.empty(self.meshShape)
        varTmpY=np.empty(self.meshShape)
        varTmpZ=np.empty(self.meshShape)
        for i in range(self.meshShape[0]):
          for j in range(self.meshShape[1]):
            for k in range(self.meshShape[2]):
              r=self.rectVars[self.getVarID("r")][i][0][0]
              theta=self.rectVars[self.getVarID("theta")][0][j][0]
              phi=self.rectVars[self.getVarID("phi")][0][0][k]
              varTmpX[i][j][k]=r
              varTmpY[i][j][k]=theta
              varTmpZ[i][j][k]=phi
        self.rectVars.append(varTmpX)
        self.rectVars.append(varTmpY)
        self.rectVars.append(varTmpZ)
  def calZoneCenteredConvVelocity(self,curvature=True):
    """Zone centered velocity minus the grid velocity
    
    """
    
    if curvature:
      
      #register the variables
      self._addVarID("vx_con_cen","rect")
      self._addVarID("vy_con_cen","rect")
      self._addVarID("vz_con_cen","rect")
      if self.numDims==1:
        
        #add vx_cen
        varTmpX=np.empty(self.cellCenteredShape)
        for i in range(self.cellCenteredShape[0]):
          for j in range(self.cellCenteredShape[1]):
            for k in range(self.cellCenteredShape[2]):
              varTmpX[i][j][k]=(self.rectVars[self.getVarID("u")][i][j][k]
                -self.rectVars[self.getVarID("u_0")][i][0][0]
                +self.rectVars[self.getVarID("u")][i+1][j][k]
                -self.rectVars[self.getVarID("u_0")][i+1][0][0])*0.5
        self.rectVars.append(varTmpX)
        
        #add vy_cen
        self.rectVars.append(np.zeros(self.cellCenteredShape))
        
        #add vz_cen
        self.rectVars.append(np.zeros(self.cellCenteredShape))
      elif self.numDims==2:
        
        #add x_cen and y_cen
        varTmpX=np.empty(self.cellCenteredShape)
        varTmpY=np.empty(self.cellCenteredShape)
        for i in range(self.cellCenteredShape[0]):
          for j in range(self.cellCenteredShape[1]):
            for k in range(self.cellCenteredShape[2]):
              vr=(self.rectVars[self.getVarID("u")][i][j][k]
                -self.rectVars[self.getVarID("u_0")][i][0][0]
                +self.rectVars[self.getVarID("u")][i+1][j][k]
                -self.rectVars[self.getVarID("u_0")][i+1][0][0])*0.5
              vt=(self.rectVars[self.getVarID("v")][i][j][k]
                +self.rectVars[self.getVarID("v")][i][j+1][k])*0.5
              theta=self.rectVars[self.getVarID("theta_cen")][0][j][0]
              varTmpX[i][j][k]=vr*math.sin(theta)+vt*math.cos(theta)
              varTmpY[i][j][k]=vr*math.cos(theta)+vt*math.sin(theta)
        self.rectVars.append(varTmpX)
        self.rectVars.append(varTmpY)
        
        #add z_cen
        self.rectVars.append(np.zeros(self.cellCenteredShape))
      elif self.numDims==3:
        
        #add x_cen, y_cen and z_cen
        varTmpX=np.empty(self.cellCenteredShape)
        varTmpY=np.empty(self.cellCenteredShape)
        varTmpZ=np.empty(self.cellCenteredShape)
        sinTheta=np.empty((self.cellCenteredShape[1]))
        cosTheta=np.empty((self.cellCenteredShape[1]))
        sinPhi=np.empty((self.cellCenteredShape[2]))
        cosPhi=np.empty((self.cellCenteredShape[2]))
        for j in range((self.cellCenteredShape[1])):
          sinTheta[j]=math.sin(self.rectVars[self.getVarID("theta_cen")][0][j][0])
          cosTheta[j]=math.cos(self.rectVars[self.getVarID("theta_cen")][0][j][0])
        for k in range((self.cellCenteredShape[2])):
          sinPhi[k]=math.sin(self.rectVars[self.getVarID("phi_cen")][0][0][k])
          cosPhi[k]=math.cos(self.rectVars[self.getVarID("phi_cen")][0][0][k])
        for i in range(self.cellCenteredShape[0]):
          for j in range(self.cellCenteredShape[1]):
            for k in range(self.cellCenteredShape[2]):
              vr=(self.rectVars[self.getVarID("u")][i][j][k]
                -self.rectVars[self.getVarID("u_0")][i][0][0]
                +self.rectVars[self.getVarID("u")][i+1][j][k]
                -self.rectVars[self.getVarID("u_0")][i+1][0][0])*0.5
              vt=(self.rectVars[self.getVarID("v")][i][j][k]
                +self.rectVars[self.getVarID("v")][i][j+1][k])*0.5
              vp=(self.rectVars[self.getVarID("w")][i][j][k]
                +self.rectVars[self.getVarID("w")][i][j][k+1])*0.5
              varTmpX[i][j][k]=(
                vr*sinTheta[j]*cosPhi[k]
                +vt*cosTheta[j]*cosPhi[k]
                +vp*sinPhi[k])
              varTmpY[i][j][k]=(
                vr*sinTheta[j]*sinPhi[k]
                +vt*cosTheta[j]*sinPhi[k]
                +vp*cosTheta[j])
              varTmpZ[i][j][k]=(
                vr*cosTheta[j]
                -vt*sinTheta[j])
        self.rectVars.append(varTmpX)
        self.rectVars.append(varTmpY)
        self.rectVars.append(varTmpZ)
    else:
      
      #register the variables
      self._addVarID("vr_con_cen","rect")
      self._addVarID("vt_con_cen","rect")
      self._addVarID("vp_con_cen","rect")
      if self.numDims==1:
        
        #add x_cen
        varTmpX=np.empty(self.cellCenteredShape)
        for i in range(self.cellCenteredShape[0]):
          for j in range(self.cellCenteredShape[1]):
            for k in range(self.cellCenteredShape[2]):
              varTmpX[i][j][k]=(self.rectVars[self.getVarID("r")][i][j][k]
                +self.rectVars[self.getVarID("r")][i+1][j][k])*0.5
        self.rectVars.append(varTmpX)
        
        #add y_cen
        self.rectVars.append(np.zeros(self.cellCenteredShape))
        
        #add z_cen
        self.rectVars.append(np.zeros(self.cellCenteredShape))
        
      elif self.numDims==2:
        
        #add x_cen and y_cen
        varTmpX=np.empty(self.cellCenteredShape)
        varTmpY=np.empty(self.cellCenteredShape)
        for i in range(self.cellCenteredShape[0]):
          for j in range(self.cellCenteredShape[1]):
            for k in range(self.cellCenteredShape[2]):
              vr=(self.rectVars[self.getVarID("u")][i][j][k]
                -self.rectVars[self.getVarID("u_0")][i][0][0]
                +self.rectVars[self.getVarID("u")][i+1][j][k]
                -self.rectVars[self.getVarID("u_0")][i+1][0][0])*0.5
              vt=(self.rectVars[self.getVarID("v")][i][j][k]
                +self.rectVars[self.getVarID("v")][i][j+1][k])*0.5
              varTmpX[i][j][k]=vr
              varTmpY[i][j][k]=vt
        self.rectVars.append(varTmpX)
        self.rectVars.append(varTmpY)
        
        #add z_cen
        self.rectVars.append(np.zeros(self.cellCenteredShape))
      elif self.numDims==3:
        
        #add x_cen, y_cen and z_cen
        varTmpX=np.empty(self.cellCenteredShape)
        varTmpY=np.empty(self.cellCenteredShape)
        varTmpZ=np.empty(self.cellCenteredShape)
        for i in range(self.cellCenteredShape[0]):
          for j in range(self.cellCenteredShape[1]):
            for k in range(self.cellCenteredShape[2]):
              vr=(self.rectVars[self.getVarID("u")][i][j][k]
                -self.rectVars[self.getVarID("u_0")][i][0][0]
                +self.rectVars[self.getVarID("u")][i+1][j][k]
                -self.rectVars[self.getVarID("u_0")][i+1][0][0])*0.5
              vt=(self.rectVars[self.getVarID("v")][i][j][k]
                +self.rectVars[self.getVarID("v")][i][j+1][k])*0.5
              vp=(self.rectVars[self.getVarID("w")][i][j][k]
                +self.rectVars[self.getVarID("w")][i][j][k+1])*0.5
              varTmpX[i][j][k]=vr
              varTmpY[i][j][k]=vt
              varTmpZ[i][j][k]=vp
        self.rectVars.append(varTmpX)
        self.rectVars.append(varTmpY)
        self.rectVars.append(varTmpZ)
  def calZoneCenteredVelocity(self,curvature=True):
    
    if curvature:
      
      #register the variables
      self._addVarID("vx_cen","rect")
      self._addVarID("vy_cen","rect")
      self._addVarID("vz_cen","rect")
      if self.numDims==1:
        
        #add vx_cen
        varTmpX=np.empty(self.cellCenteredShape)
        for i in range(self.cellCenteredShape[0]):
          for j in range(self.cellCenteredShape[1]):
            for k in range(self.cellCenteredShape[2]):
              varTmpX[i][j][k]=(self.rectVars[self.getVarID("u")][i][j][k]
                +self.rectVars[self.getVarID("u")][i+1][j][k])*0.5
        self.rectVars.append(varTmpX)
        
        #add vy_cen
        self.rectVars.append(np.zeros(self.cellCenteredShape))
        
        #add vz_cen
        self.rectVars.append(np.zeros(self.cellCenteredShape))
      elif self.numDims==2:
        
        #add x_cen and y_cen
        varTmpX=np.empty(self.cellCenteredShape)
        varTmpY=np.empty(self.cellCenteredShape)
        for i in range(self.cellCenteredShape[0]):
          for j in range(self.cellCenteredShape[1]):
            for k in range(self.cellCenteredShape[2]):
              vr=(self.rectVars[self.getVarID("u")][i][j][k]
                +self.rectVars[self.getVarID("u")][i+1][j][k])*0.5
              vt=(self.rectVars[self.getVarID("v")][i][j][k]
                +self.rectVars[self.getVarID("v")][i][j+1][k])*0.5
              theta=self.rectVars[self.getVarID("theta_cen")][0][j][0]
              varTmpX[i][j][k]=vr*math.sin(theta)+vt*math.cos(theta)
              varTmpY[i][j][k]=vr*math.cos(theta)+vt*math.sin(theta)
        self.rectVars.append(varTmpX)
        self.rectVars.append(varTmpY)
        
        #add z_cen
        self.rectVars.append(np.zeros(self.cellCenteredShape))
      elif self.numDims==3:
        
        #add x_cen, y_cen and z_cen
        varTmpX=np.empty(self.cellCenteredShape)
        varTmpY=np.empty(self.cellCenteredShape)
        varTmpZ=np.empty(self.cellCenteredShape)
        sinTheta=np.empty((self.cellCenteredShape[1]))
        cosTheta=np.empty((self.cellCenteredShape[1]))
        sinPhi=np.empty((self.cellCenteredShape[2]))
        cosPhi=np.empty((self.cellCenteredShape[2]))
        for j in range((self.cellCenteredShape[1])):
          sinTheta[j]=math.sin(self.rectVars[self.getVarID("theta_cen")][0][j][0])
          cosTheta[j]=math.cos(self.rectVars[self.getVarID("theta_cen")][0][j][0])
        for k in range((self.cellCenteredShape[2])):
          sinPhi[k]=math.sin(self.rectVars[self.getVarID("phi_cen")][0][0][k])
          cosPhi[k]=math.cos(self.rectVars[self.getVarID("phi_cen")][0][0][k])
        for i in range(self.cellCenteredShape[0]):
          for j in range(self.cellCenteredShape[1]):
            for k in range(self.cellCenteredShape[2]):
              vr=(self.rectVars[self.getVarID("u")][i][j][k]
                +self.rectVars[self.getVarID("u")][i+1][j][k])*0.5
              vt=(self.rectVars[self.getVarID("v")][i][j][k]
                +self.rectVars[self.getVarID("v")][i][j+1][k])*0.5
              vp=(self.rectVars[self.getVarID("w")][i][j][k]
                +self.rectVars[self.getVarID("w")][i][j][k+1])*0.5
              varTmpX[i][j][k]=(
                vr*sinTheta[j]*cosPhi[k]
                +vt*cosTheta[j]*cosPhi[k]
                +vp*sinPhi[k])
              varTmpY[i][j][k]=(
                vr*sinTheta[j]*sinPhi[k]
                +vt*cosTheta[j]*sinPhi[k]
                +vp*cosTheta[j])
              varTmpZ[i][j][k]=(
                vr*cosTheta[j]
                -vt*sinTheta[j])
        self.rectVars.append(varTmpX)
        self.rectVars.append(varTmpY)
        self.rectVars.append(varTmpZ)
    else:
      
      #register the variables
      self._addVarID("vr_cen","rect")
      self._addVarID("vt_cen","rect")
      self._addVarID("vp_cen","rect")
      if self.numDims==1:
        
        #add x_cen
        varTmpX=np.empty(self.cellCenteredShape)
        for i in range(self.cellCenteredShape[0]):
          for j in range(self.cellCenteredShape[1]):
            for k in range(self.cellCenteredShape[2]):
              varTmpX[i][j][k]=(self.rectVars[self.getVarID("r")][i][j][k]
                +self.rectVars[self.getVarID("r")][i+1][j][k])*0.5
        self.rectVars.append(varTmpX)
        
        #add y_cen
        self.rectVars.append(np.zeros(self.cellCenteredShape))
        
        #add z_cen
        self.rectVars.append(np.zeros(self.cellCenteredShape))
        
      elif self.numDims==2:
        
        #add x_cen and y_cen
        varTmpX=np.empty(self.cellCenteredShape)
        varTmpY=np.empty(self.cellCenteredShape)
        for i in range(self.cellCenteredShape[0]):
          for j in range(self.cellCenteredShape[1]):
            for k in range(self.cellCenteredShape[2]):
              vr=(self.rectVars[self.getVarID("u")][i][j][k]
                +self.rectVars[self.getVarID("u")][i+1][j][k])*0.5
              vt=(self.rectVars[self.getVarID("v")][i][j][k]
                +self.rectVars[self.getVarID("v")][i][j+1][k])*0.5
              varTmpX[i][j][k]=vr
              varTmpY[i][j][k]=vt
        self.rectVars.append(varTmpX)
        self.rectVars.append(varTmpY)
        
        #add z_cen
        self.rectVars.append(np.zeros(self.cellCenteredShape))
      elif self.numDims==3:
        
        #add x_cen, y_cen and z_cen
        varTmpX=np.empty(self.cellCenteredShape)
        varTmpY=np.empty(self.cellCenteredShape)
        varTmpZ=np.empty(self.cellCenteredShape)
        for i in range(self.cellCenteredShape[0]):
          for j in range(self.cellCenteredShape[1]):
            for k in range(self.cellCenteredShape[2]):
              vr=(self.rectVars[self.getVarID("u")][i][j][k]
                +self.rectVars[self.getVarID("u")][i+1][j][k])*0.5
              vt=(self.rectVars[self.getVarID("v")][i][j][k]
                +self.rectVars[self.getVarID("v")][i][j+1][k])*0.5
              vp=(self.rectVars[self.getVarID("w")][i][j][k]
                +self.rectVars[self.getVarID("w")][i][j][k+1])*0.5
              varTmpX[i][j][k]=vr
              varTmpY[i][j][k]=vt
              varTmpZ[i][j][k]=vp
        self.rectVars.append(varTmpX)
        self.rectVars.append(varTmpY)
        self.rectVars.append(varTmpZ)
  def setAdditionalVars(self,varsToSet=None):
    """Set additional variables
    
    Sets variuous equation of state variables specified in the keyword 
    "varsToSet".
    
    varsToSet: a list of variable names to set. Variables which can be set are
      "dlnPdlnT","dlnPdlnD","dEdT","p","e","kappa","gamma","cp","c","gamma1"
      ,"DelAd","cv". By default all variables are set.
    """
    
    settableVars=["dlnPdlnT","dlnPdlnD","dEdT","p","e","kappa","gamma","cp","c"
      ,"gamma1","DelAd","cv"]
    
    #see if specific variables are selected
    if varsToSet!=None:
      pass
      #lets allow variables we can't set to be silently ignored, perhaps this 
      #should be a warning
      #check that we can set all the varsToSet
      #for var in varsToSet:
      #  if var not in settableVars:
      #    raise Exception("variable \""+var+"\" in keyword argument varsToSet "
      #      +"not in list of settable variables:"+str(settableVars))
    else:#otherwise set all variables
      varsToSet=settableVars
    
    #see if we need to load a new eos table
    global eosTable
    global eosFileName
    
    #check to see if it is a path relative to the install directory
    if self.eosString[0]!="/" and self.eosString[0:1]!="./":
      self.eosString=os.path.join(paths.BasePath,self.eosString)
      
    if eosFileName!=self.eosString:#if not the same table, load a new one
      eosTable=eos.Eos()
      eosTable.readBin(self.eosString)
      eosFileName=self.eosString
    
    temp=self.rectVars[self.getVarID("T")]
    rho=self.rectVars[self.getVarID("rho")]
    
    #variables to set
    if "dlnPdlnT" in varsToSet:
      dlnPdlnT=np.empty(self.cellCenteredShape)
    if "dlnPdlnD" in varsToSet:
      dlnPdlnD=np.empty(self.cellCenteredShape)
    if "dEdT" in varsToSet:
      dEdT=np.empty(self.cellCenteredShape)
    if "p" in varsToSet:
      p=np.empty(self.cellCenteredShape)
    if "e" in varsToSet:
      e=np.empty(self.cellCenteredShape)
    if "kappa" in varsToSet:
      kappa=np.empty(self.cellCenteredShape)
    if "gamma" in varsToSet:
      gamma=np.empty(self.cellCenteredShape)
    if "cp" in varsToSet:
      cp=np.empty(self.cellCenteredShape)
    if "c" in varsToSet:
      c=np.empty(self.cellCenteredShape)
    if "gamma1" in varsToSet:
      gamma1=np.empty(self.cellCenteredShape)
    if "DelAd" in varsToSet:
      DelAd=np.empty(self.cellCenteredShape)
    if "cv" in varsToSet:
      cv=np.empty(self.cellCenteredShape)
    
    for i in range(self.cellCenteredShape[0]):
      for j in range(self.cellCenteredShape[1]):
        for k in range(self.cellCenteredShape[2]):
          
          #get extra quantities from equation of state
          if ("dlnPdlnT" in varsToSet
            or "dlnPdlnD" in varsToSet
            or "dEdT" in varsToSet):
            [dlnPdlnTtmp,dlnPdlnDtmp,dEdTtmp]=eosTable.getDlnPDlnTDlnPDlnPDEDT(
              temp[i][j][k],rho[i][j][k])
          if ("p" in varsToSet
            or "e" in varsToSet
            or "kappa" in varsToSet
            or "gamma" in varsToSet
            or "cp" in varsToSet):
            [ptmp,etmp,kappatmp,gammatmp,cptmp]=eosTable.getPEKappaGammaCp(
              temp[i][j][k],rho[i][j][k])
          if ("gamma1" in varsToSet
            or "DelAd" in varsToSet
            or "cv" in varsToSet):
            [gamma1tmp,DelAdtmp,cvtmp]=eosTable.getGamma1DelAdC_v(
              temp[i][j][k],rho[i][j][k])
          if "c" in varsToSet:
            ctmp=eosTable.getSoundSpeed(
              temp[i][j][k],rho[i][j][k])
          if "dRhodP" in varsToSet:
            dRhodPtmp=eosTable.getDRhoDP(
              temp[i][j][k],rho[i][j][k])
          
          #set values
          if "dlnPdlnT" in varsToSet:
            dlnPdlnT[i][j][k]=dlnPdlnTtmp
          if "dlnPdlnD" in varsToSet:
            dlnPdlnD[i][j][k]=dlnPdlnDtmp
          if "dEdT" in varsToSet:
            dEdT[i][j][k]=dEdTtmp
          if "p" in varsToSet:
            p[i][j][k]=ptmp
          if "e" in varsToSet:
            e[i][j][k]=etmp
          if "kappa" in varsToSet:
            kappa[i][j][k]=kappatmp
          if "gamma" in varsToSet:
            gamma[i][j][k]=gammatmp
          if "cp" in varsToSet:
            cp[i][j][k]=cptmp
          if "c" in varsToSet:
            c[i][j][k]=ctmp
          if "gamma1" in varsToSet:
            gamma1[i][j][k]=gamma1tmp
          if "DelAd" in varsToSet:
            DelAd[i][j][k]=DelAdtmp
          if "cv" in varsToSet:
            cv[i][j][k]=cvtmp
    
    #Add newly set variables
    if "dlnPdlnT" in varsToSet:
      self._addVarID("dlnPdlnT","rect")
      self.rectVars.append(dlnPdlnT)
    if "dlnPdlnD" in varsToSet:
      self._addVarID("dlnPdlnD","rect")
      self.rectVars.append(dlnPdlnD)
    if "dEdT" in varsToSet:
      self._addVarID("dEdT","rect")
      self.rectVars.append(dEdT)
    if "p" in varsToSet:
      self._addVarID("p","rect")
      self.rectVars.append(p)
    if "e" in varsToSet:
      self._addVarID("e","rect")
      self.rectVars.append(e)
    if "kappa" in varsToSet:
      self._addVarID("kappa","rect")
      self.rectVars.append(kappa)
    if "gamma" in varsToSet:
      self._addVarID("gamma","rect")
      self.rectVars.append(gamma)
    if "cp" in varsToSet:
      self._addVarID("cp","rect")
      self.rectVars.append(cp)
    if "c" in varsToSet:
      self._addVarID("c","rect")
      self.rectVars.append(c)
    if "gamma1" in varsToSet:
      self._addVarID("gamma1","rect")
      self.rectVars.append(gamma1)
    if "DelAd" in varsToSet:
      self._addVarID("DelAd","rect")
      self.rectVars.append(DelAd)
    if "cv" in varsToSet:
      self._addVarID("cv","rect")
      self.rectVars.append(cv)
def main():
  """Performs some basic tests of the class Dump.
  
  This script should be run from the command line only as a means of testing
  Dump class methods.
  """
  
  import optparse as op
  
  #make command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] DUMPFILENAME"
    ,version="%prog 1.0"
    ,description="Performs some basic tests on the Dump class. Provide a test"
      +"dump file with DUMPFILENAME")
  
  #parse command line options
  (options,args)=parser.parse_args()
  
  
  #test out some simple uses for this class
  #print "reading a dump file ..."
  #dumpFile1=Dump(args[0])
  
  #print "printing the header ..."
  #dumpFile1.printHeader(sys.stdout)
  
  #print "printing variable 1 ..."
  #dumpFile1.printVar(1,sys.stdout)
  
  #print "printing variable 2 ..."
  #dumpFile1.printVar(2,sys.stdout)
  
  #print dumpFile1.getVarNames()
  #dumpFile1.printVarToSTDOut("v")
  
  #a=dumpFile1.getVarSlice("v",rIndexMin=100)
  #TFile=open("T_0.txt",'w')
  #dumpFile1.printVarSliceToOutInRadCol("T",TFile,tIndexMin=2,tIndexMax=20)
  #TFile.close()
  #TFile=open("T_ave_0.txt",'w')
  #dumpFile1.printHorAveVarSliceToOutInRadCol("T",TFile,tIndexMin=2,tIndexMax=20)
  #TFile.close()
  #dumpFile1.printVarSliceToOutInRadCol("v",sys.stdout,tIndexMin=2,tIndexMax=4)
if __name__ == "__main__":
  main()