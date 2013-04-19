#!/usr/bin/env python
import struct
import sys
import numpy as np

class Dump:
  """Allows manipulation of SPHERLS binary and ascii dump files.
  
  \todo should probalby add methods to write out binary/ascii dump files in the
  version that SPHERLS knows how to read.
  """
  
  def __init__(self,fileName):
    """Initilizes the dump by reading in a binary file."""
    
    self.read(fileName)
  def _setVarIDs(self):
    """Sets names for the interger values of the grid varibles"""
    
    self._varIDs={}
    self.varNames=[]
    if self.gamma!=None: #using gamma law gas
      if self.numDims==1:
        self._varIDs["M_r"]=0
        self._varIDs["Delta_M_r"]=1
        self._varIDs["r"]=2
        self._varIDs["rho"]=3
        self._varIDs["u"]=4
        self._varIDs["u_0"]=5
        self._varIDs["e"]=6
        self._varIDs["theta"]=None
        self._varIDs["phi"]=None
        self._varIDs["T"]=None
        self._varIDs["v"]=None
        self._varIDs["w"]=None
        
        self.varNames.append("M_r")
        self.varNames.append("Delta_M_r")
        self.varNames.append("r")
        self.varNames.append("rho")
        self.varNames.append("u")
        self.varNames.append("u_0")
        self.varNames.append("e")
        self.varNames.append(None)
        self.varNames.append(None)
        self.varNames.append(None)
        self.varNames.append(None)
        self.varNames.append(None)
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
        self._varIDs["phi"]=None
        self._varIDs["T"]=None
        self._varIDs["w"]=None
        
        self.varNames.append("M_r")
        self.varNames.append("theta")
        self.varNames.append("Delta_M_r")
        self.varNames.append("r")
        self.varNames.append("rho")
        self.varNames.append("u")
        self.varNames.append("u_0")
        self.varNames.append("v")
        self.varNames.append("e")
        self.varNames.append(None)
        self.varNames.append(None)
        self.varNames.append(None)
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
        self._varIDs["T"]=None
        
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
        self.varNames.append(None)
    else: #using a tabulated equaiton of state
      if self.numDims==1:
        self._varIDs["M_r"]=0
        self._varIDs["Delta_M_r"]=1
        self._varIDs["r"]=2
        self._varIDs["rho"]=3
        self._varIDs["u"]=4
        self._varIDs["u_0"]=5
        self._varIDs["T"]=6
        self._varIDs["theta"]=None
        self._varIDs["phi"]=None
        self._varIDs["e"]=None
        self._varIDs["v"]=None
        self._varIDs["w"]=None
        
        self.varNames.append("M_r")
        self.varNames.append("Delta_M_r")
        self.varNames.append("r")
        self.varNames.append("rho")
        self.varNames.append("u")
        self.varNames.append("u_0")
        self.varNames.append("T")
        self.varNames.append(None)
        self.varNames.append(None)
        self.varNames.append(None)
        self.varNames.append(None)
        self.varNames.append(None)
      elif self.numDims==2:
        self._varIDs["M_r"]=0
        self._varIDs["theta"]=1
        self._varIDs["Delta_M_r"]=2
        self._varIDs["r"]=3
        self._varIDs["rho"]=4
        self._varIDs["u"]=5
        self._varIDs["u_0"]=6
        self._varIDs["v"]=7
        self._varIDs["T"]=8
        self._varIDs["phi"]=None
        self._varIDs["e"]=None
        self._varIDs["w"]=None
        
        self.varNames.append("M_r")
        self.varNames.append("theta")
        self.varNames.append("Delta_M_r")
        self.varNames.append("r")
        self.varNames.append("rho")
        self.varNames.append("u")
        self.varNames.append("u_0")
        self.varNames.append("v")
        self.varNames.append("T")
        self.varNames.append(None)
        self.varNames.append(None)
        self.varNames.append(None)
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
        self._varIDs["e"]=None
        
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
        self.varNames.append(None)
  def _readHeaderBinary(self):
    """Reads a header from a binary file, after the type has been read in."""
    
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
  def _readHeaderAscii(self):
    """Reads a header from a ascii file, after the type has been read in."""
    
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
    """Writes a variable to out."""
    
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
  def readHeader(self):
    """Reads header information from binary dump file.
    
    This version calls either the _readHeaderAscii or the _readHeaderBinary
    
    """
    
    #read header
    self.type=struct.unpack('c',self.f.read(1))[0]#file type, either a or b
    if self.type=='a':
      self._readHeaderAscii()
    else:
      self._readHeaderBinary()
  def read(self,fileName):
    """Reads in a binary dump file."""
    
    self.fileName=fileName
    self.f=open(fileName,'rb')
    self.readHeader()
    if self.type=='b':
      for i in range(self.numVars):
        self._readBinaryVar(i)
    elif self.type=='a':
      for i in range(self.numVars):
        self._readAsciiVar(i)
    self._setVarIDs()
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
  def getVarNames(self):
    """Returns a list of variable names that are availble."""
    
    #make a list of all variable keys that aren't None
    varNames=[]
    for key in self._varIDs.keys():
      if self._varIDs[key]!=None:
        varNames.append(key)
    return varNames
  def getVarID(self,var):
    """Returns the array index (ID) of a variable given by name"""
    
    return self._varIDs[var]
  def printVarToOut(self,var,out):
    """Print var to out"""
    
    varNames=self.getVarNames()
    if var not in varNames:
      raise Exception(str(var)+" not in avaible variable names "+str(varNames))
    
    self._printVarByID(self._varIDs[var],out)
  def printVarToSTDOut(self,var):
    """Prints varible to standard output"""
    
    self.printVarToOut(var,sys.stdout)
  def getRecVar(self,var):
    """Returns a rectangular numpy array version of a varible.
    
    Variables are stored as a 1D part plus a 2D or 3D part. This function
    returns a variable that has the 1D part copied to match the 2D or 3D part.
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
  print "reading a dump file ..."
  dumpFile1=Dump(args[0])
  
  #print "printing the header ..."
  #dumpFile1.printHeader(sys.stdout)
  
  #print "printing variable 1 ..."
  #dumpFile1.printVar(1,sys.stdout)
  
  #print "printing variable 2 ..."
  #dumpFile1.printVar(2,sys.stdout)
  
  print dumpFile1.getVarNames()
  #dumpFile1.printVarToSTDOut("v")
  
  #a=dumpFile1.getVarSlice("v",rIndexMin=100)
  TFile=open("T_0.txt",'w')
  dumpFile1.printVarSliceToOutInRadCol("T",TFile,tIndexMin=2,tIndexMax=20)
  TFile.close()
  TFile=open("T_ave_0.txt",'w')
  dumpFile1.printHorAveVarSliceToOutInRadCol("T",TFile,tIndexMin=2,tIndexMax=20)
  TFile.close()
  #dumpFile1.printVarSliceToOutInRadCol("v",sys.stdout,tIndexMin=2,tIndexMax=4)
if __name__ == "__main__":
  main()