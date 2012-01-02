#!/usr/bin/env python
import struct
class dump:
  def __init__(self,fileName):
    '''Initilizes the dump by reading in a binary file.'''
    
    self.read(fileName)
  def read(self,fileName):
    '''Reads in a binary dump file.'''
    self.fileName=fileName
    self.f=open(fileName,'rb')
    self.readHeader()
    if self.type=='b':
      for i in range(self.numVars):
        self.readBinaryVar(i)
    elif self.type=='a':
      for i in range(self.numVars):
        self.readAsciiVar(i)
    self.setVarIDs()
  def setVarIDs(self):
    '''Sets names for the interger values of the grid varibles'''
    self.varIDs=[]
    if self.gamma!=None: #using gamma law gas
      if self.numDims==1:
        self.varIDs.append("M_r")
        self.varIDs.append("Delta M_r")
        self.varIDs.append("r")
        self.varIDs.append("rho")
        self.varIDs.append("u")
        self.varIDs.append("u_0")
        self.varIDs.append("e")
      elif self.numDims==2:
        self.varIDs.append("M_r")
        self.varIDs.append("theta")
        self.varIDs.append("Delta M_r")
        self.varIDs.append("r")
        self.varIDs.append("rho")
        self.varIDs.append("u")
        self.varIDs.append("u_0")
        self.varIDs.append("v")
        self.varIDs.append("e")
      elif self.numDims==3:
        self.varIDs.append("M_r")
        self.varIDs.append("theta")
        self.varIDs.append("phi")
        self.varIDs.append("Delta M_r")
        self.varIDs.append("r")
        self.varIDs.append("rho")
        self.varIDs.append("u")
        self.varIDs.append("u_0")
        self.varIDs.append("v")
        self.varIDs.append("w")
        self.varIDs.append("e")
    else: #using a tabulated equaiton of state
      if self.numDims==1:
        self.varIDs.append("M_r")
        self.varIDs.append("Delta M_r")
        self.varIDs.append("r")
        self.varIDs.append("rho")
        self.varIDs.append("u")
        self.varIDs.append("u_0")
        self.varIDs.append("T")
      elif self.numDims==2:
        self.varIDs.append("M_r")
        self.varIDs.append("theta")
        self.varIDs.append("Delta M_r")
        self.varIDs.append("r")
        self.varIDs.append("rho")
        self.varIDs.append("u")
        self.varIDs.append("u_0")
        self.varIDs.append("v")
        self.varIDs.append("T")
      elif self.numDims==3:
        self.varIDs.append("M_r")
        self.varIDs.append("theta")
        self.varIDs.append("phi")
        self.varIDs.append("Delta M_r")
        self.varIDs.append("r")
        self.varIDs.append("rho")
        self.varIDs.append("u")
        self.varIDs.append("u_0")
        self.varIDs.append("v")
        self.varIDs.append("w")
        self.varIDs.append("T")
  def readHeader(self):
    '''Reads header information from binary dump file.'''
    
    #read header
    self.type           =struct.unpack('c',self.f.read(1))[0]#file type, either a or b
    if self.type=='a':
      self.readHeaderAscii()
    else:
      self.readHeaderBinary()
  def readHeaderBinary(self):
    '''Reads a header from a binary file, after the type has been read in.'''
    self.version        =struct.unpack('i',self.f.read(4))[0]#file version integer
    self.time           =struct.unpack('d',self.f.read(8))[0]#time
    self.timeStepIndex  =struct.unpack('i',self.f.read(4))[0]#time step index
    self.delta_t_nm1half=struct.unpack('d',self.f.read(8))[0]#delta t nm1half
    self.delta_t_np1half=struct.unpack('d',self.f.read(8))[0]#delta t np1half
    self.alpha          =struct.unpack('d',self.f.read(8))[0]#alpha
    self.eosStringLen   =struct.unpack('i',self.f.read(4))[0]#size of equation of state string
    if self.eosStringLen>0:
      self.eosString=""
      self.gamma=None
      for i in range(self.eosStringLen):
        self.eosString+=struct.unpack('s',self.f.read(1))[0]#get a character of the eos string
    else:
      self.eosString=None
      self.gamma=struct.unpack('d',self.f.read(8))[0]#get value of gamma
    self.av             =struct.unpack('d',self.f.read(8))[0]#artificial viscosity
    self.avthreshold    =struct.unpack('d',self.f.read(8))[0]#artificial viscosity threshold
    self.globalDims=[]
    self.globalDims.append(struct.unpack('i',self.f.read(4))[0])#x0 dim
    self.globalDims.append(struct.unpack('i',self.f.read(4))[0])#x1 dim
    self.globalDims.append(struct.unpack('i',self.f.read(4))[0])#x2 dim
    self.boundaryConditions=[]
    self.boundaryConditions.append(struct.unpack('i',self.f.read(4))[0])#x0 BC
    self.boundaryConditions.append(struct.unpack('i',self.f.read(4))[0])#x1 BC
    self.boundaryConditions.append(struct.unpack('i',self.f.read(4))[0])#x2 BC
    self.num1DZones=struct.unpack('i',self.f.read(4))[0]# number of 1D zones
    self.numGhostCells=struct.unpack('i',self.f.read(4))[0]# number of ghost cells at boundaries
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
  def readHeaderAscii(self):
    '''Reads a header from a ascii file, after the type has been read in.'''
    
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
  def readBinaryVar(self,var):
    '''Read in a variable from a binary dump file. Must be called with var increasing from 0 to 
    self.numVars.'''
    
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
  def readAsciiVar(self,var):
    '''Read in a variable from an ascii dump file. Must be called with var increasing from 0 to 
    self.numVars.'''
    
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
  def printHeader(self):
    '''Prints the header of a binary dump file to the standard output.'''
    
    #print header
    print "type=",self.type
    print "version=",self.version
    print "time=",self.time
    print "timeStepIndex=",self.timeStepIndex
    print "delta_t_nm1half=",self.delta_t_nm1half
    print "delta_t_np1half=",self.delta_t_np1half
    print "alpha=",self.alpha
    print "eosStringLen=",self.eosStringLen
    if self.eosStringLen>0:
      print "eosString=",self.eosString
    else:
      print "gamma=",self.gamma
    print "av=",self.av
    print "avthreshold=",self.avthreshold
    print "globalDims=",self.globalDims
    print "boundaryConditions=",self.boundaryConditions
    print "num1DZones=",self.num1DZones
    print "numGhostCells=",self.numGhostCells
    print "numVars=",self.numVars
    print "varInfo=",self.varInfo
    print "varSize=",self.varSize
    print "numDims=",self.numDims
  def printVar(self,var):
    '''Prints a variable to the standard output.'''
    
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
          print "var["+str(i)+"]["+str(j)+"]["+str(k)+"]=",self.vars[var][i][j][k]
    
    #print out MD part
    sizeX02=self.varSize[var][0]+ghostCellsInX0*2*self.numGhostCells
    sizeX1=self.varSize[var][1]+ghostCellsInX1*2*self.numGhostCells
    sizeX2=self.varSize[var][2]+ghostCellsInX2*2*self.numGhostCells
    for i in range(sizeX01,sizeX02):
      for j in range(sizeX1):
        for k in range(sizeX2):
          print "var["+str(i)+"]["+str(j)+"]["+str(k)+"]=",self.vars[var][i][j][k]
  def printDump(self):
    self.printHeader()
    for i in range(self.numVars):
      self.printVar(i)
def main():
  
  #test out some simple uses for this class
  dumpFile1=dump("3DNARestartTest3_t00000002")
  dumpFile1.printHeader()
  dumpFile1.printVar(1)
  dumpFile1.printVar(2)
if __name__ == "__main__":
  main()