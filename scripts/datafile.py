import os
import numpy as np

class DataFile:
  '''A generic class for holding a file consisting of a header and columns of floats'''
  sColumnNames=None
  fColumnValues=None
  sHeader=None
  
  def setFileSize(self,nNumRows,nNumColumns):
    self.fColumnValues=np.zeros( (nNumRows,nNumColumns) )
  def readFile(self,sFileName):
    ''' a wrapper to determine which readFile function should be used'''
    if type(self.fColumnValues)==type(None):
      self.readFileUnFixed(sFileName)
    else:
      self.readFileFixed(sFileName)
  def readFileFixed(self,sFileName):
    '''Reads in a file when the size has already been set using \ref setFileSize, or by a previous
    file read using \ref readFileUnFixed.'''
    
    #open file
    if not os.access(sFileName,os.F_OK|os.R_OK):
      print "error opening \"",sFileName,"\" for reading"
      return
    f=open(sFileName,'r')
    
    #read in first line
    line1=f.readline()
    
    #read in second line
    line2=f.readline()
    
    #split them at spaces
    line1Split=line1.split()
    line2Split=line2.split()
    
    #if both have the same number of columns, then the first line is assumed to be column names
    nLine=0
    if len(line1Split)==len(line2Split):
      self.sHeader=None
      self.sColumnNames=line1Split
      
      #convert columns to floats
      fLineValues=[]
      
      for i in range(len(self.sColumnNames)):
        if line2Split[i]=='-' :
          self.fColumnValues[nLine][i]=None
        else:
          self.fColumnValues[nLine][i]=float(line2Split[i])
      nLine=nLine+1
    else:# if not the same number of columns, then the first line is assumed to be a header
      self.sHeader=line1
      self.sColumnNames=line2Split
    
    #read in columns
    for line in f:
      
      #split up line into columns
      rowTemp=line.split()
      
      #convert columns to floats
      for i in range(len(self.sColumnNames)-1):
        if rowTemp[i]=='-' :
          self.fColumnValues[nLine][i]=None
        else:
          self.fColumnValues[nLine][i]=float(rowTemp[i])
      nLine=nLine+1
  def readFileUnFixed(self,sFileName):
    '''Reads in a file when the size is not fixed and needs to be determined from the input file 
    being read in'''
    
    self.sColumnNames=[]
    self.fColumnValues=[]
    self.sHeader=""
    
    #open file
    if not os.access(sFileName,os.F_OK|os.R_OK):
      print "error opening \"",sFileName,"\" for reading"
      return
    f=open(sFileName,'r')
    
    #read in first line
    line1=f.readline()
    
    #read in second line
    line2=f.readline()
    
    #split them at spaces
    line1Split=line1.split()
    line2Split=line2.split()
    
    #if both have the same number of columns, then the first line is assumed to be column names
    if len(line1Split)==len(line2Split):
      self.sHeader=None
      self.sColumnNames=line1Split
      
      #convert columns to floats
      fLineValues=[]
      
      for i in range(len(self.sColumnNames)):
        if line2Split[i]=='-' :
          fLineValues.append(None)
        else:
          fLineValues.append(float(line2Split[i]))
      
      #add line of values to list
      self.fColumnValues.append(fLineValues)
      
    else:# if not the same number of columns, then the first line is assumed to be a header
      self.sHeader=line1
      self.sColumnNames=line2Split
    
    #read in columns
    nLine=0
    for line in f:
      
      #split up line into columns
      rowTemp=line.split()
      
      #convert columns to floats
      fLineValues=[]
      for i in range(len(self.sColumnNames)):
        if rowTemp[i]=='-' :
          fLineValues.append(None)
        else:
          fLineValues.append(float(rowTemp[i]))
      
      #add line of values to list
      self.fColumnValues.append(fLineValues)
      
    #convert to a numpy array
    self.fColumnValues=np.array(self.fColumnValues)
      