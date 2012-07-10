import parser
from math import *
def getFormula(text):
  
  formulaOrig=None
  formula=None
  nColumn=None
  nRowShift=None
  code=None
  
  if text.isdigit():
    nColumn=int(text)-1
  else:
    #if there is a foluma get list of columns
    formulaOrig=text
    columnsTemp=formulaOrig.split('$')
    formula=columnsTemp[0]
    if len(columnsTemp)>1:#turn into a list
      nColumn=[]
      nRowShift=[]
      for i in range(1,len(columnsTemp)):
        columnAndFormula=splitFirstInt(columnsTemp[i])
        if columnAndFormula[1][0:4] =="_ip1":
          columnAndFormula[1]=columnAndFormula[1][4:]
          nRowShift.append(1)
        elif columnAndFormula[1][0:4] =="_im1":
          columnAndFormula[1]=columnAndFormula[1][4:]
          self.nRowShift.append(-1)
        elif columnAndFormula[1][0:2] =="_i":
          columnAndFormula[1]=columnAndFormula[1][2:]
          nRowShift.append(0)
        else:
          nRowShift.append(0)
        formula+="a["+str(i-1)+"]"+columnAndFormula[1]
        nColumn.append(columnAndFormula[0]-1)
      code=parser.expr(formula).compile()
    else:
      formula=formulaOrig
      code=parser.expr(formula).compile()
      if text == None or text=="":
        raise Exception("No column number given for curve, must have curve text be either an "
          +"integer, or a mathatmical formula containing column references prefixed with a \"$\".")
  return [formulaOrig,formula,nColumn,nRowShift,code]
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
def getY(nRowShift,nColumn,fileData,code,i):
  if nColumn!=None:
    if isinstance(nColumn,int):
      #normal simple 1 column data
      return fileData.fColumnValues[i][nColumn]
    else:
      
      #add elements to list a for columns used in the "code"
      a=[]
      for j in range(len(nColumn)):
        if i+nRowShift[j]<0:#referencing data too far in
          return None
        elif i+nRowShift[j]>=len(fileData.fColumnValues):#referencing data too far out
          return None
        elif fileData.fColumnValues[i+nRowShift[j]][nColumn[j]]==None:#referencing empty data
          return None
        else:
          a.append(fileData.fColumnValues[i+nRowShift[j]][nColumn[j]])
  return eval(code)
