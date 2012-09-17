#!/usr/bin/env python
import xml.etree.ElementTree as xml
import mywarnings
import warnings
import optparse as op

def parseOptions():
  
  #setup command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] XMLFILE"
    ,version="%prog 1.0",description=r"Reads in the XMLFILE which specifies the stars to phase.")
    
  #parse command line options
  return parser.parse_args()
def parseXML(xmlConfig):
  starList={}
  tree=xml.parse(xmlConfig)
  root=tree.getroot()
  
  #get inputFile
  inputFileElements=root.findall("inputFile")
  if len(inputFileElements)>1:
    warn.warning("more than one \"inputFile\" node ignoring all but first node")
  
  inputFile=inputFileElements[0].text
  
  #get stars node
  starsElements=root.findall("stars")
  if len(starsElements)>1:
    warn.warning("more than one \"stars\" node ignoring all but first node")
  
  #get all stars
  starElements=starsElements[0].findall("star")
  for starElement in starElements:
    name=starElement.text.strip()
    starList[name]={}
    starList[name]['period']=float(starElement.get("period"))
    starList[name]['filter']=starElement.get("filter")
  return [starList,inputFile]
def main():
  
  (options,args)=parseOptions()
  
  if len(args)>1 or len(args)==0:
    raise Exception("must have an xml file")
  
  [starList,inputFile]=parseXML(args[0])
  
  #read data
  stars={}
  f=open(inputFile)
  
  #dump first set of line
  for i in range(13):
    line=f.readline()
  
  #parse first line
  line=f.readline()
  nameTmp=line.split()[0]
  timeTmp=float(line.split()[1])
  filterTmp=line.split()[2]
  if len(line.split())>3:
    mTmp=float(line.split()[3])
  else:
    mTmp=None
  nameTmpOld=nameTmp
  filterTmpOld=filterTmp
  
  #parse the rest of the lines
  time=[]
  m=[]
  for line in f:
    
    #parse line
    nameTmp=line.split()[0]
    timeTmp=float(line.split()[1])
    filterTmp=line.split()[2]
    if len(line.split())>3:
      mTmp=float(line.split()[3])
    else:
      mTmp=None
    
    #save data
    if nameTmp!=nameTmpOld or filterTmp!=filterTmpOld:
      
      #save data collected for last star
      if not(nameTmpOld in stars.keys()):#don't have a dictionary yet
        stars[nameTmpOld]={}#add dictionary for that star
      stars[nameTmpOld][filterTmpOld]=[time,m]
      
      #clean slate for next star
      time=[]
      m=[]
      
      time.append(timeTmp)
      m.append(mTmp)
    else:
      
      #add current time/magintude
      time.append(timeTmp)
      m.append(mTmp)
      
    #save for later
    nameTmpOld=nameTmp
    filterTmpOld=filterTmp
  f.close()
  
  #phase data
  for starName in starList:
    
    #phase data
    t0=stars[starName][starList[starName]['filter']][0][0]
    period=starList[starName]['period']
    phasedData=[]
    for i in range(len(stars[starName][starList[starName]['filter']][0])):
      t=stars[starName][starList[starName]['filter']][0][i]-t0
      phase=(t-int(t/period)*period)/period
      mag=stars[starName][starList[starName]['filter']][1][i]
      phasedData.append([phase,mag])
    
    minData=[]
    minDataTmp=min(phasedData,key=lambda x: x[1])
    minData.append(minDataTmp[0])
    minData.append(minDataTmp[1])
    
    for i in range(len(phasedData)):
      phasedData[i][0]=phasedData[i][0]-minData[0]
      if phasedData[i][0]<0:
        phasedData[i][0]=phasedData[i][0]+1.0
    
    phasedData.sort(key=lambda x: x[0])
    f=open(starName+"_phased.txt",'w')
    f.write("phase mag\n")
    print "saving phased light curve to \""+starName+"_phased.txt\""
    for data in phasedData:
      f.write(str(data[0])+" "+str(data[1])+"\n")
      
    #print out twice
    for data in phasedData:
      f.write(str(data[0]+1.0)+" "+str(data[1])+"\n")
    f.close()
if __name__ == "__main__":
  main()