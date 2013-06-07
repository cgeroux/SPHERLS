#!/usr/bin/env python
import os
import sys

BasePath=os.path.dirname(os.path.dirname(sys.argv[0]))#remove script name, remove scripts directory
BasePath=os.path.join(os.getcwd(),BasePath)#fix for relative paths
SPHERLSanalPath=BasePath+"/SPHERLSanal"
SPHERLSgenPath=BasePath+"/SPHERLSgen"
SPHERLSPath=BasePath+"/SPHERLS"
scriptPath=BasePath+"/scripts"
SPHERLSDataPath=BasePath+"/data"
EOSPath=SPHERLSDataPath+"/eos"
velocityProfilePath=SPHERLSDataPath+"/velocity_pro"
ref_calcsPath=SPHERLSDataPath+"/ref_calcs"
srcPath=BasePath+"/src/"
debug=False

def check_paths():
  global BasePath
  global SPHERLSPath
  global SPHERLSgenPath
  global SPHERLSanalPath
  global scriptPath
  global SPHERLSDataPath
  global EOSPath
  global velocityProfilePath
  global ref_calcsPath
  global srcPath
  
  #check that BasePath exists
  if not os.access(BasePath,os.F_OK):
    raise Exception("\""+BasePath+"\" not found!")
  
  #check for SPHERLS in SPHERLSPath
  if not os.access(SPHERLSPath,os.F_OK):#source location
    tempPath=BasePath+"/bin/SPHERLS"
    if not os.access(tempPath,os.F_OK):#installed location
      raise Exception("SPHERLS not found in\""+SPHERLSPath+"\" or \""+tempPath
        +"\"!")
    else:
      SPHERLSPath=tempPath
  
  #check for SPHERLSgen in SPHERLSgenPath
  if not os.access(SPHERLSgenPath,os.F_OK):#source location
    tempPath=BasePath+"/bin/SPHERLSgen"
    if not os.access(tempPath,os.F_OK):#installed location
      raise Exception("SPHERLSgen not found in\""+SPHERLSgenPath+"\" or \""
        +tempPath+"\"!")
    else:
      SPHERLSgenPath=tempPath
  
  #check for SPHERLSanal in SPHERLSanalPath
  if not os.access(SPHERLSanalPath,os.F_OK):#source location
    tempPath=BasePath+"/bin/SPHERLSanal"
    if not os.access(tempPath,os.F_OK):#installed location
      raise Exception("SPHERLSanal not found in \""+SPHERLSanalPath+"\" or \""
        +tempPath+"\"!")
    else:
      SPHERLSanalPath=tempPath
  
  #check that script path exists
  if not os.access(scriptPath,os.F_OK):#source location
    tempPath=BasePath+"/bin"
    if not os.access(tempPath,os.F_OK):#installed location
      raise Exception("script path \""+scriptPath+"\" and \""
        +tempPath+"\" not found!")
    else:
      scriptPath=tempPath
  
  #check that data path exists
  if not os.access(SPHERLSDataPath,os.F_OK):#source location
    raise Exception("data path \""+SPHERLSDataPath+"\" not found!")
  
  #check that equation of state path exists
  if not os.access(EOSPath,os.F_OK):#source location
    raise Exception("eos path \""+EOSPath+"\" not found!")
  
  #check that velocityProfile path exists
  if not os.access(velocityProfilePath,os.F_OK):#source location
    raise Exception("velocity profile path \""+velocityProfilePath
      +"\" not found!")
  
  #check that ref_calcs path exists
  if not os.access(ref_calcsPath,os.F_OK):#source location
    raise Exception("reference calculations path \""+ref_calcsPath
      +"\" not found!")
  
  #check that src path exists, will not exsist in installation
  if not os.access(srcPath,os.F_OK):#source location
    #print ("source  path \""+srcPath+"\" not found!")
    srcPath=""
  
  if debug:
    print "BasePath=\""+BasePath+"\""
    print "SPHERLSPath=\""+SPHERLSPath+"\""
    print "SPHERLSgenPath=\""+SPHERLSgenPath+"\""
    print "SPHERLSanalPath=\""+SPHERLSanalPath+"\""
    print "scriptPath=\""+scriptPath+"\""
    print "SPHERLSDataPath=\""+SPHERLSDataPath+"\""
    print "EOSPath=\""+EOSPath+"\""
    print "velocityProfilePath=\""+velocityProfilePath+"\""
    print "ref_calcsPath=\""+ref_calcsPath+"\""
    print "srcPath=\""+srcPath+"\""
  
  return True

#check paths whenever this script is imported
check_paths()