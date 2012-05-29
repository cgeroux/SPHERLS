#!/usr/bin/env python
import os
import sys

BasePath=os.path.dirname(os.path.dirname(sys.argv[0]))#remove script name, remove scripts directory
SPHERLSanalPath=BasePath+"/SPHERLSanal"
SPHERLSgenPath=BasePath+"/SPHERLSgen"
SPHERLSPath=BasePath+"/SPHERLS"
scriptPaths=BasePath+"/scripts/"
SPHERLSDATA=BasePath+"/data/"
EOSPath=BasePath+"/data/eos/eosY240Z002"
velocityProfilePath=BasePath+"/data/velocity_pro/"
ref_calcs=BasePath+"/data/ref_calcs/"

def check_paths():
  #check paths
  if not os.access(SPHERLSPath,os.F_OK):
    print "\""+SPHERLSPath+"\" not found, stopping!"
    return False
  if not os.access(SPHERLSgenPath,os.F_OK):
    print "\""+SPHERLSgenPath+"\" not found, stopping!"
    return False
  if not os.access(SPHERLSanalPath,os.F_OK):
    print "\""+SPHERLSanalPath+"\" not found, stopping!"
    return False
  if not os.access(EOSPath,os.F_OK):
    print "\""+EOSPath+"\" not found, stopping!"
    return False
  if not os.access(velocityProfilePath,os.F_OK):
    print "\""+velocityProfilePath+"\" not found, stopping!"
    return False
  if not os.access(ref_calcs,os.F_OK):
    print "\""+ref_calcs+"\" not found, stopping!"
    return False
  return True
