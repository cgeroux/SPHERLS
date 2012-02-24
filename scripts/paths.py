#!/usr/bin/env python
import os

SPHERLSanalPath="/home/cgeroux/SPHERLS/SPHERLSanal"
SPHERLSgenPath="/home/cgeroux/SPHERLS/SPHERLSgen"
#SPHERLSPath="/home/cgeroux/SPHERLS/SPHERLS"
SPHERLSPath="/home/cgeroux/SPHERLS_new_energy_ev/SPHERLS"
scriptPaths="/home/cgeroux/SPHERLS/scripts/"
EOSPath="/home/cgeroux/SPHERLS/data/eos/eosY240Z002"
velocityProfilePath="/home/cgeroux/SPHERLS/data/velocity_pro/T5700_L50_M575_fu.dat"
ref_calcs="/home/cgeroux/SPHERLS/data/ref_calcs/"

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
