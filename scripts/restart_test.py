#!/usr/bin/env python
'''
  Test that restarts produce the same results as calculating straight through the restart with no
 restart.
'''
useArgparse=True
try:
  import argparse
except ImportError:
  useArgparse=False
  import optparse as op
import subprocess
import os
import shutil
import diffDumps

def main():
  
  #if argparse is availble
  if useArgparse:
    #set parser options
    parser=argparse.ArgumentParser(description="Tests that model restarts produce the same output"\
      +" as if the code had been run straight through.")
    parser.add_argument('-f',action="store_true",default=False,help="Force removal of"\
      +" pre-existing temporary data automatically.")
    parser.add_argument('-k',action="store_true",default=False,help="Keep temporary directory.")
    parser.add_argument('-p',default=5e-14,type=float,help="Sets the amount of relative"\
      +" absolute difference allowed between numerical values of the two dump files."\
      +" [default: 5e-14]")
    parser.add_argument('-t',default=5e-14,type=float,help="Sets the smallest number to"\
      +" consider for numerical comparison, differences in numbers which are smaller than this "\
      +"threshold will be ignored for the comparison. [default: 5e-14]")
  
    #parse arguments
    options=parser.parse_args()
  
  #if not fall back to older optparse
  else:
    parser=op.OptionParser(usage="Usage: %prog [options]"
      ,version="%prog 1.0"
      ,description="Tests that model restarts produce the same output as if the code had been run "\
      +"straight through.")
    parser.add_option("-f",action="store_true",dest="f"
      ,help="Force removal of pre-existing temporary data automatically."
      ,default=False)
    parser.add_option('-p',default=15,type="int",dest="p",help="Sets the amount of percision to"\
      +" use in ascii files. Effectively setting the precision of the comparison. "\
      +"[default: %default]")
    
    #parse command line options
    (options,args)=parser.parse_args()
  
  #set paths
  SPHERLSPath="/home/cgeroux/SPHERLS/trunk/SPHERLS"
  SPHERLSgenPath="/home/cgeroux/SPHERLSgen/trunk/SPHERLSgen"
  SPHERLSanalPath="/home/cgeroux/SPHERLSanal/trunk/SPHERLSanal"
  EOSFile="/home/cgeroux/SPHERLS/trunk/eos/eosY240Z002"
  velocityProfile="/home/cgeroux/SPHERLSgen/trunk/velocity_pro/T5700_L50_M575_fu.dat"
  
  #check paths
  if not os.access(SPHERLSPath,os.F_OK):
    print "\""+SPHERLSPath+"\" not found, stopping!"
    return
  if not os.access(SPHERLSgenPath,os.F_OK):
    print "\""+SPHERLSgenPath+"\" not found, stopping!"
    return
  if not os.access(SPHERLSanalPath,os.F_OK):
    print "\""+SPHERLSanalPath+"\" not found, stopping!"
    return
  if not os.access(EOSFile,os.F_OK):
    print "\""+EOSFile+"\" not found, stopping!"
    return
  if not os.access(velocityProfile,os.F_OK):
    print "\""+velocityProfile+"\" not found, stopping!"
    return
  numProcs=4
  
  test3DNARestarts([SPHERLSPath,SPHERLSgenPath,SPHERLSanalPath],EOSFile,velocityProfile,numProcs
    ,options)
def test3DNARestarts(exePaths,EOSFile,velocityProfile,numProcs,options):
  '''
  input:     path to SPHERLS and SPHERLSgen executables, number of processors to run test with
  output:    sucess of the test
  algorithm:
    1) make starting model
      a) gernerate a SPHERLSgen.xml file
      b) run SPHERLSgen to make starting model
    2) run code for 2 time steps
      a) generate a SPHERLS.xml file
      b) run SPHERLS for 3 timesteps with given number of processors
    4) restart code at first time step after initial model
    5) diff last models to see if they are the same
  '''
  SPHERLSgen_xml=\
    '''
    <data>
      <model type="stellar">
        
        <output>
          <timeStepFactor>0.25</timeStepFactor>
          <fileName>3DNARestartTest_t00000000</fileName>
          <binary>true</binary>
          <writeToScreen>false</writeToScreen>
        </output>
        <EOS type="table">
          <T-eff>6.1e3</T-eff>
          <L>50.0</L>
          <eosTable>'''+EOSFile+'''</eosTable>
          <tolerance>5e-14</tolerance>
        </EOS>
        <dimensions>
          <radIndepVar>
            <M-total>5.75E-01</M-total>
            <M-delta-init>1.2e-9</M-delta-init>
            <M-delta-delta stopType="T" stopValue="1e4">2.5e-2</M-delta-delta>
            <M-delta-delta stopType="T" stopValue="6e6">5e-2</M-delta-delta>
            <alpha>0.2</alpha>
            <num-1D>100</num-1D>
          </radIndepVar>
          <num-ghost-cells>2</num-ghost-cells><!-- number of ghost cells -->
          <num-theta>3</num-theta><!-- number of theta zones before number of ghost cells -->
          <delta-theta>1.0</delta-theta><!-- in degrees -->
          <num-phi>3</num-phi><!-- number of phi zones before number of ghost cells -->
          <delta-phi>1.0</delta-phi><!-- in degrees -->
        </dimensions>
        <velocityDist type="PRO">
          <fileName>'''+velocityProfile+'''</fileName>
          <uSurf>-2.0e5</uSurf>
          <perturb type="torus">
            <r_cen_off>2e11</r_cen_off>
            <theta_cen_off>0.0</theta_cen_off>
            <phi_cen_off>0.0</phi_cen_off>
            <radius_cen>1e9</radius_cen>
            <radius_outter>5e8</radius_outter>
            <width_guassian>2.5e8</width_guassian>
            <amplitude>1.0e20</amplitude>
          </perturb>
        </velocityDist>
      </model>
    </data>
    '''
  SPHERLS_xml=\
    '''
    <data>
      <procDims>
        <x0>'''+str(numProcs)+'''</x0>
        <x1>1</x1>
        <x2>1</x2>
      </procDims>
      <startModel>3DNARestartTest_t00000000</startModel>
      <outputName>3DNARestartTest2</outputName>
      <peakKE>false</peakKE>
      <prints type="normal">
        <frequency type="timeSteps">1</frequency>
      </prints>
      <dumps>
        <frequency type="timeSteps">1</frequency>
      </dumps>
      <eos>
        <tolerance>5e-14</tolerance>
        <max-iterations>50</max-iterations>
      </eos>
      <av>1.4</av>
      <av-threshold>0.01</av-threshold>
      <time>
        <endTime>13.0</endTime>
        <timeStep>7.0</timeStep>
      </time>
      <adiabatic>false</adiabatic>
      <turbMod>
        <type>smagorinsky</type>
        <eddyVisc>0.17</eddyVisc>
      </turbMod>
      <implicit>
        <numImplicitZones>150</numImplicitZones>
        <derivativeStepFraction>5e-7</derivativeStepFraction>
        <tolerance>5.0e-14</tolerance>
        <max-iterations>100</max-iterations>
      </implicit>
    </data>
    '''
  SPHERLS2_xml=\
    '''
    <data>
      <procDims>
        <x0>'''+str(numProcs)+'''</x0>
        <x1>1</x1>
        <x2>1</x2>
      </procDims>
      <startModel>3DNARestartTest2_t00000001</startModel>
      <outputName>3DNARestartTest3</outputName>
      <peakKE>false</peakKE>
      <prints type="normal">
        <frequency type="timeSteps">1</frequency>
      </prints>
      <dumps>
        <frequency type="timeSteps">1</frequency>
      </dumps>
      <eos>
        <tolerance>5e-14</tolerance>
        <max-iterations>50</max-iterations>
      </eos>
      <av>1.4</av>
      <av-threshold>0.01</av-threshold>
      <time>
        <endTime>13.0</endTime>
        <timeStep>7.0</timeStep>
      </time>
      <adiabatic>false</adiabatic>
      <turbMod>
        <type>smagorinsky</type>
        <eddyVisc>0.17</eddyVisc>
      </turbMod>
      <implicit>
        <numImplicitZones>150</numImplicitZones>
        <derivativeStepFraction>5e-7</derivativeStepFraction>
        <tolerance>5.0e-14</tolerance>
        <max-iterations>100</max-iterations>
      </implicit>
    </data>
    '''
  
  #check to see if we can make a tmp directory, and creating it if we can make tmp directory 
  #to perform test
  tmpDir="./restart_test_tmp"
  print "making a temporary directory, \""+tmpDir+"\" to store test data ...",
  
  #check for read and write permission in cwd
  if not os.access("./",os.R_OK) :
    print "FAILED"
    print "\n  Do not have read access in current directory \"",os.getcwd()\
      ,"\", which is need for this test."
    return False
  if not os.access("./",os.W_OK) :
    print "FAILED"
    print "  Do not have write access in current directory \"",os.getcwd()\
      ,"\", which is need for this test."
    return False
  
  #check to see if the temporary directory already exists, if so remove it and keep going
  if os.access(tmpDir,os.F_OK):
    if options.f:
      shutil.rmtree(tmpDir)
      os.mkdir(tmpDir)
      print "SUCCESS"
      print "  \""+tmpDir+"\" already existed, removed it"
    else:
      print "FAILED"
      print "  \""+tmpDir+"\" already exists, not removing it and stopping! Use \"-f\" to force"\
        +" removal."
      return False
  else:
    os.mkdir(tmpDir)
    print "SUCCESS"
  
  #make the directory
  print "changing into directory \""+tmpDir+"\" ..."
  os.chdir(tmpDir)
  
  #make SPHERLSgen.xml file
  print "making \"SPHERLSgen.xml\" ..."
  f=open("SPHERLSgen.xml",'w')
  f.write(SPHERLSgen_xml)
  f.close()
  
  #make SPHERLS.xml file
  print "making \"SPHERLS.xml\" ..."
  f=open("SPHERLS.xml",'w')
  f.write(SPHERLS_xml)
  f.close()
  
  #run SPHERLSgen
  print "running \""+exePaths[1]+"\" ...",
  log=open("log.txt",'w')
  log.write("GENERATING STARTING MODEL ...\n")
  log.close()
  log=open("log.txt",'a')
  result=subprocess.call(exePaths[1],stdout=log,stderr=log)
  if result!=0:
    print "FAILED"
    return False
  else:
    print "SUCCESS"
  
  #run 2 steps with SPHERLS
  print "running \""+exePaths[0]+"\" for 2 time steps ...",
  log.write("\nEVOLVING FOR 2 TIME STEPS ...\n")
  log.close()
  log=open("log.txt",'a')
  result=subprocess.call(["mpirun", "-np",str(numProcs),exePaths[0]],stdout=log,stderr=log)
  if result!=0:
    print "FAILED"
    return False
  else:
    print "SUCCESS"
  
  #remake SPHERLS.xml file with new configuration
  print "remaking \"SPHERLS.xml\" ..."
  f=open("SPHERLS.xml",'w')
  f.write(SPHERLS2_xml)
  f.close()
  
  #combine binary files
  print "combining binary dump \"./3DNARestartTest2_t00000001\" for restart ...",
  result=subprocess.call([exePaths[2],"-c","dbcb","./3DNARestartTest2_t00000001"])
  if result!=0:
    print "FAILED"
    return False
  else:
    print "SUCCESS"
  
  #run 1 step with SPHERLS
  print "restarting \""+exePaths[0]+"\" for 1 time step from dump \"3DNARestartTest2_t00000001\""\
    +" ...",
  log.write("\nEVOLVING FOR 1 TIME STEP ...\n")
  log.close()
  log=open("log.txt",'a')
  result=subprocess.call(["mpirun", "-np",str(numProcs),exePaths[0]],stdout=log,stderr=log)
  if result!=0:
    print "FAILED"
    return False
  else:
    print "SUCCESS"
  
  #combine final binary files and convert to ascii
  print "combining binary files for dump \"./3DNARestartTest2_t00000002\" ...",
  result=subprocess.call([exePaths[2],"-p",str(15),"-c","dbca"
    ,"./3DNARestartTest2_t00000002"])
  if result!=0:
    print "FAILED"
    return False
  else:
    print "SUCCESS"
  print "combining binary files for dump \"./3DNARestartTest3_t00000002\" ...",
  result=subprocess.call([exePaths[2],"-p",str(15),"-c","dbca"
    ,"./3DNARestartTest3_t00000002"])
  if result!=0:
    print "FAILED"
    return False
  else:
    print "SUCCESS"
  
  #compare final binary files
  print "diffing last dumps ...",
  log.write("\nDIFFING LAST MODEL DUMPS ...\n")
  log.close()
  log=open("log.txt",'a')
  result=diffDumps.diffDumps("./3DNARestartTest2_t00000002.txt","./3DNARestartTest3_t00000002.txt"
    ,options.p,options.t,log)
  log.close()
  if not result:
    print "FAILED"
    print "  model files \"./3DNARestartTest2_t00000002.txt\" and "\
      +"\"./3DNARestartTest3_t00000002.txt\" differ"
    
    #leave temporary directory so the failure can be inspected
    return False
  else:
    print "SUCCESS"
  
  #moving out of temporary directory
  os.chdir("../")
  
  #remove temporary directory if flag k not setting
  if not options.k:
    print "removing temporary directory ...",
    result=subprocess.call(["rm","-rf",tmpDir])
    if result!=0:
      print "FAILED"
      print "  unable to remove temporary directory \""+tmpDir+"\""
      return False
    else:
      print "SUCCESS"
  
  return True
if __name__ == "__main__":
  main()