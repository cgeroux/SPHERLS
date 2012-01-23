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
import paths
import datetime
def main():
  
  #if argparse is availble
  if useArgparse:
    #set parser options
    parser=argparse.ArgumentParser(description="Tests that model restarts produce the same output"\
      +" as if the code had been run straight through.")
    parser.add_argument('-f',action="store_true",default=False,help="Force removal of"\
      +" pre-existing temporary data automatically.")
    parser.add_argument('-m',action="store_true",default=False,help="Force remaking of"\
      +" reference calculations.")
    parser.add_argument('-k',action="store_true",default=False,help="Keep temporary directories.")
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
    parser.add_option('-k',action="store_true",dest="k",default=False,help="Keep temporary directories.")
    parser.add_option('-p',default=5e-14,type=float,dest="p",help="Sets the amount of relative"\
      +" absolute difference allowed between numerical values of the two dump files."\
      +" [default: 5e-14]")
    parser.add_option('-t',default=5e-14,type=float,dest="t",help="Sets the smallest number to"\
      +" consider for numerical comparison, differences in numbers which are smaller than this "\
      +"threshold will be ignored for the comparison. [default: 5e-14]")
    
    #parse command line options
    (options,args)=parser.parse_args()
  
  #check paths
  paths.check_paths()
  
  #check that we have reference calculations, if not create them
  haveRefCals=checkForRefCalcsAndRemake(options)
  
  #perform calculation tests
  failedTests=[]
  failedTestDirs=[]
  
  #if not recreating reference calculations
  if not options.m:
    for haveRef in haveRefCals:
      print "checking calculation against \""+paths.ref_calcs+haveRef+"\" calculation ...",
      if haveRefCals[haveRef]:
        if checkCalAgainstRef(haveRef,options):
          
          #remove temporary directory
          if not options.k:
            tmpDir="./"+haveRef+"CalculationTest_tmp"
            result=subprocess.call(["rm","-rf",tmpDir])
            if result!=0:
              print "  unable to remove temporary directory \""+tmpDir+"\""
              return False
          print "SUCCESS"
        else:
          print "  see \"./"+haveRef+"CalculationTest_tmp/log.txt\" for details of failure ..."
      else:
        print "FAILED"
        print "  \""+haveRef+"\" didn't have pre-calculated reference calculations, can't check"\
          +" against them."
    if len(failedTests)>0:
      print "The following tests failed:"
      i=0
      for failedTest in failedTests:
        print "  "+failedTest+" : see \""+failedTestDirs[i]+"/log.txt\" for details on why the"\
          +" test failed"
        i=i+1
  else:
    print "can not do comparison as reference calculations were just remade and will nessessarily "\
      +"be identical to test calculations"
def checkCalAgainstRef(subDir,options):
  
  #create Test calculation to compare to reference calculations
  numProcs=4
  if not createTestCalcNA(subDir,numProcs,options):
    return False
  os.chdir("./"+subDir+"CalculationTest_tmp")
  
  #compare model calculations
  log=open("./log.txt",'a')
  log.write("\nDIFFING MODEL DUMPS ...\n")
  for i in range(0,10):
    file1="./"+subDir+"Test_t0000000"+str(i)
    file2=paths.ref_calcs+subDir+"/"+subDir+"Ref_t0000000"+str(i)
    log.close()
    log=open("./log.txt",'a')
    log.write("diffing model files \""+file1+"\" and "+"\""+file2+"\" ... \n  ")
    result=diffDumps.diffDumps(file1,file2,options.p,options.t,log)
    log.close()
    log=open("./log.txt",'a')
    if not result:
      #leave temporary directory so the failure can be inspected
      os.chdir("../")
      return False
  os.chdir("../")
  return True
def createTestCalcNA(subDir,numProcs,options):
  '''
  1) create starting model
    a) create directories as needed
    b) make SPHERLSgen configuration file
    c) run SPHERLSgen
  2) run 10 time steps, keeping models at t00000000, t00000001, and t00000010
    a) create SPHERLS configurations file
    b) call SPHERLS_run.py to run the code for 10 time steps
    c) remove any models not needed for comparison
  '''
  
  #check to see if we can make a tmp directory, and creating it if we can make tmp directory 
  #to perform test
  tmpDir="./"+subDir+"CalculationTest_tmp"
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
  
  #check to see if the temporary directory already exists, if so remove it and make a new one
  if os.access(tmpDir,os.F_OK):
    if options.f:
      shutil.rmtree(tmpDir)
      os.mkdir(tmpDir)
    else:
      print "FAILED"
      print "  \""+tmpDir+"\" already exists, not removing it and stopping! Use \"-f\" to force"\
        +" removal."
      return False
  else:
    os.mkdir(tmpDir)
  
  os.chdir(tmpDir)
  log=open("log.txt",'w')
  now=datetime.datetime.now()
  subDirsAllowed=["1DNA","2DNA","3DNA"]
  if subDir not in subDirsAllowed:
    log.write(createRefCalcDNA.__name__+": doesn't know how to handle calculation of subDir \""+subDir\
      +"\" only calcuation of directories "+str(subDirsAllowed)+" is supported")
    log.close()
    return False
  log.write("GENERATING STARTING MODEL ...\n")
  SPHERLSgen_xml=\
    '''
    <data>
      <model type="stellar">
        
        <output>
          <timeStepFactor>0.25</timeStepFactor>
          <fileName>'''+subDir+'''Test_t00000000</fileName>
          <binary>true</binary>
          <writeToScreen>false</writeToScreen>
        </output>
        <EOS type="table">
          <T-eff>6.1e3</T-eff>
          <L>50.0</L>
          <eosTable>'''+paths.EOSPath+'''</eosTable>
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
          <num-ghost-cells>2</num-ghost-cells><!-- number of ghost cells -->'''
  if subDir=="3DNAD":
    SPHERLSgen_xml=SPHERLSgen_xml+'''
          <num-theta>3</num-theta><!-- number of theta zones before number of ghost cells -->
          <delta-theta>1.0</delta-theta><!-- in degrees -->
          <num-phi>3</num-phi><!-- number of phi zones before number of ghost cells -->
          <delta-phi>1.0</delta-phi><!-- in degrees -->'''
  elif subDir=="2DNAD":
    SPHERLSgen_xml=SPHERLSgen_xml+'''
          <num-theta>3</num-theta><!-- number of theta zones before number of ghost cells -->
          <delta-theta>1.0</delta-theta><!-- in degrees -->
          <num-phi>1</num-phi><!-- number of phi zones before number of ghost cells -->
          <delta-phi>1.0</delta-phi><!-- in degrees -->'''
  else:
    SPHERLSgen_xml=SPHERLSgen_xml+'''
          <num-theta>1</num-theta><!-- number of theta zones before number of ghost cells -->
          <delta-theta>1.0</delta-theta><!-- in degrees -->
          <num-phi>1</num-phi><!-- number of phi zones before number of ghost cells -->
          <delta-phi>1.0</delta-phi><!-- in degrees -->'''
  SPHERLSgen_xml=SPHERLSgen_xml+'''
      </dimensions>
      <velocityDist type="PRO">
        <fileName>'''+paths.velocityProfilePath+'''</fileName>
        <uSurf>-2.0e5</uSurf>'''
  if subDir=="3DNAD":
    SPHERLSgen_xml=SPHERLSgen_xml+'''
          <perturb type="torus">
            <r_cen_off>2e11</r_cen_off>
            <theta_cen_off>0.0</theta_cen_off>
            <phi_cen_off>0.0</phi_cen_off>
            <radius_cen>1e9</radius_cen>
            <radius_outter>5e8</radius_outter>
            <width_guassian>2.5e8</width_guassian>
            <amplitude>1.0e20</amplitude>
          </perturb>'''
  elif subDir=="2DNAD":
    SPHERLSgen_xml=SPHERLSgen_xml+''''''
  else :
    SPHERLSgen_xml=SPHERLSgen_xml+''''''
  SPHERLSgen_xml=SPHERLSgen_xml+'''
        </velocityDist>
      </model>
    </data>
    '''
  SPHERLS_xml=\
    '''
    <data>
      <job>
        <que>false</que>
      </job>
      <procDims>
        <x0>'''+str(numProcs)+'''</x0>
        <x1>1</x1>
        <x2>1</x2>
      </procDims>
      <startModel>'''+subDir+'''Test_t00000000</startModel>
      <outputName>'''+subDir+'''Test</outputName>
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
        <endTime>69.0</endTime>
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
  
  #make SPHERLSgen.xml file
  log.write("making \"SPHERLSgen.xml\" ... ")
  f=open("SPHERLSgen.xml",'w')
  f.write(SPHERLSgen_xml)
  f.close()
  log.write("SUCCESS\n")
  
  #run SPHERLSgen
  log.write("running \""+paths.SPHERLSgenPath+"\" ...\n",)
  log.close()
  log=open("log.txt",'a')
  result=subprocess.call(paths.SPHERLSgenPath,stdout=log,stderr=log)
  if result!=0:
    os.chdir("../")
    return False
  
  #run 10 steps with SPHERLS
  log.write("\nEVOLVING FOR 10 TIME STEPS ...\n")
  #make SPHERLS.xml file
  log.write("making \"SPHERLS.xml\" ...")
  f=open("SPHERLS.xml",'w')
  f.write(SPHERLS_xml)
  f.close()
  log.write("SUCCESS\n")
  log.write( "running \""+paths.scriptPaths+"SPHERLS_run.py"+"\" for 10 time steps ...\n",)
  log.close()
  log=open("log.txt",'a')
  result=subprocess.call(paths.scriptPaths+"SPHERLS_run.py",stdout=log,stderr=log)
  if result!=0:
    os.chdir("../")
    return False
  
  #combine bins and remove distributed bins
  log.write("\ncombining binary dumps \"./"+subDir+"Test_t[0-*]\" ...\n")
  log.close()
  log=open("log.txt",'a')
  result=subprocess.call([paths.scriptPaths+"combine_bins.py","./"+subDir+"Test_t[0-*]"],stdout=log,stderr=log)
  if result!=0:
    os.chdir("../")
    return False
  os.chdir("../")
  return True
def checkForRefCalcsAndRemake(options):
  
  #set number of processors to use
  numProcs=4
  
  #set paths
  refCalcs={
    '1DNA':
    [paths.ref_calcs+"/1DNA/1DNARef_t00000000"
    ,paths.ref_calcs+"/1DNA/1DNARef_t00000001"
    ,paths.ref_calcs+"/1DNA/1DNARef_t00000002"
    ,paths.ref_calcs+"/1DNA/1DNARef_t00000003"
    ,paths.ref_calcs+"/1DNA/1DNARef_t00000004"
    ,paths.ref_calcs+"/1DNA/1DNARef_t00000005"
    ,paths.ref_calcs+"/1DNA/1DNARef_t00000006"
    ,paths.ref_calcs+"/1DNA/1DNARef_t00000007"
    ,paths.ref_calcs+"/1DNA/1DNARef_t00000008"
    ,paths.ref_calcs+"/1DNA/1DNARef_t00000009"
    ,paths.ref_calcs+"/1DNA/1DNARef_t00000010"]
    ,'2DNA':
    [paths.ref_calcs+"/2DNA/2DNARef_t00000000"
    ,paths.ref_calcs+"/2DNA/2DNARef_t00000001"
    ,paths.ref_calcs+"/2DNA/2DNARef_t00000002"
    ,paths.ref_calcs+"/2DNA/2DNARef_t00000003"
    ,paths.ref_calcs+"/2DNA/2DNARef_t00000004"
    ,paths.ref_calcs+"/2DNA/2DNARef_t00000005"
    ,paths.ref_calcs+"/2DNA/2DNARef_t00000006"
    ,paths.ref_calcs+"/2DNA/2DNARef_t00000007"
    ,paths.ref_calcs+"/2DNA/2DNARef_t00000008"
    ,paths.ref_calcs+"/2DNA/2DNARef_t00000009"
    ,paths.ref_calcs+"/2DNA/2DNARef_t00000010"]
    ,'3DNA':
    [paths.ref_calcs+"/3DNA/3DNARef_t00000000"
    ,paths.ref_calcs+"/3DNA/3DNARef_t00000001"
    ,paths.ref_calcs+"/3DNA/3DNARef_t00000002"
    ,paths.ref_calcs+"/3DNA/3DNARef_t00000003"
    ,paths.ref_calcs+"/3DNA/3DNARef_t00000004"
    ,paths.ref_calcs+"/3DNA/3DNARef_t00000005"
    ,paths.ref_calcs+"/3DNA/3DNARef_t00000006"
    ,paths.ref_calcs+"/3DNA/3DNARef_t00000007"
    ,paths.ref_calcs+"/3DNA/3DNARef_t00000008"
    ,paths.ref_calcs+"/3DNA/3DNARef_t00000009"
    ,paths.ref_calcs+"/3DNA/3DNARef_t00000010"]
    }
  
  #check for calculations
  haveRefCalcs={}
  for key in refCalcs.keys():
    
    #check that all files are there
    haveRefCalc=True
    if not options.m:
      print "checking for reference calculations for \""+key+"\" calculation ... ",
      for j in range(len(refCalcs[key])):
        if not os.access(refCalcs[key][j],os.F_OK):
          haveRefCalc=False
    else:
      haveRefCalc=False
    haveRefCalcs[key]=haveRefCalc
    
    #if a file is missing recreate reference calculations
    if not haveRefCalc:
      if not options.m:
        print "FAILED"
        print "  remaking \""+key+"\" reference calculations ... ",
      else:
        print "remaking \""+key+"\" reference calculations ... ",
      if key=='1DNA':
        if createRefCalcNA("1DNA",numProcs,options):
          print "SUCCESS"
        else:
          print "FAILED"
      elif key=='2DNA':
        if createRefCalcNA("2DNA",numProcs,options):
          print "SUCCESS"
        else:
          print "FAILED"
      elif key=='3DNA':
        if createRefCalcNA("3DNA",numProcs,options):
          print "SUCCESS"
        else:
          print "FAILED"
      else:
        print "unknown reference calculation \""+key+"\", can not remake"
    else:
      print "SUCCESS"
  return haveRefCalcs
def createRefCalcNA(subDir,numProcs,options):
  '''
  1) create starting model
    a) create directories as needed
    b) make SPHERLSgen configuration file
    c) run SPHERLSgen
  2) run 10 time steps, keeping models at t00000000, t00000001, and t00000010
    a) create SPHERLS configurations file
    b) call SPHERLS_run.py to run the code for 10 time steps
    c) remove any models not needed for comparison
  '''
  tmpDir=paths.ref_calcs+subDir
  if not os.access(tmpDir,os.F_OK):
    os.mkdir(tmpDir)
  else:
    shutil.rmtree(tmpDir)
    os.mkdir(tmpDir)
  
  #change into directory
  os.chdir(tmpDir)
  
  log=open("log.txt",'w')
  now=datetime.datetime.now()
  log.write("These reference models were created on "+str(now)+"\n\n")
  subDirsAllowed=["1DNA","2DNA","3DNA"]
  if subDir not in subDirsAllowed:
    log.write(createRefCalcDNA.__name__+": doesn't know how to handle calculation of subDir \""+subDir\
      +"\" only calcuation of directories "+str(subDirsAllowed)+" is supported")
    log.close()
    return False
  log.write("GENERATING STARTING MODEL ...\n")
  SPHERLSgen_xml=\
    '''
    <data>
      <model type="stellar">
        
        <output>
          <timeStepFactor>0.25</timeStepFactor>
          <fileName>'''+subDir+'''Ref_t00000000</fileName>
          <binary>true</binary>
          <writeToScreen>false</writeToScreen>
        </output>
        <EOS type="table">
          <T-eff>6.1e3</T-eff>
          <L>50.0</L>
          <eosTable>'''+paths.EOSPath+'''</eosTable>
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
          <num-ghost-cells>2</num-ghost-cells><!-- number of ghost cells -->'''
  if subDir=="3DNAD":
    SPHERLSgen_xml=SPHERLSgen_xml+'''
          <num-theta>3</num-theta><!-- number of theta zones before number of ghost cells -->
          <delta-theta>1.0</delta-theta><!-- in degrees -->
          <num-phi>3</num-phi><!-- number of phi zones before number of ghost cells -->
          <delta-phi>1.0</delta-phi><!-- in degrees -->'''
  elif subDir=="2DNAD":
    SPHERLSgen_xml=SPHERLSgen_xml+'''
          <num-theta>3</num-theta><!-- number of theta zones before number of ghost cells -->
          <delta-theta>1.0</delta-theta><!-- in degrees -->
          <num-phi>1</num-phi><!-- number of phi zones before number of ghost cells -->
          <delta-phi>1.0</delta-phi><!-- in degrees -->'''
  else:
    SPHERLSgen_xml=SPHERLSgen_xml+'''
          <num-theta>1</num-theta><!-- number of theta zones before number of ghost cells -->
          <delta-theta>1.0</delta-theta><!-- in degrees -->
          <num-phi>1</num-phi><!-- number of phi zones before number of ghost cells -->
          <delta-phi>1.0</delta-phi><!-- in degrees -->'''
  SPHERLSgen_xml=SPHERLSgen_xml+'''
      </dimensions>
      <velocityDist type="PRO">
        <fileName>'''+paths.velocityProfilePath+'''</fileName>
        <uSurf>-2.0e5</uSurf>'''
  if subDir=="3DNAD":
    SPHERLSgen_xml=SPHERLSgen_xml+'''
          <perturb type="torus">
            <r_cen_off>2e11</r_cen_off>
            <theta_cen_off>0.0</theta_cen_off>
            <phi_cen_off>0.0</phi_cen_off>
            <radius_cen>1e9</radius_cen>
            <radius_outter>5e8</radius_outter>
            <width_guassian>2.5e8</width_guassian>
            <amplitude>1.0e20</amplitude>
          </perturb>'''
  elif subDir=="2DNAD":
    SPHERLSgen_xml=SPHERLSgen_xml+''''''
  else :
    SPHERLSgen_xml=SPHERLSgen_xml+''''''
  SPHERLSgen_xml=SPHERLSgen_xml+'''
        </velocityDist>
      </model>
    </data>
    '''
  SPHERLS_xml=\
    '''
    <data>
      <job>
        <que>false</que>
      </job>
      <procDims>
        <x0>'''+str(numProcs)+'''</x0>
        <x1>1</x1>
        <x2>1</x2>
      </procDims>
      <startModel>'''+subDir+'''Ref_t00000000</startModel>
      <outputName>'''+subDir+'''Ref</outputName>
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
        <endTime>69.0</endTime>
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
  
  #make SPHERLSgen.xml file
  log.write("making \"SPHERLSgen.xml\" ... ")
  f=open("SPHERLSgen.xml",'w')
  f.write(SPHERLSgen_xml)
  f.close()
  log.write("SUCCESS\n")
  
  #run SPHERLSgen
  log.write("running \""+paths.SPHERLSgenPath+"\" ...\n",)
  log.close()
  log=open("log.txt",'a')
  result=subprocess.call(paths.SPHERLSgenPath,stdout=log,stderr=log)
  if result!=0:
    os.chdir("../")
    return False
  
  #run 10 steps with SPHERLS
  log.write("\nEVOLVING FOR 10 TIME STEPS ...\n")
  #make SPHERLS.xml file
  log.write("making \"SPHERLS.xml\" ...")
  f=open("SPHERLS.xml",'w')
  f.write(SPHERLS_xml)
  f.close()
  log.write("SUCCESS\n")
  log.write( "running \""+paths.scriptPaths+"SPHERLS_run.py"+"\" for 10 time steps ...\n",)
  log.close()
  log=open("log.txt",'a')
  result=subprocess.call(paths.scriptPaths+"SPHERLS_run.py",stdout=log,stderr=log)
  if result!=0:
    os.chdir("../")
    return False
  
  #combine bins and remove distributed bins
  log.write("\ncombining binary dumps \"./3DNARef_t[0-*]\" and removing distributed binary files ...\n")
  log.close()
  log=open("log.txt",'a')
  result=subprocess.call([paths.scriptPaths+"combine_bins.py","-r","./"+subDir+"Ref_t[0-*]"],stdout=log,stderr=log)
  if result!=0:
    os.chdir("../")
    return False
  return True
if __name__ == "__main__":
  main()