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
import xml.etree.ElementTree as xml

#set number of processors to use
numProcs=4

#set paths to reference calculations
import ref_calcs

def main():
  
  #if argparse is availble
  if useArgparse:
    #set parser options
    parser=argparse.ArgumentParser(description="Tests that model restarts produce the same output"\
      +" as if the code had been run straight through.")
    parser.add_argument('-f',action="store_true",default=False,help="Force removal of"\
      +" pre-existing temporary data automatically.")
    parser.add_argument('-m',action="store_true",default=False,help="Force remaking of"\
      +" all reference calculations.")
    parser.add_argument('-r',action="store_true",default=False,help="Remake if missing"\
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
    parser.add_option('-m',action="store_true",dest="m",default=False,help="Force remaking of"\
      +" all reference calculations.")
    parser.add_option('-r',action="store_true",dest="r",default=False,help="Remake if missing"\
      +" reference calculations.")
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
            tmpDir="./test"+haveRef+"Calculation"
            result=subprocess.call(["rm","-rf",tmpDir])
            if result!=0:
              print "  unable to remove temporary directory \""+tmpDir+"\""
              return False
          
          print "SUCCESS"
        else:
          #print "FAILED"
          print "\n  see \"./test"+haveRef+"Calculation/log.txt\" for details of failure"
      else:
        print "FAILED"
        print "  \""+haveRef+"\" didn't have pre-calculated reference calculations, can't check"\
          +" against them."
    if len(failedTests)>0:#will never execute
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
  cwd=os.getcwd()
  if not createTestCalcNA(subDir,numProcs,options):
    return False
  os.chdir("./test"+subDir+"Calculation")
  
  #compare model calculations
  log=open("./log.txt",'a')
  log.write("\nDIFFING MODEL DUMPS ...\n")
  calculationsAreTheSame=True
  for i in range(1,len(ref_calcs.refCalcs[subDir])):
    refFileNameParts=ref_calcs.refCalcs[subDir][i].split('_t')
    file1="./"+subDir+"Test_t"+str(refFileNameParts[1])
    file2=ref_calcs.refCalcs[subDir][i]
    log.close()
    log=open("./log.txt",'a')
    log.write("diffing model files \""+file1+"\" and "+"\""+file2+"\" ... \n  ")
    result=diffDumps.diffDumps(file1,file2,options.p,options.t,log)
    log.close()
    log=open("./log.txt",'a')
    if not result:
      calculationsAreTheSame=False
  os.chdir(cwd)
  return calculationsAreTheSame
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
  tmpDir="./test"+subDir+"Calculation"
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
      print "  \""+tmpDir+"\" already exists, not removing it and skipping this comparison! Use"\
        +" \"-f\" to force removal."
      return False
  else:
    os.mkdir(tmpDir)
  cwd=os.getcwd()
  os.chdir(tmpDir)
  log=open("log.txt",'w')
  now=datetime.datetime.now()
  subDirsAllowed=["1DNA","2DNA","3DNA"]
  if subDir not in subDirsAllowed:
    log.write(createRefCalcDNA.__name__+": doesn't know how to handle calculation of subDir \""+subDir\
      +"\" only calcuation of directories "+str(subDirsAllowed)+" is supported")
    log.close()
    return False
  
  #make SPHERLS.xml file
  log.write("making \"SPHERLS.xml\" ...")
  configFile=os.path.dirname(ref_calcs.refCalcs[subDir][0])+"/SPHERLS.xml"
  cmd=["cp",configFile,"./SPHERLS.xml"]
  result=subprocess.call(cmd,stdout=log,stderr=log)
  setSPHERLSStartAndOutputModel("./SPHERLS.xml",ref_calcs.refCalcs[subDir][0],subDir+"Test")
  log.write("SUCCESS\n")
  
  #run 10 steps with SPHERLS
  log.write("\nEVOLVING FOR 10 TIME STEPS ...\n")
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
  os.chdir(cwd)
  return True
def checkForRefCalcsAndRemake(options):
  
  #check for calculations
  haveRefCalcs={}
  for key in ref_calcs.refCalcs.keys():
    
    #check that all model files are there
    haveRefCalc=True
    if not options.m:
      print "checking for reference calculations for \""+key+"\" calculation ... ",
      for j in range(len(ref_calcs.refCalcs[key])):
        if not os.access(ref_calcs.refCalcs[key][j],os.F_OK):
          haveRefCalc=False
    else:
      haveRefCalc=False
    
    #check to see if there is a SPHERLS.xml file in the test directory
    configFile=os.path.dirname(ref_calcs.refCalcs[key][0])+"/SPHERLS.xml"
    if not os.access(configFile,os.F_OK):
      haveRefCalc=False
    haveRefCalcs[key]=haveRefCalc
    
    #if a file is missing recreate reference calculations if possible
    if not haveRefCalc:
      if options.r or options.m:
        if not options.m :
          print "FAILED"
          print "  remaking \""+key+"\" reference calculations ... ",
        else:
          print "remaking \""+key+"\" reference calculations ... ",
        if key=='1DNA':
          if createRefCalcNA("1DNA",numProcs,options):
            print "SUCCESS"
        elif key=='2DNA':
          if createRefCalcNA("2DNA",numProcs,options):
            print "SUCCESS"
        elif key=='3DNA':
          if createRefCalcNA("3DNA",numProcs,options):
            print "SUCCESS"
        else:
          print "FAILED"
          print "unknown reference calculation \""+key+"\", can not remake"
      else:
        print "FAILED"
        print "  not making \""+key+"\" reference calculation, use -r option to remake reference"\
          +" calculations automatically as needed."
    else:
      print "SUCCESS"
  return haveRefCalcs
def getSPHERLSStartModel(pathToSPHERLSXML):
  '''Returns the text in the \"startModel\" node found in the \"SPHERLS.xml\" file pointed to by
  pathToSPHERLSXML.'''
  
  #open xml file, and get root node
  tree=xml.parse(pathToSPHERLSXML)
  root=tree.getroot()
  
  #get output node
  startModelElement=root.find("startModel")
  if startModelElement==None:
    return None
  elif startModelElement.text!="":
    return startModelElement.text
  else :
    return None
def getSPHERLSgenOutputModel(pathToSPHERLSgenXML):
  '''Returns the text in the \"startModel\" node found in the \"SPHERLS.xml\" file pointed to by
  pathToSPHERLSXML.'''
  
  #open xml file, and get root node
  tree=xml.parse(pathToSPHERLSgenXML)
  root=tree.getroot()
  
  #get output node
  outputElement=root.find("output")
  if outputElement==None:
    return None
  
  #get output fileName
  fileNameElement=output.find("fileName")
  if fileNameElement==None:
    return None
  elif fileNameElement.text!="":
    return fileNameElement.text
  else :
    return None
def setSPHERLSStartAndOutputModel(pathToSPHERLSxml,startModel,outputModel):
  '''modifies the file pathToSPHERLSxml to have a new output model, outputModel'''
  
  #open xml file, and get root node
  tree=xml.parse(pathToSPHERLSxml)
  root=tree.getroot()
  
  #get and set outputName node
  outputElement=root.find("outputName")
  if outputElement==None:
    return False
  outputElement.text=outputModel
  
  #get and set startModel node
  startModelElement=root.find("startModel")
  if startModelElement==None:
    return False
  startModelElement.text=startModel
  
  #write out new SPHERLS.xml file
  out=open(pathToSPHERLSxml,'w')
  out.write(xml.tostring(root))
  out.close()
  return True
def remakeSPHERLSXMLWithNewStartModel(pathToSPHERLSXML):
  
  #open xml file, and get root node
  tree=xml.parse(pathToSPHERLSXML)
  root=tree.getroot()
  
  xml.dump(tree)
  #include reference to new starting model
  #adjust end time/timestep as needed for new starting model
  return True
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
  cwd=os.getcwd()#keep cwd to return to it afterwards
  
  #establish that all required files to run SPHERLS are present, if not remake as needed keep as 
  #many pre-existing settings as is feasible
  if os.access(tmpDir,os.F_OK):#if subDir already there
    
    #change into directory
    os.chdir(tmpDir)
    log=open("log.txt",'w')
    log.write(str(datetime.datetime.now())+"\n")
    log.write("directory \""+tmpDir+"\" found\n")
    log.write("checking for needed files ...\n")
    if os.access(tmpDir+"/SPHERLS.xml",os.F_OK):#if SPHERLS.xml is already there
      
      log.write("found \""+tmpDir+"/SPHERLS.xml\"\n")
      
      #get model referenced by SPHERLS.xml
      modelPath=getSPHERLSStartModel(tmpDir+"/SPHERLS.xml")
      
      if os.access(modelPath,os.F_OK):#if model referenced by SPHERLS.xml already there
        log.write("found start model \""+modelPath+"\"\n")
        pass#nothing to be done here, except run SPHERLS
      elif os.access(tmpDir+"/SPHERLSgen.xml",os.F_OK):#elif SPHERLSgen.xml is already there
        print "FAILED"
        print "    Building reference calculations from scratch isn't yet supported, and probalby "
        print "    not as good as manually setting them up."
        return False
        '''print "start model not found, but \""+tmpDir+"/SPHERLSgen.xml\" found, using it to generate starting model"
        #remake starting model by running SPHERLSgen, and get output file name
        result=subprocess.call(paths.SPHERLSgenPath,stdout=log,stderr=log)
        startModel=getSPHERLSgenOutputModel(tmpDir+"/SPHERLSgen.xml")
        
        #remake SPHERLS.xml but keep as many settings as possible
        remakeSPHERLSXMLWithNewStartModel(tmpDir+"/SPHERLS.xml",startModel)'''
        
      else:#else no starting model and no SPHERLSgen.xml
        #remake SPHERLSgen.xml, starting model
        #remake SPHERLS.xml but keep as many settings as possible
          #include reference to new starting model
          #adjust end time/timestep as needed for new starting model
        print "FAILED"
        print "    Building reference calculations from scratch isn't yet supported, and probalby "
        print "    not as good as manually setting them up."
        return False
    else:#else no SPHERLS.xml
      print "FAILED"
      message="    file \"SPHERLS.xml\" not found in directory \""+tmpDir+"\""
      log.write(message)
      print message
  else:#else
    #make SPHERLS.xml, SPHERLSgen.xml, and starting model
    print "FAILED"
    message="    directory \""+tmpDir+"\" not found."
    print message
    return False
    
    '''
    #below is a starting of the implementation required here, much more is needed though.
    #make directory
    os.mkdir(tmpDir)
    
    #change into directory
    os.chdir(tmpDir)
    log=open("log.txt",'w')'''
  
  #run SPHERLS_run.py
  log.write("\nRUNNING SPHERLS\n")
  log.close()
  log=open("log.txt",'a')
  result=subprocess.call(paths.scriptPaths+"SPHERLS_run.py",stdout=log,stderr=log)
  if result!=0:
    os.chdir(cwd)
    return False
  
  #combine bins and remove distributed bins
  log.write("\ncombining binary dumps \"./3DNARef_t[0-*]\" and removing distributed binary files ...\n")
  log.close()
  log=open("log.txt",'a')
  modelPath=getSPHERLSStartModel(tmpDir+"/SPHERLS.xml")
  fileName=os.path.basename(modelPath)
  fileNameParts=fileName.split('_t')
  cmd=[paths.scriptPaths+"combine_bins.py","-r","-m","./"+fileNameParts[0]+"_t["\
    +str(int(fileNameParts[1]))+"-*]"]
  result=subprocess.call(cmd,stdout=log,stderr=log)
  if result!=0:
    os.chdir(cwd)
    return False
  return True
if __name__ == "__main__":
  main()