#!/usr/bin/env python
import optparse as op
import os
import xml.etree.ElementTree as xml
import paths

def addParserOptions(parser):
  """Adds options to parser
  """
  return parser
def main():
  
  #create the parser
  parser=op.OptionParser(usage="Usage: %prog [options] DIRTOARCHIVE ARCHIVENAME"
    ,version="%prog 1.0"
    ,description="Searches subdirectories of the directory DIRTOARCHIVE for "
      +"\"SPHERLS.xml\" files. When one is found, it creates the same "
      +"directory structure in the directory ARCHIVENAME and copies over a "
      +"number of files usually associated with a run (e.g. SPHERLS.xml "
      +"SPHERLSgen.xml etc.) which are desirable to archive. It then submits a "
      +"job to the sun grid engine to move combined binary files form the "
      +"output directory to a new output directory in the same directory as "
      +"the other moved files. This process allows one to then archive the "
      +"new directory, which contains all the basic information from a "
      +"calculation using tar -czf ARCHIVENAME.tar.gz ARCHIVENAME.")

  #add parser options
  parser=addParserOptions(parser)
  
  #parse command line options
  (options,args)=parser.parse_args()
  
  if len(args)!=2:
    raise Exception("Must have two arguments")
  
  #get output file names and paths for running average_PKE.py in
  outputFilePaths=[]
  outputFileNames=[]
  jobNames=[]
  startModels=[]
  count=0
  srcdir=args[0]
  destdir=args[1]
  cwd=os.path.abspath('./')
  
  #if there isn't a destination directory make items
  if not os.path.isdir(destdir):
    os.makedirs(destdir)
  for root,dirs,files in os.walk(srcdir):
    #print root,dirs,files
    if "README.txt" in files:
      path=os.path.relpath(root,srcdir)
      file=os.path.join(path,"README.txt")
      os.renames(os.path.join(srcdir,file),os.path.join(destdir,file))
    if "SPHERLS.xml" in files:
      
      (outputPath,fileName,startModel,jobName)=getConfigOutputDIR(os.path.join(root,"SPHERLS.xml"),srcdir)
      runPath=os.path.dirname(outputPath)
      spherslConfig=os.path.join(runPath,"SPHERLS.xml")
      spherslGenConfig=os.path.join(runPath,"SPHERLSgen.xml")
      stdOut=os.path.join(runPath,jobName+".out")
      stdErr=os.path.join(runPath,jobName+".err")
      #print "outputPath=",outputPath
      #print "runPath=",runPath
      #print "spherslConfig=",spherslConfig
      #print "stdOut=",stdOut
      #print "stdErr=",stdErr
      
      #make directory in destination for output (and all other higher directories)
      if(not os.path.exists(os.path.join(destdir,outputPath))):
        os.makedirs(os.path.join(destdir,outputPath))
      
      #move over SPHERLS.xml, and SPHERLSgen.xml
      if(os.path.exists(os.path.join(srcdir,spherslConfig))):
        os.renames(os.path.join(srcdir,spherslConfig),os.path.join(destdir,spherslConfig))
      if(os.path.exists(os.path.join(srcdir,spherslGenConfig))):
        os.renames(os.path.join(srcdir,spherslGenConfig),os.path.join(destdir,spherslGenConfig))
      
      #move over start model
      #print "src start model=",os.path.join(srcdir,startModel)
      #print "dest start model=",os.path.join(destdir,startModel)
      if(os.path.exists(os.path.join(srcdir,startModel))):
        os.renames(os.path.join(srcdir,startModel),os.path.join(destdir,startModel))
      
      #move over std output and std err
      if(os.path.exists(os.path.join(srcdir,stdErr))):
        os.renames(os.path.join(srcdir,stdErr),os.path.join(destdir,stdErr))
      if(os.path.exists(os.path.join(srcdir,stdOut))):
        os.renames(os.path.join(srcdir,stdOut),os.path.join(destdir,stdOut))
      
      outputFilePaths.append(outputPath)
      outputFileNames.append(fileName)
    
    #remove output directories from search, potentially lots of files in there 
    #that we know we won't need to search through
    for dir in dirs:
      testDir=os.path.join(root,dir)
      if testDir in outputFilePaths:
        dirs.remove(dir)
  
  #create new directory structure, submit scripts, and run the job
  settings={}
  settings['shell']="/bin/bash"
  for i in range(len(outputFilePaths)):
    
    jobName=os.path.relpath(outputFilePaths[i])
    jobName=jobName.replace("/","_")
    jobName=jobName.replace("\\","_")
    print jobName,outputFileNames[i],outputFilePaths[i]
    
    settings['exe']=os.path.join(paths.scriptPath,"mv_files.py")+" --cb-only"
    settings['jobName']=jobName
    settings['shell']="/bin/bash"
    fileSet=os.path.join(outputFilePaths[i],outputFileNames[i])
    print "fileSet=",fileSet
    settings['src']=os.path.join(srcdir,fileSet)+"_t[0-*]"
    settings['dest']=os.path.join(destdir,fileSet)
    print "src fileSet=",settings['src']
    print "dest fileSet=",settings['dest']
    script=makeSubScript(settings,options)
    
    #run job
    cmd="qsub "+script
    os.system(cmd)
def getConfigOutputDIR(fileName,srcdir):
  
  #get some directories
  cwd=os.path.abspath("./")
  dir=os.path.dirname(fileName)
  srcdirabs=os.path.abspath(srcdir)
  
  #open xml configuration file
  tree=xml.parse(fileName)
  root=tree.getroot()
  
  #get output path and file name
  outputNameElement=root.find("outputName")
  outputPath=None
  outputFileName=None
  if outputNameElement!=None:
    text=outputNameElement.text
    if text!="":
      
      #change to directory of config file (paths are relative to it's location)
      os.chdir(dir)
      
      #get output file name from spherls config file
      outputFileName=os.path.basename(text)
      
      #Get output path from spherls config file
      outputPath=os.path.relpath(os.path.dirname(text),srcdirabs)
      
      #switch back to cwd of script
      os.chdir(cwd)
  
  #get start model
  startModelElement=root.find("startModel")
  startModel=None
  if startModelElement!=None:
    text=startModelElement.text
    if text!="":
      os.chdir(dir)
      
      startModelPath=os.path.relpath(os.path.dirname(text),srcdirabs)
      startModelFileName=os.path.basename(text)
      startModel=os.path.join(startModelPath,startModelFileName)
      
      os.chdir(cwd)
  
  #get job name to identify std output and err
  jobElement=root.find("job")
  if jobElement!=None:
    nameElement=jobElement.find("name")
    if nameElement!=None:
      jobName=nameElement.text
      if jobName=="":
        jobName=None
  
  return (outputPath,outputFileName,startModel,jobName)
def makeSubScript(settings,options):
  """
  Creates a submit script for the sun grid engine based on settings, and returns 
  the name of the script.
  """
  
  #additional hard coded settings
  scriptName=settings['jobName']+"_que.sh"
  f=open(scriptName,'wb')
  script="#!"+settings['shell']+"\n"\
    +"##\n"\
    +"##Shell to run job in\n"\
    +"#$ -S "+settings['shell']+"\n"\
    +"##\n"\
    +"## Set job name\n"\
    +"#$ -N J"+settings['jobName']+"\n"\
    +"##\n"\
    +"## Run job from current working directory\n"\
    +"#$ -cwd\n"\
    +"##\n"\
    +"#$ -j y\n"\
    +"## Output file name\n"\
    +"#$ -o "+settings['jobName']+".out\n"\
    +"##\n"\
    +"## Transfer all environment variables when job is submitted\n"\
    +"#$ -V\n"\
    +settings['exe']+" "+settings['src']+" "+settings['dest']
  
  f.write(script)
  f.close()
  return scriptName
if __name__ == "__main__":
  main()