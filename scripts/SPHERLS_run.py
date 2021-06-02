#!/usr/bin/env python
import os
import xml.etree.ElementTree as xml
import paths
import optparse as op

def parseXML(fileName):
  '''Parses the SPHERLS.xml file to figure out settings for submit script and how to run SPHERLS
  returns settings, a dictionary with entries for:
  
  numProcs    = an integer number of process to submit the job with
  jobName     = a string containing the name of the job
  jobDuration = a string containg the duration of the job
  totalview   = logical containing either True if running with totalview debugging other wise False
  que         = a string that is either "test" for the test que or "main" for the main que, if the
                string is None, then run on the command line and don't submit to a que
  email       = set to the email to email when job is completed, if False or not given no email will
                be sent.
  '''
  
  #get root element of fileName
  settings={}
  tree=xml.parse(fileName)
  root=tree.getroot()
  
  #get job element
  jobElement=root.find("job")
  settings["runInQueue"]=True
  if jobElement==None:
    print "WARNING: No \"job\" element found under \"data\" node in file \""+fileName+"\", running directly and not submitting to the queue!"
    settings["runInQueue"]=False
  
  #get procDims element
  procDimsElement=root.find("procDims")
  if procDimsElement==None:
    print "No \"procDims\" element found under \"data\" node in file \""+fileName+"\", quitting!"
    quit()
  
  #get x0 element
  x0Element=procDimsElement.find("x0")
  if x0Element==None:
    print "No \"x0\" element found under \"procDims\" node in file \""+fileName+"\", quitting!"
    quit()
  try:
    x0=int(x0Element.text)
  except ValueERror:
    print "\"x0\" element found under \"procDims\" node in file \""+fileName+"\" is not an integer, quitting!"
    quit()
  
  #get x1 element
  x1Element=procDimsElement.find("x1")
  if x1Element==None:
    print "No \"x1\" element found under \"procDims\" node in file \""+fileName+"\", quitting!"
    quit()
  try:
    x1=int(x1Element.text)
  except ValueERror:
    print "\"x1\" element found under \"procDims\" node in file \""+fileName+"\" is not an integer, quitting!"
    quit()
  
  #get x2 element
  x2Element=procDimsElement.find("x2")
  if x2Element==None:
    print "No \"x2\" element found under \"procDims\" node in file \""+fileName+"\", quitting!"
    quit()
  try:
    x2=int(x2Element.text)
  except ValueERror:
    print "\"x2\" element found under \"procDims\" node in file \""+fileName+"\" is not an integer, quitting!"
    quit()
  settings['numProcs']=str((x0-1)*x1*x2+1)
  
  if jobElement!=None:
    
    #get cmd
    cmdElement=jobElement.find("cmd")
    if cmdElement==None:
      settings['cmd']=None
    else:
      settings['cmd']=cmdElement.text
    
    #get template-file
    templateFileElement=jobElement.find("template-file")
    if templateFileElement==None:
      settings['templateFile']=None
    else:
      if templateFileElement.text=="":
        print "\"template-file\" element is empty. Must provide a file name and path to use as a template job script file."
        quit()
      settings['templateFile']=templateFileElement.text
    
    #get generated-file
    generatedFileElement=jobElement.find("generated-file")
    if generatedFileElement==None:
      settings['generatedFile']=None
    else:
      if generatedFileElement.text=="" or generatedFileElement.text==None:
        print "\"generated-file\" element is empty. Must provide a file name and path to use as a template job script file."
        quit()
      settings['generatedFile']=generatedFileElement.text
    
    #get string replacements
    replacementsElement=jobElement.find("replacements")
    settings['replacements']=[('NUMPROCS',str(settings['numProcs']))]
    if replacementsElement!=None:
      replacementElements=replacementsElement.findall("replacement")
      for replacementElement in replacementElements:
        searchStr=replacementElement.find("search-str").text
        if searchStr=="" or searchStr==None:
          print "found empty \"search-str\" element but this element may not be empty"
          quit()
        sub=replacementElement.find("substitution").text
      settings['replacements'].append((searchStr,sub))
  
  return settings
def makeSubScript(settings):
  '''Make submission script using the specified template file and string 
    replacements
    '''
  
  file=open(settings['templateFile'],'r')
  fileText=file.read()
  file.close()
  for replacement in settings['replacements']:
    fileText=fileText.replace(replacement[0],replacement[1])
  file=open(settings['generatedFile'],'w')
  file.write(fileText)
  file.close()
def main():
  
  #make command line parser
  parser=op.OptionParser(usage="Usage: %prog [options]"
    ,version="%prog 1.0"
    ,description="Starts a SPHERLS run by parsing the SPHERLS.xml configuration file")
  parser.add_option("-d",action="store_true",dest="dryRun"
    ,help="If set it print out the run command but not execute it [not default]."
    ,default=False)
  
  #parse command line options
  (options,args)=parser.parse_args()
  
  settings=parseXML("SPHERLS.xml")
  
  if settings["runInQueue"]:
    
    makeSubScript(settings)
    
    #additional hard coded settings
    cmd=settings["cmd"]+" "+settings["generatedFile"]
    if options.dryRun:
      print cmd
    else:
      os.system(cmd)
  else:
    cmd="mpirun -np "+str(settings['numProcs'])+" SPHERLS"
    print(cmd)
    os.system(cmd)
if __name__ == "__main__":
  main()
  