#!/usr/bin/env python

import datafile
import optparse as op
import make_profiles
import glob
import math
import numpy as np
import sys
import os
import disect_filename
import average_PKE

def parseOptions():
  #setup command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START]"
    ,version="%prog 1.0",description="")
  parser.add_option("-k","--keep",action="store_true",dest="keep"
    ,help="Keeps distributed binary files [default].",default=True)
  parser.add_option("-r","--remove",action="store_false",dest="keep"
    ,help="Removes distributed binary files")
  parser.add_option("-m","--remake",action="store_true",dest="remake"
    ,help="Will remake profiles if they already exist. [not default].",default=False)
  
  parser.add_option("-o","--outputFile",dest="outputFile",default="out"
    ,help="Specifies the OUTPUTFILE. [default: %default]"
    ,metavar="OUTPUTFILE", type="string")
  parser.add_option("-f","--format",dest="format",default="png", type="choice",
    help="Sets the FMT of the OUTPUTFILE, availble formats are 'png', 'pdf', 'ps', 'eps', and 'svg'."
    +"[default: %default]", metavar="FMT",choices=('png','pdf','ps','eps','svg'))
  parser.add_option("-s","--show",dest="show",action="store_true",default=False
    ,help="Display plot to x11 display rather than saving to a file.")
  parser.add_option("--no-grid",action="store_true",dest="noGrid"
    ,help="Turns off grid [default].",default=False)
  parser.add_option("--points",action="store_true",dest="points",help="If set, will use points when"
    +" plotting in addition to lines [default: %default]",default=False)
  parser.add_option("--without-mass",action="store_false",dest="withMass",help="If set, will "
    +"include the zone mass in the calculations to give PKE in [ergs] rather than [ergs/g] [default]"
    ,default=True)
  parser.add_option("-a",action="store_true",dest="withAV",help="If set, will include the "
    +"artificial viscosity in the pressure [default: %default]",default=False)
  parser.add_option("--no-lines",action="store_true",dest="noLines",help="If set, will not use "
    +"lines when plotting, and only points [default: %default]",default=False)
  parser.add_option("--plot-curves",action="store_true",dest="plotCurves",help="If set, plot PdV "
    +"curves for all radial zones [default: %default]",default=False)
  parser.add_option("-n", "--zone",dest="startZone",
    help="Specifies the zone to start creating work plots from. Counted from the center of the "
    +"model.",type="int",default=2)
  parser.add_option("--re-sum",action="store_true",dest="resum"
    ,help="Will re-sum all model profiles, usefull when files have problems being made from "
    +"corruption and have to be re-made. Other wise should not be used as it takes more time"
    +" [not default].",default=False)
    
  #parse command line options
  return parser.parse_args()
def main():
  deltaMColumn=2#the column number for delta M in the profile file, start at 0 for the first column
  pColumn=58#the column number for P_ave in the profile file, starting at 0 for first column
  rhoColumn=6#the column number for Rho_ave in the profile file
  tColumn=45#the column number for T_ave in the profile file
  QColumn=35
  
  #parse command line options
  (options,args)=parseOptions()
  
  import matplotlib
  if not options.show:
    matplotlib.use('Agg')#needed for saving figures
  
  import matplotlib.pyplot as plt
  from matplotlib.gridspec import GridSpec
  
  #get base file name
  [start,end,baseFileName]=disect_filename.disectFileName(args[0])
  
  #make sure that all the combined binary files have profiles made
  failedFiles=make_profiles.make_profiles(options.keep,args[0],options.remake,False)
  
  #check that the averagePKE.txt file is there, and make sure it has all entries
  averagePKEFile=os.path.dirname(baseFileName)+"/averagePKE.txt"
  average_PKE.averagePKE(start,sys.maxint,baseFileName,options)
  
  #get periods in range indicated
  averagePKEFileData=datafile.DataFile()
  print __name__+": reading file \"",averagePKEFile,"\" ..."
  averagePKEFileData.readFile(averagePKEFile)
  periodCount=0
  periodRange=[]
  nSet=0
  time=[]
  for i in range(len(averagePKEFileData.fColumnValues)-1):
    if averagePKEFileData.fColumnValues[i][0]>=start and averagePKEFileData.fColumnValues[i][0]<end:
      if averagePKEFileData.fColumnValues[i][3]!=None and nSet==0:
        periodRange.append([averagePKEFileData.fColumnValues[i][0],None])
        nSet=1
      elif averagePKEFileData.fColumnValues[i][3]!=None and nSet==1:
        nSet=2
      elif averagePKEFileData.fColumnValues[i][3]!=None and nSet==2:
        if periodCount>0:
          periodRange.append([periodRange[periodCount-1][1],None])
        time.append(averagePKEFileData.fColumnValues[i][1])
        periodRange[periodCount][1]=averagePKEFileData.fColumnValues[i][0]
        nSet=1
        periodCount+=1
  dWSum=[]
  
  #get all possible file names
  extension="_pro"+".txt"
  filesExistProfiles=glob.glob(baseFileName+"*"+extension)
  filesExistProfiles.sort()
  
  #check for the t=0 model
  firstFile=baseFileName+"00000000"
  failedFiles2=make_profiles.make_profiles(options.keep,firstFile,options.remake,False)
  if not os.path.isfile(firstFile+extension):
    print "didn't find static model at ",firstFile+extension\
      ," using first model found in range given, ",filesExistProfiles[0]," instead"
    firstFile=filesExistProfiles[0]
  else:
    firstFile+=extension
  
  fileData=datafile.DataFile()
  fileData.readFile(firstFile)
  Log10T=fileData.fColumnValues[0:len(fileData.fColumnValues)-1,tColumn]
  Log10T=np.array(Log10T,dtype=np.float64)#force double precision
  Log10T=np.log10(Log10T)
  
  for n in range(len(periodRange)):
    
    #get and sort files
    files=[]
    for file in filesExistProfiles:
      intOfFile=int(file[len(baseFileName):len(file)-len(extension)])
      if intOfFile>=periodRange[n][0] and intOfFile<periodRange[n][1]:
        files.append(file)
    
    if len(files)<3:
      print "need more than 3 files to compute the work, likely alot more"
      return False
    
    #for first model dump
    fileData=datafile.DataFile()
    print "reading file ",files[0]," ..."
    fileData.readFile(files[0])
    
    #get temperature profile from first model, may want to change this
    Log10T=fileData.fColumnValues[0:len(fileData.fColumnValues)-1,tColumn]
    Log10T=np.array(Log10T,dtype=np.float64)#force double precision
    Log10T=np.log10(Log10T)
    
    #read in p, and 1/rho for first file
    p=np.zeros( (len(fileData.fColumnValues)-1,len(files)+1) )
    rhoInvert=np.zeros( (len(fileData.fColumnValues)-1,len(files)+1) )
    deltaM=np.zeros( (len(fileData.fColumnValues)-1,len(files)+1) )
    for i in range(len(fileData.fColumnValues)-1):
      p[i][0]=fileData.fColumnValues[i][pColumn]
      p[i][len(files)]=fileData.fColumnValues[i][pColumn]+fileData.fColumnValues[i][QColumn]
      if options.withAV:
        p[i][0]+=fileData.fColumnValues[i][QColumn]
        p[i][len(files)]+=fileData.fColumnValues[i][QColumn]
      rhoInvert[i][0]=1.0/fileData.fColumnValues[i][rhoColumn]
      rhoInvert[i][len(files)]=1.0/fileData.fColumnValues[i][rhoColumn]
      if options.withMass:
        deltaM[i][0]=fileData.fColumnValues[i][deltaMColumn]
        deltaM[i][len(files)]=fileData.fColumnValues[i][deltaMColumn]
      else:
        deltaM[i][0]=1.0
        deltaM[i][len(files)]=1.0
    
    #read all files in
    for i in range(1,len(files)):#for each dump
      print "reading file ",files[i]," ..."
      fileData.readFile(files[i])
      for j in range(len(fileData.fColumnValues)-1):#for each zone
        p[j][i]=fileData.fColumnValues[j][pColumn]#+fileData.fColumnValues[j][QColumn]
        if options.withAV:
          p[j][i]+=fileData.fColumnValues[j][QColumn]
        if options.withMass:
          deltaM[j][i]=fileData.fColumnValues[j][deltaMColumn]
        else:
          deltaM[j][i]=1.0
        rhoInvert[j][i]=1.0/fileData.fColumnValues[j][rhoColumn]
      
    #compute work
    dW=np.zeros(len(fileData.fColumnValues)-1)
    dWSum.append(0.0)
    for i in range(1,len(files)+1):
      for j in range(len(fileData.fColumnValues)-1):
        dW[j]+=(0.5*(rhoInvert[j][i]-rhoInvert[j][i-1])*(p[j][i]+p[j][i-1]))*deltaM[j][i]
    for j in range(len(fileData.fColumnValues)-1):
      dWSum[n]+=dW[j]
    
    #make plots of P and 1/rho for each zone
    fig=plt.figure(figsize=(13,8))
    ax1 = plt.subplot2grid((3,3), (0,0),colspan=3,rowspan=3)
    if options.plotCurves:
      for j in range(options.startZone,len(fileData.fColumnValues)-1):
        make_plot(rhoInvert[j],p[j],options,j,fig,ax1)
    
    #make work plot
    make_plot2(Log10T,dW,options,fig,ax1,"work_plot_"+str(n),"Log10(T)")
    make_plot2(range(len(fileData.fColumnValues)-1),dW,options,fig,ax1,"work_plot2_"+str(n),"zone #")
    
  #print out total work done by model
  f=open(options.outputFile+"_work_per_period.txt",'w')
  f.write("time[s]      work[ergs]\n")
  for i in range(len(dWSum)):
    f.write(str(time[i])+" "+str(dWSum[i])+"\n")
  f.close()
def make_plot(x,y,options,i,fig,ax1):
  if not options.noGrid:
    ax1.grid()
  plotString='-'
  if options.noLines:
    plotString='o'
  if options.points:
    plotString=plotString+'o'
  line1=ax1.plot(x,y,plotString,markersize=3)
  line1=ax1.plot(x[0],y[0],'x',markersize=15)
  ax1.set_xlabel("<1/rho>")
  ax1.set_ylabel("<P>")
  if options.show:
    plt.show()
  else:
    sOutFileName=options.outputFile+"_"+str(i)+"."+options.format
    print __name__+":"+main.__name__+": creating plot \""+sOutFileName+"\" ..."
    fig.savefig(sOutFileName,format=options.format,transparent=False)#save to file
  
  ax1.cla()#clear plot
  return True
def make_plot2(x,y,options,fig,ax1,name,xlabel):
  if not options.noGrid:
    ax1.grid()
  plotString='-'
  #if options.noLines:
  #  plotString='o'
  #if options.points:
  #  plotString=plotString+'o'
  #line1=ax1.plot(x,y,plotString,markersize=5)
  line1=ax1.plot(x,y,plotString,markersize=5)
  ax1.set_xlabel(xlabel)
  if options.withMass:
    ax1.set_ylabel("Work [ergs]")
  else:
    ax1.set_ylabel("Work [ergs/g]")
  if options.show:
    plt.show()
  else:
    sOutFileName=options.outputFile+"_"+name+"."+options.format
    print __name__+":"+main.__name__+": creating plot \""+sOutFileName+" ..."
    fig.savefig(sOutFileName,format=options.format,transparent=False)#save to file
  
  ax1.cla()#clear plot
  return True
if __name__ == "__main__":
  main()