#!/usr/bin/env python

import datafile
import optparse as op
import make_profiles
import glob
import math
import numpy as np
import sys

def parseOptions():
  #setup command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END]"
    ,version="%prog 1.0",description="Creates a plot of the many luminosities as a function of time.")
  parser.add_option("-o","--outputFile",dest="outputFile",default="out"
    ,help="Specifies the OUTPUTFILE. [default: %default]"
    ,metavar="OUTPUTFILE", type="string")
  parser.add_option("-f","--format",dest="format",default="png", type="choice",
    help="Sets the FMT of the OUTPUTFILE, availble formats are 'png', 'pdf', 'ps', 'eps', and 'svg'."
    +"[default: %default]", metavar="FMT",choices=('png','pdf','ps','eps','svg'))
  parser.add_option("-s","--show",dest="show",action="store_true",default=False
    ,help="Display plot to x11 display rather than saving to a file.")
  parser.add_option("--xmin",type="float",dest="xmin",default=None,help="Sets the minimum x value "
    +"to be plotted. [default: %default]")
  parser.add_option("--xmax",type="float",dest="xmax",default=None,help="Sets the maximum x value "
    +"to be plotted. [default: %default]")
  parser.add_option("--ymin",type="float",dest="ymin",default=None,help="Sets the minimum y value "
    +"to be plotted. [default: %default]")
  parser.add_option("--ymax",type="float",dest="ymax",default=None,help="Sets the maximum y value "
    +"to be plotted. [default: %default]")
  parser.add_option("--no-grid",action="store_true",dest="noGrid"
    ,help="Turns off grid [default].",default=False)
  parser.add_option("-k","--keep",action="store_true",dest="keep"
    ,help="Keeps distributed binary files [default].",default=True)
  parser.add_option("-r","--remove",action="store_false",dest="keep"
    ,help="Removes distributed binary files")
  parser.add_option("--points",action="store_true",dest="points",help="If set, will use points when"
    +" plotting in addition to lines [default: %default]",default=False)
  parser.add_option("--no-lines",action="store_true",dest="noLines",help="If set, will not use "
    +"lines when plotting, and only points [default: %default]",default=False)
  #parse command line options
  return parser.parse_args()
  
def main():
  #parse command line options
  (options,args)=parseOptions()
  
  import matplotlib
  if not options.show:
    matplotlib.use('Agg')#needed for saving figures
  
  import matplotlib.pyplot as plt
  from matplotlib.gridspec import GridSpec
  
  #get base file name
  fileName=args[0]
  parts0=fileName.partition('[')
  baseFileName=parts0[0]
  if(parts0[1]=='' or parts0[2]==''):
    print __name__+":"+main.__name__+": error parsing file range from argument \""+fileName+"\""
    print "  expecting something like BASEFILENAME_t[STARTINDEX-ENDINDEX] ."
    return False
  
  #get interation range
  parts1=parts0[2].partition(']')
  if(parts1[1]==''):
    print __name__+":"+main.__name__+": error parsing file range from argument \""+fileName+"\""
    print "  expecting something like BASEFILENAME_t[STARTINDEX-ENDINDEX] ."
    return False
  parts2=parts1[0].partition('-')
  if(parts2[1]=='' or parts2[2]==''):
    print __name__+":"+main.__name__+": error parsing file range from argument \""+fileName+"\""
    print "  expecting something like BASEFILENAME_t[STARTINDEX-ENDINDEX] ."
    return False
  start=int(parts2[0])
  if parts2[2]=="*":
    end=sys.maxint
  else:
    end=int(parts2[2])
  
  #make sure that all the combined binary files have profiles made
  make_profiles.make_profiles(options.keep,fileName)
  
  #get and sort files
  extension="_pro"+".txt"
  filesExistProfiles=glob.glob(baseFileName+"*"+extension)
  files=[]
  for file in filesExistProfiles:
    intOfFile=int(file[len(baseFileName):len(file)-len(extension)])
    if intOfFile>=start and intOfFile<end:
      files.append(file)
  
  files.sort()
  if len(files)==0:
    print __name__+":"+main.__name__+": no files found matching "+baseFileName+"*"+extension
    return False
    
  #get surface velocity for all files in range, and find max/min of value to plotting
  L_LogTeff5mLogTeff45=[]
  L_LogTeff5mLogTeff5_t0=[]
  L_LogTeff45m3ZonesIn=[]
  fileTimes=[]
  nCount=0
  max=-1e300
  min=1e300
  L_at_LogTeff50=None
  U0_3zonesin=[]
  for file in files:
    fileData=datafile.DataFile()
    fileData.readFile(file)
    
    L_at_LogTeff5=None
    L_at_LogTeff45=None
    L_at_3ZonesIn=None
    for i in range(len(fileData.fColumnValues)):
      if fileData.fColumnValues[i][25]!=None:
        
        #find L at Log(Teff)=5
        logTTest=math.log10(fileData.fColumnValues[i][25])
        if logTTest<5.0 and L_at_LogTeff5==None:
          L_at_LogTeff5=fileData.fColumnValues[i][31]
          if L_at_LogTeff50==None:
            L_at_LogTeff50=fileData.fColumnValues[i][31]
          
        #find L at Log(Teff)=4.5
        if logTTest<4.5 and L_at_LogTeff45==None:
          L_at_LogTeff45=fileData.fColumnValues[i][31]
    
    L_at_3ZonesIn=fileData.fColumnValues[len(fileData.fColumnValues)-4][31]
    L_LogTeff5mLogTeff45.append(L_at_LogTeff5-L_at_LogTeff45)
    L_LogTeff45m3ZonesIn.append(L_at_LogTeff45-L_at_3ZonesIn)
    L_LogTeff5mLogTeff5_t0.append(L_at_LogTeff5-50.0)
    U0_3zonesin.append(fileData.fColumnValues[len(fileData.fColumnValues)-4][12])
    
    #get file times
    fileHeader=fileData.sHeader.split()
    fileTimes.append(float(fileHeader[1]))
  
  #make plot
  fig=plt.figure(figsize=(13,8))
  if True:
    ax1 = plt.subplot2grid((3,3), (0,0),colspan=3,rowspan=2)
    ax2 = plt.subplot2grid((3,3), (2,0),colspan=3)
  else:
    ax1 = plt.subplot2grid((3,3), (0,0),colspan=3,rowspan=3)
    #ax2 = plt.subplot2grid((3,3), (2,0),colspan=3)
  
  #turn on grid?
  if not options.noGrid:
    ax1.grid()
    ax2.grid()
    
  #creat plot string, (format of plot)
  plotString='-'
  if options.noLines:
    plotString='o'
  if options.points:
    plotString=plotString+'o'
  
  #set plot limits
  if options.ymin!=None:
    min=options.ymin
  if options.ymax!=None:
    max=options.ymax
  print "min=",min," max=",max
  
  #make plot
  ax1.plot(fileTimes,L_LogTeff5mLogTeff45,plotString,label="L(LogTeff=5)-L(LogTeff=4.5)",markersize=6)
  ax1.plot(fileTimes,L_LogTeff45m3ZonesIn,plotString,label="L(LogTeff=45)-L(3 zones in)",markersize=4)
  ax1.plot(fileTimes,L_LogTeff5mLogTeff5_t0,plotString,label="L(LogTeff=5)-L(LogTeff=5)_t=0",markersize=2)
  ax1.legend()
  #plt.ylim(options.ymin,options.ymax)
  ax2.set_xlabel("t [s]")
  ax1.set_ylabel("L [L_sun]")
  ax2.set_ylabel("U0 3 zones in")
  ax2.plot(fileTimes,U0_3zonesin,'-')
  
  if options.show:
    plt.show()
  else:
    sOutFileName=options.outputFile+"."+options.format
    print __name__+":"+main.__name__+": creating plot \""+sOutFileName+"\" from file "+file+"..."
    fig.savefig(sOutFileName,format=options.format,transparent=False)#save to file
  
  ax1.cla()#clear plot
  
if __name__ == "__main__":
  main()