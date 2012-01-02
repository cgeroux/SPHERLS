#!/usr/bin/env python

import datafile
import optparse as op
import make_profiles
import glob
import math
import numpy as np
import sys
import disect_filename

def parseOptions():
  #setup command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END]"
    ,version="%prog 1.0",description="Creates a plot of maximum variance of temperature and theta"
    +" velocity along side the radial grid velocity as a function of time. The input files should "
    +"be created using \"make_profiles.py\" or \"SPHERLSanal\". The resulting plot will be"
    +" saved as max_variance.FMT where. OUTPUTFILE can be specified with"
    +" the \"-o\" option and defaults to \"tmp\".")
  parser.add_option("-o","--outputFile",dest="outputFile",default="out"
    ,help="Specifies the OUTPUTFILE. [default: %default]"
    ,metavar="OUTPUTFILE", type="string")
  parser.add_option("-f","--format",dest="format",default="png", type="choice",
    help="Sets the FMT of the OUTPUTFILE, availble formats are 'png', 'pdf', 'ps', 'eps', and 'svg'."
    +"[default: %default]", metavar="FMT",choices=('png','pdf','ps','eps','svg'))
  parser.add_option("-s","--show",dest="show",action="store_true",default=False
    ,help="Display plot to x11 display rather than saving to a file.")
  parser.add_option("-t","--title",dest="title",default="Title"
    ,help="Sets the title of the plot [deafult: %default]",metavar="VARIABLE",type="string")
  parser.add_option("-k","--keep",action="store_true",dest="keep"
    ,help="Keeps distributed binary files [default].",default=True)
  parser.add_option("-r","--remove",action="store_false",dest="keep"
    ,help="Removes distributed binary files")
  parser.add_option("-m","--remake",action="store_true",dest="remake"
    ,help="Will remake profiles if they already exist. [not default].",default=False)
  parser.add_option("-p","--points",action="store_true",dest="points",help="If set, will use points when"
    +" plotting in addition to lines [default: %default]",default=False)
  parser.add_option("--no-grid",action="store_true",dest="noGrid"
    ,help="Turns off grid [default].",default=False)
  parser.add_option("--no-lines",action="store_true",dest="noLines",help="If set, will not use "
    +"lines when plotting, and only points [default: %default]",default=False)
  parser.add_option("--unscale-vel",action="store_true",dest="unscaleVel",help="If set, will scale "
    +"velocity to cm/s instead of km/s [default: %default]",default=False)
  parser.add_option("--xmin",type="float",dest="xMin",help="Sets the minimum x value of the"
    +" x-axis. If not set the minimum of the data range is used. [default: %default]"
    ,default=None)
  parser.add_option("--xmax",type="float",dest="xMax",help="Sets the maximum x value of the"
    +" x-axis. If not set the maximum of the data range is used. [default:%default]"
    ,default=None)
  parser.add_option("--period",type="float",dest="period",help="Sets the period to use for "
    +"calculating the phase of plots. If not set it plots time instead of phase [default:%default]"
    ,default=None)
    
  #parse command line options
  return parser.parse_args()
  
def main():

  #parse command line options
  (options,args)=parseOptions()
  
  import matplotlib
  if not options.show:
    matplotlib.use('Agg')#needed for saving figures
  
  import matplotlib.pyplot as plt
  
  #need a file name
  if len(args)!=1:
    print "need one and only one argument"
  
  fileName=args[0]
  
  #get base file name start and end indices of iteration range
  [start,end,baseFileName]=disect_filename.disectFileName(fileName)
  
  #make sure that all the combined binary files in range have profiles
  make_profiles.make_profiles(options.keep,fileName,options.remake,False)
  
  #get and sort files
  extension="_pro.txt"
  filesExistSlices=glob.glob(baseFileName+"*"+extension)
  files=[]
  for file in filesExistSlices:
    intOfFile=int(file[len(baseFileName):len(file)-len(extension)])
    if intOfFile>=start and intOfFile<end:
      files.append(file)
  
  files.sort()
  
  if len(files)==0:
    print "no files matching \""+baseFileName+"\" in the range specified"
    return False
  
  #lists to hold data
  times=[]
  phase=[]
  vThetaMaxs=[]
  vU_U0Maxs=[]
  maxTempVariations=[]
  surfaceGridVelocities=[]
  
  radialProfile=datafile.DataFile()
  for file in files:
    
    print "reading in profile \""+file+"\" ..."
    
    #read in profile
    radialProfile.readFile(file)
    
    #get time
    parts=radialProfile.sHeader.split()
    times.append(float(parts[1]))
    if options.period!=None:
      phase.append((float(parts[1])-times[0])/options.period)
    
    #get maximum theta velocity of entire model
    vThetaMaxTemp=math.fabs(radialProfile.fColumnValues[0][22])#used first value 
    for i in range(len(radialProfile.fColumnValues)-1):
      vThetaMaxTemp1=max(math.fabs(radialProfile.fColumnValues[i][22])
        ,math.fabs(radialProfile.fColumnValues[i][23]))
      vThetaMaxTemp=max(vThetaMaxTemp,vThetaMaxTemp1)
    if options.unscaleVel:
      vThetaMaxs.append(vThetaMaxTemp)
    else:
      vThetaMaxs.append(vThetaMaxTemp/1.0e5)
      
    
    #get maximum U-U0 velocity of entire model
    vU_U0MaxsTemp=math.fabs(radialProfile.fColumnValues[0][14]-radialProfile.fColumnValues[0][20])#used first value 
    for i in range(len(radialProfile.fColumnValues)-1):
      vU_U0MaxsTemp1=max(math.fabs(radialProfile.fColumnValues[i][14]-radialProfile.fColumnValues[i][20])
        ,math.fabs(radialProfile.fColumnValues[i][17]-radialProfile.fColumnValues[i][20]))
      vU_U0MaxsTemp=max(vU_U0MaxsTemp,vU_U0MaxsTemp1)
    if options.unscaleVel:
      vU_U0Maxs.append(vU_U0MaxsTemp)
    else:
      vU_U0Maxs.append(vU_U0MaxsTemp/1.0e5)
    
    #get maximum temperature vairation of entire model
    tempVariance=(radialProfile.fColumnValues[0][46]-radialProfile.fColumnValues[0][49]) 
      #/radialProfile.fColumnValues[0][25]
    for i in range(len(radialProfile.fColumnValues)-1):

      tempVariance2=(radialProfile.fColumnValues[i][46]-radialProfile.fColumnValues[i][49])
        #/radialProfile.fColumnValues[i][25]
      
      tempVariance=max(tempVariance,tempVariance2)
    maxTempVariations.append(tempVariance)
    
    #get surface grid velocity
    if options.unscaleVel:
      surfaceGridVelocities.append(radialProfile.fColumnValues[len(radialProfile.fColumnValues)-4][20])
    else:
      surfaceGridVelocities.append(radialProfile.fColumnValues[len(radialProfile.fColumnValues)-4][20]
        /1.0e5)
  
  '''
  #write to a file
  f=open('variation.txt','w')
  f.write("time[s] MaxTheta[cm/s] MaxRelHorTempVar U0[cm/s]\n")
  for i in range(len(times)):
    line= str(times[i])+" "+str(vThetaMaxs[i])+" "+str(maxTempVariations[i])+" "+str(surfaceGridVelocities[i])+"\n"
    f.write(line)
  f.close()'''
  
  #plt.subplots_adjust(hspace=0.001)
  #fig=plt.figure(figsize=(6,8))
  #ax=plt.subplot(1,1,1)
  #ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=6,steps=None,trim=True,integer=False,symmetric=False,prune='both'))
  #ax.yaxis.major.formatter.set_powerlimits((0,0))
  #y_max=max(watchZones[1][index])
  #y_min=min(watchZones[1][index])
  #y_max_abs=max([abs(y_max),abs(y_min)])
  #exponent=int(math.log10(y_max_abs))
  #if exponent>5:
  #  exponent=5
  #scale=10.0**float(-1.0*exponent)
  #y1=list(y*scale for y in watchZones[1][index])
  #line1=plt.plot(phase1,y1,'-',label="period "+str(period1))
  #line2=plt.plot(phase2,y1,'--',label="period "+str(period2))
  #line3=plt.plot(phase3,y1,':',label="period "+str(period3))
  #plt.title(options.title)
  #plt.annotate(shells[1],(0.88,0.2),xycoords='axes fraction',va="top",ha="center")
  #plt.setp(ax.get_xticklabels(), visible=False)
  #plt.ylabel(" 1e"+str(exponent))
  
  print "generating plot ..."
  plotString='-'
  if options.noLines:
    plotString='o'
  if options.points and not options.noLines:
    plotString=plotString+'o'
  pointSize=3
  fig=plt.figure(figsize=(13,8))
  fig.subplots_adjust(hspace=0.0000)
  ax1=plt.subplot(4,1,1)
  if not options.noGrid:
    ax1.grid()
  if options.period!=None:
    plt.plot(phase,surfaceGridVelocities,plotString,markersize=pointSize)
  else:
    plt.plot(times,surfaceGridVelocities,plotString,markersize=pointSize)
  if options.unscaleVel:
    ax1.set_ylabel("u_0 [cm/s]")
  else:
    ax1.set_ylabel("u_0 [km/s]")
  ax2=plt.subplot(4,1,2)
  if not options.noGrid:
    ax2.grid()
  if options.period!=None:
    plt.plot(phase,vThetaMaxs,plotString,markersize=pointSize)
  else:
    plt.plot(times,vThetaMaxs,plotString,markersize=pointSize)
  
  if options.unscaleVel:
    ax2.set_ylabel("|v| max [cm/s]")
  else:
    ax2.set_ylabel("|v| max [km/s]")
  ax3=plt.subplot(4,1,3)
  if not options.noGrid:
    ax3.grid()
  if options.period!=None:
    plt.plot(phase,maxTempVariations,plotString,markersize=pointSize)
  else:
    plt.plot(times,maxTempVariations,plotString,markersize=pointSize)
  ax3.set_ylabel("|Tmax-Tmin| [K]")
  ax4=plt.subplot(4,1,4)
  if not options.noGrid:
    ax4.grid()
  if options.period!=None:
    plt.plot(phase,vU_U0Maxs,plotString,markersize=pointSize)
  else:
    plt.plot(times,vU_U0Maxs,plotString,markersize=pointSize)
  if options.unscaleVel:
    ax4.set_ylabel("|u-u_0| max [cm/s]")
  else:
    ax4.set_ylabel("|u-u_0| max [km/s]")
  if options.period!=None:
    ax4.set_xlabel("phase")
  else:
    ax4.set_xlabel("t [s]")
  ax1.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(prune='both'))
  ax2.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(prune='both'))
  ax3.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(prune='both'))
  ax4.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(prune='both'))
  xticklables=ax1.get_xticklabels()+ax2.get_xticklabels()+ax3.get_xticklabels()
  plt.setp(xticklables,visible=False)
  ax1.set_title(options.title)
  ax1.set_xlim(options.xMin,options.xMax)
  ax2.set_xlim(options.xMin,options.xMax)
  ax3.set_xlim(options.xMin,options.xMax)
  ax4.set_xlim(options.xMin,options.xMax)
  if options.show:
    plt.show()
  else:
    print __name__+":"+main.__name__+": saving figure to file \""+options.outputFile+"."+options.format
    fig.savefig(options.outputFile+"."+options.format,format=options.format,transparent=False)#save to file
  
if __name__ == "__main__":
  main()
