#!/usr/bin/env python

import datafile
import optparse as op
import make_profiles
import glob
import math
import numpy as np
import sys
import disect_filename
def addParserOptions(parser):
  
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
  parser.add_option("--points",action="store_true",dest="points",help="If set, will use points when"
    +" plotting in addition to lines [default: %default]",default=False)
  parser.add_option("--no-lines",action="store_true",dest="noLines",help="If set, will not use "
    +"lines when plotting, and only points [default: %default]",default=False)
  
  make_profiles.addParserOptions(parser)
  
  #parser.add_option("-k","--keep",action="store_true",dest="keep"
  #  ,help="Keeps distributed binary files [default].",default=True)
  #parser.add_option("-r","--remove",action="store_false",dest="keep"
  #  ,help="Removes distributed binary files")
  #parser.add_option("-m","--remake",action="store_true",dest="remake"
  #  ,help="Will remake profiles if they already exist. [not default].",default=False)
  #parser.add_option("--remake-bins",action="store_true",dest="remakeBins"
  #  ,help="Will remake binaries even if they already exist. [not default].",default=False)

def parseOptions():
  #setup command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END]"
    ,version="%prog 1.0",description="Creates a plot of the surface luminosity as a function of time.")
  addParserOptions(parser)
  
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
  [start,end,baseFileName]=disect_filename.disectFileName(fileName)
  
  #make sure that all the combined binary files have profiles made
  #make_profiles.make_profiles(options.keep,fileName,options.remake,options.remakeBins)
  #create profile files, and save list of files
  failedFiles=make_profiles.make_fileSet(args[0],options)
  
  if __name__=="__main__":#keeps from redundantly reporting errors
    #report profiles failed files
    for file in failedFiles:
      print file
  
  #get and sort files
  extension="_pro"+".txt"
  filesExistProfiles=glob.glob(baseFileName+"*"+extension)
  files=[]
  for file in filesExistProfiles:
    intOfFile=int(file[len(baseFileName):len(file)-len(extension)])
    if intOfFile>=start and intOfFile<end:
      files.append(file)
  
  files.sort()
  
  #get surface velocity for all files in range, and find max/min of value to plotting
  surfaceLuminosity=[]
  fileTimes=[]
  nCount=0
  max=-1e60
  min=1e60
  nColumn=55
  nZonesFromSurface=3
  for file in files:
    
    print "reading in profile ",file," ..."
    fileData=datafile.DataFile()
    fileData.readFile(file)
    surfaceLuminosity.append(fileData.fColumnValues[len(fileData.fColumnValues)-nZonesFromSurface-1][nColumn])
    fileHeader=fileData.sHeader.split()
    fileTimes.append(float(fileHeader[1]))
    nCount=nCount+1
    
    if max<fileData.fColumnValues[len(fileData.fColumnValues)-nZonesFromSurface-1][nColumn] and fileData.fColumnValues[len(fileData.fColumnValues)-nZonesFromSurface-1][nColumn] !=None:
      max=fileData.fColumnValues[len(fileData.fColumnValues)-nZonesFromSurface-1][nColumn]
    if min>fileData.fColumnValues[len(fileData.fColumnValues)-nZonesFromSurface-1][nColumn] and fileData.fColumnValues[len(fileData.fColumnValues)-nZonesFromSurface-1][nColumn] !=None:
      min=fileData.fColumnValues[len(fileData.fColumnValues)-nZonesFromSurface-1][nColumn]

  #make plots
  nCount=0
  fig=plt.figure(figsize=(13,8))
  ax1 = plt.subplot2grid((3,3), (0,0),colspan=3,rowspan=3)
  
  if not options.noGrid:
    ax1.grid()
  plotString='-'
  if options.noLines:
    plotString='o'
  if options.points:
    plotString=plotString+'o'
  
  if options.ymin!=None:
    min=options.ymin
  if options.ymax!=None:
    max=options.ymax
  print "min=",min," max=",max
  ax1.plot(fileTimes,surfaceLuminosity,plotString)
  plt.ylim(min,max)
  ax1.set_xlabel("t [s]")
  ax1.set_ylabel("L [L_sun]")
  ax1.set_title("Ligh Curve")
  
  if options.show:
    plt.show()
  else:
    sOutFileName=options.outputFile+"."+options.format
    print __name__+":"+main.__name__+": creating plot \""+sOutFileName+"\" from file "+file+"..."
    fig.savefig(sOutFileName,format=options.format,transparent=False)#save to file
  
  ax1.cla()#clear plot
  
if __name__ == "__main__":
  main()