#!/usr/bin/env python

#import matplotlib
import optparse as op
import os
import math
import numpy as np
  
def main():
  
  #parse command line options
  (options,args)=parseOptions()
  import matplotlib
  if not options.show:
    matplotlib.use('Agg')#needed for saving figures
  
  import matplotlib.pyplot as plt
  
  ############################################
  #velocity as a function of theta
  '''#read in input file
  print "reading in ",args[0]," ..."
  slice=read2DSlice(args[0])
  
  x=[]
  y=[]
  y2=[]
  count=0
  for i in range(len(slice['coord1'])-1):
    for j in range(len(slice['coord2'])-1):
      if (j==5):
        #y.append(slice['vars'][0][count]-slice['vars'][1][count])
        y.append(slice['vars'][0][count])
        #y2.append(slice['vars'][0][count])
        #x.append((slice['coord1'][i+1]+slice['coord1'][i])*0.5)
        x.append(i)
      count=count+1
  
  #plot the figure
  fig=plt.figure()
  ax=fig.add_subplot(111)
  l1=ax.plot(x,y,'x:b',label="u")
  ax.set_xlabel("r [cm]")
  ax.set_ylabel("u [cm/s]")
  ax.set_title("u, v, velocity corelation at t="+str(round(slice['time'],2))+"[s]")
  #ax2=ax.twinx()
  #l2=ax2.plot(x,y2,'x-r',label="v")
  #ax2.set_ylabel("v [cm/s]")
  #ax.legend((l1,l2),('u','v'),'upper left')
  plt.show()#x11'''
  
  #####################################
  #phase plots
  watchZones=[]
  for file in args:
    print "reading in ",file," ..."
    watchZones.append(readWatchzone(file))
  
  #calculate the phase from watchZone
  print "calculating phase ..."
  phase=calcPhase(watchZones[0])
  period1=200
  period2=203
  period3=250
  phase1=list(phi-float(period1) for phi in phase)
  phase2=list(phi-float(period2) for phi in phase)
  phase3=list(phi-float(period3) for phi in phase)
  
  print "making plot ..."
  xLabelSt="Phase"
  yLabelSt="None Set"
  if options.variable=="U0":
    index=4
    yLabelSt="u_0"
    yUnits="[cm/s]"
  if options.variable=="R":
    index=13
    yLabelSt="r"
    yUnits="[cm]"
  if options.variable=="DA":
    index=16
    yLabelSt="H. A. D."
    yUnits="[g/cm^3]"
  if options.variable=="D":
    index=15
    yLabelSt="D"
    yUnits="[g/cm^3]"
  if options.variable=="E":
    index=17
    yLabelSt="E"
    yUnits="[ergs]"
  if options.variable=="P":
    index=18
    yLabelSt="P"
    yUnits="[dyn/cm^2]"
  if options.variable=="T":
    index=19
    yLabelSt="T"
    yUnits="[K]"
  
  #move plots close together
  print "i=1"
  shells=["Shell 77","Shell 125","Shell 85","Shell 65","Shell 45","Shell 25","Shell 5"]
  plt.subplots_adjust(hspace=0.001)
  fig=plt.figure(figsize=(6,8))
  ax=plt.subplot(len(watchZones)-1,1,1)
  ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=6,steps=None,trim=True,integer=False,symmetric=False,prune='both'))
  ax.yaxis.major.formatter.set_powerlimits((0,0))
  y_max=max(watchZones[1][index])
  y_min=min(watchZones[1][index])
  y_max_abs=max([abs(y_max),abs(y_min)])
  exponent=int(math.log10(y_max_abs))
  if exponent>5:
    exponent=5
  scale=10.0**float(-1.0*exponent)
  y1=list(y*scale for y in watchZones[1][index])
  line1=plt.plot(phase1,y1,'-',label="period "+str(period1))
  line2=plt.plot(phase2,y1,'--',label="period "+str(period2))
  line3=plt.plot(phase3,y1,':',label="period "+str(period3))
  #plt.title(options.title)
  plt.annotate(shells[1],(0.88,0.2),xycoords='axes fraction',va="top",ha="center")
  plt.setp(ax.get_xticklabels(), visible=False)
  plt.ylabel(" 1e"+str(exponent))
  for i in range(2,len(watchZones)):
    print "i=",i
    y_max=max(watchZones[i][index])
    y_min=min(watchZones[i][index])
    y_max_abs=max([abs(y_max),abs(y_min)])
    exponent=int(math.log10(y_max_abs))
    if exponent>5:
      exponent=5
    scale=10.0**float(-1.0*exponent)
    y1=list(y*scale for y in watchZones[i][index])
    ax1=plt.subplot(len(watchZones)-1,1,i,sharex=ax)
    ax1.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=6,steps=None,trim=True,integer=False,symmetric=False,prune='both'))
    ax1.yaxis.major.formatter.set_powerlimits((0,0))
    plt.plot(phase1,y1,'-',label="period "+str(period1))
    plt.plot(phase2,y1,'--',label="period "+str(period2))
    plt.plot(phase3,y1,':',label="period "+str(period3))
    plt.annotate(shells[i],(0.88,0.2),xycoords='axes fraction',va="top",ha="center")
    plt.ylabel(" 1e"+str(exponent))
    if i !=len(watchZones)-1:
      plt.setp(ax1.get_xticklabels(), visible=False)
  plt.figlegend( (line1,line2,line3), ("period "+str(period1),"period "+str(period2),"period "+str(period3)),'upper right')
  plt.xlabel(xLabelSt)
  ax.set_xlim([0.0,2.0])
  
  if options.show:
    plt.show()
  else :
    fig.savefig(options.outputFile+"."+options.format,format=options.format,transparent=False)#save to file
  
def read2DSlice(fileName):
  
  slice={}
  
  #open file for reading
  if not os.access(fileName,os.F_OK|os.R_OK):
    print "error opening ",fileName," for reading"
    return
  f=open(fileName,'r')
  
  #read in time
  line=f.readline()
  row=[]
  rowTemp=line.split()
  for num in rowTemp:
    row.append(float(num))
  time=row[0]
  
  slice['time']=time
  
  #skip first line
  line=f.readline()
  
  #read in r
  line=f.readline()
  row=[]
  rowTemp=line.split()
  for num in rowTemp:
    row.append(float(num))
  r=np.array(row)
  slice['coord1']=r
  
  #read in theta
  line=f.readline()
  row=[]
  rowTemp=line.split()
  for num in rowTemp:
    row.append(float(num))
  theta=np.array(row)
  slice['coord2']=theta
  
  #read in first line and add a list for each variable
  vars=[]
  line=f.readline()
  rowTemp=line.split()
  i=0
  for num in rowTemp:
    vars.append([])
    vars[i].append(float(num))
    i=i+1
  
  #read in the rest of the file
  for line in f:
    i=0
    rowTemp=line.split()
    for num in rowTemp:
      vars[i].append(float(num))
      i=i+1
  slice['vars']=vars
  return slice
  
def readWatchzone(fileName):
  
  #open file for reading
  if not os.access(fileName,os.F_OK|os.R_OK):
    print "error opening ",fileName," for reading"
    return
  f=open(fileName,'r')
  
  #skip first two lines
  line=f.readline()
  line=f.readline()
  
  watchZone=[]
  
  #add columns, and add first number to column
  line=f.readline()
  lineTemp=line.split()
  colNum=0
  for num in lineTemp:
    watchZone.append([])
    if num=="-":
      num=0
    watchZone[colNum].append(float(num))
    colNum=colNum+1
  
  #fill up columns
  i=1
  for line in f:
    lineTemp=line.split()
    colNum=0
    #if(i>1022674 and i<1022682):
    if(lineTemp[0]!=str(i)):
      print "line i=",i," time step index, is not i+1 from last time step index"
      print line
      quit()
    for num in lineTemp:
      if num=="-":
        num=0
      watchZone[colNum].append(float(num))
      colNum=colNum+1
    i=i+1
  return watchZone
  
def calcPhase(watchZone):
  phase=[]
  i=0
  indexGridVelocity=4
  indexTime=1
  #find ia1
  ia1=None
  while ia1==None:
    if (watchZone[indexGridVelocity][i+1]>0.0 and watchZone[indexGridVelocity][i]<0.0):
      ia1=i+1
    i=i+1
  
  #find ib1
  ib1=None
  while ib1==None:
    if (watchZone[indexGridVelocity][i+1]<0.0 and watchZone[indexGridVelocity][i]>0.0):
      ib1=i+1
    i=i+1
    
  #find max1
  imax1=ia1
  for j in range(ia1,ib1):
    if (watchZone[indexGridVelocity][j]>watchZone[indexGridVelocity][imax1]):
      imax1=j
      
  #calculate tmax1
  p1=(watchZone[indexTime][imax1-1],watchZone[indexGridVelocity][imax1-1])
  p2=(watchZone[indexTime][imax1  ],watchZone[indexGridVelocity][imax1  ])
  p3=(watchZone[indexTime][imax1+1],watchZone[indexGridVelocity][imax1+1])
  tmax1=calPMax(p1,p2,p3)
  
  phiCalIntStart=0
  bContinue=True
  nCount=0
  while bContinue:
    
    #find ia2
    ia2=None
    while ia2==None and bContinue:
      if (watchZone[indexGridVelocity][i+1]>0.0 and watchZone[indexGridVelocity][i]<0.0):
        ia2=i+1
      i=i+1
      if(i>=len(watchZone[indexGridVelocity])-1):#past the end of the file
        bContinue=False;
    
    #find ib2
    ib2=None
    while ib2==None and bContinue:
      if (watchZone[indexGridVelocity][i+1]<0.0 and watchZone[indexGridVelocity][i]>0.0):
        ib2=i+1
      i=i+1
      if(i>=len(watchZone[indexGridVelocity])-1):#past the end of the file
        bContinue=False;
    
    if bContinue:
      #find max2
      imax2=ia2
      for j in range(ia2,ib2):
        if (watchZone[indexGridVelocity][j]>watchZone[indexGridVelocity][imax2]):
          imax2=j
      
      #calculate tmax2
      p1=(watchZone[indexTime][imax2-1],watchZone[indexGridVelocity][imax2-1])
      p2=(watchZone[indexTime][imax2  ],watchZone[indexGridVelocity][imax2  ])
      p3=(watchZone[indexTime][imax2+1],watchZone[indexGridVelocity][imax2+1])
      tmax2=calPMax(p1,p2,p3)
      
      #calculate slope in linear relation between phase, and time
      m=1.0/(tmax2-tmax1)
      
      #calculate t_0
      t_0=tmax1 #phase zero is defined to be at the first maximum
      
      #calculate phase for region
      for k in range(phiCalIntStart,imax2+1):
        phase.append(m*(watchZone[indexTime][k]-t_0)+1.0*float(nCount))
      
      #set starting values for next loop
      nCount=nCount+1
      phiCalIntStart=imax2+1
      tmax1=tmax2
    else:
      for k in range(phiCalIntStart,len(watchZone[1])):
        phase.append(m*(watchZone[indexTime][k]-t_0)+1.0*float(nCount))
  return phase
def calPMax(p1,p2,p3):
  """This function returns the point at which the derivative of the porabola defined by points 
    p1,p2, and p3 is zero.
    
    If the second derivative is positive then it is a maximum, if the second derivative is negative
    it is a minimum.
  """
  a=((p1[1]-p2[1])*(p2[0]-p3[0])-(p2[1]-p3[1])*(p1[0]
    -p2[0]))/((p1[0]*p1[0]-p2[0]*p2[0])*(p2[0]-p3[0])
    -(p2[0]*p2[0]-p3[0]*p3[0])*(p1[0]-p2[0]))
  b=(p1[1]-p2[1]-a*(p1[0]*p1[0]-p2[0]*p2[0]))/(p1[0]-p2[0])
  return -1.0*b/(2.0*a)
  
def parseOptions():
  #setup command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] INPUT1 INPUT2",version="1.0"
    ,description="Creastes plots of phased variables. The phase is determined using INPUT1 "
    +"(a watchzone file) and applied to INPUT2 (a second watchZone file. The watchZone file used "
    +"for phasing should be near the node of the first overtone so that it will be assured that "
    +"the phase is of the fundamental mode.")
  parser.add_option("-o","--outputFile",dest="outputFile",default="out"
    ,help="Specifies the OUTPUTFILE. [default: %default]"
    ,metavar="OUTPUTFILE", type="string")
  parser.add_option("-f","--format",dest="format",default="png", type="choice",
    help="Sets the FMT of the OUTPUTFILE, availble formats are 'png', 'pdf', 'ps', 'eps', and 'svg'."
    +"[default: %default]", metavar="FMT",choices=('png','pdf','ps','eps','svg'))
  parser.add_option("-s","--show",dest="show",action="store_true",default=False
    ,help="Display plot to x11 display rather than saving to a file.")
  parser.add_option("-v","--variable",dest="variable",default="U0"
    ,help="Sets the variable to be plotted. Possible values of VARIABLE are 'U0' "
    +"(Radial grid velocity), 'R' (radius), 'DA' (average density), 'D' (density)"
    +", 'E' (energy), 'P' (pressure), 'T' (temperature) [deafult: %default]",metavar="VARIABLE",type="choice"
    ,choices=('U0','R','DA','D','E','P','T'))
  parser.add_option("-t","--title",dest="title",default="Title"
    ,help="Sets the title of the plot [deafult: %default]",metavar="VARIABLE",type="string")
  #parse command line options
  return parser.parse_args()
if __name__ == "__main__":
  main()
