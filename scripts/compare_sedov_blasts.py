#!/usr/bin/env python

import os
import numpy
import math
import datafile
from scipy.interpolate import interp1d


def main():
  fColumnAnYScale=[1.0e2,1e8,1.0,1e16]  
  nColumnsAnY=[1,2,3]
  nColumnsPrY=[8,5,24]
  sVariables=["u","d","p"]
  radialProfile=dataFile.DataFile()
  analyticSolution=dataFile.DataFile()
  profileFiles=[
    "./30m/run1_t00020886_pro_t100.txt",
    "./30m/run1_t00020539_pro_t095.txt",
    "./30m/run1_t00020177_pro_t090.txt",
    "./30m/run1_t00019796_pro_t085.txt",
    "./30m/run1_t00019393_pro_t080.txt",
    "./30m/run1_t00018974_pro_t075.txt",
    "./30m/run1_t00018528_pro_t070.txt",
    "./30m/run1_t00018051_pro_t065.txt",
    "./30m/run1_t00017550_pro_t060.txt",
    "./30m/run1_t00017013_pro_t055.txt",
    "./30m/run1_t00016460_pro_t050.txt",
    "./30m/run1_t00015859_pro_t045.txt",
    "./30m/run1_t00015203_pro_t040.txt",
    "./30m/run1_t00014475_pro_t035.txt",
    "./30m/run1_t00013663_pro_t030.txt",
    "./30m/run1_t00012735_pro_t025.txt",
    "./30m/run1_t00011637_pro_t020.txt",
    "./30m/run1_t00010361_pro_t015.txt",
    "./30m/run1_t00008760_pro_t010.txt",
    "./30m/run1_t00006334_pro_t005.txt"
    ]
  analyticFiles=[
    "./100.txt",
    "./095.txt",
    "./090.txt",
    "./085.txt",
    "./080.txt",
    "./075.txt",
    "./070.txt",
    "./065.txt",
    "./060.txt",
    "./055.txt",
    "./050.txt",
    "./045.txt",
    "./040.txt",
    "./035.txt",
    "./030.txt",
    "./025.txt",
    "./020.txt",
    "./015.txt",
    "./010.txt",
    "./005.txt"
    ]
  times=[
    "100",
    "095",
    "090",
    "085",
    "080",
    "075",
    "070",
    "065",
    "060",
    "055",
    "050",
    "045",
    "040",
    "035",
    "030",
    "025",
    "020",
    "015",
    "010",
    "005"
    ]
  for p in range(len(profileFiles)):
    radialProfile.readFile(profileFiles[p])
    analyticSolution.readFile(analyticFiles[p])
    for k in range(3):
      fileOutput=open("./interpolate_t"+times[p]+"_"+sVariables[k]+".txt",'w')
      analyticSolution.initInterpolate(0,'cubic')
      
      #loop over radial profile in range of analytic solution and calculate stadard deviation
      fMaxAnalytic=analyticSolution.fColumnValues[0].max()*fColumnAnYScale[0]
      fMinAnalytic=analyticSolution.fColumnValues[0].min()*fColumnAnYScale[0]
      fileOutput.write(radialProfile.sColumnNames[3]+" "+radialProfile.sColumnNames[nColumnsPrY[k]]
        +"(Comp) "+radialProfile.sColumnNames[nColumnsPrY[k]]+"(Ana)\n")
      fSum1=0.0
      fSum2=0.0
      nCount=0
      fX=[]
      fY1=[]
      fY2=[]
      fDiff=[]
      fMaxError=0.0
      for i in range(30,len(radialProfile.fColumnValues[4])):
        
        #if zone centered quantity, center the radius
        if nColumnsPrY[k]==6 or nColumnsPrY[k]==19 or nColumnsPrY[k]==25 :
          fXTemp=(radialProfile.fColumnValues[4][i]+radialProfile.fColumnValues[3][i])*0.5
        else:
          fXTemp=radialProfile.fColumnValues[4][i]
        
        #if it is in the range of the analytic solution, interpolate analytic solution to this value
        
        if fXTemp>=fMinAnalytic and fXTemp<=fMaxAnalytic: 
          fX.append(fXTemp)
          
          if sVariables[k] =="B":
            fE1=radialProfile.fColumnValues[24][i]/(radialProfile.fColumnValues[5][i]*0.6)
            fD=analyticSolution.interpf[1](fXTemp/fColumnAnYScale[0])*fColumnAnYScale[2]
            fP=analyticSolution.interpf[2](fXTemp/fColumnAnYScale[0])*fColumnAnYScale[3]
            fE2=fP/(fD*0.6)
            fY1.append(fE1)
            fY2.append(fE2)
          else:
            fY1.append(radialProfile.fColumnValues[nColumnsPrY[k]][i])
            fY2.append(analyticSolution.interpf[nColumnsAnY[k]-1](fXTemp/fColumnAnYScale[0])*fColumnAnYScale[nColumnsAnY[k]])
          fDiff.append((fY2[nCount]-fY1[nCount])/fY2[nCount])
          fDiff1=(fY2[nCount]-fY1[nCount])
          fDiff2=(fY2[nCount]-fY1[nCount])/fY2[nCount]
          fSum1=fSum1+fDiff1*fDiff1
          fSum2=fSum2+fDiff2*fDiff2
          nCount=nCount+1
          if fDiff2>fMaxError:
            fMaxError=fDiff2
      fSTD=math.sqrt(fSum1/float(nCount))
      fSTDRel=math.sqrt(fSum2/float(nCount))
      print sVariables[k]+" "+times[p]+" "+str(fSTDRel)+" "+str(fMaxError)
      
      #print out info
      for i in range(len(fX)):
        s=str(fX[i])+" "+str(fY1[i])+" "+str(fY2[i])+" "+" "+str(fDiff[i])+" "+str(fSTD)+"\n"
        fileOutput.write(s)
      fileOutput.close()
  
if __name__ == "__main__":
  main()