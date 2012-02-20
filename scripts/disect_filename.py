#!/usr/bin/env python
import sys
def disectFileName(fileName):

  #get base file name
  parts0=fileName.split('[',1)
  if len(parts0)<2:# maybe a single file?
    parts3=fileName.rsplit('_t',1)
    if parts3[1].isdigit():
      start=int(parts3[1])
      end=start+1
      baseFileName=parts3[0]+'_t'
    else:
      print disectFileName.__name__+": unrecognized file type "+fileName\
        +" expecting something with an _tXXXXXXXX suffix, where X's denote digits."
      quit()
  
  else:# a range of files
    baseFileName=parts0[0]
    
    #get interation range
    parts1=parts0[1].split(']')
    parts2=parts1[0].split('-')
    start=int(parts2[0])
    if parts2[1]=="*":
      end=sys.maxint
    else:
      end=int(parts2[1])
  
  return [start,end,baseFileName]
