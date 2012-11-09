#!/usr/bin/env python
import os
import getopt
import sys
import glob
import optparse as op
import combine_bins
import paths
import disect_filename

def main():
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME[START-END]",version="%prog 1.0"
    , description="Slices the 3D dump files of \"SPHERLS\" into"
    +" 2D slices. This is a front end for the \"SPHERLSanal\" program allowing multiple files "
    +"to be sliced at once from START upto but not including END. Ex. there are a number of files "
    +"in the current directory from run1_t00000000, to run1_t00012300, one could create 2D slice "
    +"of all these files with the command %prog run1_t[0-12300] or %prog run1_t[0-*]. The \"*\" "
    +"wild character is used to include all files upto the maximum integer starting from START.")
  #make plane group options
  plane_group=op.OptionGroup(parser,"Plane Indices","Only one of these three options should be "
    +"specified at once. [default: -k 2]")
  plane_group.add_option("-i", "--radial-index",dest="radialIndex",
    help="Specifies the RINDEX in the r-direction at which a theta-phi slice should be taken.",
    metavar="RINDEX",type="int",default=-1)
  plane_group.add_option("-j", "--theta-index",dest="thetaIndex",
    help="Specifies the THETAINDEX in the theta-direction at which an r-phi slice should be taken.",
    metavar="THETAINDEX",type="int",default=-1)
  plane_group.add_option("-k", "--phi-index",dest="phiIndex",
    help="Specifies the PHIINDEX in the phi-direction at which an r-theta slice should be taken.",
    metavar="PHIINDEX",type="int",default=-1)
  parser.add_option_group(plane_group)
  parser.add_option("--keep",action="store_true",dest="keep"
    ,help="Keeps distributed binary files [default].",default=True)
  parser.add_option("--remove",action="store_false",dest="keep"
    ,help="Removes distributed binary files")
  parser.add_option("-m","--remake",action="store_true",dest="remake"
    ,help="Remake 2D slice even if it already exists.")
  
  #parse arguments
  (options,args)=parser.parse_args()
  
  #make sure only one plane is spceified
  bPlaneSet=False
  bMoreThanOnePlaneSet=False
  if options.radialIndex!=-1 :
    bPlaneSet=True
  if options.thetaIndex!=-1 :
    if bPlaneSet :
      bMoreThanOnePalneSet=True
    bPlaneSet=True
  if options.phiIndex!=-1 :
    if bPlaneSet:
      bMoreThanOnePlaneSet=True
    bPlaneSet=True
      
  if bMoreThanOnePlaneSet :
    parser.error("Only specify one of -i, -j, and -k")
  if not bPlaneSet :#if no plane is specified use default of "-k 2"
    options.phiIndex=2
  
  #need a file name
  if len(args)!=1:
    parser.error("need one and only one argument")
  
  planeID=None
  planeIndex=None
  planeID=None
  planeIndex=None
  if options.radialIndex!=-1:
    planeID=1
    planeIndex=options.radialIndex
  if options.thetaIndex!=-1:
    planeID=2
    planeIndex=options.thetaIndex
  if options.phiIndex!=-1:
    planeID=0
    planeIndex=options.phiIndex
  
  #should have at least one argument
  if len(args) < 1:
    print "no file name"
    usage()
    quit()
  if planeID==None:
    print "no plane specified"
    usage()
    quit()
  if planeIndex==None:
    print "no plane index specified"
    usage()
    quit()
  
  #create profile files, and save list of files
  make_2DSlices(options.keep,args[0],planeID,planeIndex,options.remake)
def make_2DSlices(keep,fileName,nPlane,nPlaneIndex,remake):
  
  #get base file name
  [start,end,baseFileName]=disect_filename.disectFileName(fileName)
  
  #make sure that all the distributed binary files in this range have been converted to combined binary files
  combine_bins.combine_bin_files(keep,fileName,False)
  
  #check for combined binary files in interation range starting with baseFileName
  filesExistCombBin=glob.glob(baseFileName+"[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]")
  files=[]
  for file in filesExistCombBin:
    intOfFile=int(file[len(baseFileName):len(file)])
    if intOfFile>=start and intOfFile<end:
      files.append(file)
  
  extension="_2D"
  if nPlane==0:
    extension=extension+"k="+str(nPlaneIndex)+".txt"
  if nPlane==1:
    extension=extension+"i="+str(nPlaneIndex)+".txt"
  if nPlane==2:
    extension=extension+"j="+str(nPlaneIndex)+".txt"
  count=1
  for file in files:
    #if this particular 2D slice doesn't already exsist for this binary file create it
    if not (os.path.exists(file+extension)) or remake:
      
      #make 2D slice
      print __name__+":"+make_2DSlices.__name__+": creating 2D slice from \""+file+"\" "\
        +str(count)+"/"+str(len(files))+" ..."
      success=os.system(paths.SPHERLSanalPath+' -s cb '+str(nPlane)+' '+str(nPlaneIndex)+' '+file)
      if success==0:
        pass
      else :
        
        #say there was an error and return
        print __name__+":"+make_2DSlices.__name__+": error making 2D slice "+file+extension
        return False
    else:
      print __name__+":"+make_2DSlices.__name__+": 2D slice \""+file+extension+"\" already exists"
    count+=1
  return True
if __name__ == "__main__":
  main()
  