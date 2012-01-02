#!/usr/bin/env python
import combine_bins
import optparse as op
import time

def main():
  
  #make command line parser
  parser=op.OptionParser(usage="Usage: %prog [options] BASEFILENAME",version="%prog 1.0"
    ,description="Used to combine binary files while a program is running. It does this by"
    +" continuously checking to see if new files exists and if they do combining them. It will"
    +" continue doing this until \"ctrl+c\" is pressed.")
  
  parser.add_option("-k","--keep",action="store_true",dest="keep"
    ,help="Keeps distributed binary files [default].",default=True)
  parser.add_option("-r","--remove",action="store_false",dest="keep"
    ,help="Removes distributed binary files")
  parser.add_option("-w","--wait",type="float",dest="wait",default=2.0,help="Sets the time to wait "
    +"before attempting to look for more files in seconds. [default: %default s]")
    
  #parse command line options
  (options,args)=parser.parse_args()
    
  #should have at least one argument
  if len(args)!=1 :
    parser.error(" need one and one argument only.")
    
  #create profile files, and save list of files
  bContinue=True
  while(bContinue):
    print "combining all current files ..."
    bContinue=combine_bins.combine_bin_files(options.keep,args[0]+"_t[0-*]")
    while(not bContinue):#keep looking for files
      print "waiting for more files. press \"ctrl+c\" to stop ..."
      bContinue=combine_bins.combine_bin_files(options.keep,args[0]+"_t[0-*]")
      time.sleep(options.wait)
  
if __name__=="__main__":
  main()