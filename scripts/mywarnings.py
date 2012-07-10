import warnings
import sys
def myshowwarning(message,category,filename,lineno,file=None,line=None):
  """Custom format warning messages to not show code line"""
  
  sys.stderr.write(warnings.formatwarning(message,category,filename,lineno,line=""))
  
#use a custom warning function
warnings.showwarning=myshowwarning