#this version works when all files are in CWD need to set include paths some how
#I can't get include paths to work no matter how hard I try
#my solution is to just include in the build method a copy of the needed files 
#into this directory is is clunky cause they will always be rebuilt but that
#is the only way I can get this thing to work.
from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize(
           "eos.pyx",                 # our Cython source
           sources=["./eos_tmp.cpp","./exception2.cpp"],  # additional source file(s)
           language="c++",             # generate C++ code
      )
      )
