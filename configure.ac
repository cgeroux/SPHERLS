AC_INIT([SPHERLS], [1.0.1], [chris.m.geroux@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_CONFIG_SRCDIR([src/SPHERLS/physEquations.cpp])#safty check to make sure src directory is what it should be
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_LANG([C++])

AC_CHECK_PROG(MPICPP_CHECK,mpic++,yes)
if test x"$MPICPP_CHECK" != x"yes" ; then
    AC_MSG_ERROR([mpic++ not found, aborting! Please make sure that the bin directory where OpenMPI is installed is in your path.])
fi

# Checks for header files.
AC_CHECK_HEADERS([fenv.h memory.h stddef.h stdlib.h string.h unistd.h wchar.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memmove memset pow sqrt strcasecmp strncasecmp strstr])

#Should optional 3rd party libraries be included in the distribution
INCLUDE_OPTIONAL_LIBS_IN_DIST=no
AC_ARG_ENABLE([include-optional-libs],
  [AS_HELP_STRING([--enable-include-optional-libs],
  [Optional third party libraries will be included in distribution package created from "make dist"])],
  [INCLUDE_OPTIONAL_LIBS_IN_DIST="$enableval"],
  [])
AM_CONDITIONAL([INCLUDE_OPTIONAL_LIBS_IN_DIST],[test "$INCLUDE_OPTIONAL_LIBS_IN_DIST" = "yes"])

#Should crucial 3rd party libraries be included in the distribution
INCLUDE_CRUCIAL_LIBS_IN_DIST=yes
AC_ARG_ENABLE([include-crucial-libs],
  [AS_HELP_STRING([--disable-include-crucial-libs],
  [Crucial third party libraries will not be included in distribution package created from "make dist"])],
  [INCLUDE_CRUCIAL_LIBS_IN_DIST="$enableval"],
  [])
AM_CONDITIONAL([INCLUDE_CRUCIAL_LIBS_IN_DIST],[test "$INCLUDE_CRUCIAL_LIBS_IN_DIST" = "yes"])

#Should documents be made with "make" or "make all"
MAKE_DOCS=no
AC_ARG_ENABLE([make-docs],
  [AS_HELP_STRING([--enable-make-docs],
  [If enabled it will remake the documentation from the source files when running "make" or "make all". This requires the installation of pdflatex and is not recommended unless you are actively making changes to the source files.])],
  [MAKE_DOCS="$enableval"],
  [])
AM_CONDITIONAL([MAKE_DOCS],[test "$MAKE_DOCS" = "yes"])
AS_IF(
  [test "$MAKE_DOCS" = "yes"],
  [#
  
  #check for pdflatex
  AC_CHECK_PROG(PS2PDF_EXISTS,pdflatex,"yes","no")
  AS_IF(
  [test "$PS2PDF_EXISTS" = "no"],
  [
  AC_MSG_ERROR([creating documentation requires pdflatex, remove the option --enable-make-docs or install pdflatex and put it in your path.])
  
  ])
  ])


#################################################################
## Check for FFTW3 include and library
#################################################################
#
#check to see if user disabled fftw
FFTW_ENABLE=yes
AC_ARG_ENABLE([fftw],
  [AS_HELP_STRING([--disable-fftw],
  [Disable fftw features, such as performing frequency analysis of time varying quantities.])],
  [FFTW_ENABLE="$enableval"],
  [])

#
#check to make sure we have a library and include path, if we are using fftw
AS_IF(
  [test "$FFTW_ENABLE" = "yes"],
  [#
  AC_CHECK_HEADERS(fftw3.h,[],[
    AC_MSG_ERROR([
---------------------------------------------------------------------
  Unable to find the fftw3.h header file.
 
  Try adding a path to the header file to the CPPFLAGS environment
  variable e.g. export CPPFLAGS="-I<include dir> \${CPPFLAGS}", or 
  using the --disable-fftw option to disable FFTW functionality if
  it isn't needed.
---------------------------------------------------------------------
    ])
  ])
  for func in fftw_execute fftw_plan_dft_1d fftw_destroy_plan fftw_free; do
    AC_SEARCH_LIBS([${func}], fftw3,[],[
      AC_MSG_ERROR([
---------------------------------------------------------------------
  Unable to find an FFTW library containing the ${func} function.

  If you know the path to the library try adding it to the LDFLAGS
  environment variable. e.g. export LDFLAGS="-L<lib dir> \${LDFLAGS}"
  or use the --disable-fftw option to disable FFTW functionality if
  it isn't needed.
---------------------------------------------------------------------
      ])
    ])
  done
  #define FFTW_ENABLE in include file
  AC_DEFINE([FFTW_ENABLE],[],[Defined if fftw is enabled])
  ])
AM_CONDITIONAL([FFTW_ENABLE],[test "$FFTW_ENABLE" = "yes"])


#################################################################
## Check for HDF4 include and library
#################################################################

#check to see if user disabled hdf
HDF_ENABLE=yes
AC_ARG_ENABLE([hdf],
  [AS_HELP_STRING([--disable-hdf],
  [Disable hdf features. This includes not being able to create HDF4 files from model dumps.])],
  [HDF_ENABLE="$enableval"],
  [])

#check to make sure we have a library and include path, if we are using hdf
AS_IF([test "$HDF_ENABLE" = "yes"],[
  AC_CHECK_HEADERS(hdf.h,[],[
    AC_MSG_ERROR([
---------------------------------------------------------------------
  Unable to find the hdf.h header file.
 
  Try adding a path to the header file to the CPPFLAGS environment
  variable e.g. export CPPFLAGS="-I<include dir> \${CPPFLAGS}", or 
  using the --disable-hdf option to disable HDF functionality if
  it isn't needed.
---------------------------------------------------------------------
    ])
  ])
  AC_SEARCH_LIBS([cos],[m],[],[
      AC_MSG_ERROR([
---------------------------------------------------------------------
  Unable to find an m library containing the cos function.

  If you know the path to the library try adding it to the LDFLAGS
  environment variable. e.g. export LDFLAGS="-L<lib dir> \${LDFLAGS}"
  or use the --disable-hdf option to disable HDF functionality if
  it isn't needed.
---------------------------------------------------------------------
      ])
  ])
  
  LIBS="$LIBS -ldf -ljpeg -lz -lsz"
  for func in SDstart SDend SDcreate SDendaccess SDwritedata; do
    AC_SEARCH_LIBS([$func],[mfhdf],[],[
      AC_MSG_ERROR([
---------------------------------------------------------------------
  Unable to find an HDF library containing the $func function.

  If you know the path to the library try adding it to the LDFLAGS
  environment variable. e.g. export LDFLAGS="-L<lib dir> \${LDFLAGS}"
  or use the --disable-hdf option to disable HDF functionality if
  it isn't needed.
---------------------------------------------------------------------
      ])
    ],[])
  done
  AC_DEFINE([HDF_ENABLE],[],[Defined if hdf is enabled])
])
AM_CONDITIONAL([HDF_ENABLE],[test "$HDF_ENABLE" = "yes"])

#################################################################
#
#################################################################
## Check for petsc include and library
#################################################################
#
#check to see if user disabled petsc
PETSC_ENABLE=yes

#
#check to make sure we have a library and include path, if we are using petsc
AS_IF([test "$PETSC_ENABLE" = "yes"],[
  AC_CHECK_HEADERS(petscksp.h,[],[
    AC_MSG_ERROR([
---------------------------------------------------------------------
  Unable to find the petscksp.h header file.
 
  Try adding a path to the header file to the CPPFLAGS environment
  variable e.g. export CPPFLAGS="-I<include dir> \${CPPFLAGS}".
---------------------------------------------------------------------
    ])
  ])
  PETSCFUNCS="KSPSolve MatSetValues MatAssemblyBegin VecSetValues
    VecAssemblyBegin VecAssemblyEnd MatAssemblyEnd KSPSetOperators
    VecScatterBegin VecScatterEnd"
  #VecGetArray is also used, but for some reason the tests fail, however when 
  #building/linking with the main code that uses it there is no problem
  for func in $PETSCFUNCS; do
    AC_SEARCH_LIBS([$func], petsc,[],[
      AC_MSG_ERROR([
---------------------------------------------------------------------
  Unable to find an PETSc library containing the $func function.

  If you know the path to the library try adding it to the LDFLAGS
  environment variable. e.g. export LDFLAGS="-L<lib dir> \${LDFLAGS}".
---------------------------------------------------------------------
      ])
    ],[])
  done
  #
  #define PETSC_ENABLE in include file
  AC_DEFINE([PETSC_ENABLE],[],[Defined if PETSc is enabled])
  ])

#tell automake about PETSC_ENABLE variable
AM_CONDITIONAL([PETSC_ENABLE],[test "$PETSC_ENABLE" = "yes"])
#################################################################


#
#################################################################
## Check for CYTHON
#################################################################
#
#check to see if user disabled cython
CYTHON_ENABLE=yes
AC_ARG_ENABLE([cython],
  [AS_HELP_STRING([--disable-cython],
  [Disable cython dependent features, such as making vtk files for 
 visualization. Cython install should be added to your PYTHONPATH.])],
  [CYTHON_ENABLE="$enableval"],
  [])

#If we want to include Cython make sure it works
AS_IF(
  [test "$CYTHON_ENABLE" = "yes"],[
  
  #Check if Cython.Distutils.build_ext is available, needed for building evtk
  AC_MSG_CHECKING(for module Cython.Distutils.build_ext in python)
  echo "import Cython.Distutils.build_ext" | python -
  if test $? -ne 0 ; then
  AC_MSG_RESULT(not found)
  AC_MSG_ERROR(You need the module Cython.Distutils.build_ext available to 
  python to include cython dependent functionality. This can be disabled with 
  the --disable-cython option.)
  fi
  AC_MSG_RESULT(found)
  
  #Check if Cython.Build.cythonize is available, needed for building eos python
  #module
  AC_MSG_CHECKING(for module Cython.Build import of cythonize in python)
  echo "from Cython.Build import cythonize" | python -
  if test $? -ne 0 ; then
  AC_MSG_RESULT(not found)
  AC_MSG_ERROR(You need the module Cython.Build and the ability to import 
  cythonize available to python to include cython dependent functionality. This 
  can be disabled with the --disable-cython option.)
  fi
  AC_MSG_RESULT(found)
  
  #Check that numpy is installed
  AC_MSG_CHECKING(for module numpy in python)
  echo "import numpy as np" | python -
  if test $? -ne 0 ; then
  AC_MSG_RESULT(not found)
  AC_MSG_ERROR(You need the module numpy available to python to include cython dependent functionality. This 
  can be disabled with the --disable-cython option, but disabling cython support is not recommended.)
  fi
  AC_MSG_RESULT(found)
  ])
 
 #check to see if builds put things in lib or lib64
 AC_SUBST(PYTHONLIBDIR,lib)

AM_CONDITIONAL([CYTHON_ENABLE],[test "$CYTHON_ENABLE" = "yes"])
#################################################################
#AC_SUBST([prefix])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/pythonextensions/hdf/setup.py])
#scripts/basePath.py
AC_OUTPUT
