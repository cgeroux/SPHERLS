AC_INIT([SPHERLS], [1.0], [chris.m.geroux@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CXX
AC_LANG([C++])
#
#common places for fftw libraries and includes
FFTW_INCLUDE_PATHS="/home/cgeroux/lib/fftw-3.2.2/include  /share/apps/fftw/include "
FFTW_LIB_PATHS="/home/cgeroux/lib/fftw-3.2.2/lib /share/apps/fftw/lib /usr/local/FFTW3/lib /usr/lib /share/apps/fftw/lib"
#
#common places for fftw libraries and includes
HDF_INCLUDE_PATHS="/home/cgeroux/lib/hdf4/include /share/apps/hdf4/include"
HDF_LIB_PATHS="/home/cgeroux/lib/hdf4/lib /share/apps/hdf4/lib"
#
#common places for PETSc libraries and includes
PETSC_INCLUDE_PATHS="/home/cgeroux/lib/petsc/include"
PETSC_LIB_PATHS="/home/cgeroux/lib/petsc/lib"
#
#################################################################
## Check for FFTW3 include and library
#################################################################
#
#check to see if user disabled fftw
FFTW_ENABLE=yes
AC_ARG_ENABLE([fftw],
  [AS_HELP_STRING([--disable-fftw],
  [Disable fftw features])],
  [FFTW_ENABLE="$enableval"],
  [])
#
#check to see if user set an overide for include path
AC_ARG_ENABLE([fftw-include-path],
 [AS_HELP_STRING([--enable-fftw-include-path],
 [set location of the fftw3 headers])],
 [FFTW_INCLUDE_PATHS="$enableval"],
 [])
#
#check to see if user set an overide for library path
AC_ARG_ENABLE([fftw-lib-path],
  [AS_HELP_STRING([--enable-fftw-lib-path],
  [set location of the fftw3 library])],
  [FFTW_LIB_PATHS="$enableval"],
  [])
#
#check to make sure we have a library and include path, if we are using fftw
AS_IF([test "$FFTW_ENABLE" = "yes"],[
  #
  #search include paths for fftw3.h
  FFTW_INCLUDE_PATH=""
  for TEST_PATH in $FFTW_INCLUDE_PATHS; do
    AC_CHECK_FILE([${TEST_PATH}/fftw3.h],[FFTW_INCLUDE_PATH=${TEST_PATH};break],[])
  done
  #
  #seach lib paths for libfftw3.so and libfftw3.a
  FFTW_LIB_PATH=""
  for TEST_PATH in $FFTW_LIB_PATHS; do
    AC_CHECK_FILES([${TEST_PATH}/libfftw3.so ${TEST_PATH}/libfftw3.a ]\
      ,[FFTW_LIB_PATH=${TEST_PATH};break],[])
  done
  #
  #check that a lib path was found
  AS_IF([test "$FFTW_LIB_PATH" = ""],[
    AC_MSG_RESULT(["No fftw3 library found in search paths \"$FFTW_LIB_PATHS\". Either specify the\
    path to libfftw3.a/.so with --enable-fftw-lib-path or disable fftw features with --disable-fftw\
    ."])
    exit
    ])
  #
  #check that an include path was found
  AS_IF([test "$FFTW_INCLUDE_PATH" = ""],[
    AC_MSG_RESULT(["No fftw3 include file found in search paths \"$FFTW_INCLUDE_PATHS\". Either \
    specify the path to fftw3.h with --enable-fftw-include-path or disable fftw features with \
    --disable-fftw ."])
    exit
    ])
  #
  #sub strings for lib and include paths
  AC_SUBST(FFTW_LIB_PATH)
  AC_SUBST(FFTW_INCLUDE_PATH)
  #
  #define FFTW_ENABLE in include file
  AC_DEFINE([FFTW_ENABLE],[],[Defined if fftw is enabled])
  ])
#################################################################

#################################################################
## Check for HDF4 include and library and dependent libraries
#################################################################

#check to see if user disabled hdf
HDF_ENABLE=yes
AC_ARG_ENABLE([hdf],
  [AS_HELP_STRING([--disable-hdf],
  [Disable hdf features])],
  [HDF_ENABLE="$enableval"],
  [])

#check to see if user set an overide for include path
AC_ARG_ENABLE([hdf-include-path],
 [AS_HELP_STRING([--enable-hdf-include-path],
 [set location of the hdf headers])],
 [HDF_INCLUDE_PATHS="$enableval"],
 [])

#check to see if user set an overide for library path
AC_ARG_ENABLE([hdf-lib-path],
  [AS_HELP_STRING([--enable-hdf-lib-path],
  [set location of the hdf library])],
  [HDF_LIB_PATHS="$enableval"],
  [])

#check to make sure we have a library and include path, if we are using hdf
AS_IF([test "$HDF_ENABLE" = "yes"],[
  # Get HDF library and include locations
  
  # #search include paths for hdf.h
  HDF_INCLUDE_PATH=""
  for TEST_PATH in $HDF_INCLUDE_PATHS; do
    AC_CHECK_FILE([${TEST_PATH}/hdf.h],[HDF_INCLUDE_PATH=${TEST_PATH};break],[])
  done
  
  #seach lib paths for libmfhdf.so and libmfhdf.a
  HDF_LIB_PATH=""
  for TEST_PATH in $HDF_LIB_PATHS; do
    AC_CHECK_FILES([${TEST_PATH}/libmfhdf.so ${TEST_PATH}/libmfhdf.a ]
      ,[HDF_LIB_PATH=${TEST_PATH};break],[])
  done

  AS_IF([test "$HDF_LIB_PATH" = ""],[
    AC_MSG_RESULT(["No hdf library found in search paths \"$HDF_LIB_PATHS\". Either specify the\
    path to libmfhdf.a/.so with --enable-hdf-lib-path or disable hdf features with --disable-hdf\
    ."])
    exit
    ])
  AS_IF([test "$HDF_INCLUDE_PATH" = ""],[
    AC_MSG_RESULT(["No hdf include file found in search paths \"$HDF_INCLUDE_PATHS\". Either \
    specify the path to hdf.h with --enable-hdf-include-path or disable hdf features with \
    --disable-hdf."])
    exit
    ])
  AC_DEFINE([HDF_ENABLE],[],[Defined if hdf is enabled])
  AC_SUBST([HDF_INCLUDE_PATH])
  AC_SUBST([HDF_LIB_PATH])
  ])
#################################################################
#
#################################################################
## Check for petsc include and library
#################################################################
#
#check to see if user disabled petsc
PETSC_ENABLE=yes
AC_ARG_ENABLE([petsc],
  [AS_HELP_STRING([--disable-petsc],
  [Disable petsc features])],
  [PETSC_ENABLE="$enableval"],
  [])
#
#check to see if user set an overide for include path
AC_ARG_ENABLE([petsc-include-path],
 [AS_HELP_STRING([--enable-petsc-include-path],
 [set location of the petsc headers])],
 [PETSC_INCLUDE_PATHS="$enableval"],
 [])
#
#check to see if user set an overide for library path
AC_ARG_ENABLE([petsc-lib-path],
  [AS_HELP_STRING([--enable-petsc-lib-path],
  [set location of the petsc library])],
  [PETSC_LIB_PATHS="$enableval"],
  [])
#
#check to make sure we have a library and include path, if we are using petsc
AS_IF([test "$PETSC_ENABLE" = "yes"],[
  #
  #search include paths for petscksp.h
  PETSC_INCLUDE_PATH=""
  for TEST_PATH in $PETSC_INCLUDE_PATHS; do
    AC_CHECK_FILE([${TEST_PATH}/petscksp.h],[PETSC_INCLUDE_PATH=${TEST_PATH};break],[])
  done
  #
  #seach lib paths for libfftw3.so and libfftw3.a
  PETSC_LIB_PATH=""
  for TEST_PATH in $PETSC_LIB_PATHS; do
    AC_CHECK_FILES([${TEST_PATH}/libpetsc.so ${TEST_PATH}/libpetsc.a ]\
      ,[PETSC_LIB_PATH=${TEST_PATH};break],[])
  done
  #
  #check that a lib path was found
  AS_IF([test "$PETSC_LIB_PATH" = ""],[
    AC_MSG_RESULT(["No petsc library found in search paths \"$PETSC_LIB_PATHS\". Either specify the\
    path to libpetsc.a/.so with --enable-petsc-lib-path or disable petsc features with \
    --disable-petsc."])
    exit
    ])
  #
  #check that an include path was found
  AS_IF([test "$PETSC_INCLUDE_PATH" = ""],[
    AC_MSG_RESULT(["No petsc include file found in search paths \"$PETSC_INCLUDE_PATHS\". Either \
    specify the path to petscksp.h with --enable-petsc-include-path or disable petsc features with \
    --disable-petsc ."])
    exit
    ])
  #
  #sub strings for lib and include paths
  AC_SUBST(PETSC_LIB_PATH)
  AC_SUBST(PETSC_INCLUDE_PATH)
  #
  #define PETSC_ENABLE in include file
  AC_DEFINE([PETSC_ENABLE],[],[Defined if PETSc is enabled])
  ])
#################################################################

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
   Makefile\
   src/Makefile\
   src/SPHERLS/Makefile\
   src/SPHERLSgen/Makefile\
   src/SPHERLSanal/Makefile\
  ])
AC_OUTPUT
