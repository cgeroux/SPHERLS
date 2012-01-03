AC_INIT([SPHERLS], [1.0], [chris.m.geroux@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CXX
AC_LANG([C++])

#################################################################
## Check for FFTW3 include and library 
#################################################################

#common places for fftw libraries and includes
FFTW_INCLUDE_PATHS="/home/cgeroux/lib/fftw-3.2.2/include /usr/local/FFTW3/include \
  /usr/include /share/apps/fftw/include"
FFTW_LIB_PATHS="/home/cgeroux/lib/fftw-3.2.2/lib /usr/local/FFTW3/lib /usr/lib /share/apps/fftw/lib"

#check to see if user set an overide for include path
AC_ARG_ENABLE([fftw-include-path],
  [AS_HELP_STRING([--enable-fftw-include-path],
    [set location of the fftw3 headers])],
  [FFTW_INCLUDE_PATHS="$enableval"],
  [])

#check to see if user set an overide for library path
AC_ARG_ENABLE([fftw-lib-path],
  [AS_HELP_STRING([--enable-fftw-lib-path],
    [set location of the fftw3 library])],
  [FFTW_LIB_PATHS="$enableval"],
  [])

#check to see if user disabled fftw
FFTW_ENABLE=yes
AC_ARG_ENABLE([fftw],
  [AS_HELP_STRING([--disable-fftw],
    [Disable fftw features])],
  [FFTW_ENABLE="$enableval"],
  [])

#check to make sure we have a library and include path, if we are using fftw
AS_IF([test "$FFTW_ENABLE" = "yes"],[
  # Get FFTW library and include locations
  
  #search include paths for fftw3.h
  FFTW_INCLUDE_PATH=""
  for TEST_PATH in $FFTW_INCLUDE_PATHS; do
    AC_CHECK_FILE([${TEST_PATH}/fftw3.h],[FFTW_INCLUDE_PATH=${TEST_PATH};break],[])
  done
  
  #seach lib paths for libfftw3.so and libfftw3.a
  FFTW_LIB_PATH=""
  for TEST_PATH in $FFTW_LIB_PATHS; do
    AC_CHECK_FILES([${TEST_PATH}/libfftw3.so ${TEST_PATH}/libfftw3.a ]
      ,[FFTW_LIB_PATH=${TEST_PATH};break],[])
  done

  AS_IF([test "$FFTW_LIB_PATH" = ""],[
    AC_MSG_RESULT(["No fftw3 library found in search paths \"$FFTW_LIB_PATHS\". Either specify the\
    path to libfftw3.a/.so with --enable-fftw-lib-path or disable fftw features with --disable-fftw\
    ."])
    exit
    ])
  AS_IF([test "$FFTW_INCLUDE_PATH" = ""],[
    AC_MSG_RESULT(["No fftw3 include file found in search paths \"$FFTW_INCLUDE_PATHS\". Either \
    specify the path to fftw3.h with --enable-fftw-include-path or disable fftw features with \
    --disable-fftw ."])
    exit
    ])
  AC_DEFINE([FFTW_ENABLE],[],[Defined if fftw is enabled])
  AC_SUBST([FFTW_INCLUDE_PATH])
  AC_SUBST([FFTW_LIB_PATH])
  ])
#################################################################

#################################################################
## Check for HDF4 include and library 
#################################################################

#common places for HDF libraries and includes
HDF4_INCLUDE_PATHS="/home/cgeroux/lib/hdf4/include"
HDF4_LIB_PATHS="/home/cgeroux/lib/hdf4/lib"

#check to see if user set an overide for include path
AC_ARG_ENABLE([hdf-include-path],
  [AS_HELP_STRING([--enable-hdf-include-path],
    [set location of the hdf4 headers])],
  [HDF4_INCLUDE_PATHS="$enableval"],
  [])

#check to see if user set an overide for library path
AC_ARG_ENABLE([hdf-lib-path],
  [AS_HELP_STRING([--enable-hdf-lib-path],
    [set location of the hdf4 library])],
  [HDF4_LIB_PATHS="$enableval"],
  [])

#check to see if user disabled HDF
HDF_ENABLE=yes
AC_ARG_ENABLE([hdf],
  [AS_HELP_STRING([--disable-hdf],
    [Disable hdf features])],
  [HDF_ENABLE="$enableval"],
  [])

#check to make sure we have a library and include path, if we are using HDF
AS_IF([test "$HDF_ENABLE" = "yes"],[
  # Get HDF library and include locations
  
  #search include paths for hdf.h
  HDF4_INCLUDE_PATH=""
  for TEST_PATH in $HDF4_INCLUDE_PATHS; do
    AC_CHECK_FILE([${TEST_PATH}/hdf.h],[HDF4_INCLUDE_PATH=${TEST_PATH};break],[])
  done
  
  #seach lib paths for libmfhdf.so and libmfhdf.a
  HDF4_LIB_PATH=""
  for TEST_PATH in $HDF4_LIB_PATHS; do
    AC_CHECK_FILES([${TEST_PATH}/libmfhdf.so ${TEST_PATH}/libmfhdf.a ]
      ,[HDF4_LIB_PATH=${TEST_PATH};break],[])
  done

  AS_IF([test "$HDF4_LIB_PATH" = ""],[
    AC_MSG_RESULT(["No hdf4 library found in search paths \"$HDF4_LIB_PATHS\". Either specify the\
    path to libmhdf.a/.so with --enable-hdf-lib-path or disable hdf features with --disable-hdf\
    ."])
    exit
    ])
  AS_IF([test "$HDF4_INCLUDE_PATH" = ""],[
    AC_MSG_RESULT(["No hdf4 include file found in search paths \"$HDF4_INCLUDE_PATHS\". Either \
    specify the path to hdf.h with --enable-hdf-include-path or disable hdf features with \
    --disable-hdf ."])
    exit
    ])
  AC_DEFINE([HDF_ENABLE],[],[Defined if hdf is enabled])
  AC_SUBST([HDF_INCLUDE_PATH])
  AC_SUBST([HDF_LIB_PATH])
  ])
#################################################################


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
   Makefile
   src/Makefile
   src/SPHERLS/Makefile
   src/SPHERLSgen/Makefile
   src/SPHERLSanal/Makefile
  ])
AC_OUTPUT