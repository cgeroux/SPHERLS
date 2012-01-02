AC_INIT([SPHERLS], [1.0], [chris.m.geroux@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CXX
AC_LANG([C++])

# Get FFTW library and include locations

#common places for fftw libraries and includes
FFTW_INCLUDE_PATHS="/usr/local/FFTW3/include /usr/include /share/apps/fftw/include"
# /share/apps/fftw/include => lachesis  ( CentOS release 5.5 (Final), Rocks release 5.4 (Maverick) )
# /usr/local/FFTW3/include => ACEnet    ( Red Hat Enterprise Linux AS release 4 (Nahant Update 9) )
# /usr/include             => Guenevere ( Ubuntu 10.04.3 LTS )
FFTW_LIB_PATHS="/usr/local/FFTW3/lib /usr/lib /share/apps/fftw/lib"
# /share/apps/fftw/lib => lachesis  ( CentOS release 5.5 (Final), Rocks release 5.4 (Maverick) )
# /usr/local/FFTW3/lib => ACEnet    ( Red Hat Enterprise Linux AS release 4 (Nahant Update 9) )
# /usr/lib             => Guenevere ( Ubuntu 10.04.3 LTS )

#search include paths for fftw3.h
FFTW_INCLUDE_PATH=""
for TEST_PATH in $FFTW_INCLUDE_PATHS; do
  AC_CHECK_FILE([${TEST_PATH}/fftw3.h],[FFTW_INCLUDE_PATH=${TEST_PATH};break],[])
done

#seach lib paths for libfftw3.so and libfftw3.a
FFTW_LIB_PATH=""
for TEST_PATH in $FFTW_LIB_PATHS; do
  AC_CHECK_FILES([${TEST_PATH}/libfftw3.so ${TEST_PATH}/libfftw3.a ],[FFTW_LIB_PATH=${TEST_PATH};break],[])
done

#check to see if user set an overide for include path
AC_ARG_ENABLE([fftw-include-path],
  [AS_HELP_STRING([--enable-fftw-include-path],
    [set location of the fftw3 headers])],
  [FFTW_INCLUDE_PATH="$enableval"],
  [])
AC_SUBST([FFTW_INCLUDE_PATH])

#check to see if user set an overide for library path
AC_ARG_ENABLE([fftw-lib-path],
  [AS_HELP_STRING([--enable-fftw-lib-path],
    [set location of the fftw3 library])],
  [FFTW_LIB_PATH="$enableval"],
  [])
AC_SUBST([FFTW_LIB_PATH])

#check to make sure we have a library and include path, if not let user know and quit


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT