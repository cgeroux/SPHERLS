######################
#### MPI MAKEFILE ####
######################

#program name
TARGET = SPHERLSgen

#files needed to make target
LIBDIR=/home/cgeroux/lib
SPHERLSDIR=/home/cgeroux/SPHERLS/trunk
EXTOBJS =\
  ${LIBDIR}/exception/exception2.o\
  ${LIBDIR}/xmlParser/xmlParser.o\
  ${LIBDIR}/xmlFunctions/xmlFunctions.o\
  ${SPHERLSDIR}/eos.o
  
INTOBJS=\
  main.o
  
INTSOURCE=\
  main.cpp
  
INTHEADERS=\
  main.h

#macros to set compiler name/options and file naming conventions
COMPILER  = g++
LINKER    = g++
LINKFLAGS = -fsignaling-nans -Wall -g -o
COMPBEFFLAGS = -fsignaling-nans -Wall -g -c
COMPAFTFLAGS = -o
OBJECTEXT = .o
SOURCEEXT = .cpp
HEADEREXT = .h
LIBS      = 
INCLUDE   =\
  -I${LIBDIR}/exception\
  -I${LIBDIR}/xmlParser\
  -I${LIBDIR}/xmlFunctions\
  -I${SPHERLSDIR}\

##################
#### make bin ####
##################

#link program to make binary
${TARGET}: ${EXTOBJS} ${INTOBJS} ${INTHEADERS}
	${LINKER} ${LINKFLAGS} ${TARGET} ${EXTOBJS} ${INTOBJS} ${LIBS}

#compiles source files, but doesn't link
%${OBJECTEXT}:%${SOURCEEXT}
	${COMPILER} ${COMPBEFFLAGS} $< ${COMPAFTFLAGS} $@ ${INCLUDE}

####################
##### clean up #####
####################

#cleans up object files, and executable
clean:cleano cleanx cleancore

#cleans up object files
cleano:
  #remove object files if they exsist
	@(for i in ${EXTOBJS} ${TARGET}.o ${INTOBJS}; do\
	  if [ -f $$i ]; then rm $$i;echo rm $$i;\
	  fi;\
	  done;)

#cleans up executable
cleanx:
  #remove executable
	@(if [ -f ${TARGET} ]; then rm ${TARGET};echo rm ${TARGET};fi)

cleancore:
#cleans up any core dumps
	@(for i in core*; do\
	  if [ -f $$i ]; then rm $$i;echo rm $$i;\
	  fi;\
	  done;)

#delete documentation
cleandocs:
	rm -r -f ./docs

#make documentation
docs:${INTSOURCE} Doxyfile ${INTHEADERS}
	@echo ${INTSOURCE}
	@doxygen
	@(cd ./docs/latex; make pdf)

