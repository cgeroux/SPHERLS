######################
#### MPI MAKEFILE ####
######################

#program name
TARGET = SPHERLS

#path
TARGETPATH = /home/cgeroux/SPHERLS/branches/LES_New

#number of processors
NUMPROCS = 4

#time in hours minutes seconds
RUNTIME=35:00:00

#job name
#JOBNAME=########
JOBNAME=NewTurb

#email to mail job notifications to
EMAIL=chris.m.geroux@gmail.com

#virtual memory
VMEM=0.1G

#stack size
STACK=0.1G

#files needed to make target
LIBDIR=/home/cgeroux/lib
EXTOBJS =\
  ${LIBDIR}/exception/exception2.o\
  ${LIBDIR}/xmlParser/xmlParser.o\
  ${LIBDIR}/xmlFunctions/xmlFunctions.o

INTOBJS=\
  main.o\
  physEquations.o\
  dataManipulation.o\
  watchzone.o\
  dataMonitoring.o\
  global.o\
  eos.o
  
INTSOURCE=\
  main.cpp\
  physEquations.cpp\
  dataManipulation.cpp\
  watchzone.cpp\
  dataMonitoring.cpp\
  eos.cpp\
  global.cpp
  
INTHEADERS=\
  main.h\
  physEquations.h\
  dataManipulation.h\
  watchzone.h\
  dataMonitoring.h\
  eos.h\
  global.h\
  userguide.h

#macros to set compiler name/options and file naming conventions
COMPILER  = mpic++
LINKER    = mpic++
LINKBEFFLAGS = -fsignaling-nans -Wall
debug: LINKBEFFLAGS += -g

LINKAFTFLAGS=-o

COMPBEFFLAGS = -fsignaling-nans -Wall
debug: COMPBEFFLAGS += -g
COMPBEFFLAGS+=-c
COMPAFTFLAGS = -o
OBJECTEXT = .o
SOURCEEXT = .cpp
HEADEREXT = .h
STDLIBS   = 
LIBS      = -L${LIBDIR}/petsc/lib -lpetsc
INCLUDE   =\
  -I${LIBDIR}/exception\
  -I${LIBDIR}/xmlParser\
  -I${LIBDIR}/xmlFunctions\
  -I${LIBDIR}/petsc/include

##################
#### make bin ####
##################


#link program to make binary
${TARGET}: ${EXTOBJS} ${INTOBJS} ${INTHEADERS}
	@echo "linking ..."
	${LINKER} ${LINKBEFFLAGS} ${LINKAFTFLAGS} ${TARGET} ${EXTOBJS} ${INTOBJS} ${LIBS}

#make bin with debugging flags
debug: ${TARGET}

#compiles source files that have header files, but doesn't link, i.e. rebuilds if the file,
#  if the header has changed
%${OBJECTEXT}:%${SOURCEEXT} %${HEADEREXT}
	@echo "compiling object that has header ..."
	${COMPILER} ${COMPBEFFLAGS} $< ${COMPAFTFLAGS} $@ ${INCLUDE}

#compiles source files that don't have header files, but doesn't link, i.e. allows building of files
#  without header files
%${OBJECTEXT}:%${SOURCEEXT}
	@echo "compiling object that doesn't have header ..."
	${COMPILER} ${COMPBEFFLAGS} $< ${COMPAFTFLAGS} $@ ${INCLUDE}

####################
##### clean up #####
####################

#cleans up object files, and executable
clean:cleano cleanout cleanx cleansh cleancore

#cleans up object files
cleano:
  #remove object files if they exsist
	@(for i in ${EXTOBJS} ${TARGET}.o ${INTOBJS}; do\
	  if [ -f $$i ]; then rm $$i;echo rm $$i;\
	  fi;\
	  done;)

#cleans up output files of run
cleanout:
  #remove outputfiles
	@(for i in ${TARGET}*.out; do\
	  if [ -f $$i ]; then rm $$i;echo rm $$i;\
	  fi;\
	  done;)
	@(if [ -f ${TARGET}.err ]; then rm ${TARGET}.err;echo rm ${TARGET}.err;fi)

#cleans up executable
cleanx:
  #remove executable
	@(if [ -f ${TARGET} ]; then rm ${TARGET};echo rm ${TARGET};fi)

#cleans up submit script
cleansh:
  #remove submit scripts
	@(if [ -f ${TARGET}_test_que.sh ]; then rm ${TARGET}_test_que.sh;echo rm ${TARGET}_test_que.sh;fi)
	@(if [ -f ${TARGET}.sh ]; then rm ${TARGET}.sh;echo rm ${TARGET}.sh;fi)

cleancore:
#cleans up any core dumps
	@(for i in core*; do\
	  if [ -f $$i ]; then rm $$i;echo rm $$i;\
	  fi;\
	  done;)

#delete documentation
cleandocs:
	rm -r -f ./docs

#################
#### run job ####
#################

#production job
run:submit_script
	@qsub ${JOBNAME}_que.sh

#runs on test que
runt:submit_script_testq
	@qsub ${JOBNAME}_test_que.sh

#run on current machine
rund:
	mpirun -np ${NUMPROCS} ${TARGETPATH}/${TARGET}

#run mpi code with totalview debugger
rundd:
	mpirun --debug -np ${NUMPROCS} ${TARGETPATH}/${TARGET}

#creates the submit script with correct names for the test que
submit_script_testq:
	@echo '#!/bin/bash'           >${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Shell to run job in'>>${JOBNAME}_test_que.sh
	@echo '#$$ -S /bin/bash'      >>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Set job name'>>${JOBNAME}_test_que.sh
	@echo '#$$ -N ${JOBNAME}'>>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Run job form Current Working Directory'>>${JOBNAME}_test_que.sh
	@echo '#$$ -cwd'            >>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Output file name'            >>${JOBNAME}_test_que.sh
	@echo '#$$ -o ${JOBNAME}.out'>>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Error output file'            >>${JOBNAME}_test_que.sh
	@echo '#$$ -e ${JOBNAME}.err'>>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Error is not merged with standard out'            >>${JOBNAME}_test_que.sh
	@echo '#$$ -j n'            >>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Number or processors to run job with'>>${JOBNAME}_test_que.sh
	@echo '#$$ -pe ompi* ${NUMPROCS}'>>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Run time in seconds'>>${JOBNAME}_test_que.sh
	@echo '#$$ -l h_rt=${RUNTIME}'>>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Amount of memory used'>>${JOBNAME}_test_que.sh
	@echo '##$$ -l h_vmem=${VMEM}'>>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Submit to the test que'>>${JOBNAME}_test_que.sh
	@echo '#$$ -l test=true'>>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Transfer all environment variables when job is submitted'>>${JOBNAME}_test_que.sh
	@echo '#$$ -V'>>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Stack size??'>>${JOBNAME}_test_que.sh
	@echo '##$$ -l h_stack=${STACK}'>>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Email on end, abort or suspend'>>${JOBNAME}_test_que.sh
	@echo '#$$ -m eas'>>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
	@echo '## Who to email'>>${JOBNAME}_test_que.sh
	@echo '#$$ -M ${EMAIL}'>>${JOBNAME}_test_que.sh
	@echo '##'>>${JOBNAME}_test_que.sh
#	@echo '## Checkpointing via grid engine'>>${JOBNAME}_test_que.sh
#	@echo '#$$ -ckpt transparent'>>${JOBNAME}_test_que.sh
#	@echo '#$$ -c xs'>>${JOBNAME}_test_que.sh
#	@echo '#$$ -notify'>>${JOBNAME}_test_que.sh
#	@echo "trap '' usr2">>${JOBNAME}_test_que.sh
		@echo 'mpirun ${TARGETPATH}/${TARGET}'>>${JOBNAME}_test_que.sh

#creates the submit script with correct names for the normal que
submit_script:
	@echo '#!/bin/bash'>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
	@echo '## Shell to run job in'>>${JOBNAME}_que.sh
	@echo '#$$ -S /bin/bash'>>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
	@echo '## Set job name'>>${JOBNAME}_que.sh
	@echo '#$$ -N ${JOBNAME}'>>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
	@echo '## Run job form Current Working Directory'>>${JOBNAME}_que.sh
	@echo '#$$ -cwd'>>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
	@echo '## Output file name'>>${JOBNAME}_que.sh
	@echo '#$$ -o ${JOBNAME}.out'>>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
	@echo '## Error output file'>>${JOBNAME}_que.sh
	@echo '#$$ -e ${JOBNAME}.err'>>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
	@echo '## Error is not merged with standard out'>>${JOBNAME}_que.sh
	@echo '#$$ -j n'>>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
	@echo '## Number or processors to run job with'>>${JOBNAME}_que.sh
	@echo '#$$ -pe ompi* ${NUMPROCS}'>>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
	@echo '## Run time in seconds'>>${JOBNAME}_que.sh
	@echo '#$$ -l h_rt=${RUNTIME}'>>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
	@echo '## Amount of memory used'>>${JOBNAME}_que.sh
	@echo '##$$ -l h_vmem=${VMEM}'>>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
	@echo '## Transfer all environment variables when job is submitted'>>${JOBNAME}_que.sh
	@echo '#$$ -V'>>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
	@echo '## Stack size??'>>${JOBNAME}_que.sh
	@echo '##$$ -l h_stack=${STACK}'>>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
	@echo '## Email on end, abort or suspend'>>${JOBNAME}_que.sh
	@echo '#$$ -m eas'>>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
	@echo '## Who to email'>>${JOBNAME}_que.sh
	@echo '#$$ -M ${EMAIL}'>>${JOBNAME}_que.sh
	@echo '##'>>${JOBNAME}_que.sh
#	@echo '## Checkpointing via grid engine'>>${JOBNAME}_que.sh
#	@echo '#$$ -ckpt transparent'>>${JOBNAME}_que.sh
#	@echo '#$$ -c xs'>>${JOBNAME}_que.sh
#	@echo '#$$ -notify'>>${JOBNAME}_que.sh
#	@echo "trap '' usr2">>${JOBNAME}_que.sh
	@echo 'mpirun ${TARGETPATH}/${TARGET}'>>${JOBNAME}_que.sh

#make documentation
docs:${INTSOURCE} Doxyfile ${INTHEADERS}
	@echo ${INTSOURCE}
	@doxygen
	@(cd ./docs/latex; make pdf)

