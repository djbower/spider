# executable name
EXNAME = main

# Required objects (rules included below build these from .c files)
OBJ = main.o ctx.o rhs.o ic.o mesh.o lookup.o bc.o util.o twophase.o monitor.o aug.o energy.o matprop.o atmosphere.o output.o

ALL: ${EXNAME}

clean ::
	rm -f ${EXNAME} ${OBJ}

.PHONY: ALL clean

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# This is to ensure that the debugging commands work in LLDB
# Don't use for optimized code!
CFLAGS_EXTRA=-g -O0
# You can specify this from the command line, e.g.
#   make clean; make -j CFLAGS_EXTRA="-O0"
#CFLAGS_EXTRA=-g -O3
CFLAGS+=${CFLAGS_EXTRA}

# Provide the current directory so that absolute
# paths to data files can be constructed. Note that this assumes that CURDIR
# is correct (you can break this with make -f /somewhere/else/to/here/Makefile)
CFLAGS+=-DMAGMA_ROOT_DIR=${CURDIR}

${EXNAME} : ${OBJ}
	-${CLINKER}  -o $@ $^ ${PETSC_TS_LIB}
	#${RM} $^

${OBJ} : global_defs.h

### Output ####################################################################

# We do not by default clear any of the output data (so beware of stale files)
#  use "make clear_output" to delete the output
clear_output :
	rm -f output/*.m

.PHONY : clear_output

### Tests #####################################################################
test :
	cd tests && ./runTests.py && cd ..
.PHONY : test
