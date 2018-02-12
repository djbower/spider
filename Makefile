# Executable Name
EXNAME = spider

# Source files, each corresponding to a .o object file and .d dependency file
SRC_C = main.c ctx.c rhs.c ic.c mesh.c lookup.c bc.c util.c \
				twophase.c monitor.c aug.c energy.c matprop.c output.c parameters.c \
				scalablefield.c
SRC_O = ${SRC_C:%.c=%.o}
SRC_D = ${SRC_C:%.c=%.d}

all :: ${EXNAME}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# This is to ensure that the debugging commands work in LLDB
# Don't use for optimized code!
#CFLAGS_EXTRA=-O0

# You can specify this from the command line, e.g.
#   make clean; make -j CFLAGS_EXTRA="-O0"
CFLAGS+=${CFLAGS_EXTRA}

# Generate dependencies
CFLAGS+=${C_DEPFLAGS}

# Provide the current directory so that absolute
# paths to data files can be constructed. Note that this assumes that CURDIR
# is correct (you can break this with make -f /somewhere/else/to/here/Makefile)
CFLAGS+=-DSPIDER_ROOT_DIR=${CURDIR}

${EXNAME} : ${SRC_O}
	-${CLINKER} -o $@ $^ ${PETSC_TS_LIB}
	#${RM} $^

clean ::
	rm -f ${EXNAME} ${SRC_O} ${SRC_D}

.PHONY: ALL clean

### Tests ######################################################################
PYTHON=python
test :
	cd tests && ${PYTHON} ./runTests.py -t bottomUpLiquid && cd ..

testall :
	cd tests && ${PYTHON} ./runTests.py && cd ..

.PHONY : testall

### Dependencies ###############################################################
# Indicate that SRC_D is up to date. Prevents the include from having quadratic complexity.
$(SRC_D) : ;

# Include dependency files
-include $(SRC_D)

.DELETE_ON_ERROR:
