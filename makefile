# executable name
EXNAME = main

# Required objects (rules included below build these from .c files)
OBJ = main.o ctx.o rhs.o ic.o mesh.o lookup.o bc.o util.o twophase.o monitor.o

ALL: ${EXNAME}

clean ::
	rm -f ${EXNAME} ${OBJ}

.PHONY: ALL clean

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

PETSC_CC_INCLUDES+=-I/opt/local/include

GSL_LIB=-L/opt/local/lib -lgsl -lgslcblas

${EXNAME} : ${OBJ}
	-${CLINKER}  -o $@ $^  ${GSL_LIB} ${PETSC_TS_LIB} 
	#${RM} $^

${OBJ} : global_defs.h

### Output ####################################################################

# We do not by default clear any of the output data (so beware of stale files)
#  use "make clear_output" to delete the output
clear_output :
	rm -rf output/TIMESTEPPER/S_s
	mkdir output/TIMESTEPPER/S_s
	touch output/TIMESTEPPER/S_s/.gitignore
	rm -rf output/TIMESTEPPER/rhs
	mkdir output/TIMESTEPPER/rhs
	touch output/TIMESTEPPER/rhs/.gitignore

.PHONY : clear_output

### Tests #####################################################################

# Note: to update the reference for a test, you can
#  1. Comment out the @rm -f testX.tmp line
#  2. Run the test
#  3. cp testX.tmp testref/testX.ref

test : test1 test2 test3 test4

TEST_OPTIONS = -test_view -ts_view 

# SINIT = 1600
test1: ${EXNAME}
	@rm -f test1.tmp
	@echo "\033[34mRunning Test 1\033[0m"
	@${MPIEXEC} -n 1 ./${EXNAME} ${TEST_OPTIONS} -nstepsmacro 0 -sinit 1600 \
    2>&1 > test1.tmp
	@diff test1.tmp testref/test1.ref && \
    echo "\033[32mSuccess\033[0m" || \
    echo "\033[31mFailure: output does not match reference (see diff above)\033[0m"
	@rm -f test1.tmp

# SINIT = 2500
test2: ${EXNAME}
	@rm -f test2.tmp
	@echo "\033[34mRunning Test 2\033[0m"
	@${MPIEXEC} -n 1 ./${EXNAME} ${TEST_OPTIONS} -nstepsmacro 0  -sinit 2500 \
    2>&1 > test2.tmp
	@diff test2.tmp testref/test2.ref && \
    echo "\033[32mSuccess\033[0m" || \
    echo "\033[31mFailure: output does not match reference (see diff above)\033[0m"
	@rm -f test2.tmp

# SINIT = 3000
test3: ${EXNAME}
	@rm -f test3.tmp
	@echo "\033[34mRunning Test 3\033[0m"
	@${MPIEXEC} -n 1 ./${EXNAME} ${TEST_OPTIONS} -nstepsmacro 0 -sinit 3000 \
    2>&1 > test3.tmp
	@diff test3.tmp testref/test3.ref && \
    echo "\033[32mSuccess\033[0m" || \
    echo "\033[31mFailure: output does not match reference (see diff above)\033[0m"
	@rm -f test3.tmp

# SINIT = 2500, parallel
test4: ${EXNAME}
	@rm -f test4.tmp
	@echo "\033[34mRunning Test 4\033[0m"
	@${MPIEXEC} -n 5 ./${EXNAME} ${TEST_OPTIONS} -nstepsmacro 0 -sinit 2500 \
    2>&1 > test4.tmp
	@diff test4.tmp testref/test4.ref && \
    echo "\033[32mSuccess\033[0m" || \
    echo "\033[31mFailure: output does not match reference (see diff above)\033[0m"
	@rm -f test4.tmp

.PHONY: test test1 test2 test3 test4

include ${PETSC_DIR}/lib/petsc/conf/test
