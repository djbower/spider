################################################################################
# SPIDER Makefile                                                              #
################################################################################

# Executable Name
EXNAME = spider

# Source files, each corresponding to a .o object file and .d dependency file
SRC_C = main.c ctx.c rhs.c ic.c mesh.c lookup.c bc.c util.c twophase.c \
        monitor.c energy.c matprop.c parameters.c \
        dimensionalisablefield.c \
				cJSON.c

# Main Target
all :: ${EXNAME}

### PETSc ######################################################################
# Include PETSc variables and rules
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

### Flags ######################################################################
# Extra flags
# Use -O0 to turn off optimization, to debug with LLDB/GDB
# You can specify this from the command line, e.g.
#   make clean; make -j CFLAGS_EXTRA="-O0"
CFLAGS+=${CFLAGS_EXTRA}

# Generate dependency (.d) files as we compile
CFLAGS+=${C_DEPFLAGS}

# Provide the current directory so that absolute
# paths to data files can be constructed. Note that this assumes that CURDIR
# is correct (you can break this with make -f /somewhere/else/to/here/Makefile)
CFLAGS+=-DSPIDER_ROOT_DIR=${CURDIR}

### Compiling/Linking  #########################################################

# Objects (PETSc rules provides recipe)
SRC_O = ${SRC_C:%.c=%.o}

# Main executable
${EXNAME} : ${SRC_O}
	-${CLINKER} -o $@ $^ ${PETSC_TS_LIB}
	#${RM} $^

### Tests ######################################################################
test :
	cd tests && ./runTests.py -t bottomUpLiquid && cd ..
	@printf "If on a batch system, wait until jobs complete and then\n"
	@printf "  make testcheck\n"

testcheck :
	cd tests && ./runTests.py -v -t bottomUpLiquid && cd ..

testall :
	cd tests && ./runTests.py && cd ..
	@printf "If on a batch system, wait until jobs complete and then\n"
	@printf "  make testallcheck\n"

testallcheck :
	cd tests && ./runTests.py -v && cd ..

.PHONY: test testall

### Dependencies ###############################################################
SRC_D = ${SRC_C:%.c=%.d}

# Indicate that SRC_D is up to date. Prevents the include from having quadratic complexity.
$(SRC_D) : ;

# Include dependency files
-include $(SRC_D)

### Helper Targets #############################################################
clean ::
	rm -f ${EXNAME} ${SRC_O} ${SRC_D}

.PHONY: clean

### Misc #######################################################################

# Remove results of partial, failed builds
.DELETE_ON_ERROR:
