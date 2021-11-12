################################################################################
# SPIDER Makefile                                                              #
################################################################################

# Executable Name
EXNAME = spider

# Source files, each corresponding to a .o object file and .d dependency file
SRC_C = \
        atmosphere.c \
        bc.c \
        cJSON.c \
	constants.c \
        ctx.c \
        dimensionalisablefield.c \
        energy.c \
        eos.c \
        eos_adamswilliamson.c \
        eos_composite.c \
        eos_lookup.c \
        eos_output.c \
        ic.c \
        interp.c \
        main.c \
        matprop.c \
        mesh.c \
        monitor.c \
        parameters.c \
        poststep.c \
        reaction.c \
        rheologicalfront.c \
        rhs.c \
        rollback.c \
        twophase.c \
        util.c \

# Main Target
all :: ${EXNAME}

### SPIDER Root Directory #####################################################
# Placement of this line matters. It will only work before any "include"s
SPIDER_ROOT_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))

### PETSc ######################################################################
# Include PETSc variables and rules
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

### Flags ######################################################################
# Extra flags
# Use -O0 to turn off optimization, to debug with LLDB/GDB
# You can specify this from the command line, e.g.
# Available sanitizer flags with OSX
# With debugging
#  make clean; make -j CFLAGS_EXTRA="-O0"
# with Address Sanitizer
#   make clean; make -j CFLAGS_EXTRA="-O0 -fsanitize=address"
# with UndefinedBehaviorSanitizer
#   make clean; make -j CFLAGS_EXTRA="-O0 -fsanitize=undefined"
CFLAGS+=${CFLAGS_EXTRA}

# Generate dependency (.d) files as we compile
CFLAGS+=${C_DEPFLAGS}

# Provide the current directory so that absolute paths to data files can be constructed.
CFLAGS+=-DSPIDER_ROOT_DIR=${SPIDER_ROOT_DIR}

### Compiling/Linking  #########################################################

# Objects (PETSc rules provides recipe)
SRC_O = ${SRC_C:%.c=%.o}

# Main executable
${EXNAME} : ${SRC_O}
	-${CLINKER} -o $@ $^ ${PETSC_TS_LIB}
	#${RM} $^

### Tests ######################################################################
SPIDER_TEST_DIR=${SPIDER_ROOT_DIR}/test_dir
SPIDER_TEST_SCRIPT=PYTHONPATH=${PYTHONPATH}:${SPIDER_ROOT_DIR}/tests/sciath python -m sciath ${SPIDER_ROOT_DIR}/tests/tests.yml
SPIDER_TEST_CONF=${SPIDER_TEST_DIR}/pth.conf

check_sciath:
	PYTHONPATH=${PYTHONPATH}:${PWD}/tests/sciath ./tests/check_sciath.sh

test_create_output_dir :
	mkdir -p ${SPIDER_TEST_DIR}

test : test_create_output_dir check_sciath
	cd ${SPIDER_TEST_DIR} && ${SPIDER_TEST_SCRIPT} -w ${SPIDER_TEST_CONF} && cd -
	@printf "Test output lives in ${SPIDER_TEST_DIR}\n"
	@printf "If on a batch system, wait until jobs complete and then\n"
	@printf "  make test_check\n"

test_check : test_create_output_dir check_sciath
	cd ${SPIDER_TEST_DIR} && ${SPIDER_TEST_SCRIPT} -w ${SPIDER_TEST_CONF} -v && cd -

.PHONY: test test_create_output_dir

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
