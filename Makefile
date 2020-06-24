################################################################################
# SPIDER Makefile                                                              #
################################################################################

# Executable Name
EXNAME = spider

# Source files, each corresponding to a .o object file and .d dependency file
SRC_C = main.c ctx.c rhs.c ic.c mesh.c bc.c util.c twophase.c \
        monitor.c energy.c matprop.c dimensionalisablefield.c rheologicalfront.c \
        cJSON.c rollback.c poststep.c parameters.c atmosphere.c \
        reaction.c eos.c \

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
SPIDER_TEST_SCRIPT=python -m sciath ${SPIDER_ROOT_DIR}/tests/tests.yml
SPIDER_TEST_CONF=${SPIDER_TEST_DIR}/pth.conf

test_create_output_dir :
	mkdir -p ${SPIDER_TEST_DIR}

# Basic Tests
test : test_create_output_dir
	cd ${SPIDER_TEST_DIR} && ${SPIDER_TEST_SCRIPT} -w ${SPIDER_TEST_CONF} -g basic  && cd -
	@printf "Test output lives in ${SPIDER_TEST_DIR}\n"
	@printf "To run more tests\n"
	@printf "  make test_all\n"
	@printf "If on a batch system, wait until jobs complete and then\n"
	@printf "  make test_check\n"

test_check : test_create_output_dir
	cd ${SPIDER_TEST_DIR} && ${SPIDER_TEST_SCRIPT} -w ${SPIDER_TEST_CONF} -v -g basic && cd -

# Atmosphere tests
test_atmos : test_create_output_dir
	cd ${SPIDER_TEST_DIR} && ${SPIDER_TEST_SCRIPT} -w ${SPIDER_TEST_CONF} -g atmos && cd -
	@printf "Test output lives in ${SPIDER_TEST_DIR}\n"
	@printf "To run more tests\n"
	@printf "  make test_all\n"
	@printf "If on a batch system, wait until jobs complete and then\n"
	@printf "  make testatmoscheck\n"

test_atmos_check : test_create_output_dir
	cd ${SPIDER_TEST_DIR} && ${SPIDER_TEST_SCRIPT} -w ${SPIDER_TEST_CONF} -v -g atmos && cd -

# All Tests
test_all : test_create_output_dir
	cd ${SPIDER_TEST_DIR} && ${SPIDER_TEST_SCRIPT} -w ${SPIDER_TEST_CONF} && cd -
	@printf "Test output lives in ${SPIDER_TEST_DIR}\n"
	@printf "If on a batch system, wait until jobs complete and then\n"
	@printf "  make test_all_check\n"

test_all_check : test_create_output_dir
	cd ${SPIDER_TEST_DIR} && ${SPIDER_TEST_SCRIPT} -w ${SPIDER_TEST_CONF} -v && cd -

.PHONY: test testatmos testall test_create_output_dir

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
