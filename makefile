ALL: main

.PHONY: ALL

# TODO: add simple test against a ref output

clean ::
	rm -f *.o main

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

PETSC_CC_INCLUDES+=-I/opt/local/include

#CFLAGS+= -ferror-limit=3

GSL_LIB=-L/opt/local/lib -lgsl -lgslcblas

main: main.o ctx.o rhs.o
	-${CLINKER}  -o $@ $^  ${GSL_LIB} ${PETSC_TS_LIB} 
	#${RM} $<

include ${PETSC_DIR}/lib/petsc/conf/test
