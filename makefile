ALL: main

.PHONY: ALL

clean ::
	rm -f *.o main

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

CFLAGS+= -ferror-limit=3

main: main.o ctx.o 
	-${CLINKER} - -o main main.o ctx.o  ${PETSC_TS_LIB}
	#${RM} main.o ctx.o

include ${PETSC_DIR}/lib/petsc/conf/test
