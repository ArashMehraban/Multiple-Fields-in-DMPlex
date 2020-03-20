export PETSC_DIR=/home/arash/math/petsc-3.14/petsc
export PETCS_ARCH=arch-linux2-c-debug

CFLAGS           = -std=c99
FFLAGS           =
CPPFLAGS         =
FPPFLAGS         =

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

.DEFAULT_GOAL:= main

main: main.o 
	-${CLINKER} -o main main.o ${PETSC_KSP_LIB} ${PETSC_LIB} ${PETSC_SNES_LIB}
	${RM} main.o 





