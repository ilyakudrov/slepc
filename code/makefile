#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-2019, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

CFLAGS     =
FFLAGS     =
CPPFLAGS   =
FPPFLAGS   =
LOCDIR     = src/eps/examples/tutorials/
EXAMPLESC  = ex1.c ex2.c ex3.c ex4.c ex5.c ex7.c ex9.c ex10.c ex11.c ex12.c ex13.c \
             ex18.c ex19.c ex24.c ex25.c ex29.c ex30.c ex31.c ex34.c ex35.c ex36.c ex41.c
EXAMPLESF  = ex1f.F ex1f90.F90 ex6f90.F90
MANSEC     = EPS

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

ex1: ex1.o chkopts
	-${CLINKER} -o ex1 ex1.o ${SLEPC_EPS_LIB}
	${RM} ex1.o
ex3: ex3.o matrix.o data.o link.o chkopts
	g++ -o ex3 ex3.o matrix.o data.o link.o ${SLEPC_EPS_LIB}
	${RM} ex3.o matrix.o data.o link.o
ex4: ex4.o chkopts
	-${CLINKER} -o ex4 ex4.o ${SLEPC_EPS_LIB}
	${RM} ex4.o
eigen_test: eigen_test.o matrix.o data.o link.o eigen.o
	g++ -o eigen_test eigen_test.o matrix.o data.o link.o eigen.o ${SLEPC_EPS_LIB}
	#${RM} eigen_test.o matrix.o data.o link.o
	#mpic++ -o eigen_test eigen_test.o matrix.o data.o link.o eigen.o -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -g3   -I/home/ilya/soft/slepc-3.12.0/include -I/home/ilya/soft/slepc-3.12.0/arch-linux-c-debug/include  -I/home/ilya/soft/petsc/include -I/home/ilya/soft/petsc/arch-linux-c-debug/include
	#mpicc -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -g3  -o eigen_test eigen_test.o -Wl,-rpath,/home/ilya/soft/slepc-3.12.0/arch-linux-c-debug/lib -L/home/ilya/soft/slepc-3.12.0/arch-linux-c-debug/lib -lslepc       -Wl,-rpath,/home/ilya/soft/petsc/arch-linux-c-debug/lib -L/home/ilya/soft/petsc/arch-linux-c-debug/lib -Wl,-rpath,/home/ilya/soft/petsc/arch-linux-c-debug/lib -L/home/ilya/soft/petsc/arch-linux-c-debug/lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/x86_64-linux-gnu/openmpi/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/8 -L/usr/lib/gcc/x86_64-linux-gnu/8 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lpetsc -lflapack -lfblas -lm -lX11 -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lquadmath -lstdc++ -ldl
	#${RM} eigen_test.o matrix.o data.o link.o
	${RM} eigen_test.o matrix.o data.o link.o eigen.o
mat_mult_test_parallel: mat_mult_test_parallel.o matrix.o data.o link.o eigen.o
	g++ -o mat_mult_test_parallel mat_mult_test_parallel.o matrix.o data.o link.o eigen.o ${SLEPC_EPS_LIB}
	${RM} mat_mult_test_parallel.o matrix.o data.o link.o eigen.o
ex3.o: ex3.c
	mpic++ -o ex3.o -c -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -g3   -I/home/ilya/soft/slepc-3.12.0/include -I/home/ilya/soft/slepc-3.12.0/arch-linux-c-debug/include  -I/home/ilya/soft/petsc/include -I/home/ilya/soft/petsc/arch-linux-c-debug/include -I/usr/local/cuda/include    `pwd`/ex3.c
eigen_test.o: eigen_test.c
	mpic++ -o eigen_test.o -c -w -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -g3   -I/home/ilya/soft/slepc-3.12.0/include -I/home/ilya/soft/slepc-3.12.0/arch-linux-c-debug/include  -I/home/ilya/soft/petsc/include -I/home/ilya/soft/petsc/arch-linux-c-debug/include -I/usr/local/cuda/include    `pwd`/eigen_test.c
mat_mult_test_parallel.o: mat_mult_test_parallel.c
	mpic++ -o mat_mult_test_parallel.o -c -w -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -g3   -I/home/ilya/soft/slepc-3.12.0/include -I/home/ilya/soft/slepc-3.12.0/arch-linux-c-debug/include  -I/home/ilya/soft/petsc/include -I/home/ilya/soft/petsc/arch-linux-c-debug/include -I/usr/local/cuda/include    `pwd`/mat_mult_test_parallel.c
matrix.o: matrix.cpp
	mpic++ -c matrix.cpp
data.o: data.cpp
	mpic++ -c data.cpp
link.o: link.cpp
	mpic++ -c link.cpp
eigen.o: eigen.cpp
	mpic++ -c eigen.cpp