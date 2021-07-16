INCLUDE_PATH=./
SRC_PATH=./
OBJ_PATH=./obj

CUDA_PATH=/usr/local/cuda-11.0
NVCC=${CUDA_PATH}/bin/nvcc

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

OBJFILES=$(OBJ_PATH)/link.o $(OBJ_PATH)/data.o $(OBJ_PATH)/eigen.o $(OBJ_PATH)/matrix.o

FLAGS1=-Xcompiler -fPIC -Xptxas -O3  -gencode arch=compute_61,code=sm_61  -Wno-deprecated-gpu-targets -c \
--compiler-options="-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -I${CUDA_PATH}/include    \
-L${SLEPC_DIR}/${PETSC_ARCH}/lib -lslepc -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include \
-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -I${CUDA_PATH}/include -fPIC -Wall \
-Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -O3 -march=native -mtune=native -ffast-math  \
 -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -O3 -march=native -mtune=native -ffast-math "

FLAGS2=-fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -O3 -march=native -mtune=native -ffast-math  \
-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include  \
-I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include \

FLAGS3=-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include  \
-L${SLEPC_DIR}/${PETSC_ARCH}/lib -lslepc -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include \
-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -I/usr/lib/cuda/include -fPIC -Wall \
-Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -O3 -march=native -mtune=native -ffast-math  \
 -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -O3 -march=native -mtune=native -ffast-math

eigen_double1: $(OBJ_PATH)/eigen_double.o ${OBJFILES}
	echo ${SLEPC_SYS_LIB}
	g++   $(OBJ_PATH)/eigen_double.o ${OBJFILES} -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib \
	-L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib \
	-L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/usr/local/cuda-11.0/lib64 -L/usr/local/cuda-11.0/lib64 \
	-Wl,-rpath,${SLEPC_DIR}/${PETSC_ARCH}/lib -L${SLEPC_DIR}/${PETSC_ARCH}/lib -lslepc \
	-Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/gcc/x86_64-linux-gnu/9 -Wl,-rpath,/usr/lib/x86_64-linux-gnu \
	-L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu \
	${SLEPC_SYS_LIB} -o eigen_double
#	-lpetsc -lf2clapack -lf2cblas -lm -lcufft -lcublas -lcudart -lcusparse -lcusolver -lcurand -lX11 -lquadmath -lstdc++ -ldl -o eigen_double

eigen_double: $(OBJ_PATH)/eigen_double.o ${OBJFILES}
	g++ $(OBJ_PATH)/eigen_double.o ${OBJFILES} ${SLEPC_SYS_LIB} -o eigen_double

eigen_double_cpu: $(OBJ_PATH)/eigen_double_cpu.o ${OBJFILES}
	g++ $(OBJ_PATH)/eigen_double_cpu.o ${OBJFILES} -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib \
	-L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib \
	-L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/usr/lib/cuda/lib64 -L/usr/lib/cuda/lib64 \
	-Wl,-rpath,${SLEPC_DIR}/${PETSC_ARCH}/lib -L${SLEPC_DIR}/${PETSC_ARCH}/lib -lslepc \
	-Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/gcc/x86_64-linux-gnu/9 -Wl,-rpath,/usr/lib/x86_64-linux-gnu \
	-L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lpetsc -lflapack -lfblas -lpthread \
	-lm -lcufft -lcublas -lcudart -lcusparse -lcusolver -lcurand -lX11 -lstdc++ -ldl -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -ldl -o eigen_double_cpu

test_multiplication: $(OBJ_PATH)/test_multiplication.o ${OBJFILES}
	g++ $(OBJ_PATH)/test_multiplication.o ${OBJFILES} -o test_multiplication

$(OBJ_PATH)/matrix.o: $(SRC_PATH)/matrix.cpp
	g++ -O3 -c $(SRC_PATH)/matrix.cpp -o $(OBJ_PATH)/matrix.o
$(OBJ_PATH)/eigen.o: $(SRC_PATH)/eigen.cpp
	g++ -O3 -c $(SRC_PATH)/eigen.cpp -o $(OBJ_PATH)/eigen.o
$(OBJ_PATH)/link.o: $(SRC_PATH)/link.cpp
	g++ -O3 -c $(SRC_PATH)/link.cpp -o $(OBJ_PATH)/link.o
$(OBJ_PATH)/data.o: $(SRC_PATH)/data.cpp
	g++ -O3 -c $(SRC_PATH)/data.cpp -o $(OBJ_PATH)/data.o

$(OBJ_PATH)/eigen_double.o: eigen_double.cpp
	$(NVCC) ${FLAGS1} eigen_double.cpp -o $(OBJ_PATH)/eigen_double.o

$(OBJ_PATH)/eigen_double_cpu.o: eigen_double.cpp
	nvcc ${FLAGS1} eigen_double.cpp -o $(OBJ_PATH)/eigen_double_cpu.o

$(OBJ_PATH)/test_multiplication.o: $(SRC_PATH)/test_multiplication.cpp
	g++ -O3 -c $(SRC_PATH)/test_multiplication.cpp -o $(OBJ_PATH)/test_multiplication.o

clean_local:
	rm -rf $(OBJ_PATH)/*.o
