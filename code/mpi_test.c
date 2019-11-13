#include <mpi.h>
#include <stdio.h>

void print_rank(int* a){
    int rank, size;
    printf("adding 1\n");
    *a += 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("I'm %d of %d\n", rank, size);
}

void func1(){
    printf("executing function1\n");

}

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);
    func1();
    MPI_Finalize();
}