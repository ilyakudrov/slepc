#include <slepceps.h>
#include <iostream>
#include "data.h"
#include "link.h"
#include "matrix.h"
#include "eigen.h"

static char help[] = "Testing creation of AIJMPI matrix";

int x_size = 32;
int y_size = 32;
int z_size = 32;
int t_size = 32;

void mat_set_index(PetscInt* d_nnz, PetscInt* o_nnz, int low, int high);
void mat_insert(Mat A, int low, int high, data& conf, double mu_q, double mass);
PetscErrorCode create_test_vec(Vec vecx);

int main(int argc,char **argv){
    data conf;
    conf.read_float("/home/ilya/lattice/slepc/conf/nosmeared/time_32/mu0.00/conf_0001.fl");
    double mass = 0.0075;
    double mu_q = 0.1;
    Mat A;
    PetscErrorCode ierr;

    ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
    int rank, mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    PetscInt vec_size = x_size*y_size*z_size*t_size*2;
    PetscInt vec_size_local;
    if(mpi_size == 1) vec_size_local = vec_size;
    else if(rank != mpi_size - 1) vec_size_local = vec_size/mpi_size/2*2;
    else if(rank == mpi_size - 1) vec_size_local = vec_size - vec_size/mpi_size/2*2*(mpi_size - 1);
    int low, high;
    if(rank != mpi_size - 1){
        low = rank * vec_size_local;
        high = low + vec_size_local;
    }
    if(rank == mpi_size - 1){
        low = vec_size - vec_size_local;
        high = vec_size;
    }
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetType(A, MATMPIAIJ);
    MatSetSizes(A, vec_size_local, vec_size_local, vec_size, vec_size);
    PetscInt* d_nnz;
    PetscInt* o_nnz;
    if(!(d_nnz = (PetscInt*) malloc(vec_size_local * sizeof(PetscInt)))) PetscPrintf(PETSC_COMM_WORLD, "err malloc d_nnz");
    if(!(o_nnz = (PetscInt*) malloc(vec_size_local * sizeof(PetscInt)))) PetscPrintf(PETSC_COMM_WORLD, "err malloc o_nnz");
    for(int i = 0;i < vec_size_local;i++){
        d_nnz[i] = 0;
        o_nnz[i] = 0;
    }
    mat_set_index(d_nnz, o_nnz, low, high);

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank is %d; low = %d; high = %d\n", rank, low, high);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank is %d; d_nnz = %d; o_nnz = %d\n", rank, d_nnz[0], o_nnz[0]);
    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

    MatMPIAIJSetPreallocation(A, PETSC_DECIDE, d_nnz, PETSC_DECIDE, o_nnz);
    MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
    mat_insert(A, low, high, conf, mu_q, mass);

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    // MatView(A,PETSC_VIEWER_STDOUT_WORLD);

    Vec vecx, vecy;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, vec_size_local, vec_size, &vecx); CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, vec_size_local, vec_size, &vecy); CHKERRQ(ierr);

    create_test_vec(vecx);
    // VecView(vecx, PETSC_VIEWER_STDOUT_WORLD);
    MatMult(A, vecx, vecy);

    PetscScalar *px;
    PetscScalar *py;
    ierr = VecGetArray(vecx,&px);
    ierr = VecGetArray(vecy,&py);
    // for(int i = vec_size - 10;i < vec_size;i++){
    //     PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank is %d; px = %g+%gi; py = %g+%gi\n",
    //     rank, (double)PetscRealPart(px[i]), (double)PetscImaginaryPart(px[i]), (double)PetscRealPart(py[i]), (double)PetscImaginaryPart(py[i]));
    // }
    // PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
    VecView(vecy, PETSC_VIEWER_STDOUT_WORLD);
    ierr = VecRestoreArray(vecy,&py);CHKERRQ(ierr);

    VecDestroy(&vecx);
    VecDestroy(&vecy);
    free(d_nnz);
    free(o_nnz);
    MatDestroy(&A);
    PetscFinalize();
}

void mat_set_index(PetscInt* d_nnz, PetscInt* o_nnz, int low, int high){
    int x, y, z, t;
    link1 link_ferm(x_size, y_size, z_size, t_size);
    for(int place = low;place < high;place+=2){
        int t = place/(x_size * y_size * z_size*2);
        int z = (place - (x_size * y_size * z_size*2)*t)/(x_size * y_size*2);
        int y = (place - (x_size * y_size * z_size*2)*t - (2*x_size * y_size)*z)/(x_size*2);
        int x = (place - (x_size * y_size * z_size*2)*t - (2*x_size * y_size)*z - 2*x_size*y)/2;
        d_nnz[place - low]+=1;
        d_nnz[place - low + 1]+=1;

        link_ferm.go(x, y, z ,t);
        for(int mu = 1;mu <= 4;mu++){
            link_ferm.move(mu, 1);
            if(low <= complex_place(link_ferm) && complex_place(link_ferm) < high){
                d_nnz[place - low]+=2;
                d_nnz[place - low + 1]+=2;
            }
            if(low > complex_place(link_ferm) || complex_place(link_ferm) >= high){
                o_nnz[place - low]+=2;
                o_nnz[place - low + 1]+=2;
            }
            link_ferm.move(mu, -2);
            if(low <= complex_place(link_ferm) && complex_place(link_ferm) < high){
                d_nnz[place - low]+=2;
                d_nnz[place - low + 1]+=2;
            }
            if(low > complex_place(link_ferm) || complex_place(link_ferm) >= high){
                o_nnz[place - low]+=2;
                o_nnz[place - low + 1]+=2;
            }
            link_ferm.move(mu, 1);
        }
    }
}

void mat_insert(Mat A, int low, int high, data& conf, double mu_q, double mass){
    int rank, mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    int x, y, z, t;
    link1 link(x_size, y_size, z_size, t_size);
    link1 link_ferm(x_size, y_size, z_size, t_size);
    double delta_4 = 0;
    double sign;
    double border_sign;
    matrix B;
    cout<<"rank is "<<rank<<"; "<<low<<" "<<high<<endl;
    for(int place = low;place < high;place+=2){
        int t = place/(x_size * y_size * z_size*2);
        int z = (place - (x_size * y_size * z_size*2)*t)/(x_size * y_size*2);
        int y = (place - (x_size * y_size * z_size*2)*t - (2*x_size * y_size)*z)/(x_size*2);
        int x = (place - (x_size * y_size * z_size*2)*t - (2*x_size * y_size)*z - 2*x_size*y)/2;

        MatSetValue(A, place, place, mass, INSERT_VALUES);
        MatSetValue(A, place + 1, place + 1, mass, INSERT_VALUES);

        link_ferm.go(x, y, z ,t);
        link.go(x, y, z, t);
        for(int mu = 1;mu <= 4;mu++){
            if(mu == 4) delta_4 = 1;
            else delta_4 = 0;
            sign = eta_sign(mu, link_ferm);
            border_sign = link_ferm.border_sign(mu);
            link_ferm.move(mu, 1);
            link.move_dir(mu);
            B = link.get_matrix(conf.array);
            MatSetValue(A, place, complex_place(link_ferm), exp(mu_q * delta_4)/2 * border_sign * sign * (B.a0 + B.a3*PETSC_i), INSERT_VALUES);
            MatSetValue(A, place, complex_place(link_ferm) + 1, exp(mu_q * delta_4)/2 * border_sign * sign * (B.a2 + B.a1*PETSC_i), INSERT_VALUES);
            MatSetValue(A, place + 1, complex_place(link_ferm), exp(mu_q * delta_4)/2 * border_sign * sign * (-B.a2 + B.a1*PETSC_i), INSERT_VALUES);
            MatSetValue(A, place + 1, complex_place(link_ferm) + 1, exp(mu_q * delta_4)/2 * border_sign * sign * (B.a0 - B.a3*PETSC_i), INSERT_VALUES);
            link.move_dir(-mu);
            link_ferm.move(-mu, 1);
            border_sign = link_ferm.border_sign(-mu);
            link_ferm.move(-mu, 1);
            B = link.get_matrix(conf.array);
            MatSetValue(A, place, complex_place(link_ferm), -exp(-mu_q * delta_4)/2 * border_sign * sign * (B.a0 + B.a3*PETSC_i), INSERT_VALUES);
            MatSetValue(A, place, complex_place(link_ferm) + 1, -exp(-mu_q * delta_4)/2 * border_sign * sign * (B.a2 + B.a1*PETSC_i), INSERT_VALUES);
            MatSetValue(A, place + 1, complex_place(link_ferm), -exp(-mu_q * delta_4)/2 * border_sign * sign * (-B.a2 + B.a1*PETSC_i), INSERT_VALUES);
            MatSetValue(A, place + 1, complex_place(link_ferm) + 1, -exp(-mu_q * delta_4)/2 * border_sign * sign * (B.a0 - B.a3*PETSC_i), INSERT_VALUES);
            link_ferm.move(mu, 1);
        }
    }
}

PetscErrorCode create_test_vec(Vec vecx){
    PetscScalar *px;
    PetscErrorCode ierr;
    ierr = VecGetArray(vecx,&px);CHKERRQ(ierr);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PetscInt local_size;
    PetscInt low, high;
    VecGetLocalSize(vecx, &local_size);
    VecGetOwnershipRange(vecx, &low, &high);
    for(int i = 0;i < local_size/2;i++){
        int place = low + i * 2;
        int t = place/(x_size * y_size * z_size*2);
        int z = (place - (x_size * y_size * z_size*2)*t)/(x_size * y_size*2);
        int y = (place - (x_size * y_size * z_size*2)*t - (2*x_size * y_size)*z)/(x_size*2);
        int x = (place - (x_size * y_size * z_size*2)*t - (2*x_size * y_size)*z - 2*x_size*y)/2;
        px[i * 2] = 0.1 * (x+1) + 0.2 * (y+1) + 0.4 * (z+1) + 0.5 * (t+1) + (1 + 0.6 * (x+1) + 0.7 * (y+1) + 0.8 * (z+1) + 0.9 * (t+1)) * PETSC_i;
        px[i * 2 + 1] = 1 + 0.1 * (x+1) + 0.2 * (y+1) + 0.4 * (z+1) + 0.5 * (t+1) + (0.6 * (x+1) + 0.7 * (y+1) + 0.8 * (z+1) + 0.9 * (t+1)) * PETSC_i;
    }

    ierr = VecRestoreArray(vecx,&px);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}