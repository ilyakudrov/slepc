#include <slepceps.h>
#include <iostream>
#include "data.h"
#include "link.h"
#include "matrix.h"
#include "eigen.h"

static char help[] = "Testing multiplication with MPI";

typedef struct {
   data conf;
   int x_size;
   int y_size;
   int z_size;
   int t_size;
   double mass;
   double mu_q;
} mat_data;

PetscErrorCode MatMult_test_parallel(Vec x,Vec y, const mat_data ctx);
PetscErrorCode create_test_vec(Vec vecx);

int x_size = 32;
int y_size = 32;
int z_size = 32;
int t_size = 32;

int main(int argc,char **argv){
    MPI_Comm comm = MPI_COMM_WORLD;
    PetscErrorCode ierr;
    mat_data my_data;
    my_data.x_size = x_size;
    my_data.y_size = y_size;
    my_data.z_size = z_size;
    my_data.t_size = t_size;
    my_data.conf.read_float("/home/ilya/lattice/slepc/conf/nosmeared/time_32/mu0.10/conf_0001.fl"/*argv[5]*/);
    my_data.mass = 0.0075/*atof(argv[6])*/;
    my_data.mu_q = 10/*atof(argv[7])*/;

    ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
    int rank, mpi_size;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &rank);
    Vec vecx, vecx1, vecy;
    PetscInt vec_size = x_size*y_size*z_size*t_size*2;
    PetscInt vec_size_local;
    if(rank != mpi_size - 1) vec_size_local = vec_size/mpi_size/2*2;
    if(rank == mpi_size - 1) vec_size_local = vec_size - vec_size/mpi_size/2*2*(mpi_size - 1);
    ierr = VecCreateMPI(comm, vec_size, vec_size, &vecx); CHKERRQ(ierr);
    ierr = VecCreateMPI(comm, vec_size_local, vec_size, &vecy); CHKERRQ(ierr);

    PetscInt local_size;
    PetscInt low, high;
    VecGetLocalSize(vecx, &local_size);
    VecGetOwnershipRange(vecx, &low, &high);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"I'm the rank %d process; low is %d; high is %d\n",rank, low, high);
    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

    create_test_vec(vecx);



    MatMult_test_parallel(vecx, vecy, my_data);

    VecDestroy(&vecx);
    VecDestroy(&vecy);
}

void MatVecMult(const matrix& A, const PetscScalar* x, PetscScalar* y, int border_sign){
    PetscScalar z1[2], z2[2];
    z1[0] = A.a0 + A.a3*PETSC_i;
    z2[0] = A.a2 + A.a1*PETSC_i;
    z1[1] = -A.a2 + A.a1*PETSC_i;
    z2[1] = A.a0 - A.a3*PETSC_i;
    y[0] = border_sign * (x[0] * z1[0] + x[1] * z2[0]);
    y[1] = border_sign * (x[0] * z1[1] + x[1] * z2[1]);
}

PetscErrorCode MatMult_test_parallel(const Vec vecx,Vec vecy, const mat_data ctx){
    const PetscScalar *px;
    PetscScalar       *py;
    PetscErrorCode    ierr;
    //PetscFunctionBeginUser;
    ierr = VecGetArrayRead(vecx,&px);CHKERRQ(ierr);
    ierr = VecGetArray(vecy,&py);CHKERRQ(ierr);
    double mass = ctx.mass;
    double mu_q = ctx.mu_q;
    int x_size = ctx.x_size;
    int y_size = ctx.y_size;
    int z_size = ctx.z_size;
    int t_size = ctx.t_size;
    data conf = ctx.conf;
    int lat_size = x_size*y_size*z_size*t_size;
    // PetscSynchronizedPrintf(PETSC_COMM_WORLD,"ok1\n");
    // PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PetscInt local_size;
    PetscInt low, high;
    VecGetLocalSize(vecx, &local_size);
    VecGetOwnershipRange(vecx, &low, &high);
    // PetscSynchronizedPrintf(PETSC_COMM_WORLD,"ok2\n");
    // PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
    double delta_4 = 0;
    double sign;
    double border_sign;
    int place;
    PetscScalar vec[2];
    PetscScalar res[2];
    link1 link(x_size, y_size, z_size, t_size);
    link1 link_ferm(x_size, y_size, z_size, t_size);
    // PetscSynchronizedPrintf(PETSC_COMM_WORLD,"ok3\n");
    // PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
    for(int i = 0;i < local_size/2;i++){
        int place = low + i * 2;
        int t = place/(x_size * y_size * z_size*2);
        int z = (place - (x_size * y_size * z_size*2)*t)/(x_size * y_size*2);
        int y = (place - (x_size * y_size * z_size*2)*t - (2*x_size * y_size)*z)/(x_size*2);
        int x = (place - (x_size * y_size * z_size*2)*t - (2*x_size * y_size)*z - 2*x_size*y)/2;

        link.go(x, y, z, t);
        link_ferm.go(x, y, z, t);
        place = complex_place(link_ferm);
        for(int mu = 1;mu <= 4;mu++){
            if(mu == 4) delta_4 = 1;
            else delta_4 = 0;
            link.move_dir(mu);
            sign = eta_sign(mu, link_ferm);
            border_sign = link_ferm.border_sign(mu);
            link_ferm.move(mu, 1);
            // matrix_mult_complex1(exp(mu_q * delta_4)/2 * sign * link.get_matrix1(data_conf), &src[complex_place(link_ferm)], &vec, src_index, border_sign);
            MatVecMult(exp(mu_q * delta_4)/2 * sign * link.get_matrix(conf.array), &px[complex_place(link_ferm)], vec, border_sign);
            res[0] = res[0] + vec[0];
            res[1] = res[1] + vec[1];
            link.move_dir(-mu);
            link_ferm.move(-mu, 1);
            border_sign = link_ferm.border_sign(-mu);
            link_ferm.move(-mu, 1);
            // matrix_mult_complex1(exp(-mu_q * delta_4)/2 * sign * link.get_matrix1(data_conf), &src[complex_place(link_ferm)], &vec, src_index, border_sign);
            MatVecMult(exp(-mu_q * delta_4)/2 * sign * link.get_matrix(conf.array), &px[complex_place(link_ferm)], vec, border_sign);
            res[0] = res[0] - vec[0];
            res[1] = res[1] - vec[1];
            link_ferm.move(mu, 1);
        }
        py[i * 2] = res[0] + mass*px[complex_place(link_ferm)];
        py[i * 2 + 1] = res[1] + mass*px[complex_place(link_ferm) + 1];
        res[0] = 0;
        res[1] = 0;

        ierr = VecRestoreArrayRead(vecx,&px);CHKERRQ(ierr);
        ierr = VecRestoreArray(vecy,&py);CHKERRQ(ierr);
        PetscFunctionReturn(0);
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
        int t = place/(x_size * y_size * z_size*2) + 1;
        int z = (place - (x_size * y_size * z_size*2)*(t - 1))/(x_size * y_size*2) + 1;
        int y = (place - (x_size * y_size * z_size*2)*(t - 1) - (2*x_size * y_size)*(z - 1))/(x_size*2) + 1;
        int x = (place - (x_size * y_size * z_size*2)*(t - 1) - (2*x_size * y_size)*(z - 1) - 2*x_size*(y - 1))/2 + 1;
        px[i * 2] = 0.1 * (x+1) + 0.2 * (y+1) + 0.4 * (z+1) + 0.5 * (t+1) + (1 + 0.6 * (x+1) + 0.7 * (y+1) + 0.8 * (z+1) + 0.9 * (t+1)) * PETSC_i;
        px[i * 2 + 1] = 1 + 0.1 * (x+1) + 0.2 * (y+1) + 0.4 * (z+1) + 0.5 * (t+1) + (0.6 * (x+1) + 0.7 * (y+1) + 0.8 * (z+1) + 0.9 * (t+1)) * PETSC_i;
    }

    ierr = VecRestoreArray(vecx,&px);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}