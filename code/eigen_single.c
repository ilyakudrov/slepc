/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2019, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Solves eigenproblem. "
                     "The problem is a standard symmetric eigenproblem "
                     "corresponding to the 2-D Laplacian operator.\n\n"
                     "The command line options are:\n"
                     "  -n <n>, where <n> = number of grid subdivisions in "
                     "both x and y dimensions.\n\n";

#include "data.h"
#include "eigen.h"
#include "link.h"
#include "matrix.h"
#include <iostream>
#include <slepceps.h>

/*
   User-defined routines
*/
typedef struct {
  data conf;
  int x_size;
  int y_size;
  int z_size;
  int t_size;
  double mass;
  double mu_q;
} mat_data;

void mat_set_index(PetscInt *d_nnz, PetscInt *o_nnz, int low, int high);
void mat_insert(Mat A, int low, int high, data &conf, double mu_q, double mass);
PetscErrorCode MatMult_eigen_sequential(Mat A, Vec x, Vec y);
PetscErrorCode MatGetDiagonal_Laplacian2D(Mat A, Vec diag);
PetscErrorCode MatMult_simple(Mat A, Vec vecx, Vec vecy);
PetscErrorCode TestMatMul(mat_data &my_data, const PetscScalar *px,
                          PetscScalar *py);
void CheckMatMult(mat_data &my_data);

int x_size = 32;
int y_size = 32;
int z_size = 32;
int t_size;

int main(int argc, char **argv) {
  // int x_size = 32/*atof(argv[1])*/;
  // int y_size = 32/*atof(argv[2])*/;
  // int z_size = 32/*atof(argv[3])*/;
  t_size = atof(argv[2]);
  mat_data my_data;
  my_data.x_size = x_size /*atof(argv[1])*/;
  my_data.y_size = y_size /*atof(argv[2])*/;
  my_data.z_size = z_size /*atof(argv[3])*/;
  my_data.t_size = t_size /*atof(argv[4])*/;
  my_data.conf.read_double(
      /*"/home/ilya/lattice/slepc/conf/nosmeared/time_32/mu0.00/conf_0001.fl"*/
      argv[1]);
  my_data.mass = 0.0075 /*atof(argv[6])*/;
  my_data.mu_q = atof(argv[3]);
  cout << "mu = " << my_data.mu_q << endl;
  cout << "path = " << argv[1] << endl;
  cout << "time = " << t_size << endl;
  Mat A;   /* operator matrix */
  EPS eps; /* eigenproblem solver context */
  EPSType type;
  PetscMPIInt size;
  PetscInt N, n = x_size * y_size * z_size * t_size * 2, nev;
  PetscBool terse;
  PetscErrorCode ierr;

  // CheckMatMult();

  ierr = SlepcInitialize(&argc, &argv, (char *)0, help);
  if (ierr)
    return ierr;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nProblem with %D processors\n", size);
  CHKERRQ(ierr);
  // if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example
  // only");

  ierr = PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL);
  CHKERRQ(ierr);
  N = n * n;
  ierr = PetscPrintf(
      PETSC_COMM_WORLD,
      "\nDirac Eigenproblem (matrix-free version), N=%D (%Dx%D grid)\n\n", N, n,
      n);
  CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  //   ierr =
  //   MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,&my_data,&A);CHKERRQ(ierr);
  //   ierr =
  //   MatShellSetOperation(A,MATOP_MULT,(void(*)(void))MatMult_eigen_sequential);CHKERRQ(ierr);
  // ierr =
  // MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_Laplacian2D);CHKERRQ(ierr);
  // ierr =
  // MatShellSetOperation(A,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Laplacian2D);CHKERRQ(ierr);

  int rank, mpi_size;
  MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscInt vec_size = x_size * y_size * z_size * t_size * 2;
  PetscInt vec_size_local;
  if (mpi_size == 1)
    vec_size_local = vec_size;
  else if (rank != mpi_size - 1)
    vec_size_local = vec_size / mpi_size / 2 * 2;
  else if (rank == mpi_size - 1)
    vec_size_local = vec_size - vec_size / mpi_size / 2 * 2 * (mpi_size - 1);
  int low, high;
  if (rank != mpi_size - 1) {
    low = rank * vec_size_local;
    high = low + vec_size_local;
  }
  if (rank == mpi_size - 1) {
    low = vec_size - vec_size_local;
    high = vec_size;
  }
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetType(A, MATMPIAIJ);
  MatSetSizes(A, vec_size_local, vec_size_local, vec_size, vec_size);
  PetscInt *d_nnz;
  PetscInt *o_nnz;
  if (!(d_nnz = (PetscInt *)malloc(vec_size_local * sizeof(PetscInt))))
    PetscPrintf(PETSC_COMM_WORLD, "err malloc d_nnz");
  if (!(o_nnz = (PetscInt *)malloc(vec_size_local * sizeof(PetscInt))))
    PetscPrintf(PETSC_COMM_WORLD, "err malloc o_nnz");
  for (int i = 0; i < vec_size_local; i++) {
    d_nnz[i] = 0;
    o_nnz[i] = 0;
  }
  mat_set_index(d_nnz, o_nnz, low, high);
  MatMPIAIJSetPreallocation(A, PETSC_DECIDE, d_nnz, PETSC_DECIDE, o_nnz);
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  mat_insert(A, low, high, my_data.conf, my_data.mu_q, my_data.mass);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  // MatView(A, PETSC_VIEWER_STDOUT_SELF);

  free(d_nnz);
  free(o_nnz);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD, &eps);
  CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps, A, NULL);
  CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_NHEP);
  CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);
  CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);
  CHKERRQ(ierr);

  // TESTING
  //   CheckMatMult(my_data);
  //   cout<<"check"<<endl;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  cout << "rank is " << rank << "; starting to solve" << endl;
  ierr = EPSSolve(eps);
  CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetType(eps, &type);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, " Solution method: %s\n\n", type);
  CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps, &nev, NULL, NULL);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of requested eigenvalues: %D\n",
                     nev);
  CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  ierr = PetscOptionsHasName(NULL, NULL, "-terse", &terse);
  CHKERRQ(ierr);
  if (terse) {
    ierr = EPSErrorView(eps, EPS_ERROR_RELATIVE, NULL);
    CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,
                                 PETSC_VIEWER_ASCII_INFO_DETAIL);
    CHKERRQ(ierr);
    ierr = EPSReasonView(eps, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);
    ierr = EPSErrorView(eps, EPS_ERROR_RELATIVE, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);
    PetscScalar eigr, eigi;
    /*for(int i = 0;i < 2;i++){
       EPSGetEigenpair(eps,i,&eigr,&eigi,NULL,NULL);
       cout<<"eigenvalue "<<i<<" : "<<eigr<<" "<<eigi<<endl;
    }*/
  }
  ierr = EPSDestroy(&eps);
  CHKERRQ(ierr);
  ierr = MatDestroy(&A);
  CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*
    Compute the matrix vector multiplication y<---T*x where T is a nx by nx
    tridiagonal matrix with DD on the diagonal, DL on the subdiagonal, and
    DU on the superdiagonal.
 */
static void tv(int nx, const PetscScalar *x, PetscScalar *y) {
  PetscScalar dd, dl, du;
  int j;

  dd = 4.0;
  dl = -1.0;
  du = -1.0;

  y[0] = dd * x[0] + du * x[1];
  for (j = 1; j < nx - 1; j++)
    y[j] = dl * x[j - 1] + dd * x[j] + du * x[j + 1];
  y[nx - 1] = dl * x[nx - 2] + dd * x[nx - 1];
}

/*
    Matrix-vector product subroutine for the 2D Laplacian.

    The matrix used is the 2 dimensional discrete Laplacian on unit square with
    zero Dirichlet boundary condition.

    Computes y <-- A*x, where A is the block tridiagonal matrix

                 | T -I          |
                 |-I  T -I       |
             A = |   -I  T       |
                 |        ...  -I|
                 |           -I T|

    The subroutine TV is called to compute y<--T*x.
 */

void mat_set_index(PetscInt *d_nnz, PetscInt *o_nnz, int low, int high) {
  int x, y, z, t;
  link1 link_ferm(x_size, y_size, z_size, t_size);
  for (int place = low; place < high; place += 2) {
    int t = place / (x_size * y_size * z_size * 2);
    int z =
        (place - (x_size * y_size * z_size * 2) * t) / (x_size * y_size * 2);
    int y = (place - (x_size * y_size * z_size * 2) * t -
             (2 * x_size * y_size) * z) /
            (x_size * 2);
    int x = (place - (x_size * y_size * z_size * 2) * t -
             (2 * x_size * y_size) * z - 2 * x_size * y) /
            2;
    d_nnz[place - low] += 1;
    d_nnz[place - low + 1] += 1;

    link_ferm.go(x, y, z, t);
    for (int mu = 1; mu <= 4; mu++) {
      link_ferm.move(mu, 1);
      if (low <= complex_place(link_ferm) && complex_place(link_ferm) < high) {
        d_nnz[place - low] += 2;
        d_nnz[place - low + 1] += 2;
      }
      if (low > complex_place(link_ferm) || complex_place(link_ferm) >= high) {
        o_nnz[place - low] += 2;
        o_nnz[place - low + 1] += 2;
      }
      link_ferm.move(mu, -2);
      if (low <= complex_place(link_ferm) && complex_place(link_ferm) < high) {
        d_nnz[place - low] += 2;
        d_nnz[place - low + 1] += 2;
      }
      if (low > complex_place(link_ferm) || complex_place(link_ferm) >= high) {
        o_nnz[place - low] += 2;
        o_nnz[place - low + 1] += 2;
      }
      link_ferm.move(mu, 1);
    }
  }
}

void mat_insert(Mat A, int low, int high, data &conf, double mu_q,
                double mass) {
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
  cout << "rank is " << rank << "; " << low << " " << high << endl;
  for (int place = low; place < high; place += 2) {
    int t = place / (x_size * y_size * z_size * 2);
    int z =
        (place - (x_size * y_size * z_size * 2) * t) / (x_size * y_size * 2);
    int y = (place - (x_size * y_size * z_size * 2) * t -
             (2 * x_size * y_size) * z) /
            (x_size * 2);
    int x = (place - (x_size * y_size * z_size * 2) * t -
             (2 * x_size * y_size) * z - 2 * x_size * y) /
            2;

    MatSetValue(A, place, place, mass, INSERT_VALUES);
    MatSetValue(A, place + 1, place + 1, mass, INSERT_VALUES);

    link_ferm.go(x, y, z, t);
    link.go(x, y, z, t);
    for (int mu = 1; mu <= 4; mu++) {
      if (mu == 4)
        delta_4 = 1;
      else
        delta_4 = 0;
      sign = eta_sign(mu, link_ferm);
      border_sign = link_ferm.border_sign(mu);
      link_ferm.move(mu, 1);
      link.move_dir(mu);
      B = link.get_matrix(conf.array);
      MatSetValue(A, place, complex_place(link_ferm),
                  (float)(exp(mu_q * delta_4) / 2 * border_sign * sign) *
                      ((float)B.a0 + (float)B.a3 * PETSC_i),
                  INSERT_VALUES);
      MatSetValue(A, place, complex_place(link_ferm) + 1,
                  (float)(exp(mu_q * delta_4) / 2 * border_sign * sign) *
                      ((float)B.a2 + (float)B.a1 * PETSC_i),
                  INSERT_VALUES);
      MatSetValue(A, place + 1, complex_place(link_ferm),
                  (float)(exp(mu_q * delta_4) / 2 * border_sign * sign) *
                      (-(float)B.a2 + (float)B.a1 * PETSC_i),
                  INSERT_VALUES);
      MatSetValue(A, place + 1, complex_place(link_ferm) + 1,
                  (float)(exp(mu_q * delta_4) / 2 * border_sign * sign) *
                      ((float)B.a0 - (float)B.a3 * PETSC_i),
                  INSERT_VALUES);
      link.move_dir(-mu);
      link_ferm.move(-mu, 1);
      border_sign = link_ferm.border_sign(-mu);
      link_ferm.move(-mu, 1);
      B = link.get_matrix(conf.array);
      MatSetValue(A, place, complex_place(link_ferm),
                  (float)(-exp(-mu_q * delta_4) / 2 * border_sign * sign) *
                      ((float)B.a0 + (float)B.a3 * PETSC_i),
                  INSERT_VALUES);
      MatSetValue(A, place, complex_place(link_ferm) + 1,
                  (float)(-exp(-mu_q * delta_4) / 2 * border_sign * sign) *
                      ((float)B.a2 + (float)B.a1 * PETSC_i),
                  INSERT_VALUES);
      MatSetValue(A, place + 1, complex_place(link_ferm),
                  (float)(-exp(-mu_q * delta_4) / 2 * border_sign * sign) *
                      (-(float)B.a2 + (float)B.a1 * PETSC_i),
                  INSERT_VALUES);
      MatSetValue(A, place + 1, complex_place(link_ferm) + 1,
                  (float)(-exp(-mu_q * delta_4) / 2 * border_sign * sign) *
                      ((float)B.a0 - (float)B.a3 * PETSC_i),
                  INSERT_VALUES);
      link_ferm.move(mu, 1);
    }
  }
}

void MatVecMult(matrix A, const PetscScalar *x, PetscScalar *y,
                int border_sign) {
  y[0] = border_sign *
         (x[0] * (A.a0 + A.a3 * PETSC_i) + x[1] * (A.a2 + A.a1 * PETSC_i));
  y[1] = border_sign *
         (x[0] * (-A.a2 + A.a1 * PETSC_i) + x[1] * (A.a0 - A.a3 * PETSC_i));
}

PetscErrorCode MatMult_eigen_sequential(Mat A, Vec vecx, Vec vecy) {
  mat_data *ctx;
  int nx, lo, i, j;
  int x_size, y_size, z_size, t_size;
  const PetscScalar *px;
  PetscScalar *py;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A, (void **)&ctx);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(vecx, &px);
  CHKERRQ(ierr);
  ierr = VecGetArray(vecy, &py);
  CHKERRQ(ierr);
  // cout<<"MatMult_Laplacian2D started"<<endl;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  PetscInt local_size;
  PetscInt low, high;
  VecGetLocalSize(vecx, &local_size);
  VecGetOwnershipRange(vecx, &low, &high);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
                          "I'm the rank %d process; low is %d; high is %d\n",
                          rank, low, high);
  PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

  double mass = ctx->mass;
  double mu_q = ctx->mu_q;
  x_size = ctx->x_size;
  y_size = ctx->y_size;
  z_size = ctx->z_size;
  t_size = ctx->t_size;
  data conf = ctx->conf;
  double delta_4 = 0;
  double sign;
  double border_sign;
  int place;
  PetscScalar vec[2];
  PetscScalar res[2];
  link1 link(x_size, y_size, z_size, t_size);
  link1 link_ferm(x_size, y_size, z_size, t_size);
  for (int x = 0; x < x_size; x++) {
    for (int y = 0; y < y_size; y++) {
      for (int z = 0; z < z_size; z++) {
        for (int t = 0; t < t_size; t++) {
          link.go(x, y, z, t);
          link_ferm.go(x, y, z, t);
          place = complex_place(link_ferm);
          for (int mu = 1; mu <= 4; mu++) {
            if (mu == 4)
              delta_4 = 1;
            else
              delta_4 = 0;
            link.move_dir(mu);
            sign = eta_sign(mu, link_ferm);
            border_sign = link_ferm.border_sign(mu);
            link_ferm.move(mu, 1);
            // matrix_mult_complex1(exp(mu_q * delta_4)/2 * sign *
            // link.get_matrix1(data_conf), &src[complex_place(link_ferm)],
            // &vec, src_index, border_sign);
            MatVecMult(exp(mu_q * delta_4) / 2 * sign *
                           link.get_matrix(conf.array),
                       &px[complex_place(link_ferm)], vec, border_sign);
            res[0] = res[0] + vec[0];
            res[1] = res[1] + vec[1];
            link.move_dir(-mu);
            link_ferm.move(-mu, 1);
            border_sign = link_ferm.border_sign(-mu);
            link_ferm.move(-mu, 1);
            // matrix_mult_complex1(exp(-mu_q * delta_4)/2 * sign *
            // link.get_matrix1(data_conf), &src[complex_place(link_ferm)],
            // &vec, src_index, border_sign);
            MatVecMult(exp(-mu_q * delta_4) / 2 * sign *
                           link.get_matrix(conf.array),
                       &px[complex_place(link_ferm)], vec, border_sign);
            res[0] = res[0] - vec[0];
            res[1] = res[1] - vec[1];
            link_ferm.move(mu, 1);
          }
          py[place] = res[0] + mass * px[complex_place(link_ferm)];
          py[place + 1] = res[1] + mass * px[complex_place(link_ferm) + 1];
          res[0] = 0;
          res[1] = 0;
        }
      }
    }
  }
  ierr = VecRestoreArrayRead(vecx, &px);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(vecy, &py);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestMatMul(mat_data &my_data, const PetscScalar *px,
                          PetscScalar *py) {
  int x_size = 32;
  int y_size = 32;
  int z_size = 32;
  int t_size = 32;
  int nx, lo, i, j;
  PetscErrorCode ierr;

  double mass = my_data.mass;
  double mu_q = my_data.mu_q;
  data conf = my_data.conf;
  double delta_4 = 0;
  double sign;
  double border_sign;
  int place;
  PetscScalar vec[2];
  PetscScalar res[2];
  link1 link(x_size, y_size, z_size, t_size);
  link1 link_ferm(x_size, y_size, z_size, t_size);
  for (int x = 0; x < x_size; x++) {
    for (int y = 0; y < y_size; y++) {
      for (int z = 0; z < z_size; z++) {
        for (int t = 0; t < t_size; t++) {
          link.go(x, y, z, t);
          link_ferm.go(x, y, z, t);
          place = complex_place(link_ferm);
          for (int mu = 1; mu <= 4; mu++) {
            if (mu == 4)
              delta_4 = 1;
            else
              delta_4 = 0;
            link.move_dir(mu);
            sign = eta_sign(mu, link_ferm);
            border_sign = link_ferm.border_sign(mu);
            link_ferm.move(mu, 1);
            MatVecMult(exp(mu_q * delta_4) / 2 * sign *
                           link.get_matrix(conf.array),
                       &px[complex_place(link_ferm)], vec, border_sign);
            res[0] = res[0] + vec[0];
            res[1] = res[1] + vec[1];
            link.move_dir(-mu);
            link_ferm.move(-mu, 1);
            border_sign = link_ferm.border_sign(-mu);
            link_ferm.move(-mu, 1);
            MatVecMult(exp(-mu_q * delta_4) / 2 * sign *
                           link.get_matrix(conf.array),
                       &px[complex_place(link_ferm)], vec, border_sign);
            res[0] = res[0] - vec[0];
            res[1] = res[1] - vec[1];
            link_ferm.move(mu, 1);
          }
          py[place] = res[0] + mass * px[complex_place(link_ferm)];
          py[place + 1] = res[1] + mass * px[complex_place(link_ferm) + 1];
          res[0] = 0;
          res[1] = 0;
        }
      }
    }
  }
}

void CheckMatMult(mat_data &my_data) {
  int x_size = 32;
  int y_size = 32;
  int z_size = 32;
  int t_size = 32;
  PetscErrorCode ierr;
  PetscInt place;
  PetscScalar a;
  PetscInt size = x_size * y_size * z_size * t_size * 2;
  PetscScalar *vecy;
  PetscScalar *vecx;
  vecy = (PetscScalar *)malloc(size * sizeof(PetscScalar));
  vecx = (PetscScalar *)malloc(size * sizeof(PetscScalar));

  for (int x = 0; x < x_size; x++) {
    for (int y = 0; y < y_size; y++) {
      for (int z = 0; z < z_size; z++) {
        for (int t = 0; t < t_size; t++) {
          place = t * 2 * x_size * y_size * z_size + z * 2 * x_size * y_size +
                  y * 2 * x_size + x * 2;
          vecx[place] = 0.1 * (x + 1) + 0.2 * (y + 1) + 0.4 * (z + 1) +
                        0.5 * (t + 1) +
                        (1 + 0.6 * (x + 1) + 0.7 * (y + 1) + 0.8 * (z + 1) +
                         0.9 * (t + 1)) *
                            PETSC_i;
          vecx[place + 1] =
              1 + 0.1 * (x + 1) + 0.2 * (y + 1) + 0.4 * (z + 1) +
              0.5 * (t + 1) +
              (0.6 * (x + 1) + 0.7 * (y + 1) + 0.8 * (z + 1) + 0.9 * (t + 1)) *
                  PETSC_i;
        }
      }
    }
  }
  TestMatMul(my_data, vecx, vecy);
  for (int i = 0; i < 10; i++) {
    cout << vecx[i] << " " << vecy[i] << endl;
  }
  // PetscFunctionReturn(0);
}

PetscErrorCode MatMult_simple(Mat A, Vec vecx, Vec vecy) {
  mat_data *ctx;
  int x_size, y_size, z_size, t_size;
  const PetscScalar *px;
  PetscScalar *py;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A, (void **)&ctx);
  CHKERRQ(ierr);
  ierr = VecGetArrayRead(vecx, &px);
  CHKERRQ(ierr);
  ierr = VecGetArray(vecy, &py);
  CHKERRQ(ierr);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  PetscInt local_size;
  PetscInt low, high;
  VecGetLocalSize(vecx, &local_size);
  VecGetOwnershipRange(vecx, &low, &high);

  // PetscSynchronizedPrintf(PETSC_COMM_WORLD,"I'm the rank %d process; low is
  // %d; high is %d\n",rank, low, high);
  // PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  for (int i = 0; i < local_size; i++) {
    py[i] = px[i] * (low + i + (low + i + 1) * PETSC_i);
  }

  ierr = VecRestoreArrayRead(vecx, &px);
  CHKERRQ(ierr);
  ierr = VecRestoreArray(vecy, &py);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatGetDiagonal_Laplacian2D(Mat A, Vec diag) {
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecSet(diag, 4.0);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   test:
      suffix: 1
      args: -n 72 -eps_nev 4 -eps_ncv 20 -terse
      requires: !single

TEST*/
