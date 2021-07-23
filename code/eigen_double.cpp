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
// PetscErrorCode MatMult_eigen_sequential(Mat A,Vec x,Vec y);
// PetscErrorCode MatGetDiagonal_Laplacian2D(Mat A,Vec diag);
// PetscErrorCode MatMult_simple(Mat A,Vec vecx,Vec vecy);
PetscErrorCode TestMatMul(mat_data &my_data, const PetscScalar *px,
                          PetscScalar *py);
void CheckMatMult(mat_data &my_data);
void output(string output_eigenvalues, string output_eigenvectors, EPS eps,
            Mat A) {
  PetscInt nconv;
  EPSGetConverged(eps, &nconv);

  PetscScalar eigenvalue;
  Vec eigenvector;
  MatCreateVecs(A, NULL, &eigenvector);

  ofstream stream_eigenvalues;
  ofstream stream_eigenvectors;

  stream_eigenvalues.open(output_eigenvalues);
  stream_eigenvectors.open(output_eigenvectors);

  stream_eigenvalues << "number,eigenvalue_real,eigenvalue_imaginary" << endl;

  PetscReal real;
  PetscReal imag;
  PetscInt vec_size;
  PetscScalar *array_eigenvector;

  for (int i = 0; i < nconv; i++) {
    EPSGetEigenpair(eps, i, &eigenvalue, NULL, eigenvector, NULL);

    real = PetscRealPart(eigenvalue);
    imag = PetscImaginaryPart(eigenvalue);
    stream_eigenvalues << i + 1 << "," << real << "," << imag << endl;

    VecGetSize(eigenvector, &vec_size);
    VecGetArray(eigenvector, &array_eigenvector);

    for (int j = 0; j < vec_size; j++) {
      real = PetscRealPart(array_eigenvector[j]);
      imag = PetscImaginaryPart(array_eigenvector[j]);
      stream_eigenvectors.write((char *)&real, sizeof(double));
      stream_eigenvectors.write((char *)&imag, sizeof(double));
      // stream_eigenvectors << real << imag << endl;
    }
  }

  VecDestroy(&eigenvector);

  stream_eigenvalues.close();
  stream_eigenvectors.close();
}

int x_size;
int y_size;
int z_size;
int t_size;

int main(int argc, char **argv) {
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  string path_conf;
  string conf_format;
  string output_eigenvalues;
  string output_eigenvectors;
  double mu_q;

  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "-path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "-output_eigenvalues") {
      output_eigenvalues = argv[++i];
    } else if (string(argv[i]) == "-output_eigenvectors") {
      output_eigenvectors = argv[++i];
    } else if (string(argv[i]) == "-mu_q") {
      mu_q = stod(string(argv[++i]));
    } else if (string(argv[i]) == "-x_size") {
      x_size = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-y_size") {
      y_size = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-z_size") {
      z_size = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-t_size") {
      t_size = stoi(string(argv[++i]));
    }
  }

  cout << "conf_format " << conf_format << endl;
  cout << "path_conf " << path_conf << endl;
  cout << "mu_q " << mu_q << endl;
  cout << "x_size " << x_size << endl;
  cout << "y_size " << y_size << endl;
  cout << "z_size " << z_size << endl;
  cout << "t_size " << t_size << endl;

  data conf_su2;

  // for ml5 configuration
  vector<float> ml5_data;
  int ml5_conf_num = 1;

  // read configuration
  if (conf_format == "double") {
    conf_su2.read_double(path_conf);
  } else if (conf_format == "double_fortran") {
    conf_su2.read_double_fortran(path_conf);
  } else if (conf_format == "double_qc2dstag") {
    conf_su2.read_double_qc2dstag(path_conf);
  } else if (conf_format == "ml5") {
    ml5_data = read_full_ml5(path_conf, ml5_conf_num);
    conf_su2.read_float_ml5(ml5_data, ml5_conf_num - 1);
  } else {
    cout << "wrong conf format: " << conf_format << endl;
    return 0;
  }

  PetscInt N, n = x_size * y_size * z_size * t_size * 2, nev;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc, &argv, (char *)0, help);
  if (ierr)
    return ierr;

  // set dimension of the problem
  N = n * n;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nDirac Eigenproblem %Dx%D grid\n\n", n,
                     n);
  CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  Mat A; /* operator matrix */
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetFromOptions(A);

  MatSetSizes(A, n, n, n, n);
  // PetscInt *d_nnz;
  // PetscInt *o_nnz;
  // if (!(d_nnz = (PetscInt *)malloc(n * sizeof(PetscInt))))
  //   PetscPrintf(PETSC_COMM_WORLD, "err malloc d_nnz");
  // if (!(o_nnz = (PetscInt *)malloc(n * sizeof(PetscInt))))
  //   PetscPrintf(PETSC_COMM_WORLD, "err malloc o_nnz");
  // for (int i = 0; i < n; i++) {
  //   d_nnz[i] = 0;
  //   o_nnz[i] = 0;
  // }
  // mat_set_index(d_nnz, o_nnz, 0, n);
  MatSeqAIJSetPreallocation(A, 16, NULL);
  // MatMPIAIJSetPreallocation(A, 16, NULL, 0, NULL);
  // MatSetValue(A, 0, 1, 0.123 - 1.2134 * PETSC_i, INSERT_VALUES);
  mat_insert(A, 0, n, conf_su2, mu_q, 0);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  // free(d_nnz);
  // free(o_nnz);

  Vec vec_test_init;
  Vec vec_test_final;

  PetscScalar *array;

  VecCreateSeq(PETSC_COMM_SELF, n, &vec_test_init);
  VecCreateSeq(PETSC_COMM_SELF, n, &vec_test_final);

  VecSetType(vec_test_init, VECSEQCUDA);
  VecSetType(vec_test_final, VECSEQCUDA);

  // VecSetSizes(vec_test_init, PETSC_DECIDE, n);
  // VecSetSizes(vec_test_final, PETSC_DECIDE, n);

  VecSetUp(vec_test_init);
  VecSetUp(vec_test_final);

  VecGetArray(vec_test_init, &array);
  int place;
  for (int x = 0; x < x_size; x++) {
    for (int y = 0; y < y_size; y++) {
      for (int z = 0; z < z_size; z++) {
        for (int t = 0; t < t_size; t++) {
          place = t * 2 * x_size * y_size * z_size + z * 2 * x_size * y_size +
                  y * 2 * x_size + x * 2;
          array[place] = 0.1 * (x + 1) + 0.2 * (y + 1) + 0.4 * (z + 1) +
                         0.5 * (t + 1) +
                         (1 + 0.6 * (x + 1) + 0.7 * (y + 1) + 0.8 * (z + 1) +
                          0.9 * (t + 1)) *
                             PETSC_i;
          array[place + 1] =
              1 + 0.1 * (x + 1) + 0.2 * (y + 1) + 0.4 * (z + 1) +
              0.5 * (t + 1) +
              (0.6 * (x + 1) + 0.7 * (y + 1) + 0.8 * (z + 1) + 0.9 * (t + 1)) *
                  PETSC_i;
        }
      }
    }
  }
  VecRestoreArray(vec_test_init, &array);

  MatMult(A, vec_test_init, vec_test_final);

  VecGetArray(vec_test_final, &array);

  for (int i = 0; i < 10; i++) {
    cout << array[i] << endl;
  }

  PetscReal vec_norm = 0;
  PetscReal real_part;
  PetscReal imaginary_part;

  for (int i = 0; i < n; i++) {
    real_part = PetscRealPart(array[i]);
    imaginary_part = PetscImaginaryPart(array[i]);
    vec_norm += real_part * real_part + imaginary_part * imaginary_part;
  }
  cout << "vec_norm " << vec_norm << endl;

  VecDestroy(&vec_test_init);
  VecDestroy(&vec_test_final);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  EPS eps; /* eigenproblem solver context */
  ierr = EPSCreate(PETSC_COMM_WORLD, &eps);
  CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps, A, NULL);
  CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_NHEP);
  CHKERRQ(ierr);
  //   ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);CHKERRQ(ierr);

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
  // cout << "rank is " << rank << "; starting to solve" << endl;
  start_time = clock();

  ierr = EPSSolve(eps);
  CHKERRQ(ierr);

  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "EPS solve time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;

  /*
     Optional: Get some information from the solver and display it
  */
  EPSType type;
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
  PetscBool terse;
  ierr = PetscOptionsHasName(NULL, NULL, "-terse", &terse);
  CHKERRQ(ierr);
  if (terse) {
    ierr = EPSErrorView(eps, EPS_ERROR_RELATIVE, NULL);
    CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,
                                 PETSC_VIEWER_ASCII_INFO_DETAIL);
    CHKERRQ(ierr);
    ierr = EPSErrorView(eps, EPS_ERROR_RELATIVE, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);
  }

  start_time = clock();
  output(output_eigenvalues, output_eigenvectors, eps, A);
  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "eigenvectors output time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;

  ierr = EPSDestroy(&eps);
  CHKERRQ(ierr);
  ierr = MatDestroy(&A);
  CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

// set indexes for sparce dirac operator matrix
void mat_set_index(PetscInt *d_nnz, PetscInt *o_nnz, int low, int high) {
  int x, y, z, t;
  link1 link_ferm(x_size, y_size, z_size, t_size);
  for (int place = low; place < high; place += 2) {
    t = place / (x_size * y_size * z_size * 2);
    z = (place - (x_size * y_size * z_size * 2) * t) / (x_size * y_size * 2);
    y = (place - (x_size * y_size * z_size * 2) * t -
         (2 * x_size * y_size) * z) /
        (x_size * 2);
    x = (place - (x_size * y_size * z_size * 2) * t -
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

// insert values into matrix
void mat_insert(Mat A, int low, int high, data &conf, double mu_q,
                double mass) {
  // int rank, mpi_size;
  // MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
  // MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  int x, y, z, t;
  link1 link(x_size, y_size, z_size, t_size);
  link1 link_ferm(x_size, y_size, z_size, t_size);
  double delta_4 = 0;
  double sign;
  double border_sign;
  matrix B;
  // cout << "rank is " << rank << "; " << low << " " << high << endl;
  for (int place = low; place < high; place += 2) {
    t = place / (x_size * y_size * z_size * 2);
    z = (place - (x_size * y_size * z_size * 2) * t) / (x_size * y_size * 2);
    y = (place - (x_size * y_size * z_size * 2) * t -
         (2 * x_size * y_size) * z) /
        (x_size * 2);
    x = (place - (x_size * y_size * z_size * 2) * t -
         (2 * x_size * y_size) * z - 2 * x_size * y) /
        2;

    // MatSetValue(A, place, place, mass, INSERT_VALUES);
    // MatSetValue(A, place + 1, place + 1, mass, INSERT_VALUES);

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
                  exp(mu_q * delta_4) / 2 * border_sign * sign *
                      (B.a0 + B.a3 * PETSC_i),
                  INSERT_VALUES);
      MatSetValue(A, place, complex_place(link_ferm) + 1,
                  exp(mu_q * delta_4) / 2 * border_sign * sign *
                      (B.a2 + B.a1 * PETSC_i),
                  INSERT_VALUES);
      MatSetValue(A, place + 1, complex_place(link_ferm),
                  exp(mu_q * delta_4) / 2 * border_sign * sign *
                      (-B.a2 + B.a1 * PETSC_i),
                  INSERT_VALUES);
      MatSetValue(A, place + 1, complex_place(link_ferm) + 1,
                  exp(mu_q * delta_4) / 2 * border_sign * sign *
                      (B.a0 - B.a3 * PETSC_i),
                  INSERT_VALUES);
      link.move_dir(-mu);
      link_ferm.move(-mu, 1);
      border_sign = link_ferm.border_sign(-mu);
      link_ferm.move(-mu, 1);
      B = link.get_matrix(conf.array);
      MatSetValue(A, place, complex_place(link_ferm),
                  -exp(-mu_q * delta_4) / 2 * border_sign * sign *
                      (B.a0 + B.a3 * PETSC_i),
                  INSERT_VALUES);
      MatSetValue(A, place, complex_place(link_ferm) + 1,
                  -exp(-mu_q * delta_4) / 2 * border_sign * sign *
                      (B.a2 + B.a1 * PETSC_i),
                  INSERT_VALUES);
      MatSetValue(A, place + 1, complex_place(link_ferm),
                  -exp(-mu_q * delta_4) / 2 * border_sign * sign *
                      (-B.a2 + B.a1 * PETSC_i),
                  INSERT_VALUES);
      MatSetValue(A, place + 1, complex_place(link_ferm) + 1,
                  -exp(-mu_q * delta_4) / 2 * border_sign * sign *
                      (B.a0 - B.a3 * PETSC_i),
                  INSERT_VALUES);
      link_ferm.move(mu, 1);
    }
  }
}

void MatVecMult(matrix A, const PetscScalar *x, PetscScalar *y,
                int border_sign) {
  y[0] = (double)border_sign *
         (x[0] * (A.a0 + A.a3 * PETSC_i) + x[1] * (A.a2 + A.a1 * PETSC_i));
  y[1] = (double)border_sign *
         (x[0] * (-A.a2 + A.a1 * PETSC_i) + x[1] * (A.a0 - A.a3 * PETSC_i));
}

PetscErrorCode TestMatMul(mat_data &my_data, const PetscScalar *px,
                          PetscScalar *py) {
  int x_size = 40;
  int y_size = 40;
  int z_size = 40;
  int t_size = 40;
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
  int x_size = 40;
  int y_size = 40;
  int z_size = 40;
  int t_size = 40;
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
}
