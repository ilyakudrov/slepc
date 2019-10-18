/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2019, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Solves the same eigenproblem as in example ex2, but using a shell matrix. "
  "The problem is a standard symmetric eigenproblem corresponding to the 2-D Laplacian operator.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions in both x and y dimensions.\n\n";

#include <slepceps.h>
#include <iostream>
#include "data.h"
#include "link.h"
#include "matrix.h"
#include "eigen.h"

/*
   User-defined routines
*/
typedef struct {
   data conf;
   /*int x_size;
   int y_size;
   int z_size;
   int t_size;*/
   double mass;
   double mu_q;
} mat_data;

PetscErrorCode MatMult_Laplacian2D(Mat A,Vec x,Vec y);
PetscErrorCode MatGetDiagonal_Laplacian2D(Mat A,Vec diag);
PetscErrorCode TestMatMul(mat_data& my_data, const PetscScalar* px, PetscScalar* py);
void CheckMatMult(mat_data& my_data);

int x_size = 32;
int y_size = 32;
int z_size = 32;
int t_size = 32;

int main(int argc,char **argv)
{
   //int x_size = 32/*atof(argv[1])*/;
   //int y_size = 32/*atof(argv[2])*/;
   //int z_size = 32/*atof(argv[3])*/;
   //int t_size = 32/*atof(argv[4])*/;
   mat_data my_data;
   /*my_data.x_size = atof(argv[1]);
   my_data.y_size = atof(argv[2]);
   my_data.z_size = atof(argv[3]);
   my_data.t_size = atof(argv[4]);*/
   my_data.conf.read_float("/home/ilya/lattice/slepc/conf/nosmeared/time_32/mu0.10/conf_0001.fl"/*argv[5]*/);
   for(int i = 0;i < 10;i++){
      cout<<"matrix "<<my_data.conf.array[i].a0<<endl;
   }
   my_data.mass = 0.0075/*atof(argv[6])*/;
   my_data.mu_q = 0.1/*atof(argv[7])*/;
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  PetscMPIInt    size;
  PetscInt       N,n=x_size*y_size*z_size*t_size*2,nev;
  PetscBool      terse;
  PetscErrorCode ierr;

   //CheckMatMult();

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only");

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  N = n*n;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nDirac Eigenproblem (matrix-free version), N=%D (%Dx%D grid)\n\n",N,n,n);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Create the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreateShell(PETSC_COMM_WORLD,N,N,N,N,&my_data,&A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void))MatMult_Laplacian2D);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMult_Laplacian2D);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Laplacian2D);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  //TESTING
  CheckMatMult(my_data);
  cout<<"check"<<endl;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* show detailed info unless -terse option is given by user */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
  if (terse) {
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*
    Compute the matrix vector multiplication y<---T*x where T is a nx by nx
    tridiagonal matrix with DD on the diagonal, DL on the subdiagonal, and
    DU on the superdiagonal.
 */
static void tv(int nx,const PetscScalar *x,PetscScalar *y)
{
  PetscScalar dd,dl,du;
  int         j;

  dd  = 4.0;
  dl  = -1.0;
  du  = -1.0;

  y[0] =  dd*x[0] + du*x[1];
  for (j=1;j<nx-1;j++)
    y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
  y[nx-1] = dl*x[nx-2] + dd*x[nx-1];
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

void MatVecMult(matrix A, const PetscScalar* x, PetscScalar* y, int border_sign){
   PetscScalar z1[2], z2[2];
   z1[0] = A.a0 + A.a3*PETSC_i;
   z2[0] = A.a2 + A.a1*PETSC_i;
   z1[1] = -A.a2 + A.a1*PETSC_i;
   z2[1] = A.a0 - A.a3*PETSC_i;
   y[0] = border_sign * (x[0] * z1[0] + x[1] * z2[0]);
   y[1] = border_sign * (x[0] * z1[1] + x[1] * z2[1]);
}

PetscErrorCode MatMult_Laplacian2D(Mat A,Vec vecx,Vec vecy)
{
  void              *ctx1;
  int               nx,lo,i,j;
  int x_size, y_size, z_size, t_size;
  const PetscScalar *px;
  PetscScalar       *py;
  PetscErrorCode    ierr;

  PetscFunctionBeginUser;
  ierr = MatShellGetContext(A,&ctx1);CHKERRQ(ierr);
  ierr = VecGetArrayRead(vecx,&px);CHKERRQ(ierr);
  ierr = VecGetArray(vecy,&py);CHKERRQ(ierr);

   mat_data *ctx = (mat_data*)ctx1;
   int mass = ctx->mass;
   int mu_q = ctx->mu_q;
   data conf = ctx->conf;
   double delta_4 = 0;
   double sign;
   double border_sign;
   int place;
   PetscScalar vec[2];
   PetscScalar res[2];
   link1 link(x_size, y_size, z_size, t_size);
   link1 link_ferm(x_size, y_size, z_size, t_size);
   for(int x = 0;x < x_size;x++){
      for(int y = 0;y < y_size;y++){
         for(int z = 0;z < z_size;z++){
            for(int t = 0;t < t_size;t++){
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
               py[place] = res[0] + mass*px[complex_place(link_ferm)];
               py[place + 1] = res[1] + mass*px[complex_place(link_ferm) + 1];
               res[0] = 0;
               res[1] = 0;
            }
         }
      }
   }
  ierr = VecRestoreArrayRead(vecx,&px);CHKERRQ(ierr);
  ierr = VecRestoreArray(vecy,&py);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestMatMul(mat_data& my_data, const PetscScalar* px, PetscScalar* py){
   int x_size = 32;
   int y_size = 32;
   int z_size = 32;
   int t_size = 32;
   //void *ctx1;
   int nx,lo,i,j;
   //int x_size, y_size, z_size, t_size;
   //const PetscScalar *px;
   //PetscScalar *py;
   PetscErrorCode ierr;

   //ierr = MatShellGetContext(A,&ctx1);CHKERRQ(ierr);
   //mat_data *ctx = (mat_data*)ctx1;
   double mass = my_data.mass;
   double mu_q = my_data.mu_q;
   cout<<"test matmul mu "<<mu_q<<endl;
   data conf = my_data.conf;
   for(int i = 0;i < 10;i++){
      cout<<"matrix "<<conf.array[i].a0<<endl;
   }
   double delta_4 = 0;
   double sign;
   double border_sign;
   int place;
   PetscScalar vec[2];
   PetscScalar res[2];
   link1 link(x_size, y_size, z_size, t_size);
   link1 link_ferm(x_size, y_size, z_size, t_size);
   cout<<"ok matmul"<<endl;
   for(int x = 0;x < x_size;x++){
      //cout<<"x "<<x<<endl;
      for(int y = 0;y < y_size;y++){
         for(int z = 0;z < z_size;z++){
            for(int t = 0;t < t_size;t++){
               //cout<<"t "<<t<<endl;
               link.go(x, y, z, t);
               link_ferm.go(x, y, z, t);
               place = complex_place(link_ferm);
               for(int mu = 1;mu <= 4;mu++){
                  if(mu == 4) delta_4 = 1;
                  else delta_4 = 0;
                  link.move_dir(mu);
                  sign = eta_sign(mu, link_ferm);
                  border_sign = link_ferm.border_sign(mu);
                  //if(place == 0) cout<<border_sign<<" "<<link_ferm.coordinate[3]<<endl;
                  link_ferm.move(mu, 1);
                  if(mu == 4 && place == 0) cout<<exp(mu_q * delta_4)/2<<endl;
                  // matrix_mult_complex1(exp(mu_q * delta_4)/2 * sign * link.get_matrix1(data_conf), &src[complex_place(link_ferm)], &vec, src_index, border_sign);
                  MatVecMult(exp(mu_q * delta_4)/2 * sign * link.get_matrix(conf.array), &px[complex_place(link_ferm)], vec, border_sign);
                  if(place == 0) cout<<"vec "<<vec[0]<<" "<<mu<<endl;
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
               py[place] = res[0] + mass*px[complex_place(link_ferm)];
               py[place + 1] = res[1] + mass*px[complex_place(link_ferm) + 1];
               //if(place == 0) cout<<py[place]<<endl;
               res[0] = 0;
               res[1] = 0;
            }
         }
      }
   }
}

void CheckMatMult(mat_data& my_data){
   int x_size = 32;
   int y_size = 32;
   int z_size = 32;
   int t_size = 32;
   cout<<"chaeck mu "<<my_data.mu_q<<endl;
   PetscErrorCode ierr;
   PetscInt place;
   PetscScalar a;
   PetscInt size = x_size*y_size*z_size*t_size*2;
   PetscScalar *vecy;
   PetscScalar *vecx;
   vecy = (PetscScalar*) malloc(size * sizeof(PetscScalar));
   vecx = (PetscScalar*) malloc(size * sizeof(PetscScalar));

   for(int x = 0;x < x_size;x++){
      for(int y = 0;y < y_size;y++){
         for(int z = 0;z < z_size;z++){
            for(int t = 0;t < t_size;t++){
               place = t * 2 *  x_size * y_size * z_size +
               z * 2 * x_size * y_size +
               y * 2 * x_size +
               x * 2;
               //a = 0.1 * x + 0.2 * y + 0.4 * z + 0.5 * t + (1 + 0.6 * x + 0.7 * y + 0.8 * z + 0.9 * t) * PETSC_i;
               //ierr = VecSetValue(vecy, place, a, INSERT_VALUES); CHKERRQ(ierr);
               vecx[place] = 0.1 * (x+1) + 0.2 * (y+1) + 0.4 * (z+1) + 0.5 * (t+1) + (1 + 0.6 * (x+1) + 0.7 * (y+1) + 0.8 * (z+1) + 0.9 * (t+1)) * PETSC_i;
               //a = 1 + 0.1 * x + 0.2 * y + 0.4 * z + 0.5 * t + (0.6 * x + 0.7 * y + 0.8 * z + 0.9 * t) * PETSC_i;
               //ierr = VecSetValue(vecy, place + 1, a, INSERT_VALUES); CHKERRQ(ierr);
               vecx[place+1] = 1 + 0.1 * (x+1) + 0.2 * (y+1) + 0.4 * (z+1) + 0.5 * (t+1) + (0.6 * (x+1) + 0.7 * (y+1) + 0.8 * (z+1) + 0.9 * (t+1)) * PETSC_i;
            }
         }
      }
   }
   TestMatMul(my_data, vecx, vecy);
   cout<<"ok TestMatMul"<<endl;
   for(int i = 0;i < 10;i++){
      cout<<vecy[i]<<endl;
   }
   //PetscFunctionReturn(0);
}

PetscErrorCode MatGetDiagonal_Laplacian2D(Mat A,Vec diag)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecSet(diag,4.0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   test:
      suffix: 1
      args: -n 72 -eps_nev 4 -eps_ncv 20 -terse
      requires: !single

TEST*/
