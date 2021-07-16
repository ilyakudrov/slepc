#include "data.h"
#include "eigen.h"
#include "link.h"
#include "matrix.h"
#include <iostream>

int x_size = 40;
int y_size = 40;
int z_size = 40;
int t_size = 40;

int main(int argc, char **argv) {
  x_size = 40;
  y_size = 40;
  z_size = 40;
  t_size = 40;
  data conf;
  string conf_path = "../confs/qc2dstag/40^4/mu0.05/s0/CONF0201";
  conf.read_double_qc2dstag(conf_path);
  double mu = 0.05;
  double mass = 0;

  int vec_size = 2 * x_size * y_size * z_size * t_size;

  vector<scomplex_t> vec_init(vec_size);
  vector<scomplex_t> vec_final(vec_size);

  // for (int i = 0; i < vec_size; i++) {
  //   vec_init[i].im = 1;
  //   vec_init[i].re = 1;
  //   vec_final[i].im = 1;
  //   vec_final[i].re = 1;
  // }

  int place;
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          place = t * 2 * x_size * y_size * z_size + z * 2 * x_size * y_size +
                  y * 2 * x_size + x * 2;
          vec_init[place].re =
              0.1 * (x + 1) + 0.2 * (y + 1) + 0.4 * (z + 1) + 0.5 * (t + 1);
          vec_init[place].im =
              1 + 0.6 * (x + 1) + 0.7 * (y + 1) + 0.8 * (z + 1) + 0.9 * (t + 1);
          vec_init[place + 1].re =
              1 + 0.1 * (x + 1) + 0.2 * (y + 1) + 0.4 * (z + 1) + 0.5 * (t + 1);
          vec_init[place + 1].im =
              0.6 * (x + 1) + 0.7 * (y + 1) + 0.8 * (z + 1) + 0.9 * (t + 1);
        }
      }
    }
  }

  scomplex_t res;

  // int sign;
  // int border_sign;
  // link1 link(x_size, y_size, z_size, t_size);
  // link1 link_ferm(x_size, y_size, z_size, t_size);
  // link.go(0, 0, 0, 0);
  // link_ferm.go(0, 0, 0, 0);
  // scomplex_t res_test;
  // res_test.re = 0;
  // res_test.im = 0;

  // link.move_dir(1);
  // sign = eta_sign(mu, link_ferm);
  // border_sign = link_ferm.border_sign(mu);

  matrix *matrix_vec;
  matrix_vec = (matrix *)calloc(vec_size * 2, sizeof(matrix));
  for (int i = 0; i < vec_size * 2; i++) {
    matrix_vec[i] = conf.array[i];
  }

  // matrix_mult_complex1(exp(mu * 0.) / 2 * sign *
  // link.get_matrix1(matrix_vec),
  //                      &vec_init[complex_place(link_ferm)], &res_test, 0,
  //                      border_sign);

  // cout << res_test.re << " " << res_test.im << endl;

  for (int i = 0; i < vec_size; i++) {
    res.re = 0;
    res.im = 0;
    dirac_mult(&res, &vec_init[0], i, matrix_vec, x_size, y_size, z_size,
               t_size, mass, mu);
    vec_final[i].re = res.re;
    vec_final[i].im = res.im;
  }

  cout.precision(17);

  for (int i = 0; i < 10; i++) {
    cout << vec_final[i].re << " + (" << vec_final[i].im << ")i" << endl;
  }
  double norm = 0;
  for (int i = 0; i < vec_size; i++) {
    norm +=
        vec_final[i].re * vec_final[i].re + vec_final[i].im * vec_final[i].im;
  }
  cout << norm << endl;
}