#include "data.h"
#include "link.h"

#define SPACE_ITER_START                                                       \
  for (int t = 0; t < t_size; t++) {                                           \
    for (int z = 0; z < z_size; z++) {                                         \
      for (int y = 0; y < y_size; y++) {                                       \
        for (int x = 0; x < x_size; x++) {
#define SPACE_ITER_END                                                         \
  }                                                                            \
  }                                                                            \
  }                                                                            \
  }

#define PLACE_QC2DSTAG                                                         \
  (dir - 1) * x_size *y_size *z_size *t_size * 4 +                             \
      (t)*x_size *y_size *z_size * 4 + (z)*x_size *y_size * 4 +                \
      (y)*x_size * 4 + (x)*4

using namespace std;

data::data() { array.reserve(4 * x_size * y_size * z_size * t_size); }
void data::read_float(string &file_name) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<float> v(data_size * 4);
  if (!stream.read((char *)&v[0], data_size * 4 * sizeof(float)))
    cout << "read_float error: " << file_name << endl;
  matrix A;
  for (int i = 0; i < data_size; i++) {
    A.a0 = (double)v[i * 4];
    A.a1 = (double)v[i * 4 + 1];
    A.a2 = (double)v[i * 4 + 2];
    A.a3 = (double)v[i * 4 + 3];
    array.push_back(A);
  }
  stream.close();
}
void data::read_float_abelian(string &file_name) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<float> v(data_size * 4 + 2);
  if (!stream.read((char *)&v[0], (data_size * 4 + 2) * sizeof(float)))
    cout << "read_float error: " << file_name << endl;
  matrix A;
  for (int i = 0; i < data_size; i++) {
    A.a0 = (double)v[i * 4 + 1];
    A.a1 = (double)v[i * 4 + 2];
    A.a2 = (double)v[i * 4 + 3];
    A.a3 = (double)v[i * 4 + 4];
    array.push_back(A);
  }
  stream.close();
}
void data::read_double(string &file_name) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<double> v(data_size * 4);
  if (!stream.read((char *)&v[0], data_size * 4 * sizeof(double)))
    cout << "read_double error: " << file_name << endl;
  matrix A;
  for (int i = 0; i < data_size; i++) {
    A.a0 = v[i * 4];
    A.a1 = v[i * 4 + 1];
    A.a2 = v[i * 4 + 2];
    A.a3 = v[i * 4 + 3];
    array.push_back(A);
  }
  stream.close();
}

vector<float> read_full_ml5(string &file_name, int conf_num) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  ifstream stream(file_name);
  vector<float> v(conf_num * data_size1 * 4);
  if (!stream.read((char *)&v[0], conf_num * data_size1 * 4 * sizeof(float)))
    cout << "read_full_ml5 error: " << file_name << endl;
  return v;
  stream.close();
}

void data::read_float_ml5(const vector<float> &array_ml5, int conf_num) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  matrix A;
  int i;
  SPACE_ITER_START
  for (int mu = 1; mu <= 4; mu++) {
    if (mu == 4) {
      i = 4 *
          (t * x_size * y_size * z_size + x + y * x_size + z * x_size * y_size);
    } else {
      i = 4 * (t * x_size * y_size * z_size + x + y * x_size +
               z * x_size * y_size) +
          mu;
    }
    A.a0 = array_ml5[data_size1 * 4 * conf_num + i * 4];
    A.a3 = array_ml5[data_size1 * 4 * conf_num + i * 4 + 1];
    A.a2 = array_ml5[data_size1 * 4 * conf_num + i * 4 + 2];
    A.a1 = array_ml5[data_size1 * 4 * conf_num + i * 4 + 3];
    array.push_back(A);
  }
  SPACE_ITER_END
}

void data::read_double_fortran(string &file_name) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<double> v(data_size * 4);
  stream.ignore(4);
  if (!stream.read((char *)&v[0], (data_size * 4) * sizeof(double)))
    cout << "read_double_fortran<su2> error: " << file_name << endl;
  matrix A;
  for (int i = 0; i < data_size; i++) {
    A.a0 = v[i * 4];
    A.a1 = v[i * 4 + 1];
    A.a2 = v[i * 4 + 2];
    A.a3 = v[i * 4 + 3];
    array.push_back(A);
  }
  stream.close();
}

void data::read_double_abelian(string &file_name) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<double> v(data_size * 4 + 2);
  if (!stream.read((char *)&v[0], (data_size * 4 + 2) * sizeof(double)))
    cout << "read_float error: " << file_name << endl;
  matrix A;
  for (int i = 0; i < data_size; i++) {
    A.a0 = v[i * 4 + 1];
    A.a1 = v[i * 4 + 2];
    A.a2 = v[i * 4 + 3];
    A.a3 = v[i * 4 + 4];
    array.push_back(A);
  }
  stream.close();
}
void data::read_double_qc2dstag(string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<double> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], data_size1 * 4 * sizeof(double)))
    cout << "read_double_qc2dstag<su2> error: " << file_name << endl;
  matrix A;
  int dir;
  SPACE_ITER_START
  for (int dir1 = 1; dir1 <= 4; dir1++) {
    if (dir1 == 4)
      dir = 1;
    else
      dir = dir1 + 1;
    A.a0 = v[PLACE_QC2DSTAG + 0];
    A.a3 = v[PLACE_QC2DSTAG + 1];
    A.a2 = v[PLACE_QC2DSTAG + 2];
    A.a1 = v[PLACE_QC2DSTAG + 3];
    array.push_back(A);
  }
  SPACE_ITER_END
  stream.close();
}
void data::write_float(string &file_name) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  ofstream stream(file_name);
  vector<float> v(data_size * 4);
  for (int i = 0; i < data_size; i++) {
    v[i * 4] = (float)array[i].a0;
    v[i * 4 + 1] = (float)array[i].a1;
    v[i * 4 + 2] = (float)array[i].a2;
    v[i * 4 + 3] = (float)array[i].a3;
  }
  if (!stream.write((char *)&v[0], data_size * 4 * sizeof(float)))
    cout << "write_float error: " << file_name << endl;
  stream.close();
}
void data::write_double(string &file_name) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  ofstream stream(file_name);
  vector<double> v(data_size * 4);
  for (int i = 0; i < data_size; i++) {
    v[i * 4] = array[i].a0;
    v[i * 4 + 1] = array[i].a1;
    v[i * 4 + 2] = array[i].a2;
    v[i * 4 + 3] = array[i].a3;
  }
  if (!stream.write((char *)&v[0], data_size * 4 * sizeof(double)))
    cout << "write_double error: " << file_name << endl;
  stream.close();
}
int data::eta(int mu, int x, int y, int z, int t) {
  int a = 0;
  int coordinate[4] = {x, y, z, t};
  for (int i = 0; i < mu - 1; i++) {
    a += coordinate[i];
  }
  return sign(a);
}
int data::sign(int x) {
  if (x % 2 == 0)
    return 1;
  if (x % 2 == 1)
    return -1;
}
