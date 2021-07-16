// Data class

#ifndef __DATA_H__
#define __DATA_H__

#include "math.h"
#include "matrix.h"
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class data {
public:
  std::vector<matrix> array;
  data();
  void read_float(string &file_name); // read conf file of floats
  void read_float_abelian(string &file_name);
  void read_double(string &file_name); // read conf file of double
  void read_float_ml5(const vector<float> &array_ml5, int conf_num);
  void read_double_abelian(string &file_name);
  void read_double_qc2dstag(string &file_name);
  void read_double_fortran(string &file_name);
  void write_double(string &file_name); // writes in file
  void write_float(string &file_name);
  int eta(int mu, int x, int y, int z, int t);
  int sign(int x);
};

vector<float> read_full_ml5(string &file_name, int conf_num);

#endif
