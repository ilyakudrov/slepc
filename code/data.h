//Data class

#ifndef __DATA_H__
#define __DATA_H__

#include <iostream>
#include <fstream>
#include <vector>
#include "matrix.h"
#include "math.h"

class data {
	public:
	std::vector<matrix> array;
	data();
	void read_float(char file_name[]);//read conf file of floats
	void read_float_abelian(char file_name[]);
	void read_double(char file_name[]);//read conf file of double
	void read_double_abelian(char file_name[]);
	void write_double(char file_name[]);//writes in file
	void write_float(char file_name[]);
	int eta(int mu, int x, int y, int z, int t);
	int sign(int x);
};
#endif
