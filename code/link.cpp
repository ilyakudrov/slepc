# define Pi 3.141592653589793238462643383279502884
#define data_size 4*lattice_size[0]*lattice_size[1]*lattice_size[2]*lattice_size[3]
#define PLACE (coordinate[3])*3*lattice_size[0]*lattice_size[1]*lattice_size[2] \
		+ (coordinate[2])*3*lattice_size[0]*lattice_size[1] \
		+ (coordinate[1])*3*lattice_size[0] \
		+ (coordinate[0])*3+abs(direction)-1
#define PLACE1_LINK (coordinate[3])*4*lattice_size[0]*lattice_size[1]*lattice_size[2] \
		+ (coordinate[2])*4*lattice_size[0]*lattice_size[1] \
		+ (coordinate[1])*4*lattice_size[0] \
		+ (coordinate[0])*4+abs(direction)-1
#define PLACE_FIELD_LINK (coordinate[3]-1)*lattice_size[0]*lattice_size[1]*lattice_size[2] \
		+ (coordinate[2])*lattice_size[0]*lattice_size[1] \
		+ (coordinate[1])*lattice_size[0] \
		+ (coordinate[0])
		
#include "link.h"
#include "data.h"

link1::link1(int lattice_size_x, int lattice_size_y, int lattice_size_z, int lattice_size_t) {
	lattice_size[0] = lattice_size_x;
	lattice_size[1] = lattice_size_y;
	lattice_size[2] = lattice_size_z;
	lattice_size[3] = lattice_size_t;
	direction = 1;
	for (int i = 0; i < 4; i++) {
		coordinate[i] = 0;
	}
}

link1::link1(int dir, int lattice_size_x, int lattice_size_y, int lattice_size_z, int lattice_size_t) {
	lattice_size[0] = lattice_size_x;
	lattice_size[1] = lattice_size_y;
	lattice_size[2] = lattice_size_z;
	lattice_size[3] = lattice_size_t;
	direction = dir;
	for (int i = 0; i < 4; i++) {
		coordinate[i] = 0;
	}
}

link1::link1(int coord[4], int dir, int lattice_size_x, int lattice_size_y, int lattice_size_z, int lattice_size_t) {
	lattice_size[0] = lattice_size_x;
	lattice_size[1] = lattice_size_y;
	lattice_size[2] = lattice_size_z;
	lattice_size[3] = lattice_size_t;
	direction = dir;
	for (int i = 0; i < 4; i++) {
		coordinate[i] = coord[i];
	}
}

link1::link1(int x, int y, int z, int t, int dir, int lattice_size_x, int lattice_size_y, int lattice_size_z, int lattice_size_t) {
	lattice_size[0] = lattice_size_x;
	lattice_size[1] = lattice_size_y;
	lattice_size[2] = lattice_size_z;
	lattice_size[3] = lattice_size_t;
	direction = dir;
	coordinate[0] = x;
	coordinate[1] = y;
	coordinate[2] = z;
	coordinate[3] = t;
}

link1::link1(int x, int y, int z, int t, int lattice_size_x, int lattice_size_y, int lattice_size_z, int lattice_size_t) {
	lattice_size[0] = lattice_size_x;
	lattice_size[1] = lattice_size_y;
	lattice_size[2] = lattice_size_z;
	lattice_size[3] = lattice_size_t;
	direction = 1;
	coordinate[0] = x;
	coordinate[1] = y;
	coordinate[2] = z;
	coordinate[3] = t;
}

link1::link1(const link1& link) {
	lattice_size[0] = link.lattice_size[0];
	lattice_size[1] = link.lattice_size[1];
	lattice_size[2] = link.lattice_size[2];
	lattice_size[3] = link.lattice_size[3];
	direction = link.direction;
	coordinate[0] = link.coordinate[0];
	coordinate[1] = link.coordinate[1];
	coordinate[2] = link.coordinate[2];
	coordinate[3] = link.coordinate[3];
}

link1::link1() {
	lattice_size[0] = x_size;
	lattice_size[1] = y_size;
	lattice_size[2] = z_size;
	lattice_size[3] = t_size;
	direction = 1;
	coordinate[0] = 0;
	coordinate[1] = 0;
	coordinate[2] = 0;
	coordinate[3] = 0;
}

void link1::print_link(){
	cout<<coordinate[0]<<" "<<coordinate[1]<<" "<<coordinate[2]<<" "<<coordinate[3]<<" "<<direction<<endl;
}

void link1::move(int dir, int step) {
	coordinate[abs(dir) - 1] += step * dir / abs(dir);
	while (coordinate[abs(dir) - 1] < 0 || coordinate[abs(dir) - 1] > lattice_size[abs(dir) - 1]-1) {
		if (coordinate[abs(dir) - 1] < 0) coordinate[abs(dir) - 1] += lattice_size[abs(dir) - 1];
		else coordinate[abs(dir) - 1] -= lattice_size[abs(dir) - 1];
	}
}

void link1::go(int x, int y, int z, int t) {
	coordinate[0] = x;
	coordinate[1] = y;
	coordinate[2] = z;
	coordinate[3] = t;
}

void link1::move_dir(int dir) {
	direction = dir;
}

int link1::get_place(){
	return PLACE1_LINK;
}

int link1::get_place1(){
	return PLACE_FIELD_LINK;
}

matrix link1::get_matrix(const vector<matrix>& vec) {
	matrix A;
	if (direction / abs(direction) == 1) {
		A = vec[PLACE1_LINK];
		return A;
	}
	if (direction / abs(direction) == -1) {
		move(direction, 1);
		A = vec[PLACE1_LINK];
		move(direction, -1);
		return A.conj();
	}
	else {
		cout << "error in get_matrix" << endl;
		return A;
	}
}

matrix link1::get_matrix1(matrix* vec) {
        matrix A;
        if (direction / abs(direction) == 1) {
                A = vec[PLACE1_LINK];
                return A;
        }
        if (direction / abs(direction) == -1) {
                move(direction, 1);
                A = vec[PLACE1_LINK];
                move(direction, -1);
                return A.conj();
        }
        else {
                return A;
        }
}

double link1::border_sign(int mu){
	if(abs(mu) == 4){
		if(0 > (coordinate[3] + mu/abs(mu)) || (coordinate[3] + mu/abs(mu)) > lattice_size[3]-1) return -1.;
		else return 1.;
	}
	else return 1.;
}

double link1::get_angle_abelian(const vector<matrix>& vec) {
	double angle;
	if (direction / abs(direction) == 1) {
		angle = atan2(vec[PLACE1_LINK].a3, vec[PLACE1_LINK].a0);
		return angle;
	}
	if (direction / abs(direction) == -1) {
		move(direction, 1);
        angle = atan2(vec[PLACE1_LINK].a3, vec[PLACE1_LINK].a0);
        move(direction, -1);
		return -angle;
	}
	else {
		cout << "error in get_matrix" << endl;
		return -1;
	}
}

matrix link1::schwinger_line(const data& conf, int d, int dir, int x) {
	int dir1 = direction;
	matrix A;
	for (int i = 0; i < d; i++) {
		A = A * get_matrix(conf.array);
		move(dir1, 1);
	}
	move_dir(dir);
	for (int i = 0; i < x; i++) {
		A = A * get_matrix(conf.array);
		move(dir, 1);
	}
	move_dir(dir1);
	move(-dir, x);
	move(-dir1, d);
	return A;
}

matrix link1::plaket(const data& conf) {
	int dir = direction;
	matrix A = get_matrix(conf.array);
	move(dir, 1);
	move_dir(-4);//the same temporal direction as polyakov loop, connected to schwinger line, has
	A = A * get_matrix(conf.array);
	move(-4, 1);
	move_dir(-dir);
	A = A * get_matrix(conf.array);
	move(-dir, 1);
	move_dir(4);
	A = A * get_matrix(conf.array);
	move(4, 1);
	move_dir(dir);
	return A;
}

matrix link1::plaket_mu(const data& conf, int mu) {
	int dir = direction;
	matrix A = get_matrix(conf.array);
	move(dir, 1);
	move_dir(-mu);//the same temporal direction as polyakov loop, connected to schwinger line, has
	A = A * get_matrix(conf.array);
	move(-mu, 1);
	move_dir(-dir);
	A = A * get_matrix(conf.array);
	move(-dir, 1);
	move_dir(mu);
	A = A * get_matrix(conf.array);
	move(mu, 1);
	move_dir(dir);
	return A;
}

double link1::plaket_abelian_mu(const data& conf, int mu) {
	int dir = direction;
	double angle = get_angle_abelian(conf.array);
	move(dir, 1);
	move_dir(-mu);//the same temporal direction as polyakov loop, connected to schwinger line, has
	angle = angle + get_angle_abelian(conf.array);
	move(-mu, 1);
	move_dir(-dir);
	angle = angle + get_angle_abelian(conf.array);
	move(-dir, 1);
	move_dir(mu);
	angle = angle + get_angle_abelian(conf.array);
	move(mu, 1);
	move_dir(dir);
	return cos(angle);
}

matrix link1::plaket_implement4(const data& conf, int mu){
	matrix A = plaket_mu(conf, mu);
	move(mu, 1);
	A = A + plaket_mu(conf, mu);
	move(-direction, 1);
	A = A + plaket_mu(conf, mu);
	move(-mu, 1);
	A = A + plaket_mu(conf, mu);
	move(direction, 1);
	return 1./4*A;
}

double link1::plaket_abelian_implement4(const data& conf, int mu){
	double a = plaket_abelian_mu(conf, mu);
	move(mu, 1);
	a = a + plaket_abelian_mu(conf, mu);
	move(-direction, 1);
	a = a + plaket_abelian_mu(conf, mu);
	move(-mu, 1);
	a = a + plaket_abelian_mu(conf, mu);
	move(direction, 1);
	return 1./4*a;
}

matrix link1::plaket_implement2(const data& conf, int mu){
	matrix A = plaket_mu(conf, mu);
	move(-direction, 1);
	A = A + plaket_mu(conf, mu);
	move(direction, 1);
	return 1./2*A;
}

matrix link1::polyakov_loop(const data& conf) {
	matrix A;
	for (int i = 0; i < lattice_size[3]; i++) {
		A = A * get_matrix(conf.array);
		move(direction, 1);
	}
	return A;
}

double link1::polyakov_loop_abelian(const data& conf) {
	double angle = 0;
	for (int i = 0; i < lattice_size[3]; i++) {
		//cout<<coordinate[0]<<" "<<coordinate[1]<<" "<<coordinate[2]<<" "<<coordinate[3]<<" "<<direction<<endl;
		angle += get_angle_abelian(conf.array);
		move(direction, 1);
	}
	return angle;
}

matrix link1::wilson_loop(const data& conf, int R, int T) {
	int dir = direction;
	matrix A;
	for (int i = 0; i < R; i++) {
		A = A * get_matrix(conf.array);
		move(dir, 1);
	}
	move_dir(4);
	for (int i = 0; i < T; i++) {
		A = A * get_matrix(conf.array);
		move(4, 1);
	}
	move_dir(-dir);
	for (int i = 0; i < R; i++) {
		A = A * get_matrix(conf.array);
		move(-dir, 1);
	}
	move_dir(-4);
	for (int i = 0; i < T; i++) {
		A = A * get_matrix(conf.array);
		move(-4, 1);
	}
	move_dir(dir);
	return A;
}

double link1::wilson_loop_abelian(const data& conf, int R, int T) {
	int dir = direction;
	double angle = 0;
	for (int i = 0; i < R; i++) {
		angle += get_angle_abelian(conf.array);
		move(dir, 1);
	}
	move_dir(4);
	for (int i = 0; i < T; i++) {
		angle += get_angle_abelian(conf.array);
		move(4, 1);
	}
	move_dir(-dir);
	for (int i = 0; i < R; i++) {
		angle += get_angle_abelian(conf.array);
		move(-dir, 1);
	}
	move_dir(-4);
	for (int i = 0; i < T; i++) {
		angle += get_angle_abelian(conf.array);
		move(-4, 1);
	}
	move_dir(dir);
	return cos(angle);
}

double link1::field1(const vector<vector<matrix> >& schwinger_line, const vector<matrix>& plaket, const vector<matrix>& polyakov_loop, int d, int D, int dir, int x) {
	int dir1 = direction;
	matrix C = schwinger_line[dir - 1][PLACE];
	move(dir1, d);
	move(dir, x);
	matrix B = plaket[PLACE];
	move(-dir, x);
	move(-dir1, d);
	move_dir(-4);
	matrix A = polyakov_loop[PLACE_FIELD_LINK];
	A = A.conj() * C * B;
	A = A * C.conj();
	move(dir1, D);
	move_dir(4);
	B = polyakov_loop[PLACE_FIELD_LINK];
	move(-dir1, D);
	move_dir(dir1);
	return A.tr() * B.tr();
}

double link1::field2(const vector<matrix>& plaket, const vector<matrix>& polyakov_loop, int d, int D, int dir, int x) {
	int dir1 = direction;
	move_dir(-4);
	matrix A = polyakov_loop[PLACE_FIELD_LINK];
	move(dir1, d);
	move(dir, x);
	move_dir(dir1);
	matrix B = plaket[PLACE];
	move(-dir, x);
	move(dir1, D - d);
	move_dir(4);
	matrix C = polyakov_loop[PLACE_FIELD_LINK];
	move(-dir1, D);
	move_dir(dir1);
	return B.tr() * C.tr() * A.conj().tr();
}

double link1::field3(const vector<matrix>& polyakov_loop, int D, int x) {
	int dir1 = direction;
	move_dir(-4);
	matrix A = polyakov_loop[PLACE_FIELD_LINK];
	move(dir1, D);
	move_dir(4);
	matrix B = polyakov_loop[PLACE_FIELD_LINK];
	move(-dir1, D);
	move_dir(dir1);
	return A.conj().tr() * B.tr();
}

matrix link1::staples_first(const vector<matrix>& vec, int eta) {
	matrix A;
	matrix B;
	int dir = direction;
	move_dir(eta);
	A = get_matrix(vec);
	move(eta, 1);
	move_dir(dir);
	A = A * get_matrix(vec);
	move(dir, 1);
	move_dir(-eta);
	A = A * get_matrix(vec);
	move(-eta, 1);
	move(-dir, 1);
	B = get_matrix(vec);
	move(-eta, 1);
	move_dir(dir);
	B = B * get_matrix(vec);
	move(dir, 1);
	move_dir(eta);
	B = B * get_matrix(vec);
	move(eta, 1);
	move(-dir, 1);
	move_dir(dir);
	return (A + B);
}

matrix link1::staples_second(const vector<vector<matrix> >& smearing_first, int eta, int nu) {
	matrix A;
	matrix B;
	int dir = direction;
	move_dir(eta);
	A = get_matrix(smearing_first[position_first(nu, dir)]);
	move(eta, 1);
	move_dir(dir);
	A = A * get_matrix(smearing_first[position_first(nu, eta)]);
	move(dir, 1);
	move_dir(-eta);
	A = A * get_matrix(smearing_first[position_first(nu, dir)]);
	move(-eta, 1);
	move(-dir, 1);
	B = get_matrix(smearing_first[position_first(nu, dir)]);
	move(-eta, 1);
	move_dir(dir);
	B = B * get_matrix(smearing_first[position_first(nu, eta)]);
	move(dir, 1);
	move_dir(eta);
	B = B * get_matrix(smearing_first[position_first(nu, dir)]);
	move(eta, 1);
	move(-dir, 1);
	move_dir(dir);
	return (A + B);
}

matrix link1::staples_second_refresh(const vector<matrix>& vec, int eta, int nu, double alpha3) {
	matrix A;
	matrix B;
	int dir = direction;
	move_dir(eta);
	A = smearing_first_refresh(vec, nu, dir, alpha3);
	move(eta, 1);
	move_dir(dir);
	A = A * smearing_first_refresh(vec, nu, eta, alpha3);
	move(dir, 1);
	move_dir(-eta);
	A = A * smearing_first_refresh(vec, nu, dir, alpha3);
	move(-eta, 1);
	move(-dir, 1);
	B = smearing_first_refresh(vec, nu, dir, alpha3);
	move(-eta, 1);
	move_dir(dir);
	B = B * smearing_first_refresh(vec, nu, eta, alpha3);
	move(dir, 1);
	move_dir(eta);
	B = B * smearing_first_refresh(vec, nu, dir, alpha3);
	move(eta, 1);
	move(-dir, 1);
	move_dir(dir);
	return (A + B);
}

matrix link1::staples_third(const vector<vector<matrix> >& smearing_second, int eta) {
	matrix A;
	matrix B;
	int dir = direction;
	move_dir(eta);
	A = get_matrix(smearing_second[dir - 1]);
	move(eta, 1);
	move_dir(dir);
	A = A * get_matrix(smearing_second[eta - 1]);
	move(dir, 1);
	move_dir(-eta);
	A = A * get_matrix(smearing_second[dir - 1]);
	move(-eta, 1);
	move(-dir, 1);
	B = get_matrix(smearing_second[dir - 1]);
	move(-eta, 1);
	move_dir(dir);
	B = B * get_matrix(smearing_second[eta - 1]);
	move(dir, 1);
	move_dir(eta);
	B = B * get_matrix(smearing_second[dir - 1]);
	move(eta, 1);
	move(-dir, 1);
	move_dir(dir);
	return (A + B);
}

matrix link1::staples_third_refresh(const vector<matrix>& vec, int eta, double alpha2, double alpha3) {
	matrix A;
	matrix B;
	int dir = direction;
	move_dir(eta);
	A = smearing_second_refresh(vec, dir, alpha2, alpha3);
	move(eta, 1);
	move_dir(dir);
	A = A * smearing_second_refresh(vec, eta, alpha2, alpha3);
	move(dir, 1);
	move_dir(-eta);
	A = A * smearing_second_refresh(vec, dir, alpha2, alpha3);
	move(-eta, 1);
	move(-dir, 1);
	B = smearing_second_refresh(vec, dir, alpha2, alpha3);
	move(-eta, 1);
	move_dir(dir);
	B = B * smearing_second_refresh(vec, eta, alpha2, alpha3);
	move(dir, 1);
	move_dir(eta);
	B = B * smearing_second_refresh(vec, dir, alpha2, alpha3);
	move(eta, 1);
	move(-dir, 1);
	move_dir(dir);
	return (A + B);
}

vector<matrix> link1::smearing_first(const data& conf, double alpha3, int nu, int rho) {
	vector<matrix> vec(data_size);
	matrix A;
	for (int t = 0; t < t_size; t++) {
		for (int z = 0; z < z_size; z++) {
			for (int y = 0; y < y_size; y++) {
				for (int x = 0; x < x_size; x++) {
					go(x, y, z, t);
					for (int i = 1; i < 5; i++) {
						if (i != abs(nu) && i != abs(rho)) {
							move_dir(i);
							vec[PLACE1_LINK] = (1 - alpha3)*get_matrix(conf.array);
							for (int d = 1; d < 5; d++) {
								if (d != abs(nu) && d != abs(rho) && d != i) {
									vec[PLACE1_LINK] = vec[PLACE1_LINK]
									+ alpha3 / 2. * staples_first(conf.array, d);
								}
							}
							vec[PLACE1_LINK].proj();
						}
					}
				}
			}
		}
	}
	return vec;
}

vector<vector<matrix> > link1::smearing_first_full(const  data& conf, double alpha3) {
	vector<vector<matrix> > smearing(6, vector<matrix>(data_size));
	for (int i = 1; i < 4; i++) {
		for (int j = i + 1; j < 5; j++) {
			smearing[position_first(i, j)] = smearing_first(conf, alpha3, i, j);
		}
	}
	return smearing;
}

vector<matrix> link1::smearing_second(const data& conf, vector<vector<matrix> >& smearing_first, double alpha2, int nu) {
	vector<matrix> vec(data_size);
	matrix A;
	for (int t = 0; t < t_size; t++) {
		for (int z = 0; z < z_size; z++) {
			for (int y = 0; y < y_size; y++) {
				for (int x = 0; x < x_size; x++) {
					go(x, y, z, t);
					for (int i = 1; i < 5; i++) {
						if (i != abs(nu)) {
							move_dir(i);
							vec[PLACE1_LINK] = (1 - alpha2)*get_matrix(conf.array);
							for (int d = 1; d < 5; d++) {
								if (d != abs(nu) && d != i) {
									vec[PLACE1_LINK] = vec[PLACE1_LINK]
										+ alpha2 / 4. * staples_second(smearing_first, d, nu);
								}
							}
							vec[PLACE1_LINK].proj();
						}
					}
				}
			}
		}
	}
	return vec;
}

vector<vector<matrix> > link1::smearing_second_full(const data& conf, vector<vector<matrix> >& smearing_first, double alpha2) {
	vector<vector<matrix> > smearing(4, vector<matrix>(data_size));
	for (int i = 1; i < 5; i++) {
		smearing[i - 1] = smearing_second(conf, smearing_first, alpha2, i);
	}
	return smearing;
}

vector<matrix> link1::smearing_HYP(data& conf, vector<vector<matrix> >& smearing_second, double alpha1) {
	vector<matrix> vec(data_size);
	matrix A;
	matrix B;
	move_dir(4);
	for (int t = 0; t < t_size; t++) {
		for (int z = 0; z < z_size; z++) {
			for (int y = 0; y < y_size; y++) {
				for (int x = 0; x < x_size; x++) {
					go(x, y, z, t);
					vec[PLACE1_LINK] = (1 - alpha1)*get_matrix(conf.array);
					for (int d = 1; d < 4; d++) {
						vec[PLACE1_LINK] = vec[PLACE1_LINK]
						+ alpha1 / 6. * staples_third(smearing_second, d);
					}
					vec[PLACE1_LINK].proj();
				}
			}
		}
	}
	for (int d = 1; d < 4; d++) {
		move_dir(d);
		for (int t = 0; t < t_size; t++) {
			for (int z = 0; z < z_size; z++) {
				for (int y = 0; y < y_size; y++) {
					for (int x = 0; x < x_size; x++) {
						go(x, y, z, t);
						vec[PLACE1_LINK] = conf.array[PLACE1_LINK];
					}
				}
			}
		}
	}
	return vec;
}

int link1::position_first(int a, int b) {
	int i;
	int j;
	if (a < b) {
		i = a;
		j = b;
	}
	if (a > b) {
		i = b;
		j = a;
	}
	if (a == b) std::cout << "wrong directions " << a << " " << b << endl;
	return (i - 1)*(8 - i) / 2 + j - i - 1;
}

vector<matrix> link1::smearing_APE(data& conf, double alpha_APE) {
	vector<matrix> vec(data_size);
	for (int t = 0; t < t_size; t++) {
		for (int z = 0; z < z_size; z++) {
			for (int y = 0; y < y_size; y++) {
				for (int x = 0; x < x_size; x++) {
					go(x, y, z, t);
					for (int i = 1; i < 4; i++) {
						move_dir(i);
						vec[PLACE1_LINK] = (1 - alpha_APE) * conf.array[PLACE1_LINK];
						for (int d = 1; d < 4; d++) {
							if (d != i) vec[PLACE1_LINK]  = 
								vec[PLACE1_LINK]
								+ (alpha_APE / 6.)  * staples_first(conf.array, d);
						}
						vec[PLACE1_LINK].proj();
					}
				}
			}
		}
	}
	move_dir(4);
	for (int t = 0; t < t_size; t++) {
		for (int z = 0; z < z_size; z++) {
			for (int y = 0; y < y_size; y++) {
				for (int x = 0; x < x_size; x++) {
					go(x, y, z, t);
					vec[PLACE1_LINK] = conf.array[PLACE1_LINK];
				}
			}
		}
	}
	return vec;
}

matrix link1::smearing_first_refresh(const vector<matrix>& vec, int nu, int rho, double alpha3) {
	matrix A;
	A = (1 - alpha3) * get_matrix(vec);
	for (int d = 1; d < 5; d++) {
		if (d != abs(nu) && d != abs(rho) && d != abs(direction)) {
			A = A + alpha3 / 2. * staples_first(vec, d);
		}
	}
	A.proj();
	return A;
}

matrix link1::smearing_second_refresh(const vector<matrix>& vec, int nu,  double alpha2, double alpha3) {
	matrix A;
	A = (1 - alpha2) * get_matrix(vec);
	for (int d = 1; d < 5; d++) {
		if (d != abs(nu) && d != abs(direction)) {
			A = A + alpha2 / 4. * staples_second_refresh(vec, d, nu, alpha3);
		}
	}
	A.proj();
	return A;
}

vector<matrix> link1::smearing_HYP_refresh(data& conf, double alpha1, double alpha2, double alpha3) {
	vector<matrix> vec(data_size);
	vec = conf.array;
	matrix A;
	move_dir(4);
	for (int t = 0; t < t_size; t++) {
		for (int z = 0; z < z_size; z++) {
			for (int y = 0; y < y_size; y++) {
				for (int x = 0; x < x_size; x++) {
					go(x, y, z, t);
					A = (1 - alpha1) * get_matrix(vec);
					for (int d = 1; d < 4; d++) {
						A = A + alpha1 / 6. * staples_third_refresh(vec, d, alpha2, alpha3);
					}
					A.proj();
					vec[PLACE1_LINK] = A;
				}
			}
		}
	}
	return vec;
}

vector<matrix> link1::smearing_APE_refresh(data& conf, double alpha_APE) {
	vector<matrix> vec(data_size);
	vec = conf.array;
	matrix A;
	for (int t = 0; t < t_size; t++) {
		for (int z = 0; z < z_size; z++) {
			for (int y = 0; y < y_size; y++) {
				for (int x = 0; x < x_size; x++) {
					go(x, y, z, t);
					for (int i = 1; i < 4; i++) {
						move_dir(i);
						A = (1 - alpha_APE) * conf.array[PLACE1_LINK];
						for (int d = 1; d < 4; d++) {
							if (d != i) A = A + (alpha_APE / 4.)  * staples_first(vec, d);
						}
						A.proj();
						vec[PLACE1_LINK] = A;			
					}
				}
			}
		}
	}
	return vec;
}

vector<matrix> link1::smearing_stout(data& conf, double rho){
	vector<matrix> vec(data_size);
	vec = conf.array;
	for(int t = 0; t < t_size; t++) {
        for (int z = 0; z < z_size; z++) {
            for (int y = 0; y < y_size; y++) {
                for (int x = 0; x < x_size; x++) {
                    go(x, y, z, t);
					for (int i = 1; i < 5; i++) {
                        move_dir(i);
						vec[PLACE1_LINK] = stout_factor(conf, rho) * conf.array[PLACE1_LINK];
					}
				}
			}
		}
	}
	return vec;
}

matrix link1::stout_factor(data& conf, double rho){
	matrix A;
	matrix B;
	matrix C;
	matrix C1;
	C = stout_omega(conf, rho);
	C1.a0 = C.a0;
	C1.a1 = -C.a1;
	C1.a2 = -C.a2;
	C1.a3 = -C.a3;
	A = (-1.) / 2 * (C1 + (-1.) * C + (-1.) / 2 * A * (C1.tr() + (-1.) * C.tr()));
	B.a0 = exp(A.a0) * cos(powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5));
	B.a1 = exp(A.a0) * A.a1 * sin(powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5)) / powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5);
	B.a2 = exp(A.a0) * A.a2 * sin(powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5)) / powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5);
	B.a3 = exp(A.a0) * A.a3 * sin(powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5)) / powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5);
	return B;
}

matrix link1::stout_omega(data& conf, double rho){
	int dir = direction;
	matrix A;
	matrix B(0., 0., 0., 0.);
	A = conf.array[PLACE1_LINK];
	A.inverse();
	for(int i = 1;i < 5;i++){
		if(i != dir){
			B = B + staples_first(conf.array, i);
		}
	}
	A = B * A * rho;
	return A;
}

void link1::gauge_transform(data& conf){
	matrix A = conf.array[0];
	matrix C = A;
	C.conj();
	double a;
	for (int t = 0; t < t_size; t++) {
		for (int z = 0; z < z_size; z++) {
			for (int y = 0; y < y_size; y++) {
				for (int x = 0; x < x_size; x++) {
					for (int i = 1; i < 5; i++) {
						a = PLACE1_LINK;
						conf.array[a] = C * conf.array[a] * A;
					}
				}
			}
		}
	}
}

//monopoles

double link1::T(data& conf, int i, int j){
	int dir = direction;
	move_dir(i);
	float angle = get_angle_abelian(conf.array);
	move(direction, 1);
	move_dir(j);
	angle += get_angle_abelian(conf.array);
	move(direction, 1);
	move_dir(-i);
	angle += get_angle_abelian(conf.array);
	move(direction, 1);
	move_dir(-j);
	angle += get_angle_abelian(conf.array);
	move(direction, 1);
	move_dir(dir);
    while((angle>Pi)||(angle<-Pi))
        {
           	if(angle>Pi) angle=angle-2*Pi;
           	if(angle<-Pi) angle=angle+2*Pi;
        }
    return angle;
}

double link1::get_current(data& conf){
     double jj=0.;
     link1 Lx(coordinate[0], coordinate[1], coordinate[2], coordinate[3], lattice_size[0], lattice_size[1], lattice_size[2], lattice_size[3]); Lx.move(1,1);
     link1 Ly(coordinate[0], coordinate[1], coordinate[2], coordinate[3], lattice_size[0], lattice_size[1], lattice_size[2], lattice_size[3]); Ly.move(2,1);
     link1 Lz(coordinate[0], coordinate[1], coordinate[2], coordinate[3], lattice_size[0], lattice_size[1], lattice_size[2], lattice_size[3]); Lz.move(3,1);
     link1 Lt(coordinate[0], coordinate[1], coordinate[2], coordinate[3], lattice_size[0], lattice_size[1], lattice_size[2], lattice_size[3]); Lt.move(4,1);
     if(direction==4) jj=Lx.T(conf,2,3)-T(conf,2,3)-(Ly.T(conf,1,3)-T(conf,1,3))+Lz.T(conf,1,2)-T(conf,1,2);
     if(direction==1) jj=-(Lt.T(conf,2,3)-T(conf,2,3))+(Ly.T(conf,4,3)-T(conf,4,3))-(Lz.T(conf,4,2)-T(conf,4,2));
     if(direction==2) jj=Lt.T(conf,1,3)-T(conf,1,3)-(Lx.T(conf,4,3)-T(conf,4,3))+Lz.T(conf,4,1)-T(conf,4,1);
     if(direction==3) jj=-(Lt.T(conf,1,2)-T(conf,1,2))+(Lx.T(conf,4,2)-T(conf,4,2))-(Ly.T(conf,4,1)-T(conf,4,1));
     return jj/2./Pi;
}

int link1::current_test(double* J){
    for(int i=1; i<=4; i++){
    	move_dir(i);
		if ( (J[get_place()]>0.3)||(J[get_place()]<-0.3) ){
			return i;
		}
	}

    for(int i=1; i<=4; i++){
    	move_dir(i);
    	move(i,-1);
        if ( (J[get_place()]>0.3)||(J[get_place()]<-0.3) ) {
        	return -i;
        }
        move(i, 1);
    }
    return 0;
}
