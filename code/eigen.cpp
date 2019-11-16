#include "eigen.h"

void dirac_mult(scomplex_t* res, const scomplex_t* src, int place, matrix* data_conf, int x_size, int y_size, int z_size, int t_size, double mass, double mu_q){
        int t = place/(x_size * y_size * z_size*2);
        int z = (place - (x_size * y_size * z_size*2)*t)/(x_size * y_size*2);
        int y = (place - (x_size * y_size * z_size*2)*t - (2*x_size * y_size)*z)/(x_size*2);
        int x = (place - (x_size * y_size * z_size*2)*t - (2*x_size * y_size)*z - 2*x_size*y)/2;
        int src_index = place - ((t - 1) * 2 * x_size * y_size * z_size
                + (z - 1) * 2 * x_size * y_size
                + (y - 1) * 2 * x_size
                + (x - 1) * 2);
        link1 link(x_size, y_size, z_size, t_size);
        link1 link_ferm(x_size, y_size, z_size, t_size);
        link.go(x, y, z, t);
        link_ferm.go(x, y, z, t);
        //matrix A(1, 0, 0, 0);
        double delta_4 = 0;
                double sign;
                double border_sign;
                scomplex_t vec;
                for(int mu = 1;mu <= 4;mu++){
                        if(mu == 4) delta_4 = 1;
                        else delta_4 = 0;
                        link.move_dir(mu);
                        sign = eta_sign(mu, link_ferm);
                        border_sign = link_ferm.border_sign(mu);
                        link_ferm.move(mu, 1);
                        matrix_mult_complex1(exp(mu_q * delta_4)/2 * sign * link.get_matrix1(data_conf), &src[complex_place(link_ferm)], &vec, src_index, border_sign);
                        res->re = res->re + vec.re;
                        res->im = res->im + vec.im;
                        link.move_dir(-mu);
                        link_ferm.move(-mu, 1);
                        border_sign = link_ferm.border_sign(-mu);
                        link_ferm.move(-mu, 1);
                        matrix_mult_complex1(exp(-mu_q * delta_4)/2 * sign * link.get_matrix1(data_conf), &src[complex_place(link_ferm)], &vec, src_index, border_sign);
                        res->re = res->re - vec.re;
                        res->im = res->im - vec.im;
                        link_ferm.move(mu, 1);
                }
                res->re = res->re + mass*src[complex_place(link_ferm) + src_index].re;
                res->im = res->im + mass*src[complex_place(link_ferm) + src_index].im;
}

/*void dirac_mult_conj(scomplex_t* res, const scomplex_t* src, int place, matrix* data_conf, int x_size, int y_size, int z_size, int t_size){
        int t = place/(x_size * y_size * z_size*2) + 1;
        int z = (place - (x_size * y_size * z_size*2)*(t - 1))/(x_size * y_size*2) + 1;
        int y = (place - (x_size * y_size * z_size*2)*(t - 1) - (2*x_size * y_size)*(z - 1))/(x_size*2) + 1;
        int x = (place - (x_size * y_size * z_size*2)*(t - 1) - (2*x_size * y_size)*(z - 1) - 2*x_size*(y - 1))/2 + 1;
        int src_index = place - ((t - 1) * 2 * x_size * y_size * z_size
                + (z - 1) * 2 * x_size * y_size
                + (y - 1) * 2 * x_size
                + (x - 1) * 2);
        link1 link(x_size, y_size, z_size, t_size);
        link1 link_ferm(x_size, y_size, z_size, t_size);
        link.go(x, y, z, t);
        link_ferm.go(x, y, z, t);
                double sign;
                double sign1;
                double sign2;

                scomplex_t vec;
                for(int mu = 1;mu <= 4;mu++){
                        link.move_dir(mu);
                        // sign = eta_sign(mu, link_ferm);
                        // sign1 = eta_sign_5(link_ferm);
                        link_ferm.move(-mu, 1);
                        sign = eta_sign(mu, link_ferm);
                        link.move(-mu, 1);
                        // sign2 = eta_sign_5(link_ferm);
                        matrix_mult_complex1(/*sign1 * sign2 **/ /*sign * link.get_matrix1(data_conf).conj(), &src[complex_place(link_ferm)], &vec, src_index);
                        res->re = res->re + vec.re;
                        res->im = res->im + vec.im;
                        link.move_dir(-mu);
                        link_ferm.move(mu, 2);
                        link.move(mu, 2);
                        sign = eta_sign(mu, link_ferm);
                        // sign2 = eta_sign_5(link_ferm);
                        matrix_mult_complex1(/*sign1 * sign2 **/ /*sign * link.get_matrix1(data_conf).conj(), &src[complex_place(link_ferm)], &vec, src_index);
                        res->re = res->re - vec.re;
                        res->im = res->im - vec.im;
                        link_ferm.move(-mu, 1);
                        link.move(-mu, 1);
                }
}*/

int complex_place(link1& link){
        return (link.coordinate[3]) * 2 * link.lattice_size[0]*link.lattice_size[1]*link.lattice_size[2]
                + (link.coordinate[2]) * 2 * link.lattice_size[0]*link.lattice_size[1]
                + (link.coordinate[1]) * 2 * link.lattice_size[0]
                + (link.coordinate[0]) * 2;
}

void matrix_mult_complex1(matrix A, const scomplex_t* a, scomplex_t* a1, int i, double border_sign){
        scomplex_t z1, z2;
        //z1[0].re = A.a0; z1[0].im = A.a3;
        //z2[0].re = A.a2; z2[0].im = A.a1;
        //z1[1].re = -A.a2; z1[1].im = A.a1;
        //z2[1].re = A.a0; z2[1].im = -A.a3;
        if(i == 0){
                z1.re = A.a0; z1.im = A.a3;
                z2.re = A.a2; z2.im = A.a1;
        }
        if(i == 1){
                z1.re = -A.a2; z1.im = A.a1;
                z2.re = A.a0; z2.im = -A.a3;
        }
        a1->re = border_sign*(a[0].re * z1.re - a[0].im * z1.im + a[1].re * z2.re - a[1].im * z2.im);
        a1->im = border_sign*(a[0].re * z1.im + a[0].im * z1.re + a[1].re * z2.im + a[1].im * z2.re);
}

double test_module(const scomplex_t* vec, int size){
        double module;
        for(int i = 0;i < size;i++){
                module += vec[i].re * vec[i].re + vec[i].im * vec[i].im;
        }
        module = sqrt(module);
        return module;
}

void test_eigenvector(const scomplex_t* eigenvector, scomplex_t eigenvalue, int size, matrix* data_conf, int x_size, int y_size, int z_size, int t_size, double mass, double mu_q, double tolerance){
        scomplex_t* vec;
        vec = (scomplex_t*) malloc(size * sizeof(scomplex_t));
        scomplex_t res;
        for(int i = 0;i < size;i++){
                res.re = 0;
                res.im = 0;
                //test_mult(&res, eigenvector, i);
                dirac_mult(&res, eigenvector, i, data_conf, x_size, y_size, z_size, t_size, mass, mu_q);
                vec[i].re = res.re;
                vec[i].im = res.im;
        }
        for(int i = 0;i < size;i++){
                vec[i].re -= eigenvector[i].re * eigenvalue.re - eigenvector[i].im * eigenvalue.im;
                vec[i].im -= eigenvector[i].re * eigenvalue.im + eigenvector[i].im * eigenvalue.re;
        }
        if(test_module(vec, size) < tolerance) printf("eigenvector is right \n");
        else printf("eigenvector is not right \n");
        free(vec);
}

double eta_sign(int mu, link1& link){
        int n = 0;
        for(int i = 0;i < mu - 1;i++){
                n += (link.coordinate[i]);
        }
        if(n%2 == 1) return (double)(-1.);
        if(n%2 == 0) return (double)1.;
}

double eta_sign_5(link1& link){
        int n = 0;
        for(int i = 0;i < 4;i++){
                n += (link.coordinate[i] - 1);
        }
        if(n%2 == 1) return (double)(-1.);
        if(n%2 == 0) return (double)1.;
}
