#!/bin/bash
t_size=2
x_size=4
y_size=4
z_size=4
mass=0.0075
mu_q=0
path="/home/ilya/lattice/slepc/conf/time_32/mu0.00"
./eigen_test $x_size $y_size $z_size $t_size $path $mass $mu_q
