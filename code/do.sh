#!/bin/bash
t_size=32
x_size=32
y_size=32
z_size=32
mass=0.0075
mu_q=0
path="/home/ilya/lattice/slepc/conf/time_32/mu0.00/conf_0001.fl"
./eigen_test $x_size $y_size $z_size $t_size $path $mass $mu_q
