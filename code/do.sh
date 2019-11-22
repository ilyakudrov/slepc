#!/bin/bash
t_size=32
x_size=32
y_size=32
z_size=32
mass=0.0075
mu=5
w=$(($mu/10))
q=$(($mu-w*10))
#path="/home/clusters/rrcmpi/nikolaev/SU_2_staggered/configurations/N_f=2/improved/32^3x$t_size/mu=0.$w$q/beta=1.8/ma=0.0075/lambda=0.00075/converted_Fort/conf_0001.fl"
path="/home/clusters/rrcmpi/kudrov/smeared_stout/time_$t_size/mu0.$w$q/conf_0001.fl"
start=$(date +%s)
/home/clusters/rrcmpi/kudrov/eigenvalues/slepc/code/eigen_test $path 0.$w$q -eps_nev 8 -eps_ncv 40 	-mat_cusparse_storage_format hyb
finish=$(date +%s)
echo $[$finish-$start]
