#! /bin/bash
cmd=calcrlb_is; echo "compiling $cmd"; make $cmd > /dev/null

dir=../../dblp
W=20 ; pt=0.22

theta=theta_W20.dat ; output=crlb_IS_W20_p${pt}.dat
./$cmd -theta $dir/$theta -output $dir/$output -mx_tc $W -p_tri $pt \
	-n $(sed -n 1p $dir/parms) -alpha 0.1
