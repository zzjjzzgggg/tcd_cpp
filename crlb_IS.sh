#! /bin/bash
cmd=calcrlb_is; echo "compiling $cmd"; make $cmd > /dev/null

dir=../../hepth; W=20; pt=0.1

n=$(sed -n 1p $dir/parms)
theta=theta_W${W}.dat ; output=crlb_IS_W${W}_p${pt}.dat
./$cmd -theta $dir/$theta -output $dir/$output -mx_tc $W -p_tri $pt -n $n -alpha 0.1
