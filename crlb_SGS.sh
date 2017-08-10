#! /bin/bash
cmd=calcrlb_fs; echo "compiling $cmd"; make $cmd > /dev/null

# dir=../../hepth
dir=../../dblp
W=20 ; pn=0.07

theta=theta_W20.dat
output=crlb_SGS_W20_p${pn}.dat

./$cmd -theta $dir/$theta -output $dir/$output -mx_tc $W -p_nd $pn \
	-n $(sed -n 1p $dir/parms)
