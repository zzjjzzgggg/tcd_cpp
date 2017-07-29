#! /bin/bash
cmd=calcrlb_sgs ; echo "compiling $cmd" ; make $cmd > /dev/null

dir=../../dblp
W=20 ; pn=0.02

theta=theta_un_W20.dat
output=crlb_un_SGS_W20_p${pn}.dat

./$cmd -theta $dir/$theta -output $dir/$output -mx_tc $W -p_nd $pn -un
