#! /bin/bash
cmd=calcrlb_its; echo "compiling $cmd"; make $cmd > /dev/null

dir=../../dblp
W=20 ; pt=0.7

theta=theta_un_W20.dat ; output=crlb_un_ITS_W20_p${pt}.dat
./$cmd -theta $dir/$theta -output $dir/$output -mx_tc $W -p_tri $pt -alpha 0.1 -un
