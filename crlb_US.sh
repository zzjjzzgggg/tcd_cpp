#! /bin/bash
cmd=calcrlb_us; echo "compiling $cmd"; make $cmd > /dev/null

dir=../../hepth; W=10; pt=0.52; pn=0.1

n=$(sed -n 1p $dir/parms)
theta=groundtruth_W${W}.dat; output=crlb_US_W${W}_pt${pt}_pn${pn}.dat

./$cmd -theta $dir/$theta -output $dir/$output -mx_tc $W -n $n -p_tri $pt -p_nd $pn
