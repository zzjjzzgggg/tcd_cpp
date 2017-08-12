#! /bin/bash
cmd=estimate_its; echo "compiling $cmd"; make $cmd > /dev/null

dir="../../hepth" ; pe=0.3 ; rpt=10 ; W=2047

cores=$(grep -c ^processor /proc/cpuinfo); total="$(($rpt * $cores))"
output="theta_ITS_binned_p${pe}_r${total}.dat"

./$cmd -graph $dir/graph_W2K.gz -theta theta_W2K_binned.dat -output $output \
	-trials $rpt -p_edge $pe -mx_tc $W
