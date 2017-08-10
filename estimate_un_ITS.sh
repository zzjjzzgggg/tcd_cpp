#! /bin/bash
cmd=estimate_un_its; echo "compiling $cmd"; make $cmd > /dev/null

dir="../../pokec"; pe=0.1 ; rpt=1; W=10000 ; mx_i=200

cores=$(grep -c ^processor /proc/cpuinfo); total="$(($rpt * $cores))"
output="theta_un_ITS_p${pe}_r${total}.dat"

./$cmd -graph $dir/graph_W10K.gz -theta $dir/theta_un_W10K.dat -output $dir/$output \
	-trials $rpt -p_edge $pe -mx_tc $W -mx_i $mx_i
