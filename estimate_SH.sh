#! /bin/bash
cmd=estimate_sh; echo "compiling $cmd"; make $cmd > /dev/null

dir="../../hepth" ; pn=0.3 ; pe=0.2 ; rpt=2 ; W=2000

cores=$(grep -c ^processor /proc/cpuinfo); total="$(($rpt * $cores))"
output="theta_SH_p${pe}_r${total}.dat"
./$cmd -graph $dir/graph_W2K.gz -theta $dir/theta_W2K.dat -output $dir/$output \
	-trials $rpt -p_edge $pe -p_node $pn -mx_tc $W
