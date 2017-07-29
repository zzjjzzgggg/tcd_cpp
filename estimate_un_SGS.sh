#! /bin/bash
cmd=estimate_un_sgs ; echo "compiling $cmd" ; make $cmd > /dev/null

p=0.05 ; rpt=10 ; W=2000
dir="../../hepth"

cores=$(grep -c ^processor /proc/cpuinfo); total="$(($rpt * $cores))"
output="theta_un_SGS_p${p}_r${total}.dat"
./$cmd -graph $dir/graph_W2K.gz -theta $dir/theta_un_W2K.dat -output $dir/$output \
	-trials $rpt -p_node $p -mx_tc $W
