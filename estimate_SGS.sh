#! /bin/bash
cmd=estimate_sgs ; echo "compiling $cmd" ; make $cmd > /dev/null

dir="../../pokec" ; p=0.02 ; rpt=3 ; W=10000 ; MX_i=200

cores=$(grep -c ^processor /proc/cpuinfo); total="$(($rpt * $cores))"
output="theta_SGS_p${p}_r${total}.dat"

./$cmd -graph $dir/graph_W10K.gz -theta $dir/theta_W10K.dat -output $dir/$output \
	-trials $rpt -p_node $p -mx_tc $W -mx_i $MX_i
