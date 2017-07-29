#! /bin/bash
cmd=estimate_un_itsc ; echo "compiling $cmd" ; make $cmd > /dev/null

dir="../../pokec" ; pe=0.2 ; rpt=1 ; W=10000 ; MX_i=500

cores=$(grep -c ^processor /proc/cpuinfo); total="$(($rpt * $cores))"
output="theta_un_ITSC_p${pe}_r${total}.dat"

./$cmd -graph $dir/graph_W10K.gz -theta $dir/theta_un_W10K.dat -output $dir/$output \
	-trials $rpt -p_edge $pe -mx_tc $W -mx_i $MX_i
