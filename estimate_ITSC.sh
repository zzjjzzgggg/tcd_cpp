#! /bin/bash
cmd=estimate_itsc ; echo "compiling $cmd" ; make $cmd > /dev/null

dir="../../hepth"
pe=0.2 ; rpt=10 ; W=2000 ; MX_i=200

cores=$(grep -c ^processor /proc/cpuinfo); total="$(($rpt * $cores))"
output="theta_ITSC_p${pe}_r${total}.dat"

./$cmd -graph $dir/graph_W2K.gz -theta $dir/theta_W2K.dat -output $dir/$output \
	-trials $rpt -p_edge $pe -mx_tc $W -mx_i $MX_i
