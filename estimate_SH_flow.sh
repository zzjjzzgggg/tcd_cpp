#! /bin/bash
cmd=estimate_sh_flow ; echo "compiling $cmd" ; make $cmd > /dev/null

dir="../../flow" ; pe=0.02 ; rpt=10 ; W=20

cores=$(grep -c ^processor /proc/cpuinfo); total="$(($rpt * $cores))"
output="theta_SH_p${pe}_r${total}.dat"
./$cmd -others $dir/flow.gz -theta $dir/groundtruth_W20.dat -output $dir/$output \
	-trials $rpt -p_edge $pe -mx_tc $W -eps_theta 0.0002 -mx_iter_theta 5000 \
	-eps_alpha 0.001 -mx_iter_alpha 100 -cores 1
