#! /bin/bash
cmd=estimate_fs; echo "compiling $cmd"; make $cmd > /dev/null

p=0.05 ; rpt=20 ; W=2000
dir="../../hepth"

cores=$(grep -c ^processor /proc/cpuinfo); total="$(($rpt * $cores))"
output="theta_FS_p${p}_r${total}.dat"
./$cmd -graph $dir/graph_W2K.gz -theta $dir/theta_W2K.dat -output $dir/$output \
	-trials $rpt -p_node $p -mx_tc $W
