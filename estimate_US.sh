#! /bin/bash
cmd=estimate_us; echo "compiling $cmd"; make $cmd > /dev/null

dir="../../hepth"; pn=0.1 ; pe=0.3 ; rpt=5; W=10

cores=$(grep -c ^processor /proc/cpuinfo); total="$(($rpt * $cores))"
output="theta_US_pe${pe}_pn${pn}_r${total}_W10.dat"
./$cmd -graph $dir/graph_W10.gz -theta $dir/groundtruth_W10.dat -output $dir/$output \
	-trials $rpt -p_edge $pe -p_node $pn -mx_tc $W
