#! /bin/bash

cmd=estimate_itsc ; echo "compiling $cmd" ; make $cmd > /dev/null

dir="../../hk/uu_daily"
pe=0.5 ; rpt=5 ; W=2000 ; MX_i=200
# graph=uu0901_0904_W2K.gz
# theta=theta_0901_0904_W2K.dat
# output="theta_ITSC_p${pe}_before.dat"
# graph=uu0928_1001_W2K.gz
# theta=theta_0928_1001_W2K.dat
# output="theta_ITSC_p${pe}_after.dat"

graph=0908_W2K.gz
theta=theta_0908_W2K.dat
output="theta_ITSC_p${pe}.dat"

./$cmd -graph $dir/$graph -theta $dir/$theta -output $dir/$output \
	-output $dir/$output -trials $rpt -p_edge $pe -mx_tc $W -mx_i $MX_i
