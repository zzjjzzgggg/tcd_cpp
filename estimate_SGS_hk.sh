#! /bin/bash

cmd=estimate_sgs ; echo "compiling $cmd" ; make $cmd > /dev/null

dir="../../hk/uu_daily"
p=0.02 ; rpt=5 ; W=2000 ; MX_i=200

# graph=uu0901_0904_W2K.gz
# theta=theta_0901_0904_W2K.dat
# output="theta_SGS_p${p}_before.dat"
# graph=uu0928_1001_W2K.gz
# theta=theta_0928_1001_W2K.dat
# output="theta_SGS_p${p}_after.dat"

graph=0908_W2K.gz
theta=theta_0908_W2K.dat
output="theta_SGS_p${p}.dat"

./$cmd -graph $dir/$graph -theta $dir/$theta -output $dir/$output \
	-trials $rpt -p_node $p -mx_tc $W -mx_i $MX_i
