/**
 * Copyright (C) by J.Z. (07/02/2017 13:52)
 * Distributed under terms of the MIT license.
 */

#include "stdafx.h"
#include "crlb.h"
#include <gflags/gflags.h>

DEFINE_int32(n, 27400, "number of nodes");
DEFINE_int32(mx_tc, 100, "maximum triadic cardinality");
DEFINE_double(p_tri, 0.001, "triangle sampling rate");
#ifdef METHOD_US
DEFINE_double(p_nd, 0.001, "node sampling rate");
#endif
DEFINE_double(alpha, 0.001, "alpha");
DEFINE_string(theta, "", "theta file name");
DEFINE_string(output, "", "output file name");

int main(int argc, char *argv[]) {
    gflags::SetUsageMessage("");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    osutils::Timer tm;

    CRLB::Config conf;
    conf.W = FLAGS_mx_tc;
    conf.n = FLAGS_n;
    conf.alpha = FLAGS_alpha;
    conf.p_tri = FLAGS_p_tri;
#ifdef METHOD_US
    conf.p_nd = FLAGS_p_nd;
#endif
    conf.theta_fnm = FLAGS_theta;
    conf.echo();

    CRLB crlb(&conf);
#ifdef METHOD_IS
    crlb.calCRLB(crlb.getISB());
#elif METHOD_US
    crlb.calCRLB(crlb.getUSB());
#endif
    crlb.save(FLAGS_output);

    printf("cost time %s\n", tm.getStr().c_str());
    gflags::ShutDownCommandLineFlags();
    return 0;
}
