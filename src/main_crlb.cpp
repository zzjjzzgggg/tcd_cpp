/**
 * Copyright (C) by J.Z. (07/02/2017 13:52)
 * Distributed under terms of the MIT license.
 */

#include "stdafx.h"

#include "crlb.h"
#include "crlb_sgs.h"

#include <gflags/gflags.h>

DEFINE_int32(n, 27400, "number of nodes");
DEFINE_int32(mx_tc, 100, "maximum triadic cardinality");
DEFINE_bool(un, false, "known graph size?");
DEFINE_double(p_tri, 0.001, "triangle sampling rate");
DEFINE_double(p_nd, 0.001, "node sampling rate");
DEFINE_double(alpha, 0.001, "alpha");
DEFINE_string(theta, "", "theta file name");
DEFINE_string(output, "", "output file name");

int main(int argc, char *argv[]) {
    gflags::SetUsageMessage("");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    osutils::Timer tm;

    CRLB::Config conf;
    conf.n = FLAGS_n;
    conf.W = FLAGS_mx_tc;
    conf.known_size = !FLAGS_un;
    conf.alpha = FLAGS_alpha;
    conf.p_tri = FLAGS_p_tri;
    conf.p_nd = FLAGS_p_nd;
    conf.theta_fnm = FLAGS_theta;
    conf.echo();

#ifdef S_ITS
    CRLB crlb(&conf);
#elif S_SGS
    CRLB_SGS crlb(&conf);
#else
    printf("Sampling method not specified\n");
    CRLB crlb(&conf);
    return -1;
#endif

    crlb.init();
    crlb.calCRLB();
    crlb.save(FLAGS_output);

    printf("cost time %s\n", tm.getStr().c_str());
    gflags::ShutDownCommandLineFlags();
    return 0;
}
