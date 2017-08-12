/**
 * Copyright (C) by J.Z. (06/30/2017 10:47)
 * Distributed under terms of the MIT license.
 */

#include "stdafx.h"

#include "ITS_sampler.h"
#include <gflags/gflags.h>

DEFINE_string(graph, "", "graph full path");
DEFINE_string(theta, "", "theta groundtruth file name");
DEFINE_string(output, "", "output file name");

DEFINE_int32(mx_tc, 2000, "maximum triadic cardinality");
DEFINE_int32(mx_i, 128, "maximum triadic cardinality");
DEFINE_int32(mx_iter_theta, 1000, "maximum iterations for estimating theta");
DEFINE_int32(mx_iter_alpha, 10, "maximum iterations for estimating alpha");
DEFINE_int32(trials, 10, "trials per core");
DEFINE_int32(cores, std::thread::hardware_concurrency(), "cores");

DEFINE_double(p_node, 0.1, "node sampling rate");
DEFINE_double(p_edge, 0.1, "edge sampling rate");
DEFINE_double(eps_theta, 1e-3, "threshold of estimating theta");
DEFINE_double(eps_alpha, 1e-3, "threshold of estimating alpha");

int main(int argc, char* argv[]) {
    printf("Built: %s %s\n", __DATE__, __TIME__);
    gflags::SetUsageMessage("usage:");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    osutils::Timer tm;

    Sampler::Config confs;
    confs.graph_fnm = FLAGS_graph;
    confs.mx_tc = FLAGS_mx_tc;
    confs.p_edge = FLAGS_p_edge;
    confs.p_node = FLAGS_p_node;
    confs.echo();

    ITSSampler sampler(&confs);
    for (int j = 0; j < 10; j++) printf("%.4e\n", sampler.bb(j, 10, 0.01));

    printf("cost time %s\n", tm.getStr().c_str());
    gflags::ShutDownCommandLineFlags();
    return 0;
}
