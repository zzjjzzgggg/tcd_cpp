/**
 * Copyright (C) by J.Z. (07/02/2017 16:11)
 * Distributed under terms of the MIT license.
 */

#include "stdafx.h"
#include "graph_sampler.h"

#include <gflags/gflags.h>

DEFINE_string(graph, "", "graph file name");
DEFINE_string(output, "", "output file name");
DEFINE_int32(mx_tc, 10000, "maximum triadic cardinality");

int main(int argc, char* argv[]) {
    gflags::SetUsageMessage("");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    osutils::Timer tm;

    Sampler::Config conf;
    conf.graph_fnm = FLAGS_graph;
    conf.mx_tc = FLAGS_mx_tc;

    Sampler sam(&conf);
    sam.getGroundtruth(FLAGS_output);

    printf("cost time %s\n", tm.getStr().c_str());
    gflags::ShutDownCommandLineFlags();
    return 0;
}
