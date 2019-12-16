/**
 * Copyright (C) by J.Z. (07/02/2017 16:11)
 * Distributed under terms of the MIT license.
 */

#include "stdafx.h"

#include <gflags/gflags.h>

DEFINE_string(graph, "", "graph full path");
DEFINE_string(output, "", "output file name");
DEFINE_bool(un, false, "graph size unknown?");

int main(int argc, char* argv[]) {
    gflags::SetUsageMessage("");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    osutils::Timer tm;
    unordered_map<int, int> tc_to_num_nodes;

    auto G = graph::loadEdgeList<graph::undir::UGraph>(FLAGS_graph,
                                                       graph::GraphType::MULTI);
    printf("N: %d, E: %d\n", G.getNodes(), G.getEdges());
    if (FLAGS_un)
        printf("unknown graph size\n");
    else
        printf("graph size known\n");

    int n_plus = 0;
    for (auto&& ni = G.beginNI(); ni != G.endNI(); ni++) {
        int tc = countNodeDirTriads(ni->first, G);
        int k = tc > 0 ? int(std::floor(std::log2(tc))) : -1;
        if (!FLAGS_un || tc > 0) {
            tc_to_num_nodes[k]++;
            n_plus++;
        }
    }

    vector<std::tuple<int, double, int>> dist;
    for (const auto& [k, n] : tc_to_num_nodes)
        dist.emplace_back(k, n / double(n_plus), n);
    // if unknown graph size, store n_+ at head
    if (FLAGS_un) dist.emplace_back(-10, n_plus, 0);

    std::sort(dist.begin(), dist.end());

    ioutils::saveTupleVec(dist,
                          strutils::subFilename(FLAGS_graph, FLAGS_output),
                          "{}\t{:.6e}\t{}\n");

    printf("cost time %s\n", tm.getStr().c_str());
    gflags::ShutDownCommandLineFlags();
    return 0;
}
