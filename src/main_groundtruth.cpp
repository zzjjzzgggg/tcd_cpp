/**
 * Copyright (C) by J.Z. (07/02/2017 16:11)
 * Distributed under terms of the MIT license.
 */

#include "stdafx.h"

#include <gflags/gflags.h>

DEFINE_string(graph, "", "graph file name");
DEFINE_string(output, "", "output file name");
DEFINE_bool(un, false, "graph size unknown?");

int main(int argc, char* argv[]) {
    gflags::SetUsageMessage("");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    osutils::Timer tm;
    unordered_map<int, int> tc_to_num_nodes;

    UGraph G = loadEdgeList<UGraph>(FLAGS_graph, GraphType::MULTI);
    printf("N: %d, E: %d\n", G.getNodes(), G.getEdges());

    int n = 0;
    for (auto&& ni = G.beginNI(); ni != G.endNI(); ni++) {
        int tc = countNodeDirTriads(ni->first, G);
        tc = tc > 0 ? int(std::floor(std::log2(tc))) : -1;
        if (!FLAGS_un || tc > 0) {
            tc_to_num_nodes[tc]++;
            n++;
        }
    }

    vector<std::tuple<int, double, int>> dist;
    for (auto&& it : tc_to_num_nodes)
        dist.emplace_back(it.first, it.second / double(n), it.second);
    // if unknown graph size, store n_+ at head
    if (FLAGS_un) dist.emplace_back(0, n, 0);

    std::sort(dist.begin(), dist.end());

    string outfnm = strutils::replaceFilename(FLAGS_graph, FLAGS_output);
    ioutils::saveTupleVec(dist, outfnm, true, "{}\t{:.6e}\t{}\n");

    printf("cost time %s\n", tm.getStr().c_str());
    gflags::ShutDownCommandLineFlags();
    return 0;
}
