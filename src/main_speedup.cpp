/**
 * Copyright (C) by J.Z. (08/16/2017 16:36)
 * Distributed under terms of the MIT license.
 */

#include "stdafx.h"

#include <gflags/gflags.h>

DEFINE_string(graph, "", "graph full path");
DEFINE_string(output, "", "output file name");

int countTriangles(const UGraph& G, const vector<int>& nodes = {}) {
    int total_triangles = 0;
    if (nodes.empty()) {
        for (auto&& ni = G.beginNI(); ni != G.endNI(); ni++)
            total_triangles += countNodeDirTriads(ni->first, G);
    } else {
        for (auto node : nodes) total_triangles += countNodeDirTriads(node, G);
    }
    return total_triangles;
}

UGraph sampleEdges(const UGraph& G, const double p) {
    rngutils::default_rng rng;
    UGraph H(GraphType::MULTI);
    for (auto&& ei = G.beginEI(); ei != G.endEI(); ei++)
        if (rng.uniform() < p) H.addEdge(ei.getSrcID(), ei.getDstID());
    H.defrag();
    return H;
}

vector<int> sampleNodes(const UGraph& G, const double p) {
    rngutils::default_rng rng;
    vector<int> nodes;
    for (auto&& ni = G.beginNI(); ni != G.endNI(); ni++) {
        if (rng.uniform() < p) nodes.push_back(ni->first);
    }
    return nodes;
}

int main(int argc, char* argv[]) {
    gflags::SetUsageMessage("usage:");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    osutils::Timer tm;

    UGraph G = loadEdgeList<UGraph>(FLAGS_graph, GraphType::MULTI);

    vector<std::pair<double, double>> p_tm_v;
    for (double p = 0.1; p < 1.01; p += 0.1) {
        // for (double p = 0.01; p < 0.11; p += 0.01) {
        tm.tick();
        countTriangles(sampleEdges(G, p));
        // countTriangles(G, sampleNodes(G, p));

        double t = tm.seconds();
        p_tm_v.emplace_back(p, t);
        printf("%g\t%.2f\n", p, t);
    }
    ioutils::savePrVec(p_tm_v, strutils::subFilename(FLAGS_graph, FLAGS_output),
                       true, "{}\t{:.3f}\n");

    printf("cost time %s\n", tm.getStr().c_str());
    gflags::ShutDownCommandLineFlags();
    return 0;
}
