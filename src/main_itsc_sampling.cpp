/**
 * Copyright (C) by J.Z. 2019-12-18 19:03
 * Distributed under terms of the MIT license.
 */

#include "stdafx.h"

#include <gflags/gflags.h>

DEFINE_string(graph, "", "graph");
DEFINE_bool(dir, false, "directed");
DEFINE_int32(r, 10, "repeat times");
DEFINE_double(p, 0.1, "nodes sampling rate");

template <class Graph>
void sampleGraph(const Graph& G) {
    int N = int(1 / FLAGS_p);
    rngutils::default_rng rng;
    double num_nodes = 0, num_edges = 0;
    for (int r = 0; r < FLAGS_r; ++r) {
        std::unordered_map<int, int> h;
        for (auto ni = G.beginNI(); ni != G.endNI(); ++ni)
            h[ni->first] = rng.randint(1, N);

        Graph S;
        for (auto&& ei = G.beginEI(); ei != G.endEI(); ++ei) {
            int src = ei.getSrcID(), dst = ei.getDstID();
            if (h[src] == h[dst]) S.addEdge(ei.getSrcID(), ei.getDstID());
        }

        num_nodes += S.getNodes();
        num_edges += S.getEdges();
    }
    printf("sampled graph nodes: %.0f, edges: %.0f\n", num_nodes / FLAGS_r,
           num_edges / FLAGS_r);
}

int main(int argc, char* argv[]) {
    gflags::SetUsageMessage("usage:");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    osutils::Timer tm;

    if (FLAGS_dir) {
        printf("directed graph\n");
        sampleGraph(graph::loadEdgeList<graph::dir::DGraph>(FLAGS_graph));
    } else {
        printf("undirected graph\n");
        sampleGraph(graph::loadEdgeList<graph::undir::UGraph>(FLAGS_graph));
    }
    printf("cost time %s\n", tm.getStr().c_str());
    gflags::ShutDownCommandLineFlags();
    return 0;
}
