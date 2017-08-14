/**
 * Copyright (C) by J.Z. (07/02/2017 15:11)
 * Distributed under terms of the MIT license.
 */

#include "stdafx.h"
#include <gflags/gflags.h>

DEFINE_int32(mx_tc, 100, "maximum triadic cardinality");
DEFINE_bool(shuffle, false, "shuffle edges?");
DEFINE_string(graph, "", "theta file name");

void truncateByNode() {
    unordered_set<int> nodes_to_del;
    auto G = loadEdgeList<DGraph>(FLAGS_graph, GraphType::MULTI);
    for (auto&& ni = G.beginNI(); ni != G.endNI(); ni++)
        if (countNodeDirTriads(ni->first, G) > FLAGS_mx_tc)
            nodes_to_del.insert(ni->first);
    DGraph H(GraphType::MULTI);
    for (auto&& ei = G.beginEI(); ei != G.endEI(); ei++) {
        int src = ei.getSrcID(), dst = ei.getDstID();
        if (nodes_to_del.find(src) == nodes_to_del.end() &&
            nodes_to_del.find(dst) == nodes_to_del.end()) {
            H.addEdge(src, dst);
        }
    }
    H.defrag();
    saveEdgelist(
        H, strutils::insertMiddle(FLAGS_graph, "W{}"_format(FLAGS_mx_tc)));
}

/**
 * If an comming edge makes a node has TC larger than the threshold, then this
 * edge is discarded.
 */
unordered_map<int, int> countNewTriadsIfAddEdge(const int u, const int v,
                                                const DGraph& G) {
    unordered_map<int, int> nbr_num, new_nd_triads;
    for (int d = 0; d < G[u].getDeg(); d++) nbr_num[G[u].getNbr(d)]++;
    for (int d = 0; d < G[v].getDeg(); d++) {
        int nbr = G[v].getNbr(d);
        if (nbr_num.find(nbr) != nbr_num.end()) {
            new_nd_triads[u] += nbr_num[nbr];
            new_nd_triads[v] += nbr_num[nbr];
            new_nd_triads[nbr] += nbr_num[nbr];
        }
    }
    return new_nd_triads;
}

void truncateByEdge() {
    auto edges = ioutils::loadPrVec<int, int>(FLAGS_graph);
    if (FLAGS_shuffle) {
        randutils::default_rng rng;
        rng.shuffle(edges);
    }

    DGraph G(GraphType::MULTI);
    unordered_map<int, int> nd_triads;
    for (auto&& edge : edges) {
        int u = edge.first, v = edge.second;
        if (!G.isNode(u) || !G.isNode(v)) {
            G.addEdge(u, v);
            continue;
        }
        auto new_nd_triads = countNewTriadsIfAddEdge(u, v, G);
        bool is_ok = true;
        for (auto&& pr : new_nd_triads)
            if (pr.second + nd_triads[pr.first] > FLAGS_mx_tc) {
                is_ok = false;
                break;
            }
        if (is_ok) {
            G.addEdge(u, v);
            for (auto&& pr : new_nd_triads) nd_triads[pr.first] += pr.second;
        }
    }
    G.defrag();
    printf("N: %d, E: %d\n", G.getNodes(), G.getEdges());
    string mid = "W{}"_format(strutils::prettyNumber(FLAGS_mx_tc));
    saveEdgelist(G, strutils::insertMiddle(FLAGS_graph, mid));
}

int main(int argc, char* argv[]) {
    gflags::SetUsageMessage("");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    osutils::Timer tm;

    truncateByEdge();

    printf("cost time %s\n", tm.getStr().c_str());
    gflags::ShutDownCommandLineFlags();
    return 0;
}
