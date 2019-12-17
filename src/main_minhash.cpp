/**
 * Copyright (C) by J.Z. 2019-12-16 22:29
 * Distributed under terms of the MIT license.
 */
#include "stdafx.h"

#include <gflags/gflags.h>

DEFINE_string(graph, "", "undirected simple graph, sorted by src");
DEFINE_int32(m, 16, "number of traversals");
int main(int argc, char* argv[]) {
    gflags::SetUsageMessage("usage:");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    osutils::Timer tm;

    // load the graph to obtain several statistics
    auto G = graph::loadEdgeList<graph::undir::UGraph>(FLAGS_graph);
    int n = G.getNodes();

    std::vector<int> deg(n);
    for (auto&& ni = G.beginNI(); ni != G.endNI(); ++ni)
        deg[ni->first] = ni->second.getDeg();

    std::vector<double> T(n);
    std::unordered_map<std::string, int> Z;

    ioutils::TSVParser ss(FLAGS_graph);

    for (int i = 0; i < FLAGS_m; ++i) {
        rngutils::default_rng rng;
        std::vector<int> min_u(n, n), hash_u(n);
        for (int k = 0; k < n; ++k) hash_u[k] = rng.randint(0, n - 1);

        ss.reset();
        while (ss.next()) {
            int src = ss.get<int>(0), dst = ss.get<int>(1);
            min_u[src] = std::min(min_u[src], hash_u[dst]);
        }

        ss.reset();
        while (ss.next()) {
            int src = ss.get<int>(0), dst = ss.get<int>(1);
            if (min_u[src] == min_u[dst]) {
                std::string key = "{}_{}"_format(src, dst);
                ++Z[key];
            }
        }
    }

    // after m iterations
    ss.reset();
    while (ss.next()) {
        int src = ss.get<int>(0), dst = ss.get<int>(1);
        std::string key = "{}_{}"_format(src, dst);
        T[src] += double(Z[key]) / (Z[key] + FLAGS_m) * (deg[src] + deg[dst]);
    }

    for (auto& t : T) t /= 2;

    std::string filename = osutils::join(strutils::getBasePath(FLAGS_graph),
                                         "minhash_m{}.dat"_format(FLAGS_m));
    ioutils::saveVec(T, filename, "{:.2f}\n");

    printf("cost time %s\n", tm.getStr().c_str());
    gflags::ShutDownCommandLineFlags();
    return 0;
}
