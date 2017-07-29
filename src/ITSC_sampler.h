/**
 * Copyright (C) by J.Z. (07/19/2017 16:03)
 * Distributed under terms of the MIT license.
 */

#ifndef __ITSC_SAMPER_H__
#define __ITSC_SAMPER_H__

#include "ITS_sampler.h"

class ITSCSampler : public ITSSampler {
public:
    ITSCSampler(const Sampler::Config* conf) : ITSSampler(conf) {
        p_tri_ = std::pow(conf_->p_edge, 2);
    }

    void info() const override {
#ifndef N_UN
        printf("ITS-color: graph size is known\n");
#else
        printf("ITS-color: graph size is unknown\n");
#endif
    }

    vector<std::pair<int, int>> sample() const override {
        randutils::default_rng rng;
        int N = int(1 / conf_->p_edge);
        unordered_map<int, int> node_color(nodes_.size());
        for (int node : nodes_) node_color[node] = rng.uniform(0, N - 1);

        UGraph G(GraphType::MULTI);
        for (auto&& edge : edges_)
            if (node_color[edge.first] == node_color[edge.second])
                G.addEdge(edge.first, edge.second);
        G.optimize();
        return statTrids(G);
    }
};

#endif /* __ITSC_SAMPER_H__ */
