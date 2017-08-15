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
        printf("ITS-color sampler is ready\n");
    }

    vector<std::pair<int, int>> sample() const override {
        rngutils::default_rng rng;
        int N = int(1 / conf_->p_edge);
        unordered_map<int, int> node_color(nodes_.size());
        for (int node : nodes_) node_color[node] = rng.uniform(0, N - 1);

        UGraph G(GraphType::MULTI);
        for (const auto & [ src, dst ] : edges_)
            if (node_color[src] == node_color[dst]) G.addEdge(src, dst);
        G.defrag();
        return statTrids(G);
    }
};

#endif /* __ITSC_SAMPER_H__ */
