/**
 * Copyright (C) by J.Z. (07/14/2017 15:54)
 * Distributed under terms of the MIT license.
 */

#ifndef __SGS_SAMPLER_H__
#define __SGS_SAMPLER_H__

#include "graph_sampler.h"

/**
 * Similar to Flow Sampling. Randomly sample a collection of nodes. Keep
 * each triangle related to a sampled node.
 */
class SGSSampler : public Sampler {
private:
    UGraph G_;

public:
    SGSSampler(const Sampler::Config* conf) : Sampler(conf) {
        G_ = loadEdgeList<UGraph>(conf_->graph_fnm, GraphType::MULTI);
    }

    void info() const override {
#ifndef N_UN
        printf("SGS: graph size is known\n");
#else
        printf("SGS: graph size is unknown\n");
#endif
    }

    vector<std::pair<int, int>> sample() const override {
        randutils::default_rng rng;
        std::vector<int> samples;  // TC -> # nodes
        for (auto node : nodes_)
            if (rng.uniform() < conf_->p_node) samples.push_back(node);
        return statTrids(G_, samples);
    }

    double bji(const int j, const int i, const double alpha) const override {
        if (i == 0 && j == 0)
            return 1;
        else if (i == j)
            return conf_->p_node;
        else if (j == 0)
            return 1 - conf_->p_node;
        return 0;
    }

    bool hasAlpha() const override { return false; }

    std::pair<double, double> getLGrad(const int i, const int j,
                                       const double alpha) const override {
        return std::pair<double, double>(0, 0);
    }
};

#endif /* __SGS_SAMPLER_H__ */
