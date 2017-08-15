/**
 * Copyright (C) by J.Z. (07/17/2017 15:20)
 * Distributed under terms of the MIT license.
 */

#ifndef __SH_FLOW_SAMPLER_H__
#define __SH_FLOW_SAMPLER_H__

#include "SH_sampler.h"

class SHFlowSampler : public SHSampler {
public:
    SHFlowSampler(const Sampler::Config* conf) : SHSampler(conf) {
        p_tri_ = conf_->p_edge;
    }

    void info() const override { printf("SH flow\n"); }

    vector<std::pair<int, int>> sample() const override {
        rngutils::default_rng rng;

        std::unordered_map<int, int> flow_size, size_freq;
        ioutils::TSVParser ss(conf_->others);
        while (ss.next()) {
            int flow = ss.get<int>(0);
            if (flow_size.find(flow) != flow_size.end() ||
                rng.uniform() < conf_->p_edge)
                flow_size[flow]++;
        }

        for (auto&& pr : flow_size) size_freq[pr.second]++;
        size_freq[0] = 10000 - flow_size.size();

        vector<std::pair<int, int>> g;
        for (auto&& pr : size_freq) g.emplace_back(pr.first, pr.second);
        std::sort(g.begin(), g.end());
        return g;
    }

    // bool hasAlpha() const override { return false; }
};

#endif /* __SH_FLOW_SAMPLER_H__ */
