/**
 * Copyright (C) by J.Z. (07/14/2017 15:51)
 * Distributed under terms of the MIT license.
 */

#ifndef __US_SAMPLER_H__
#define __US_SAMPLER_H__

#include "graph_sampler.h"

/**
 * deprecated.
 * This method has the distadvantage of both IS and FS.
 */
class USSampler : public Sampler {
private:
    // count the number of triads of a node in a probabilistic manner
    int countTriadsWithProb(const int node, const UGraph& G,
                            const unordered_set<int>& marked_nds) const {
        if (G.getNode(node).getDeg() < 2) return 0;
        randutils::default_rng rng;
        int triads = 0;
        // get unique neighbors
        std::unordered_map<int, int> nbr_cnt;
        auto&& nd_obj = G.getNode(node);
        for (int d = 0; d < nd_obj.getDeg(); d++) {
            int nbr = nd_obj.getNbr(d);
            if (nbr != node) nbr_cnt[nbr]++;
        }
        // store neighbors in the vector to allow sequential access
        std::vector<int> nbr_vec;
        nbr_vec.reserve(nbr_cnt.size());
        for (auto&& it : nbr_cnt) nbr_vec.push_back(it.first);
        // count connected neighbors
        for (auto src_it = nbr_vec.begin(); src_it != nbr_vec.end(); src_it++) {
            bool src_marked = marked_nds.find(*src_it) != marked_nds.end();
            auto& src_obj = G.getNode(*src_it);
            std::unordered_map<int, int> src_nbr_cnt;
            for (int d = 0; d < src_obj.getDeg(); d++)
                src_nbr_cnt[src_obj.getNbr(d)]++;
            for (auto dst_it = src_it + 1; dst_it != nbr_vec.end(); dst_it++) {
                if (src_nbr_cnt.find(*dst_it) != src_nbr_cnt.end()) {
                    int valid_edges = src_nbr_cnt[*dst_it];
                    if (src_marked ||
                        marked_nds.find(*dst_it) != marked_nds.end())
                        for (int k = 0; k < src_nbr_cnt[*dst_it]; k++)
                            if (rng.uniform() > conf_->p_edge) valid_edges--;
                    triads += valid_edges * nbr_cnt[*src_it] * nbr_cnt[*dst_it];
                }
            }
        }
        return triads;
    }

public:
    USSampler(const Sampler::Config* conf) : Sampler(conf) {
        p_tri_ = conf_->p_edge;
    }

    void info() const override { printf("US\n"); }

    vector<std::pair<int, int>> sample() const override {
        randutils::default_rng rng;
        // get node samples
        unordered_set<int> samples;
        for (auto&& ni = G_.beginNI(); ni != G_.endNI(); ni++)
            if (rng.uniform() < conf_->p_node) samples.insert(ni->first);

        // sample graph
        UGraph G(GraphType::MULTI);
        for (auto&& ei = G_.beginEI(); ei != G_.endEI(); ei++) {
            int src = ei.getSrcID(), dst = ei.getDstID();
            if (samples.find(src) != samples.end() ||
                samples.find(dst) != samples.end() ||
                rng.uniform() < conf_->p_edge)
                G.addEdge(src, dst);
        }
        G.optimize();

        return statTrids(G, samples);
    }

    double bji(const int j, const int i, const double alpha) const override {
        return j == 0 ? (1 - conf_->p_node + conf_->p_node * bb(0, i, alpha))
                      : conf_->p_node * bb(j, i, alpha);
    }
};

#endif /* __US_SAMPLER_H__ */
