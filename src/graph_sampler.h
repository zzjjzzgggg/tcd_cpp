/**
 * Copyright (C) by J.Z. (06/30/2017 10:31)
 * Distributed under terms of the MIT license.
 */

#ifndef __GRAPH_SAMPLER_H__
#define __GRAPH_SAMPLER_H__

#include "stdafx.h"

class Sampler {
public:
    struct Config {
        int mx_tc;         // maximum triadic cardinality
        double p_edge;     // edge sampling rate
        double p_node;     // node sampling rate
        string graph_fnm;  // multigraph file name
        string others;     // reserved string

        void echo() const {
            printf(
                "Sampling config:\n"
                "  graph: %s\n"
                "  maximum TC: %d\n"
                "  edge sample rate: %g\n"
                "  node sample rate: %g\n\n",
                graph_fnm.c_str(), mx_tc, p_edge, p_node);
        }
    };

protected:
    const Config* conf_;
    double p_tri_;  // prob. of sampling a triangle
    vector<int> nodes_;
    vector<std::pair<int, int>> edges_;

protected:
    // state triads histogram in a sampled graph G
    template <class Container = vector<int>>
    vector<std::pair<int, int>> statTrids(
        const UGraph& G, const Container& nodes = Container()) const {
        std::unordered_map<int, int> tc_to_num_nodes;

        int addition_zero_tc_nodes;
        if (nodes.empty()) {
            for (auto&& ni = G.beginNI(); ni != G.endNI(); ni++)
                tc_to_num_nodes[countNodeDirTriads(ni->first, G)]++;
            addition_zero_tc_nodes = nodes_.size() - G.getNodes();
        } else {
            for (auto node : nodes)
                tc_to_num_nodes[countNodeDirTriads(node, G)]++;
            addition_zero_tc_nodes = nodes_.size() - nodes.size();
        }
        tc_to_num_nodes[0] += addition_zero_tc_nodes;
#ifdef N_UN
        tc_to_num_nodes.erase(0);
#endif
        std::vector<std::pair<int, int>> triads_hist;
        for (auto&& it : tc_to_num_nodes)
            triads_hist.emplace_back(it.first, it.second);
        std::sort(triads_hist.begin(), triads_hist.end());
        return triads_hist;
    }

    /**
     * Beta-Binomial distribution
     */
    inline double bb(const int j, const int i, const double alpha) const {
        return SpecFun::BetaBinomial(j, i, p_tri_ / alpha,
                                     (1 - p_tri_) / alpha);
    }

public:
    Sampler(const Config* conf) : conf_(conf) {
        if (!conf_->graph_fnm.empty()) {
            printf("loading edges ...\n");
            ioutils::loadPrVec(conf_->graph_fnm, edges_);
            randutils::default_rng rng;
            // rng.shuffle(edges_);
            unordered_set<int> nodes;
            for (auto&& pr : edges_) {
                nodes.insert(pr.first);
                nodes.insert(pr.second);
            }
            nodes_.reserve(nodes.size());
            for (auto nd : nodes) nodes_.push_back(nd);
            printf("N: %lu, E: %lu\n", nodes_.size(), edges_.size());
        }
        // default triangle sampling probability
        p_tri_ = std::pow(conf_->p_edge, 3);
    }

    /**
     * P(Y=j|X=i, .)
     */
    double pji(const int j, const int i, const double alpha) const {
#ifndef N_UN
        return bji(j, i, alpha);
#else
        return bji(j, i, alpha) / (1 - bji(0, i, alpha));
#endif
    }

    virtual void info() const = 0;

    virtual vector<std::pair<int, int>> sample() const = 0;

    /**
     * P(Y=j|X=i)
     */
    virtual double bji(const int j, const int i, const double alpha) const = 0;

    virtual bool hasAlpha() const { return true; }

    /**
     * first and second order derivations of log b_{ji}(alpha) to alpha
     */
    virtual std::pair<double, double> getLGrad(const int i, const int j,
                                               const double alpha) const = 0;

    void check() const {
        for (int i = 0; i <= conf_->mx_tc; i++) {
            double sum_col = 0;
            for (int j = 0; j <= i; j++) sum_col += bji(j, i, 0.1);
            assert(std::abs(sum_col - 1) < 1e-9);
        }
    }
};

#endif /* __GRAPH_SAMPLER_H__ */
