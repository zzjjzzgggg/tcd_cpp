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
            printf("Sampling config:\n"
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

   public:
    /**
     * state triads histogram in a sampled graph G
     */
    template <class Container = vector<int>>
    vector<std::pair<int, int>> statTrids(
        const graph::undir::UGraph& G,
        const Container& nodes = Container()) const {
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
            rngutils::default_rng rng;
            rng.shuffle(edges_);
            unordered_set<int> nodes;
            for (const auto& [src, dst] : edges_) {
                nodes.insert(src);
                nodes.insert(dst);
            }
            nodes_.reserve(nodes.size());
            nodes_.insert(nodes_.end(), nodes.begin(), nodes.end());
            printf("N: %lu, E: %lu\n", nodes_.size(), edges_.size());
        }
        // default triangle sampling probability
        p_tri_ = std::pow(conf_->p_edge, 3);

#ifndef N_UN
        printf("graph size is known\n");
#else
        printf("graph size is unknown\n");
#endif
    }

    /**
     * c_{ji} represents b_{ji} or a_{ji}
     *
     *     B_k(j) = sum_{i=2^k}{2^{k+1}-1} b_ji/2^k
     *     A_k(j) = sum_{i=2^k}{2^{k+1}-1} a_ji/2^k
     */
    double pjk(const int j, const int k, const double alpha) const {
        if (j == 0 && k < 0) return 1;
        double sum = 0;
        for (int i = std::max(j, int(std::pow(2, k))); i < 2 << k; i++)
            sum += cji(j, i, alpha);
        return sum / std::pow(2, k);
    }

    virtual vector<std::pair<int, int>> sample() const = 0;

    /**
     * b_{ji} or a_{ji}, where
     *     b_{ji} = P(Y=j|X=i),
     *     a_{ji} = P(Y=j|X=i, Y>0)
     */
    double cji(const int j, const int i, const double alpha) const {
        if (i < j) return 0;
#ifndef N_UN
        return bji(j, i, alpha);
#else
        return bji(j, i, alpha) / (1 - bji(0, i, alpha));
#endif
    }

    virtual double bji(const int j, const int i, const double alpha) const = 0;

    virtual bool hasAlpha() const { return true; }

    /**
     * first and second order derivations of log B_{kj} or log A_{kj} to alpha
     */
    virtual std::pair<double, double> getLGrad(const int i, const int j,
                                               const double alpha) const {
        return std::pair<double, double>(0, 0);
    }

    void check() const {
        for (int i = 0; i <= conf_->mx_tc; i++) {
            double sum_col = 0;
            for (int j = 0; j <= i; j++) sum_col += cji(j, i, 0.1);
            assert(std::abs(sum_col - 1) < 1e-9);
        }
    }

}; /* class Sampler */

#endif /* __GRAPH_SAMPLER_H__ */
