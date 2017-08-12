/**
 * Copyright (C) by J.Z. (07/14/2017 15:56)
 * Distributed under terms of the MIT license.
 */

#ifndef __SH_SAMPLER_H__
#define __SH_SAMPLER_H__

#include "graph_sampler.h"

/**
 * Triangle Sample-and-Hold.
 * If a new triangle is formed, we randomly pick a node from this triangle,
 * and mark this node. For a marked node, its related triangles are all kept
 * since it is marked.
 * The above process is only invoked for a newly formed triangle, which must
 * consists of three unmarked nodes. Otherwise, if a node has already been
 * marked, then this triangle is already kept.
 */
class SHSampler : public Sampler {
protected:
    /**
     * Beta-Geometric distribution
     */
    inline double bg(const int j, const int i, const double alpha) const {
        // double rst = p_tri_;
        // for (int s = 0; s < i - j; s++)
        //     rst *= (1 - p_tri_ + s * alpha) / (1 + s * alpha);
        // return rst / (1 + (i - j) * alpha);

        return SpecFun::BetaGeometric(i - j + 1, p_tri_ / alpha,
                                      (1 - p_tri_) / alpha);
    }

public:
    SHSampler(const Sampler::Config* conf) : Sampler(conf) {}

    void info() const override { printf("SH: p_tri = %.3f\n", p_tri_); }

    vector<std::pair<int, int>> sample() const override {
        randutils::default_rng rng;
        unordered_set<int> candidates, held_nodes, held_nbrs;

        for (int node : nodes_)
            if (rng.uniform() < conf_->p_node) candidates.insert(node);

        UGraph G(GraphType::MULTI);  // sampled graph

        auto is_held = [&](const int node) {
            return held_nodes.find(node) != held_nodes.end() ||
                   held_nbrs.find(node) != held_nbrs.end();
        };

        // hold a node
        auto hold = [&](const int node) {
            // if the node is not a candidate, exit;
            if (candidates.find(node) == candidates.end()) return;
            // if the node is already held, exit;
            if (held_nodes.find(node) != held_nodes.end()) return;
            // otherwise, hold this node, and mark its neighbors
            held_nodes.insert(node);
            for (int d = 0; d < G[node].getDeg(); d++)
                held_nbrs.insert(G[node].getNbr(d));
        };

        for (auto&& edge : edges_) {
            int src = edge.first, dst = edge.second;
            if (is_held(src) || is_held(dst) || rng.uniform() < conf_->p_edge) {
                // if one of the two ends is held, the edge is kept; otherwise,
                // the edge is kept with probability p_edge.
                G.addEdge(src, dst);
                // If the newly added edge forms triangles, check whether these
                // triangles contain candidate nodes; if they contain, then
                // update held nodes and neighbors.
                std::unordered_set<int> src_nbrs, com_nbrs;
                for (int d = 0; d < G[src].getDeg(); d++)
                    src_nbrs.insert(G[src].getNbr(d));
                for (int d = 0; d < G[dst].getDeg(); d++) {
                    int nbr = G[dst].getNbr(d);
                    if (src_nbrs.find(nbr) != src_nbrs.end())
                        com_nbrs.insert(nbr);
                }
                if (!com_nbrs.empty()) {
                    // triangles are formed
                    hold(src);
                    hold(dst);
                    for (auto node : com_nbrs) hold(node);
                }
            }
        }
        G.defrag();
        // printf("held nodes: %lu\n", held_nodes.size());
        return statTrids(G, held_nodes);
    }

    // bool hasAlpha() const override { return false; }

    double bji(const int j, const int i, const double alpha) const override {
        return j == 0 ? 1 - conf_->p_node + conf_->p_node * bb(0, i, alpha)
                      : conf_->p_node * bg(j, i, alpha);
    }

    std::pair<double, double> getLGrad(const int i, const int j,
                                       const double alpha) const override {
        if (j == 0) {
            double pn = conf_->p_node, d1_log_bb = 0, d2_log_bb = 0;
            for (int s = 1; s < i; s++) {
                double e1 = s / (s * alpha + 1.0 - p_tri_),
                       e2 = s / (s * alpha + 1.0);
                d1_log_bb += e1 - e2;
                d2_log_bb += -e1 * e1 + e2 * e2;
            }
            if (d1_log_bb < 1e-6) return std::make_pair(0, 0);
            double bb_0i = bb(0, i, alpha), b_0i = 1 - pn + pn * bb_0i,
                   d1_log_b = pn / b_0i * bb_0i * d1_log_bb;
            double d2_bb = d1_log_bb + d2_log_bb / d1_log_bb,
                   d2_log_b = d1_log_b * (-d1_log_b + d2_bb);
            return std::make_pair(d1_log_b, d2_log_b);
        } else {
            double e0 = (i - j) / (1 + (i - j) * alpha), d1_log_bg = -e0,
                   d2_log_bg = e0 * e0;
            for (int s = 1; s < i - j; s++) {
                double e1 = s / (s * alpha + 1 - p_tri_),
                       e2 = s / (s * alpha + 1);
                d1_log_bg += e1 - e2;
                d2_log_bg += -e1 * e1 + e2 * e2;
            }
            return std::make_pair(d1_log_bg, d2_log_bg);
        }
    }

};  // end of class SHSampler

#endif /* __SH_SAMPLER_H__ */
