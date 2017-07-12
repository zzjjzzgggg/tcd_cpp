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

        void echo() const {
            printf(
                "Sampling config:\n"
                "  graph: %s\n"
                "  maximum TC: %d\n"
                "  edge sample rate: %g\n"
                "  node sample rate (US): %g\n\n",
                graph_fnm.c_str(), mx_tc, p_edge, p_node);
        }
    };

private:
    const Config* conf_;
    DGraph G_;

private:
    // state triads histogram in a sampled graph G
    vector<std::pair<int, int>> statTrids(const DGraph& G) const;

    // count the number of triads of a node in a probabilistic manner
    int countTriadsWithProb(const int node, const DGraph& G,
                            const unordered_set<int>& marked_nds) const;

public:
    Sampler(const Config* conf) : conf_(conf) {
        G_ = loadEdgeList<DGraph>(conf_->graph_fnm, GraphType::MULTI);
        printf("N: %d, E: %d\n", G_.getNodes(), G_.getEdges());
    }

    void getGroundtruth(const string& outfnm) const;

    /**
     * Identical Triangle Sampling by Edge Sampling. Each edge is idependently
     * sampled with probability p, and hence a triangle is sampled with
     * probability p^3.
     */
    vector<std::pair<int, int>> IS() const;

    /**
     * deprecated.
     * This method has the distadvantage of both IS and FS.
     */
    vector<std::pair<int, int>> US() const;

    /**
     * Similar to Flow Sampling. Randomly sample a collection of nodes. Keep
     * each triangle related to a sampled node.
     */
    vector<std::pair<int, int>> FS() const;

    /**
     * Triangle Sample-and-Hold.
     * If a new triangle is formed, we randomly pick a node from this triangle,
     * and mark this node. For a marked node, its related triangles are all kept
     * since it is marked.
     * The above process is only invoked for a newly formed triangle, which must
     * consists of three unmarked nodes. Otherwise, if a node has already been
     * marked, then this triangle is already kept.
     */
    vector<std::pair<int, int>> SH() const;
};

#endif /* __GRAPH_SAMPLER_H__ */
