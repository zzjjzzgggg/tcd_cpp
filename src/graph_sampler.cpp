/**
 * Copyright (C) by J.Z. (06/30/2017 13:19)
 * Distributed under terms of the MIT license.
 */

#include "graph_sampler.h"

vector<std::pair<int, int>> Sampler::statTrids(const DGraph& G) const {
    std::unordered_map<int, int> tc_to_num_nodes;
    for (auto&& ni = G.beginNI(); ni != G.endNI(); ni++)
        tc_to_num_nodes[countNodeDirTriads(ni->first, G)]++;
    tc_to_num_nodes[0] += G_.getNodes() - G.getNodes();

    std::vector<std::pair<int, int>> triads_hist;
    for (auto&& it : tc_to_num_nodes)
        triads_hist.emplace_back(it.first, it.second);
    std::sort(triads_hist.begin(), triads_hist.end());
    return triads_hist;
}

void Sampler::getGroundtruth(const string& outfnm) const {
    auto g_vec = statTrids(G_);
    vector<std::tuple<int, double, int>> dist;
    dist.reserve(g_vec.size());
    double n = G_.getNodes();
    for (auto& pr : g_vec)
        dist.emplace_back(pr.first, pr.second / n, pr.second);
    ioutils::saveTupleVec(dist, outfnm, true, "{}\t{:.6e}\t{}\n");
}

int Sampler::countTriadsWithProb(const int node, const DGraph& G,
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
                if (src_marked || marked_nds.find(*dst_it) != marked_nds.end())
                    for (int k = 0; k < src_nbr_cnt[*dst_it]; k++)
                        if (rng.uniform() > conf_->p_edge) valid_edges--;
                triads += valid_edges * nbr_cnt[*src_it] * nbr_cnt[*dst_it];
            }
        }
    }
    return triads;
}

vector<std::pair<int, int>> Sampler::IS() const {
    DGraph G(GraphType::MULTI);
    randutils::default_rng rng;
    for (auto&& ei = G_.beginEI(); ei != G_.endEI(); ei++)
        if (rng.uniform() <= conf_->p_edge)
            G.addEdge(ei.getSrcID(), ei.getDstID());
    G.optimize();
    return statTrids(G);
}

vector<std::pair<int, int>> Sampler::US() const {
    randutils::default_rng rng;
    // get node samples
    unordered_set<int> nd_samples;
    for (auto&& ni = G_.beginNI(); ni != G_.endNI(); ni++)
        if (rng.uniform() < conf_->p_node) nd_samples.insert(ni->first);

    // sample graph
    DGraph G(GraphType::MULTI);
    for (auto&& ei = G_.beginEI(); ei != G_.endEI(); ei++) {
        int src = ei.getSrcID(), dst = ei.getDstID();
        if (nd_samples.find(src) != nd_samples.end() ||
            nd_samples.find(dst) != nd_samples.end() ||
            rng.uniform() < conf_->p_edge)
            G.addEdge(src, dst);
    }
    G.optimize();

    // TC -> # nodes
    std::unordered_map<int, int> tc_to_num_nodes;
    for (int nd : nd_samples)
        tc_to_num_nodes[countTriadsWithProb(nd, G, nd_samples)]++;

    // calibrate
    tc_to_num_nodes[0] += G_.getNodes() - G.getNodes();

    // convert to pair vector
    std::vector<std::pair<int, int>> g_vec;
    for (auto& it : tc_to_num_nodes) g_vec.emplace_back(it.first, it.second);
    std::sort(g_vec.begin(), g_vec.end());
    return g_vec;
}

vector<std::pair<int, int>> Sampler::FS() const {
    int num_sampled_nodes = 0;
    randutils::default_rng rng;
    std::unordered_map<int, int> tc_to_num_nodes;  // TC -> # nodes
    for (auto&& ni = G_.beginNI(); ni != G_.endNI(); ni++)
        if (rng.uniform() < conf_->p_node) {
            tc_to_num_nodes[countNodeDirTriads(ni->first, G_)]++;
            num_sampled_nodes++;
        }

    // calibrate
    tc_to_num_nodes[0] += G_.getNodes() - num_sampled_nodes;

    std::vector<std::pair<int, int>> g_vec;
    for (auto& it : tc_to_num_nodes) g_vec.emplace_back(it.first, it.second);
    std::sort(g_vec.begin(), g_vec.end());
    return g_vec;
}

vector<std::pair<int, int>> Sampler::SH() const {
    randutils::default_rng rng;
    unordered_set<int> nd_samples, nbrs_of_samples;

    // sampled graph
    DGraph G(GraphType::MULTI);
    for (auto&& ei = G_.beginEI(); ei != G_.endEI(); ei++) {
        int src = ei.getSrcID(), dst = ei.getDstID();
        if (nd_samples.find(src) != nd_samples.end() ||
            nd_samples.find(dst) != nd_samples.end() ||
            rng.uniform() < conf_->p_edge)
            G.addEdge(src, dst);
    }
    G.optimize();

    // TC -> # nodes
    std::unordered_map<int, int> tc_to_num_nodes;
    for (int nd : nd_samples)
        tc_to_num_nodes[countTriadsWithProb(nd, G, nd_samples)]++;

    // calibrate
    tc_to_num_nodes[0] += G_.getNodes() - G.getNodes();

    // convert to pair vector
    std::vector<std::pair<int, int>> g_vec;
    for (auto& it : tc_to_num_nodes) g_vec.emplace_back(it.first, it.second);
    std::sort(g_vec.begin(), g_vec.end());
    return g_vec;
}
