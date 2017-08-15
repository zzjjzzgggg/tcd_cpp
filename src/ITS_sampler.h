/**
 * Copyright (C) by J.Z. (07/14/2017 15:32)
 * Distributed under terms of the MIT license.
 */

#ifndef __ITS_SAMPLER_H__
#define __ITS_SAMPLER_H__

#include "graph_sampler.h"

/**
 * Identical Triangle Sampling by Edge Sampling. Each edge is idependently
 * sampled with probability p, and hence a triangle is sampled with
 * probability p^3.
 */
class ITSSampler : public Sampler {
public:
    ITSSampler(const Sampler::Config* conf) : Sampler(conf) {}

    void info() const override {
#ifndef N_UN
        printf("ITS: graph size is known\n");
#else
        printf("ITS: graph size is unknown\n");
#endif
    }

    vector<std::pair<int, int>> sample() const override {
        UGraph G(GraphType::MULTI);
        rngutils::default_rng rng;
        for (const auto & [ src, dst ] : edges_)
            if (rng.uniform() < conf_->p_edge) G.addEdge(src, dst);
        G.defrag();
        return statTrids(G);
    }

    double bji(const int j, const int i, const double alpha) const override {
        return i >= j ? bb(j, i, alpha) : 0;
    }

    std::pair<double, double> getLGrad(const int k, const int j,
                                       const double alpha) const override {
        if (j == 0 && k < 0) return std::make_pair(0.0, 0.0);
        double d1 = 0, d2 = 0, den = 0;
        for (int i = std::max(j, int(std::pow(2, k))); i < 2 << k; i++) {
            double c = cji(j, i, alpha);
#ifndef N_UN
            auto[d1l, d2l] = dLogBji(i, j, alpha);
#else
            auto[d1l, d2l] = dLogAji(i, j, alpha);
#endif
            d1 += d1l * c;
            d2 += (d2l + d1l * d1l) * c;
            den += c;
        }
        d1 /= den;
        d2 = -d1 * d1 + d2 / den;
        return std::make_pair(d1, d2);
    }

    // bool hasAlpha() const override { return false; }

    /**
     * return [d log(b) / d alpha, d^2 log(b) /d alpha^2]
     */
    std::pair<double, double> dLogBji(const int i, const int j,
                                      const double alpha) const {
        double d1 = 0, d2 = 0;
        for (int s = 1; s < i; s++) {
            double e1 = s < j ? s / (s * alpha + p_tri_)
                              : (s - j) / ((s - j) * alpha + 1 - p_tri_),
                   e2 = s / (s * alpha + 1);
            d1 += e1 - e2;
            d2 += -e1 * e1 + e2 * e2;
        }
        return std::make_pair(d1, d2);
    }

    /**
     * return [d log(a) / d alpha, d^2 log(a) /d alpha^2]
     */
    std::pair<double, double> dLogAji(const int i, const int j,
                                      const double alpha) const {
        double d1l = 0, d2l = 0, s1 = 0, s2 = 0;
        for (int s = 1; s < i; s++) {
            double e0 = s / (s * alpha + 1 - p_tri_),
                   e1 = s < j ? s / (s * alpha + p_tri_)
                              : (s - j) / ((s - j) * alpha + 1 - p_tri_),
                   e2 = s / (s * alpha + 1);
            d1l += e1 - e2;
            d2l += -e1 * e1 + e2 * e2;
            s1 += e0 - e2;
            s2 += -e0 * e0 + e2 * e2;
        }
        double b0 = bb(0, i, alpha), r = b0 / (1 - b0), d1 = d1l + r * s1,
               d2 = d2l + r * ((1 + r) * s1 * s1 + s2);

        return std::make_pair(d1, d2);
    }

}; /* class ITSSampler */

#endif /* __ITS_SAMPLER_H__ */
