/**
 * Copyright (C) 2017 by J. Zhao
 * Distributed under terms of the MIT license.
 */

#ifndef __EMALG_H__
#define __EMALG_H__

#include "stdafx.h"
#include "graph_sampler.h"

/**
 * The general EM algorithm for estimating TCD.
 *
 * 1. when graph size is known:
 *    The sampler should provide vector g with g[0] calibrated.
 *
 * 2. when graph size is unknown:
 *    The sampler should provide vector g with g[0] removed.
 *
 */
class EM {
public:
    struct Config {
        int mx_tc;          // maximum triadic cardinality
        int mx_i;           // to speedup estimating alpha, we set another bound
                            // on tc, that is far less than mx_tc.
        int mx_iter_theta;  // maximum iterations for theta and alpha_
        int mx_iter_alpha;
        double eps_theta;  // threshold for theta and alpha_
        double eps_alpha;

        void echo() const {
            printf(
                "EM config:\n"
                "  max i: %d\n"
                "  max iterations (theta): %d\n"
                "  max iterations (alpha): %d\n"
                "  epsilon (theta): %g\n"
                "  epsilon (alpha): %g\n\n",
                mx_i, mx_iter_theta, mx_iter_alpha, eps_theta, eps_alpha);
        }
    };

private:
    const Config* conf_;
    const Sampler* sampler_;

    int K_, mx_k_;  // largest bin, k = 0, 1, ..., K_, and largest bin for
                    // estimating alpha
    double alpha_;
    vector<double> theta_;  // theta_k = p_k

    vector<std::pair<int, int>> g_vec_;

    // because sampled cardinality j is usually sparsely populated, it is better
    // to use a map, rather than a vector.
    unordered_map<int, vector<double>> z_;

    // only a small fraction of non-zero z is used for estimating alpha
    vector<std::tuple<int, int, double>> non_zero_z_;

private:
    inline int getK(const int j) const {
        return j > 0 ? int(std::floor(std::log2(j))) : -1;
    }

    void EStep();
    bool MStepTheta();
    bool MStepAlpha();
    void scale();

public:
    EM(const Config* conf, const Sampler* sampler)
        : conf_(conf), sampler_(sampler) {
        K_ = getK(conf_->mx_tc);
        mx_k_ = getK(conf_->mx_i);
    }

    void init();

    bool exec();

    double getLikelihood() const;

    double getLikelihoodTruth() const;
    /**
     * the first element is hat_alpha, remainings are hat_theta
     */
    vector<double> get();

}; /* EM */

#endif /* __EMALG_H__ */
