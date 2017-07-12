/**
 * Copyright (C) 2017 by J. Zhao
 * Distributed under terms of the MIT license.
 */

#ifndef __EMALG_H__
#define __EMALG_H__

#include "stdafx.h"

class EM {
public:
    struct Config {
        int mx_tc;          // maximum triadic cardinality
        int mx_iter_theta;  // maximum iterations for theta and alpha_
        int mx_iter_alpha;
        double eps_theta;  // threshold for theta and alpha_
        double eps_alpha;
        double p_tri;   // probability of sampling a triangle
        double p_node;  // probability of sampling a node

        void echo() const {
            printf(
                "EM config:\n"
                "  prob. triangle: %g\n"
                "  prob. node: %g\n"
                "  max iterations (theta): %d\n"
                "  max iterations (alpha): %d\n"
                "  epsilon (theta): %g\n"
                "  epsilon (alpha): %g\n\n",
                p_tri, p_node, mx_iter_theta, mx_iter_alpha, eps_theta,
                eps_alpha);
        }
    };

private:
    const Config* conf_;
    vector<std::pair<int, int>> g_vec_;

    double alpha_;
    vector<double> theta_vec_;

    // for internal use
    int M_;  // empirical maximum triadic cardinality
    vector<double> theta_pre_vec_;
    // because sampled cardinality j is usually sparsely populated, it is better
    // to use a map, rather than a vector.
    unordered_map<int, vector<double>> z_;
    vector<std::tuple<int, int, double>> nonzero_z_;

private:
    inline double bb(const int i, const int j) const {
        return SpecFun::BetaBinomial(j, i, conf_->p_tri / alpha_,
                                     (1 - conf_->p_tri) / alpha_);
    }

    double bji(const int j, const int i) const {
#ifdef METHOD_IS
        return bb(i, j);
#elif METHOD_US
        return j == 0 ? (1 - conf_->p_node + conf_->p_node * bb(i, 0))
                      : conf_->p_node * bb(i, j);
#elif METHOD_FS
        if (i == 0 && j == 0)
            return 1;
        else if (i == j)
            return conf_->p_node;
        else if (j == 0)
            return 1 - conf_->p_node;
        return 0;
#elif METHOD_SH  // triad sample and hold
        return j == 0 ? bb(i, 0) : conf_->p_tri * bb(i - j, 0);
#endif
        printf("Error! Sampling method is not specified!\n");
    }

public:
    EM(const Config* conf) : conf_(conf) {}

    void init(vector<std::pair<int, int>>&& g_vec);

    void EStep();
    bool MStepTheta();
    bool MStepAlpha();
    bool exec();

    /**
     * the first element is hat_alpha, remainings are hat_theta
     */
    vector<double> getResult() const;
};

#endif /* __EMALG_H__ */
