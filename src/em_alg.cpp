/**
 * Copyright (C) 2017 by J. Zhao
 * Distributed under terms of the MIT license.
 */

#include "em_alg.h"

void EM::init(vector<std::pair<int, int>>&& g_vec) {
    g_vec_ = g_vec;
    // find the empirical maximum tradic cardinality
    M_ = std::max_element(g_vec.begin(), g_vec.end())->first;
    // reserve space for z_ji
    for (auto&& pr : g_vec) z_[pr.first].resize(conf_->mx_tc - pr.first + 1);

    // initialize theta
    // NOTE: propally initializing theta can improve convergence speed, and
    // estimation accuracy significantly.
    theta_vec_.resize(conf_->mx_tc + 1, 0);
    double norm = 0;
    randutils::default_rng rng;
    for (int i = 0; i <= conf_->mx_tc; i++) {
        // HepTh: 0.2
        // theta_vec_[i] = std::pow(i + 1, -0.2 * rng.uniform() - 0.9);
        theta_vec_[i] = rng.uniform();
        norm += theta_vec_[i];
    }
    for (auto& theta : theta_vec_) theta /= norm;
    alpha_ = 0.001 * (rng.uniform() + 1);  // initialize alpha_
    // alpha_ = 0.01;
}

/**
 * update [z_ij]
 */
void EM::EStep() {
    nonzero_z_.clear();
    for (const auto& pr : g_vec_) {
        int j = pr.first, g = pr.second;
        double norm = 0;
        for (int i = j; i <= conf_->mx_tc; i++) {
            double tmp = theta_vec_[i] * bji(j, i);
            norm += tmp;
            z_[j][i - j] = g * tmp;  // g_j * p_{i|j}
        }
        for (int i = j; i <= conf_->mx_tc; i++) {
            z_[j][i - j] = norm > 1e-9 ? z_[j][i - j] / norm : 0;
            if (z_[j][i - j] > 1e-6)
                nonzero_z_.emplace_back(i, j, z_[j][i - j]);
        }
    }
}

bool EM::MStepTheta() {
    // store the old parameters before we update them
    theta_pre_vec_ = theta_vec_;
    std::fill(theta_vec_.begin(), theta_vec_.end(), 0);
    // now update theta
    double norm = 0;
    for (auto& pr : g_vec_) {
        int j = pr.first;
        for (int i = j; i <= conf_->mx_tc; i++) {
            double z = z_[j][i - j];
            theta_vec_[i] += z;
            norm += z;
        }
    }
    // normalize theta and calculate change
    int i = 0;
    double change = 0;
    for (auto& theta : theta_vec_) {
        theta /= norm;
        change += std::abs(theta - theta_pre_vec_[i++]);
    }
    return change <= conf_->eps_theta;
}

// TODO: remember to change the d1 and d2 for US
bool EM::MStepAlpha() {
    int iter = 0;
    double p = conf_->p_tri;
    // Newton iterations
    while (iter++ < conf_->mx_iter_alpha) {
        double d1 = 0, d2 = 0;
        for (const auto& triple : nonzero_z_) {
            int i = std::get<0>(triple), j = std::get<1>(triple);
            double z = std::get<2>(triple), e1 = 0, e2 = 0;
            for (int s = 0; s < i; s++) {
                double e = s < j ? s / (s * alpha_ + p)
                                 : (s - j) / ((s - j) * alpha_ + 1 - p);
                e1 += e - s / (s * alpha_ + 1);
                e2 += -e * e + std::pow(s / (s * alpha_ + 1), 2);
            }
            d1 += z * e1;
            d2 += z * e2;
        }
        double alpha_pre = alpha_;
        alpha_ -= d1 / d2;
        if (std::abs(alpha_pre - alpha_) <= conf_->eps_alpha) break;
    }
    return iter < conf_->mx_iter_alpha;
}

bool EM::exec() {
    for (int iter = 0; iter < conf_->mx_iter_theta; iter++) {
        EStep();
#ifdef METHOD_FS  // for FS, there is no need to estimate alpha
        if (MStepTheta()) return true;
#else
        if (MStepAlpha() && MStepTheta()) return true;
#endif
    }
    return false;
}

vector<double> EM::getResult() const {
    vector<double> rst;
    rst.reserve(theta_vec_.size() + 1);
    rst.push_back(alpha_);
    rst.insert(rst.end(), theta_vec_.begin(), theta_vec_.end());
    return rst;
}
