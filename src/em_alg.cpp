/**
 * Copyright (C) 2017 by J. Zhao
 * Distributed under terms of the MIT license.
 */

#include "em_alg.h"

void EM::init() {
    g_vec_ = sampler_->sample();
    // reserve space for z_ji
    for (auto&& pr : g_vec_) z_[pr.first].resize(conf_->mx_tc - pr.first + 1);

#ifndef N_UN
    int offset = 0;
#else
    int offset = 1;
#endif
    // initialize theta
    theta_.resize(conf_->mx_tc + 1);
    randutils::initUniform(theta_, offset);
    // initialize alpha
    randutils::default_rng rng;
    alpha_ = 0.01 * (rng.uniform() + 1);
}

/**
 * update [z_ij]
 */
void EM::EStep() {
    non_zero_z_.clear();
    for (auto&& pr : g_vec_) {
        int j = pr.first, g = pr.second;
        double norm = 0;
        for (int i = j; i <= conf_->mx_tc; i++) {
            double tmp = sampler_->pji(j, i, alpha_) * theta_[i];
            norm += tmp;
            z_[j][i - j] = g * tmp;  // g_j * p_{i|j}
        }
        if (norm > 1e-9) {
            for (int i = j; i <= conf_->mx_tc; i++) {
                z_[j][i - j] /= norm;
                double z = z_[j][i - j];
                if (i <= conf_->mx_i && z > 1e-9)
                    non_zero_z_.emplace_back(i, j, z);
            }
        } else {
            std::fill(z_[j].begin(), z_[j].end(), 0);
        }
    }
}

bool EM::MStepTheta() {
    // store the old parameters before we update them
    auto theta_pre = theta_;
    std::fill(theta_.begin(), theta_.end(), 0.0);
    // now update theta
    double norm = 0;
    for (auto&& pr : g_vec_) {
        int j = pr.first;
        for (int i = j; i <= conf_->mx_tc; i++) {
            double z = z_[j][i - j];
            theta_[i] += z;
            norm += z;
        }
    }
    // normalize theta and calculate diff
    int i = 0;
    double diff = 0;
    for (auto& theta : theta_) {
        theta /= norm;
        diff += std::abs(theta - theta_pre[i++]);
    }
    return diff < conf_->eps_theta;
}

/**
 * with L1 regularizer: max Q(alpha) - 1/2 alpha^2
 */
bool EM::MStepAlpha() {
    int iter = 0;
    while (iter++ < conf_->mx_iter_alpha) {
        double d1 = -alpha_, d2 = -1;
        for (auto&& e : non_zero_z_) {
            int i = std::get<0>(e), j = std::get<1>(e);
            double z = std::get<2>(e);
            auto grad = sampler_->getLGrad(i, j, alpha_);
            d1 += z * grad.first;
            d2 += z * grad.second;
        }
        if (std::abs(d1) < 1e-3) break;
        double alpha_pre = alpha_;
        alpha_ -= d1 / d2;
        // printf("%2d a: %2.4e d1: %+.4e d2: %+.4e\n", iter, alpha_, d1, d2);
        if (alpha_ < 0 || d2 > 0) {
            alpha_ = 0.0001;
            break;
        }
        if (std::abs(alpha_pre - alpha_) < conf_->eps_alpha) break;
    }
    return iter < conf_->mx_iter_alpha;
}

bool EM::exec() {
    // The algorithm is not sensitive to alpha. So let EM converge to an
    // approximate theta first.
    for (int iters = 0; iters < conf_->mx_iter_theta; iters++) {
        EStep();
        if (MStepTheta()) break;
    }
    // Now, we adjust alpha.
    for (int iter = 0; iter < conf_->mx_iter_theta; iter++) {
        EStep();
        bool suc_alpha = sampler_->hasAlpha() ? MStepAlpha() : true;
        bool suc_theta = MStepTheta();
        if (suc_theta && suc_alpha) return true;
    }
    return false;
}

void EM::scale() {
    theta_[0] = 0;
    for (int i = 1; i <= conf_->mx_tc; i++) {
        theta_[i] /= 1 - sampler_->bji(0, i, alpha_);
        theta_[0] += theta_[i];
    }
    double q = 0, g = 0;
    for (int i = 1; i <= conf_->mx_tc; i++) {
        theta_[i] /= theta_[0];
        q += theta_[i] * sampler_->bji(0, i, alpha_);
    }
    for (auto& pr : g_vec_) g += pr.second;
    theta_[0] = g / (1 - q);
}

vector<double> EM::get() {
#ifdef N_UN
    scale();
#endif
    vector<double> rst;
    rst.reserve(theta_.size() + 1);
    rst.push_back(alpha_);
    rst.insert(rst.end(), theta_.begin(), theta_.end());
    return rst;
}
