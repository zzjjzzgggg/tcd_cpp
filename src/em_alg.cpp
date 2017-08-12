/**
 * Copyright (C) 2017 by J. Zhao
 * Distributed under terms of the MIT license.
 */

#include "em_alg.h"

void EM::init() {
    g_vec_ = sampler_->sample();
    for (const auto& pr : g_vec_) z_[pr.first].resize(K_ + 2);

#ifndef N_UN
    int start = 0;
#else
    int start = 1;
#endif
    // initialize theta
    theta_.resize(K_ + 2);
    randutils::initUniform(theta_, start);
    // initialize alpha
    randutils::default_rng rng;
    alpha_ = 0.01 * (rng.uniform() + 1);
}

/**
 * update [z_{jk}]
 */
void EM::EStep() {
    non_zero_z_.clear();
    for (const auto & [ j, g ] : g_vec_) {
        double norm = 0;
        for (int k = getK(j); k <= K_; k++) {
            z_[j][k + 1] = sampler_->pjk(j, k, alpha_) * theta_[k + 1];
            norm += z_[j][k + 1];
        }
        if (norm < 1e-9) continue;
        for (int k = getK(j); k <= K_; k++) {
            double& z = z_[j][k + 1];
            z *= g / norm;
            if (k <= mx_k_ && z > 1e-9) non_zero_z_.emplace_back(k, j, z);
        }
    }
}

bool EM::MStepTheta() {
    // store the old parameters before we update them
    auto theta_pre = theta_;
    std::fill(theta_.begin(), theta_.end(), 0.0);
    // now update theta
    double norm = 0;
    for (const auto& pr : g_vec_) {
        const int j = pr.first;
        for (int k = getK(j); k <= K_; k++) {
            const double z = z_[j][k + 1];
            theta_[k + 1] += z;
            norm += z;
        }
    }
    // normalize theta and calculate diff
    int k = 0;
    double diff = 0;
    for (auto& theta : theta_) {
        theta /= norm;
        diff += std::abs(theta - theta_pre[k++]);
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
        for (const auto & [ k, j, z ] : non_zero_z_) {
            auto grad = sampler_->getLGrad(k, j, alpha_);
            d1 += z * grad.first;
            d2 += z * grad.second;
        }
        alpha_ -= d1 / d2;
        double diff = std::abs(d1 / d2);
        // printf("%2d a: %2.2e d1: %+.2e d2: %+.2e dif: %.4e\n", iter, alpha_,
        //        d1, d2, diff);
        if (alpha_ < 0 || d2 > 0) {
            alpha_ = 0.0001;
            break;
        }
        if (diff < conf_->eps_alpha) break;
    }
    return iter < conf_->mx_iter_alpha;
}

bool EM::exec() {
    // The algorithm is not sensitive to alpha. So let EM converge to an
    // approximate theta first.
    // for (int iters = 0; iters < conf_->mx_iter_theta; iters++) {
    //     EStep();
    //     if (MStepTheta()) break;
    // }
    // Now, we adjust alpha.
    for (int iter = 0; iter < conf_->mx_iter_theta; iter++) {
        EStep();
        if ((sampler_->hasAlpha() ? MStepAlpha() : true) && MStepTheta())
            return true;
        // if (iter % 10 == 0)
        //     printf("%.4e / %.4e\n", getLikelihood(), getLikelihoodTruth());
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

double EM::getLikelihood() const {
    double likelihood = 0;
    for (const auto & [ j, g ] : g_vec_) {
        double tmp = 0;
        for (int k = getK(j); k <= K_; k++)
            tmp += theta_[k + 1] * sampler_->pjk(j, k, alpha_);
        likelihood += g * std::log(tmp);
    }
    return likelihood;
}

double EM::getLikelihoodTruth() const {
    auto truth = ioutils::loadVec<double>(
        "/dat/workspace/tcd_journal/hepth/theta_W2K_binned.dat", 1);
    double likelihood = 0;
    for (const auto & [ j, g ] : g_vec_) {
        double tmp = 0;
        for (int k = getK(j); k <= K_; k++)
            tmp += truth[k + 1] * sampler_->pjk(j, k, alpha_);
        likelihood += g * std::log(tmp);
    }
    return likelihood;
}
