/**
 * Copyright (C) 2017 by J. Zhao
 * Distributed under terms of the MIT license.
 */

#include "em_alg.h"

void EM::init() {
    // get g
    g_vec_ = sampler_->sample();

    // initialize z
    for (const auto& pr : g_vec_) z_[pr.first].resize(K_ + 2);

    // initialize theta
    theta_.resize(K_ + 2);
#ifndef N_UN
    randutils::initUniform(theta_);
#else
    randutils::initUniform(theta_, 1);
#endif

    // initialize alpha
    randutils::default_rng rng;
    alpha_ = 0.001 * (rng.uniform() + 1);
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
    double d1 = -alpha_, d2 = -1;
    // double d1 = 0, d2 = 0;
    for (const auto & [ k, j, z ] : non_zero_z_) {
        auto grad = sampler_->getLGrad(k, j, alpha_);
        d1 += z * grad.first;
        d2 += z * grad.second;
    }
    alpha_ -= d1 / d2;

    // printf("%2d a: %2.2e d1: %+.2e d2: %+.2e\n", iter, alpha_, d1, d2);
    if (alpha_ < 0 || d2 > 0) return false;
    return std::abs(d1 / d2) < conf_->eps_alpha;
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
        bool suc_alpha = !sampler_->hasAlpha();
        if (!suc_alpha) {
            for (int iter_a = 0; iter_a < conf_->mx_iter_alpha; iter_a++)
                if (MStepAlpha()) {
                    suc_alpha = true;
                    break;
                }
        }
        if (!suc_alpha) return false;
        if (MStepTheta()) return true;

        // if (iter % 5 == 0)
        // printf("%.4e / %.4e\n", getLikelihood(), getLikelihoodTruth());
        // printf("%.4e\n", getLikelihood());
    }
    return false;
}

void EM::scale() {
    theta_[0] = 0;
    vector<double> q(K_ + 1);
    for (int k = 0; k <= K_; k++) {
        q[k] = 0;
        for (int i = int(std::pow(2, k)); i < 2 << k; i++)
            q[k] += sampler_->bji(0, i, alpha_);
        q[k] /= std::pow(2, k);
        theta_[k + 1] /= 1 - q[k];
        theta_[0] += theta_[k + 1];
    }
    double q_avg = 0, g = 0;
    for (int k = 0; k <= K_; k++) {
        theta_[k + 1] /= theta_[0];
        q_avg += theta_[k + 1] * q[k];
    }
    for (auto& pr : g_vec_) g += pr.second;
    theta_[0] = g / (1 - q_avg);
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
