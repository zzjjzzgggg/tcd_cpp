/**
 * Copyright (C) by J.Z. (07/02/2017 09:55)
 * Distributed under terms of the MIT license.
 */

#ifndef __CRLB_H__
#define __CRLB_H__

#include "stdafx.h"

#include <Eigen/Dense>

using namespace Eigen;

/**
 * ITS (ITS-color)
 */
class CRLB {
public:
    struct Config {
        int W;  // maximum triadic cardinality
        int n;  // number of nodes/sampled nodes
        bool known_size;
        double alpha;
        double p_tri;
        double p_nd;       // node sampling rate
        string theta_fnm;  // [theta_0, ..., theta_W]

        void echo() {
            printf(
                "W: %d\n"
                "n: %d\n"
                "known graph size: %d\n"
                "alpha: %g\n"
                "p_tri: %g\n"
                "theta: %s\n\n",
                W, n, known_size, alpha, p_tri, theta_fnm.c_str());
        }
    };

protected:
    const Config* conf_;
    // number of samples: (1) if we known graph size, n = graph size; (2) if
    // not, n is the number of nodes that are observed to have at least one
    // triangle.
    double n_ = 0;

    VectorXd theta_, crlb_, phi_;

public:
    CRLB(const Config* conf) : conf_(conf) {}

    void init() {
        ioutils::TSVParser ss(conf_->theta_fnm);
        if (conf_->known_size) {
            n_ = conf_->n;
            theta_ = VectorXd::Zero(conf_->W + 1);
            while (ss.next()) theta_(ss.get<int>(0)) = ss.get<double>(1);
        } else {
            theta_ = VectorXd::Zero(conf_->W);
            if (ss.next()) n_ = ss.get<double>(1);  // #nodes have triangles
            double q = 0;  // q(theta^+, alpha) = P(Y = 0 | X > 0)
            while (ss.next()) {
                int i = ss.get<int>(0);
                double theta = ss.get<double>(1);
                theta_(i - 1) = theta;
                q += bji(0, i) * theta;
            }
            n_ *= 1 - q;

            phi_ = VectorXd::Zero(conf_->W);
            for (int i = 0; i < conf_->W; i++)
                phi_(i) = theta_(i) * (1 - bji(0, i + 1));
            phi_ /= phi_.sum();
        }
        // std::cout << theta_ << std::endl;
    }

    /**
     * P(Y=j|X=i)
     */
    virtual double bji(const int j, const int i) const {
        return SpecFun::BetaBinomial(j, i, conf_->p_tri / conf_->alpha,
                                     (1 - conf_->p_tri) / conf_->alpha);
    }

    /**
     * ITS sampling matrix: sample each triangle identically
     * Known graph size
     */
    MatrixXd getB() const {
        MatrixXd B = MatrixXd::Zero(conf_->W + 1, conf_->W + 1);
        for (int j = 0; j <= conf_->W; j++)
            for (int i = j; i <= conf_->W; i++) B(j, i) = bji(j, i);
        return B;
    }

    /**
     * Unknown graph size
     */
    MatrixXd getA() const {
        MatrixXd A = MatrixXd::Zero(conf_->W, conf_->W);
        for (int j = 1; j <= conf_->W; j++)
            for (int i = j; i <= conf_->W; i++)
                A(j - 1, i - 1) = bji(j, i) / (1 - bji(0, i));
        return A;
    }

    MatrixXd getH() const {
        MatrixXd H = MatrixXd::Zero(conf_->W, conf_->W);
        for (int i = 0; i < conf_->W; i++)
            for (int k = 0; k < conf_->W; k++) {
                if (i == k)
                    H(i, i) = theta_(i) * (1 - theta_(i)) / phi_(i);
                else
                    H(i, k) = -theta_(i) * theta_(k) / phi_(i);
            }
        return H;
    }

    void calCRLBKnownSize() {
        MatrixXd B = getB();
        MatrixXd D = (B * theta_).asDiagonal().inverse();
        MatrixXd J = B.transpose() * D * B;
        MatrixXd I = J.inverse() - theta_ * theta_.transpose();
        crlb_ = (I.diagonal() / n_).cwiseSqrt();
        std::cout << "(known size) squared-crlb =\n" << crlb_ << "\n\n";
    }

    void calCRLBUnknownSize() {
        MatrixXd A = getA();
        MatrixXd D = (A * theta_).asDiagonal().inverse();
        MatrixXd J_phi = A.transpose() * D * A;
        MatrixXd J_phi_inv = J_phi.inverse() - phi_ * phi_.transpose();
        MatrixXd H = getH();
        MatrixXd I = H * J_phi_inv * H.transpose();
        crlb_ = (I.diagonal() / n_).cwiseSqrt();
        std::cout << "(unknown size) squared-crlb =\n" << crlb_ << "\n\n";
    }

    void calCRLB() {
        if (conf_->known_size)
            calCRLBKnownSize();
        else
            calCRLBUnknownSize();
    }

    void save(const string& output) {
        int offset = int(!conf_->known_size);
        vector<std::pair<int, double>> rst;
        for (int i = 0; i < crlb_.size(); i++)
            rst.emplace_back(i + offset, crlb_(i));
        ioutils::savePrVec(rst, output);
    }
};

#endif /* __CRLB_H__ */
