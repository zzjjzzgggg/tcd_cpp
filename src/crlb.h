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
 * Default to ITS (ITS-color)
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

    VectorXd theta_, crlb_;

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
        }
        // std::cout << theta_ << std::endl;
    }

    virtual double bji(const int j, const int i) const {
        return SpecFun::BetaBinomial(j, i, conf_->p_tri / conf_->alpha,
                                     (1 - conf_->p_tri) / conf_->alpha);
    }

    /**
     * ITS sampling matrix: sample each triangle identically
     */
    MatrixXd getB() const {
        if (conf_->known_size) {
            MatrixXd B = MatrixXd::Zero(conf_->W + 1, conf_->W + 1);
            for (int j = 0; j <= conf_->W; j++)
                for (int i = j; i <= conf_->W; i++) B(j, i) = bji(j, i);
            return B;
        } else {
            printf("unknown graph size\n");
            MatrixXd B = MatrixXd::Zero(conf_->W, conf_->W);
            for (int j = 1; j <= conf_->W; j++)
                for (int i = j; i <= conf_->W; i++)
                    B(j - 1, i - 1) = bji(j, i) / (1 - bji(0, i));
            return B;
        }
    }

    void calCRLB() {
        MatrixXd B = getB();
        MatrixXd D = (B * theta_).asDiagonal().inverse();
        MatrixXd J = B.transpose() * D * B;
        MatrixXd J_inv = J.inverse();
        MatrixXd I = J_inv - theta_ * theta_.transpose();

        crlb_ = (I.diagonal() / n_).cwiseSqrt();
        std::cout << "squared-crlb =\n" << crlb_ << "\n\n";
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
