/**
 * Copyright (C) by J.Z. (07/02/2017 09:55)
 * Distributed under terms of the MIT license.
 */

#ifndef __CRLB_H__
#define __CRLB_H__

#include "stdafx.h"
#include <Eigen/Dense>

using namespace Eigen;

class CRLB {
public:
    struct Config {
        int W;  // maximum triadic cardinality
        int n;  // number of nodes/sampled nodes
        double alpha;
        double p_tri;
        double p_nd;       // node sampling rate
        string theta_fnm;  // [theta_0, ..., theta_W]

        void echo() {
            printf(
                "W: %d\n"
                "n: %d\n"
                "alpha: %g\n"
                "p_tri: %g\n"
                "theta: %s\n\n",
                W, n, alpha, p_tri, theta_fnm.c_str());
        }
    };

private:
    const Config* conf_;

    VectorXd theta_, crlb_;

private:
    inline double bji(const int j, const int i) const {
        return SpecFun::BetaBinomial(j, i, conf_->p_tri / conf_->alpha,
                                     (1 - conf_->p_tri) / conf_->alpha);
    }

public:
    CRLB(const Config* conf) : conf_(conf) {
        theta_ = VectorXd::Zero(conf_->W + 1);
        ioutils::TSVParser ss(conf_->theta_fnm);
        while (ss.next()) theta_(ss.get<int>(0)) = ss.get<double>(1);
    }

    /**
     * IS sampling matrix: sample each triangle identically
     */
    MatrixXd getISB() const {
        printf("IS\n");
        MatrixXd B = MatrixXd::Zero(conf_->W + 1, conf_->W + 1);
        for (int j = 0; j <= conf_->W; j++)
            for (int i = j; i <= conf_->W; i++) B(j, i) = bji(j, i);
        return B;
    }

    /**
     * US sampling matrix: first sample nodes, then sample triangle
     */
    MatrixXd getUSB() const {
        printf("US\n");
        MatrixXd B = MatrixXd::Zero(conf_->W + 1, conf_->W + 1);
        // for j=0
        for (int i = 0; i <= conf_->W; i++)
            B(0, i) = 1 - conf_->p_nd + conf_->p_nd * bji(0, i);
        // for j>0
        for (int j = 1; j <= conf_->W; j++)
            for (int i = j; i <= conf_->W; i++)
                B(j, i) = conf_->p_nd * bji(j, i);
        return B;
    }

    void calCRLB(const MatrixXd& B);

    void save(const string& output);
};

#endif /* __CRLB_H__ */
