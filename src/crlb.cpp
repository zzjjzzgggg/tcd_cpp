/**
 * Copyright (C) by J.Z. (07/02/2017 13:08)
 * Distributed under terms of the MIT license.
 */

#include <iostream>
#include "crlb.h"

void CRLB::calCRLB(const MatrixXd& B) {
    MatrixXd D = (B * theta_).asDiagonal().inverse();
    MatrixXd J = B.transpose() * D * B;
    MatrixXd J_inv = J.inverse();
    MatrixXd I = J_inv - theta_ * theta_.transpose();

    crlb_ = (I.diagonal() / conf_->n).cwiseSqrt();
    std::cout << "crlb =\n" << crlb_ << "\n\n";
}

void CRLB::save(const string& output) {
    vector<std::pair<int, double>> rst;
    for (int i = 0; i <= conf_->W; i++) rst.emplace_back(i, crlb_(i));
    ioutils::savePrVec(rst, output);
}
