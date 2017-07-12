/**
 * Copyright (C) by J.Z. (06/30/2017 10:47)
 * Distributed under terms of the MIT license.
 */

#include "stdafx.h"

#include <Eigen/Dense>
using namespace Eigen;

int main(int argc, char *argv[]) {
    double p = 0.1, sum = 0, alpha = 0.3;
    int i = 20;
    auto BB = [&](int j, int i) {
        return SpecFun::BetaBinomial(j, i, p / alpha, (1 - p) / alpha);
    };

    sum = 0;
    for (int j = 0; j < i; j++) {
        sum += BB(1, 1) * BB(0, j);
    }
    sum += BB(0, i);
    printf("\nalpha = %.3f, sum= %.3f\n", alpha, sum);
    return 0;
}
