/**
 * Copyright (C) by J.Z. (06/30/2017 10:47)
 * Distributed under terms of the MIT license.
 */

#include "stdafx.h"

int main(int argc, char* argv[]) {
    printf("Built: %s %s\n", __DATE__, __TIME__);

    for (int i = 1; i < 20; i++)
        printf("i: %d, %d\n", i, int(std::floor(std::log2(i))));
    return 0;
}
