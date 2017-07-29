/**
 * Copyright (C) by J.Z. (07/28/2017 09:44)
 * Distributed under terms of the MIT license.
 */

#ifndef __CRLB_SGS_H__
#define __CRLB_SGS_H__

#include "crlb.h"

class CRLB_SGS : public CRLB {
public:
    CRLB_SGS(const Config* conf) : CRLB(conf) {}

    double bji(const int j, const int i) const override {
        if (i == 0 && j == 0)
            return 1;
        else if (i == j)
            return conf_->p_nd;
        else if (j == 0)
            return 1 - conf_->p_nd;
        return 0;
    }
};

#endif /* __CRLB_SGS_H__ */
