#ifndef _LBFGSB_HPP
#define _LBFGSB_HPP

int setulb_(int n, int m, double *x, double *l, double *u,
        int *nbd, double *f, double *g, double *factr,
        double *pgtol, double *wa, int *iwa,
        char *task, char *csave, int *lsave,
        int *isave, double *dsave);

#endif

