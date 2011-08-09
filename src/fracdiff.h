
// fdsim.c --------------------------------------
void fdsim(int *n, int *ip, int *iq, double *ar, double *ma,
	   double *d__, double *mu, double *y, double *s,
	   double *flmin, double *flmax, double *epmin, double *epmax);

// fdcore.c --------------------------------------

void fracdf(double *x, int *n, int *m, int *nar, int *nma,
	    double *dtol, double *drange, double *hood_etc,
	    double *d__, double *ar, double *ma, double *w,
	    int *lenw, int *iw, int *inform, // <- also use as input
	    double *flmin, double *flmax, double *epmin, double *epmax);

void fdfilt(double *x, double d,
	    /* output : */
	    double *y, double *slogvk,
	    /* using */
	    double *amk, double *ak, double *vk,
	    double *phi, double *pi);

void fdcom(int *n, int *m, int *nar, int *nma,
	   double *hood, double *flmin, double *flmax,
	   double *epmin, double *epmax);

void ajqp_(double *qp, double *a, double *ajac,
	   int lajac, int op_code, double *y);

// fdhess.c --------------------------------------

void fdhpq(double *h, int *lh, double *w);

void fdcov(double *x, double *d__, double *hh,
	   double *hd, double *cov, int *lcov, double *cor,
	   int *lcor, double *se, double *w, int *info);

// fdmin.c --------------------------------------

typedef /* Subroutine */ void (*S_fp)();

double lmder1(S_fp fcn, int m, int n,
	      double *x, double *fvec, double *fjac, int ldfjac,
	      double ftol, double xtol, double gtol, int maxfev, double *diag,
	      int mode, double factor,
	      int *info, int *nfev, int *njev,
	      int *ipvt, double *qtf,
	      double *wa1, double *wa2, double *wa3, double *wa4, double *y);

