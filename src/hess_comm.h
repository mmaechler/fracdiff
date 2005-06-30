/* included only by ./fdcore.c and ./fdhess.c : */

FD_EXTERNAL
struct { int  n,  m, p, q, pq, pq1, maxpq, maxpq1, minpq, nm; } Dims;

FD_EXTERNAL
struct { double hatmu, wnv, cllf; } filtfd_;
FD_EXTERNAL
struct { int ksvd, kcov, kcor; } hessfd_;

FD_EXTERNAL
struct { int ly, lamk, lak, lvk, lphi, lpi; } w_fil;

FD_EXTERNAL
struct { int lqp, la, lajac, ipvt, ldiag, lqtf, lwa1, lwa2, lwa3, lwa4; } w_opt;




