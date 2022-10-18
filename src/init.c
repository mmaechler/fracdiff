#include <R.h>
#include <Rinternals.h>

#include "fracdiff.h"

#include <R_ext/Rdynload.h>

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _typ)/sizeof(name ## _typ[0]), name ##_typ}


// -- ./fdsim.c --
static R_NativePrimitiveArgType fdsim_typ[13] = {
    /*n:*/ INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    /*d__:*/     REALSXP, REALSXP, REALSXP, REALSXP,
    /*flmin__:*/ REALSXP, REALSXP, REALSXP, REALSXP
};

// -- ./fdhess.c --
static R_NativePrimitiveArgType fdhpq_typ[3] = {
    REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType fdcov_typ[11] = {
    REALSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, INTSXP, REALSXP,
    INTSXP, REALSXP, REALSXP, INTSXP
};

// -- ./fdcore.c --
static R_NativePrimitiveArgType fracdf_typ[19] = {
    /* x */   REALSXP, INTSXP, INTSXP, INTSXP, INTSXP,
    /* dtol */REALSXP, REALSXP, REALSXP,
    /* d__*/  REALSXP, REALSXP, REALSXP, REALSXP,
    /* lenw */INTSXP, INTSXP, INTSXP,
    /* flmin*/REALSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType fdcom_typ[9] = {
    /* n */    INTSXP, INTSXP, INTSXP, INTSXP,
    /* hood */ REALSXP, REALSXP, REALSXP,
    /*epmin */ REALSXP, REALSXP
};

static const R_CMethodDef CEntries[]  = {
    CDEF(fdsim),
    CDEF(fdhpq),
    CDEF(fdcov),
    CDEF(fracdf),
    CDEF(fdcom),
    {NULL, NULL, 0}
};

/* static R_CallMethodDef CallEntries[] = {
 *     {NULL, NULL, 0}
 * };
 */

/* static R_FortranMethodDef FortEntries[] = {
 *     {NULL, NULL, 0}
 * };
 */

void R_init_fracdiff(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL/*CallEntries*/, NULL/*FortEntries*/, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
