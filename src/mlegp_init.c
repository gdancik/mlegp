#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void fitGPfromR(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
        {"fitGPfromR", (DL_FUNC) &fitGPfromR, 21},
            {NULL, NULL, 0}
};

void R_init_mlegp(DllInfo *dll)
{
        R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
            R_useDynamicSymbols(dll, FALSE);
}


