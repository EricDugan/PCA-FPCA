//
//  code.c
//  CLoop
//
//  Last updated by Michelle Carey 16 July 2015
//
//

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* The gateway function */
SEXP loopJuan(SEXP NLHS, SEXP PLHS,
          SEXP NRHS, SEXP PRHS, SEXP WVEC)
{
    int ncols1 = ncols(NLHS);
    int ncols2 = ncols(PLHS);
    int ncols3 = ncols(NRHS);
    int ncols4 = ncols(PRHS);
    
    int nrows1 = nrows(NLHS);
    /*PROTECT(PTR =allocVector(NLHSNEW,ncols1*nrows1));
     PROTECT(PTR =allocVector(PLHSNEW,ncols2*nrows1));
     PROTECT(PTR =allocVector(NRHSNEW,ncols3*nrows1));
     PROTECT(PTR =allocVector(PRHSNEW,ncols4*nrows1));
     PROTECT(PTR =allocVector(WVECNEW,nrows1));*/
    
    SEXP NLHSNEW;
    SEXP PLHSNEW;
    SEXP NRHSNEW;
    SEXP PRHSNEW;
    SEXP WVECNEW;
    
    PROTECT(NLHSNEW =coerceVector(NLHS,REALSXP));
    PROTECT(PLHSNEW =coerceVector(PLHS,REALSXP));
    PROTECT(NRHSNEW =coerceVector(NRHS,REALSXP));
    PROTECT(PRHSNEW =coerceVector(PRHS,REALSXP));
    PROTECT(WVECNEW =coerceVector(WVEC,REALSXP));
    
    /* variable declarations here */
    SEXP PTR;
    PROTECT(PTR =allocVector(REALSXP,ncols1*ncols2*ncols3*ncols4));
    
    for (int p = 0; p < (ncols1*ncols2*ncols3*ncols4); p++)
    {
        REAL(PTR)[p]=0;
    }
    
    for (int i = 0; i < ncols1; i++)
    {
        int mn1 = i*nrows1;
        
        for (int j = 0; j < ncols2; j++)
        {
            int mn2 = j*nrows1;
            
            for (int k = 0; k < ncols3; k++)
            {
                int mn3 = k*nrows1;
                
                for (int l = 0; l < ncols4; l++)
                {
                    int mn4 = l*nrows1;
                    
                    for (int m = 0; m < nrows1; m++)
                    {
                        double wv   = REAL(WVECNEW)[m];
                        double arri = REAL(NLHSNEW)[m + mn1];
                        double arrj = REAL(PLHSNEW)[m + mn2];
                        double arrk = REAL(NRHSNEW)[m + mn3];
                        double arrl = REAL(PRHSNEW)[m + mn4];
                        
                        int index = i*ncols4*ncols3*ncols2
                        + j*ncols4*ncols3
                        + k*ncols4
                        + l;
                        
                        REAL(PTR)[index] += arri*arrj*arrk*arrl*wv;
                        
                    }
                }
                
            }
        }
    }
    
    UNPROTECT(6);
    return PTR;
    
}
