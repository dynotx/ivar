#include "dataanalysis.h"
#include "ap.h"
#include "stdafx.h"
using namespace alglib;
/*
 * Written to extend alglib kmeans clustering, imposing 
 * and additional "compositional constraint" and fixing a noise
 * cluster at 0.03.
 */

void custom_kmeans(alglib::clusterizerstate *s,
     alglib::ae_int_t k,
     alglib::kmeansreport *rep,
     alglib_impl::ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix dummy;

    ae_frame_make(_state, &_frame_block);
    memset(&dummy, 0, sizeof(dummy));
    _kmeansreport_clear(rep);
    ae_matrix_init(&dummy, 0, 0, DT_REAL, _state, ae_true);

    ae_assert(k>=0, "ClusterizerRunKMeans: K<0", _state);

    /*
     * Incorrect distance type
     */
    if( s->disttype!=2 )
    {
        rep->npoints = s->npoints;
        rep->terminationtype = -5;
        rep->k = k;
        rep->iterationscount = 0;
        rep->energy = 0.0;
        ae_frame_leave(_state);
        return;
    }

    /*
     * No points
     */
    if( s->npoints==0 )
    {
        rep->npoints = 0;
        rep->terminationtype = 1;
        rep->k = k;
        rep->iterationscount = 0;
        rep->energy = 0.0;
        ae_frame_leave(_state);
        return;
    }

    /*
     * Normal case:
     * 1<=K<=NPoints, Euclidean distance
     */
    rep->npoints = s->npoints;
    rep->nfeatures = s->nfeatures;
    rep->k = k;
    rep->npoints = s->npoints;
    rep->nfeatures = s->nfeatures;
    kmeansgenerateinternal(&s->xy, s->npoints, s->nfeatures, k, s->kmeansinitalgo, s->seed, s->kmeansmaxits, s->kmeansrestarts, s->kmeansdbgnoits, &rep->terminationtype, &rep->iterationscount, &dummy, ae_false, &rep->c, ae_true, &rep->cidx, &rep->energy, &s->kmeanstmp, _state);
    ae_frame_leave(_state);
}

/*************************************************************************
K-means++ clusterization

INPUT PARAMETERS:
    XY          -   dataset, array [0..NPoints-1,0..NVars-1].
    NPoints     -   dataset size, NPoints>=K
    NVars       -   number of variables, NVars>=1
    K           -   desired number of clusters, K>=1
    InitAlgo    -   initialization algorithm:
                    * 0 - automatic selection of best algorithm
                    * 1 - random selection of centers
                    * 2 - k-means++
                    * 3 - fast-greedy init
                    *-1 - first K rows of dataset are used
                          (special debug algorithm)
    Seed        -   seed value for internal RNG:
                    * positive value is used to initialize RNG in order to
                      induce deterministic behavior of algorithm
                    * zero or negative value means  that random  seed   is
                      generated
    MaxIts      -   iterations limit or zero for no limit
    Restarts    -   number of restarts, Restarts>=1
    KMeansDbgNoIts- debug flag; if set, Lloyd's iteration is not performed,
                    only initialization phase.
    Buf         -   special reusable structure which stores previously allocated
                    memory, intended to avoid memory fragmentation when solving
                    multiple subsequent problems:
                    * MUST BE INITIALIZED WITH KMeansInitBuffers() CALL BEFORE
                      FIRST PASS TO THIS FUNCTION!
                    * subsequent passes must be made without re-initialization

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -3, if task is degenerate (number of distinct points is
                          less than K)
                    * -1, if incorrect NPoints/NFeatures/K/Restarts was passed
                    *  1, if subroutine finished successfully
    IterationsCount- actual number of iterations performed by clusterizer
    CCol        -   array[0..NVars-1,0..K-1].matrix whose columns store
                    cluster's centers
    NeedCCol    -   True in case caller requires to store result in CCol
    CRow        -   array[0..K-1,0..NVars-1], same as CCol, but centers are
                    stored in rows
    NeedCRow    -   True in case caller requires to store result in CCol
    XYC         -   array[NPoints], which contains cluster indexes
    Energy      -   merit function of clusterization

  -- ALGLIB --
     Copyright 21.03.2009 by Bochkanov Sergey
*************************************************************************/
void kmeansgenerateinternal(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t k,
     ae_int_t initalgo,
     ae_int_t seed,
     ae_int_t maxits,
     ae_int_t restarts,
     ae_bool kmeansdbgnoits,
     ae_int_t* info,
     ae_int_t* iterationscount,
     /* Real    */ ae_matrix* ccol,
     ae_bool needccol,
     /* Real    */ ae_matrix* crow,
     ae_bool needcrow,
     /* Integer */ ae_vector* xyc,
     double* energy,
     kmeansbuffers* buf,
     ae_state *_state)
  {
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t i1;
    double e;
    double eprev;
    double v;
    double vv;
    ae_bool waschanges;
    ae_bool zerosizeclusters;
    ae_int_t pass;
    ae_int_t itcnt;
    hqrndstate rs;

    ae_frame_make(_state, &_frame_block);
    memset(&rs, 0, sizeof(rs));
    *info = 0;
    *iterationscount = 0;
    ae_matrix_clear(ccol);
    ae_matrix_clear(crow);
    ae_vector_clear(xyc);
    *energy = 0;
    _hqrndstate_init(&rs, _state, ae_true);


    /*
     * Test parameters
     */
    if( ((npoints<k||nvars<1)||k<1)||restarts<1 )
    {
        *info = -1;
        *iterationscount = 0;
        ae_frame_leave(_state);
        return;
    }

    /*
     * TODO: special case K=1
     * TODO: special case K=NPoints
     */
    *info = 1;
    *iterationscount = 0;

    /*
     * Multiple passes of k-means++ algorithm
     */
    if( seed<=0 )
    {
        hqrndrandomize(&rs, _state);
    }
    else
    {
        hqrndseed(325355, seed, &rs, _state);
    }
    ae_vector_set_length(xyc, npoints, _state);
    rmatrixsetlengthatleast(&buf->ct, k, nvars, _state);
    rmatrixsetlengthatleast(&buf->ctbest, k, nvars, _state);
    ivectorsetlengthatleast(&buf->xycprev, npoints, _state);
    ivectorsetlengthatleast(&buf->xycbest, npoints, _state);
    rvectorsetlengthatleast(&buf->d2, npoints, _state);
    ivectorsetlengthatleast(&buf->csizes, k, _state);
    *energy = ae_maxrealnumber;
    for(pass=1; pass<=restarts; pass++)
    {

        /*
         * Select initial centers.
         *
         * Note that for performance reasons centers are stored in ROWS of CT, not
         * in columns. We'll transpose CT in the end and store it in the C.
         *
         * Also note that SelectInitialCenters() may return degenerate set of centers
         * (some of them have no corresponding points in dataset, some are non-distinct).
         * Algorithm below is robust enough to deal with such set.
         */
        clustering_selectinitialcenters(xy, npoints, nvars, initalgo, &rs, k, &buf->ct, &buf->initbuf, &buf->updatepool, _state);

        /*
         * Lloyd's iteration
         */
                if( !kmeansdbgnoits )
        {

            /*
             * Perform iteration as usual, in normal mode
             */
            for(i=0; i<=npoints-1; i++)
            {
                xyc->ptr.p_int[i] = -1;
            }
            eprev = ae_maxrealnumber;
            e = ae_maxrealnumber;
            itcnt = 0;
            while(maxits==0||itcnt<maxits)
            {

                /*
                 * Update iteration counter
                 */
                itcnt = itcnt+1;
                inc(iterationscount, _state);

                /*
                 * Call KMeansUpdateDistances(), fill XYC with center numbers,
                 * D2 with center distances.
                 */
                for(i=0; i<=npoints-1; i++)
                {
                    buf->xycprev.ptr.p_int[i] = xyc->ptr.p_int[i];
                }
                kmeansupdatedistances(xy, 0, npoints, nvars, &buf->ct, 0, k, xyc, &buf->d2, &buf->updatepool, _state);
                waschanges = ae_false;
                for(i=0; i<=npoints-1; i++)
                {
                    waschanges = waschanges||xyc->ptr.p_int[i]!=buf->xycprev.ptr.p_int[i];
                }

                /*
                 * Update centers
                 */
                for(j=0; j<=k-1; j++)
                {
                    buf->csizes.ptr.p_int[j] = 0;
                }
                for(i=0; i<=k-1; i++)
                {
                    for(j=0; j<=nvars-1; j++)
                    {
                        buf->ct.ptr.pp_double[i][j] = (double)(0);
                    }
                }
                for(i=0; i<=npoints-1; i++)
                {
                    buf->csizes.ptr.p_int[xyc->ptr.p_int[i]] = buf->csizes.ptr.p_int[xyc->ptr.p_int[i]]+1;
                    ae_v_add(&buf->ct.ptr.pp_double[xyc->ptr.p_int[i]][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
                }
                zerosizeclusters = ae_false;
                for(j=0; j<=k-1; j++)
                {
                    if( buf->csizes.ptr.p_int[j]!=0 )
                    {
                        v = (double)1/(double)buf->csizes.ptr.p_int[j];
                        ae_v_muld(&buf->ct.ptr.pp_double[j][0], 1, ae_v_len(0,nvars-1), v);
                    }
                    zerosizeclusters = zerosizeclusters||buf->csizes.ptr.p_int[j]==0;
                }
                if( zerosizeclusters )
                {

                    /*
                     * Some clusters have zero size - rare, but possible.
                     * We'll choose new centers for such clusters using k-means++ rule
                     * and restart algorithm, decrementing iteration counter
                     * in order to allow one more iteration (this one was useless
                     * and should not be counted).
                     */
                    if( !clustering_fixcenters(xy, npoints, nvars, &buf->ct, k, &buf->initbuf, &buf->updatepool, _state) )
                    {
                        *info = -3;
                        ae_frame_leave(_state);
                        return;
                    }
                    itcnt = itcnt-1;
                    continue;
                }
                                /*
                 * Stop if one of two conditions is met:
                 * 1. nothing has changed during iteration
                 * 2. energy function increased after recalculation on new centers
                 */
                e = (double)(0);
                for(i=0; i<=npoints-1; i++)
                {
                    v = 0.0;
                    i1 = xyc->ptr.p_int[i];
                    for(j=0; j<=nvars-1; j++)
                    {
                        vv = xy->ptr.pp_double[i][j]-buf->ct.ptr.pp_double[i1][j];
                        v = v+vv*vv;
                    }
                    e = e+v;
                }
                if( !waschanges||ae_fp_greater_eq(e,eprev) )
                {
                    break;
                }

                /*
                 * Update EPrev
                 */
                eprev = e;
            }
        }
        else
        {

            /*
             * Debug mode: no Lloyd's iteration.
             * We just calculate potential E.
             */
            kmeansupdatedistances(xy, 0, npoints, nvars, &buf->ct, 0, k, xyc, &buf->d2, &buf->updatepool, _state);
            e = (double)(0);
            for(i=0; i<=npoints-1; i++)
            {
                e = e+buf->d2.ptr.p_double[i];
            }
        }

        /*
         * Compare E with best centers found so far
         */
        if( ae_fp_less(e,*energy) )
        {

            /*
             * store partition.
             */
            *energy = e;
            copymatrix(&buf->ct, 0, k-1, 0, nvars-1, &buf->ctbest, 0, k-1, 0, nvars-1, _state);
            for(i=0; i<=npoints-1; i++)
            {
                buf->xycbest.ptr.p_int[i] = xyc->ptr.p_int[i];
            }
        }
    }

    /*
     * Copy and transpose
     */
    if( needccol )
    {
        ae_matrix_set_length(ccol, nvars, k, _state);
        copyandtranspose(&buf->ctbest, 0, k-1, 0, nvars-1, ccol, 0, nvars-1, 0, k-1, _state);
    }
    if( needcrow )
    {
        ae_matrix_set_length(crow, k, nvars, _state);
        rmatrixcopy(k, nvars, &buf->ctbest, 0, 0, crow, 0, 0, _state);
    }
    for(i=0; i<=npoints-1; i++)
    {
        xyc->ptr.p_int[i] = buf->xycbest.ptr.p_int[i];
    }
    ae_frame_leave(_state);
}
