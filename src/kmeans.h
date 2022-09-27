#include "dataanalysis.h"
#include "ap.h"
#include "stdafx.h"

#ifndef kmeans
#define kmeans

void custom_kmeans(alglib::clusterizerstate &s, alglib::ae_int_t k, alglib::kmeansreport &rep, alglib::xparams _xparams = alglib::xdefault);
 void custom_kmeans(alglib::clusterizerstate *s, alglib::ae_int_t k, alglib::kmeansreport *rep, alglib_impl::ae_state *_state);
#endif
