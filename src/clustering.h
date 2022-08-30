#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <sstream>
#include <vector>
#include "ap.h"
#include "dataanalysis.h"
#include "interval_tree.h"
#include "stdafx.h"
#include "primer_bed.h"

#ifndef clustering
#define clustering

void k_means(int n_clusters);
float average(std::vector<float> x);
std::vector<float> calculate_cluster_centers(alglib::real_2d_array X, alglib::kmeansreport rep, int n_clusters);
float cluster_point_distances(alglib::real_2d_array X, alglib::kmeansreport rep, float point, float center, int n_clusters);
float calculate_sil_score(alglib::real_2d_array X, alglib::kmeansreport rep,int n_clusters);
void determine_threshold(std::string bam, std::string bed, std::string pair_info, int32_t primer_offset);
#endif
