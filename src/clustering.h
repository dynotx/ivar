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

class cluster {
  public:
    int n_clusters;
    float sil_score; //average of per sample sil score
    std::vector<float> sil_scores; //per sample                                       
    std::vector<float> centers;
    std::vector<std::vector<float>> sorted_points; //all points sorted by center they belong to
    std::vector<std::vector<float>> cluster_bounds; //the upper and lower bounds of each cluster
};

//zipping and unzipping functions for haplotypes && position
void zip(const std::vector<int> &haplotypes, const std::vector<uint32_t> &positions,std::vector<std::pair<uint32_t,int>> &zipped);
void unzip(const std::vector<std::pair<uint32_t, int>> &zipped, std::vector<int> &haplotypes, std::vector<uint32_t> &positions);
void print_cluster_info(cluster cluster_results);
void k_means(int n_clusters, alglib::real_2d_array xy, cluster &cluster_results);
float average(std::vector<float> x);
std::vector<float> calculate_cluster_centers(alglib::real_2d_array X, alglib::kmeansreport rep, int n_clusters);
float cluster_point_distances(alglib::real_2d_array X, alglib::kmeansreport rep, float point, float center, int n_clusters);
void calculate_sil_score(alglib::real_2d_array X, alglib::kmeansreport rep,int n_clusters, cluster &cluster_results);
int determine_threshold(std::string bam, std::string bed, std::string pair_info, int32_t primer_offset);
std::string decoded_nucs(int tmp);
std::vector<std::vector<int>> transpose(const std::vector<std::vector<int>> data);
#endif
