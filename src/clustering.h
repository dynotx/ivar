#include<iostream>
#include "stdafx.h"
#include "dataanalysis.h"
void k_means(int n_clusters);
float average(std::vector<float> x);
std::vector<float> calculate_cluster_centers(alglib::real_2d_array X, alglib::kmeansreport rep, int n_clusters);
float cluster_point_distances(alglib::real_2d_array X, alglib::kmeansreport rep, float point, float center, int n_clusters);
float calculate_sil_score(alglib::real_2d_array X, alglib::kmeansreport rep,int n_clusters);
