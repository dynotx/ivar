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
#include "stdafx.h"
using namespace alglib;

float average(std::vector<float> x){
  float sumTotal = 0;
  for(float k=0; k < x.size(); ++k){
      sumTotal += x[k];            
  }
  return(sumTotal / x.size());
}


//calculate the cluster centers
std::vector<float> calculate_cluster_centers(alglib::real_2d_array X, alglib::kmeansreport rep, int n_clusters){
  /*
   * @params X : array containing all data points
   * @params rep : kmeans reporter object  
   * @returns centers : vector containing all center values 
   */
  std::vector<float> centers;
  float summation [n_clusters][2]; //track all points per cluster
  float point=0;
  int label=0;

  for (int i = 0; i < X.rows(); i++) {
    point = X[i][0];
    label = rep.cidx[i];
    summation[label][0] += point;
    summation[label][1] += 1;
  }
  
  for (int k = 0; k < n_clusters; k++) {
    centers.push_back(summation[k][0]/summation[k][1]); 
    //std::cout << "center " << k << " is " << summation[k][0]/summation[k][1] << std::endl;
  }
  return(centers);
}

//support function for sil score
float cluster_point_distances(alglib::real_2d_array X, alglib::kmeansreport rep, float point, float center, int n_clusters){
  /*
  * @params X : array containing all data points
  * @params rep : kemans reporter object
  * @params point : the data point of interest
  * @params center : the center of the data point of interests cluster
  * @returns ab : tuple with points a and b representing 
  */
  std::vector<float> internal_dist; //keep track of distance between point of interest and all points in cluster
  float external_dist [n_clusters][2]; //dim 1 cluster value, dim 2 distance sum points and point count
  float compare_point=0;
  int compare_label=0;
  float a=0; //sil score a term
  float b=1; //sil score b term
  float tmp=0;
  //iterate through and mind distances between points in cluster and out of clusters
  for (int i = 0; i < X.rows(); i++) {
    compare_point = X[i][0];
    compare_label = rep.cidx[i];
    //belongs to the same cluster
    if(compare_label == center){
      internal_dist.push_back(abs(point-compare_point));
    }else{ //doesn't belong to the same cluster
      external_dist[compare_label][0] += abs(point-compare_point); //sum distance between poi and cluster points
      external_dist[compare_label][1] += 1; //count number points in cluster    
    }
  } 
  a = average(internal_dist);

  //find the average distance from each external cluster
  for(int z=0; z < n_clusters; z++){
    tmp = external_dist[z][0] / external_dist[z][1];
    //find the minimum of the external cluster dists
    if(tmp < b){
      b = tmp;
    }
  }
  return((b - a) / std::max(a,b));
}

//function calculates the sil score of all points
float calculate_sil_score(alglib::real_2d_array X, alglib::kmeansreport rep, 
    int n_clusters){
  /*
  * @params X : array containing data points
  * @params rep : kmeans reported object
  * @returns sil_scores : all sil scores for all points
  */
  //std::cout << "Calculating silohuette score." << std::endl;
  int rows = X.rows();
  //int cols = X.cols();
  float point = 0;
  float center = 0;
  std::vector<float> sil_scores;
  float tmp = 0;
  for (int i = 0; i < rows; i++) {
    point = X[i][0];
    center = rep.cidx[i];
    tmp = cluster_point_distances(X, rep, point, center, n_clusters);
    sil_scores.push_back(tmp);
    //(b - a) / max(a, b) sil score formula 
    //std::cout << "i " << i << " point " << point << " center " << center << " sil " << tmp << std::endl;     
  }
  
  return(0.0);
}
//does the actual k means ++ clustering in a loop
void k_means(int n_clusters){
  /*
  * @params n_clusters : then number of clusters
  * @params amplicon : amplicon object
  */
  alglib::clusterizerstate s;
  alglib::kmeansreport rep;
  alglib::real_2d_array xy = "[[1],[1],[4],[2],[4],[5]]";
  int num_points = 6;
  clusterizercreate(s);
  //X is the data
  //the number of pointss
  //last pos is distance matrix type
  //num_features is next
  //then distance matrix
  clusterizersetpoints(s, xy, num_points, 1, 2);
  clusterizersetkmeanslimits(s, 5, 0);
  //this is the cluster size!!!
  clusterizerrunkmeans(s, n_clusters, rep);
  if (int(rep.terminationtype) != 1){
    std::cout << "Error" << std::endl;
    exit(1);
  }

  int i = 0;
  while(i < n_clusters){
    std::cout << "center "<< i << " "  << rep.c[i][0] << " " << rep.c[i][1] << std::endl;
    i++;
  }
 
  calculate_sil_score(xy, rep, n_clusters);
  calculate_cluster_centers(xy, rep, n_clusters);
}

void iterate_reads(bam1_t *r){
  //get cigar for this read
  //uint32_t *cigar = bam_get_cigar(r);
  uint32_t i = 0;  
  //iterate through each cigar operation, using only what doesn't match the ref or is SC
  //ie. we use insertions, deletions, substitutions
  while(i < r->core.n_cigar){
    break;
  }
  std::cout << "made it this far" << std::endl;
}

void determine_threshold(std::string bam){
  std::cout << "Begin threshold optimiziation" << std::endl;
  samFile *in = hts_open(bam.c_str(), "r");
  hts_idx_t *idx = sam_index_load(in, bam.c_str());
  bam_hdr_t *header = sam_hdr_read(in);
  bam1_t *aln = bam_init1();
  hts_itr_t *iter = NULL;
  std::string region_;
  //region refers to reference
  region_.assign(header->target_name[0]);
  iter  = sam_itr_querys(idx, header, region_.c_str());

  //initialize haplotype data structure
  
  //this iterates over the reads and assigns them to an amplicon
  while(sam_itr_next(in, iter, aln) >= 0) {
    iterate_reads(aln);  
  }

  //cry 
 
}
