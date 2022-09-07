#include<iostream>
#include "ap.h"
#include "stdafx.h"
#include "../src/interval_tree.h"
#include "../src/primer_bed.h"
#include "../src/clustering.h"

int _unit_test_haplotype_extraction(){
  /*
   * @return success : 0 for passing else failing 
   *
   * Looks for the correction extraction of (1) haplotypes, (2) positions, and (3) read counts
   * for an amplicon.
   */

  int success = 0;
  return(success);
}

//test the clustering function
int _unit_test_kmeans(){
  /*
   * @return success : 0 for passing else failing
   *
   * Does two checks on the kmeans function looking for (1) accurate cluster centers 
   * and (2) an accurate sil score.
   */

  int success = 0;
  //test case one
  alglib::real_2d_array xy = "[[1],[1],[1],[4],[4],[4]]";
  std::vector<float> ground_truth_centers = {1,4};
  int number_clusters = 2;

  cluster cluster_results;
  k_means(number_clusters, xy, cluster_results);
  //print_cluster_info(cluster_results);
  if(cluster_results.sil_score == 1){
    success += 1;
  }
  if(cluster_results.centers == ground_truth_centers){
    success += 1;
  }
  //std::cout << "success " << success << std::endl;
  return(success);
}

int main(){
  int success = 2; //count total number of tests
  success -= _unit_test_kmeans();
  //determine_threshold("../data/contamination_tests/test.calmd.bam", "../data/contamination_tests/sars_primers_strand.bed", "../data/contamination_tests/primer_pairs.tsv", 0);
  determine_threshold("../data/contamination_tests/simulated_alpha_beta_90_10.bam", "../data/contamination_tests/sars_primers_strand.bed", "../data/contamination_tests/primer_pairs.tsv", 0);
  return(0);
}
