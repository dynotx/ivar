#include<iostream>
#include "ap.h"
#include "stdafx.h"
#include "../src/interval_tree.h"
#include "../src/primer_bed.h"
#include "../src/clustering.h"
#include "../src/kmeans.h"

void reset_positions(std::vector<position> all_positions){

}

int _unit_test_haplotype_extraction(){
  /*
   * @return success : 0 for passing else failing 
   *
   * Looks for the correction extraction of (1) haplotypes, (2) positions, and (3) read counts
   * for an amplicon.
   */
  int success = 0;
  
  int32_t primer_offset = 0;
  std::string pair_info = "../data/contamination_tests/primer_pairs.tsv";
  std::string bed = "../data/contamination_tests/sars_primers_strand.bed";

  std::vector<double> return_frequencies;
  IntervalTree amplicons;
  std::vector<position> all_positions;
  std::vector<primer> primers; 
  //populate amplicons for sars-cov-2
  primers = populate_from_file(bed, primer_offset);
  amplicons = populate_amplicons(pair_info, primers);

  //populate all positions
  

  /* Generate some fake haplotypes / positions to test.
   * Using the amplicon between 23601-23619 23728-23751.
   *
   * Case A. Substituion in the primer region.
   *
   *
   *
   */
  

  return(success);
}

int _unit_test_frequencies(){
  /*
   * Unit test for the function create_frequency_matrix.
   */
  int success = 0;

  int32_t primer_offset = 0;
  std::string pair_info = "../data/contamination_tests/primer_pairs.tsv";
  std::string bed = "../data/contamination_tests/sars_primers_strand.bed";

  std::vector<double> return_frequencies;
  std::vector<primer> primers;
  IntervalTree amplicons;
  std::vector<position> all_positions;
  
  //populate amplicons for sars-cov-2
  primers = populate_from_file(bed, primer_offset);
  amplicons = populate_amplicons(pair_info, primers);

  //set up basic nucleotides
  std::vector<std::string> basic_nts = {"A", "C", "G", "T"};
  std::vector<allele> basic_alleles;
  for(std::string nt : basic_nts){
    allele new_allele;
    new_allele.nuc = nt;
    new_allele.depth = 0;
    basic_alleles.push_back(new_allele);
  }

  //populate all positions
  for(uint32_t i=0; i < 29903; i++){
    position new_position;
    new_position.pos = i;
    new_position.depth = 0;
    new_position.ad = basic_alleles;
    all_positions.push_back(new_position);
  }

  /* Generate some fake haplotypes / positions to test.
   * Using the amplicon between 23601-23619 23728-23751.
   *
   * In this fictional genome 23700 reference is T
   *
   * Case A. Drop position 23700 w/ mutation frequency < 0.03, return empty frequency vector
   * Case B. Drop position 23700 w/ mutation depth < 10, return empty frequency vector
   * Case C. Two mutations at 23700 in proportions 0.30, 0.20, match ref 0.50
   */

  //declare some baseline variables
  std::vector<uint32_t> positions;
  std::vector<int> haplotypes;
  std::vector<double> frequencies;
  bool reverse = false;
  uint32_t abs_start_pos = 23620;
  uint32_t abs_end_pos = 23727;
  std::vector<std::string> nucleotides;
  std::vector<uint32_t> range = {23620, 23727};

  //Case A 
  positions = {23700};
  haplotypes = {-2};
  nucleotides = {"T"};
  int i = 0;
  //add haplotypes that match the reference 1000 times
  while(i < 1000){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    update_allele_depth(all_positions, nucleotides, positions);
    i++;
  }
  nucleotides.clear();
  haplotypes.clear();
  haplotypes = {0};
  nucleotides = {"A"};
  i = 0;
  //add a mutation 10 time
  while(i < 10){
    amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range, all_positions);
    update_allele_depth(all_positions, nucleotides, positions);
    i++;
  }
  //try to create the frequency matrix
  frequencies = create_frequency_matrix(amplicons, all_positions);
  if(frequencies.size() == 0){
    success += 1;
  }

  //Case B
  nucleotides.clear();
  haplotypes.clear();

  return(0);
}

int _unit_test_sil_score(){
  /*
   * Makes sure the scoring metric is working to par.
   */
  return(0);  
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
  //if(cluster_results.centers == ground_truth_centers){
  //  success += 1;
  //}
  //std::cout << "success " << success << std::endl;
  return(success);
}

int main(){
  int success = 2; //count total number of tests
  _unit_test_frequencies();
  std::cout << "success: " << success << std::endl;
  //success -= _unit_test_kmeans();
  //success -= _unit_test_sil_score();
  //determine_threshold("../data/contamination_tests/test.calmd.bam", "../data/contamination_tests/sars_primers_strand.bed", "../data/contamination_tests/primer_pairs.tsv", 0);
  //determine_threshold("../data/contamination_tests/simulated_alpha_beta_90_10.bam", "../data/contamination_tests/sars_primers_strand.bed", "../data/contamination_tests/primer_pairs.tsv", 0);
  return(0);
}
