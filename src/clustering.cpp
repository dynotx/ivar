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
#include "primer_bed.h"
#include "stdafx.h"
#include "clustering.h"
using namespace alglib;

//get the average of a vector
float average(std::vector<float> x){
  float sumTotal = 0;
  for(float k=0; k < x.size(); ++k){
      sumTotal += x[k];            
  }
  return(sumTotal / x.size());
}

int encoded_nucs(char &tmp){
  /*
   * @param tmp : nucleotide to encode
   * @return encoded_nuc : nucleotide encoded as an int
   *
   * Encoded the char passed using the following system:
   * A:0, C:1, G:2, T:4, D:5, I:6
   */
  int encoded_nuc = -1;
  
  if (tmp == 'A') {
    encoded_nuc = 0;
  }else if(tmp == 'C'){
    encoded_nuc = 1;
  }else if (tmp == 'G'){
    encoded_nuc = 2;
  }else if(tmp == 'T'){
    encoded_nuc = 3;
  }else if(tmp == 'D'){
    encoded_nuc = 4;
  }else if(tmp == 'I'){
    encoded_nuc = 5;
  }
  return(encoded_nuc);
}

void parse_md_tag(uint8_t *aux, std::vector<int> &haplotypes, std::vector<int> &positions){
  /*
   * @param aux : the md tag
   * @param haplotypes : vector with encoded nuc haplotypes
   * @param positions : the positions of sub,ins, and dels that make the haplotype
   *
   * Parse the MD tag, populating the haplotypes and positions of the amplicon with
   * substitutions.
   * Any discrepency between the cigar string and MD tag means insertions are in the read.
   * Cigar code handles insertions and deletions elsewhere.
   */
  
  //std::cout << "aux " << aux << std::endl;
  int total = 0; //the total num bases accounted for in the MD tag, any less than cigar means insertion
  bool last_digit = false; //helping to track total nums
  bool deletion = true;
  std::string digits; //helping to track number operations 
  std::string deletion_nucs; 

  //length of this is random
  for(int i = 1; i < 100; i++){
    char tmp = aux[i];
    //if we reach the end of the unsigned char we stop iterating
    if(tmp == '\0'){
      //this makes sure if we end with a number it gets accounts for 
      if(last_digit){
        int a = std::stoi(digits);
        total += a;
        last_digit = false;
        digits.clear();
      }
      break;
    }
    if(isdigit(tmp)){ //on digit character
      if(deletion){
        deletion = false;
        positions.push_back(total);
        char del = 'D';
        haplotypes.push_back(encoded_nucs(del)); //this needs to be generalize to encode deletions that aren't a single NT
        deletion_nucs.clear();
        total += deletion_nucs.length();
      }
      last_digit = true;
      digits += tmp;
    }else if (isalpha(tmp)){ //on alpha character
      if(deletion){
        deletion_nucs += tmp;
        continue;
      }
      if(last_digit){
        int a = std::stoi(digits);
        total += a;
        last_digit = false;
        digits.clear();
      }
      positions.push_back(total);
      haplotypes.push_back(encoded_nucs(tmp));
      total += 1;
    }else if (tmp == '^'){
      deletion = true; 
    }
  }
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
void k_means(int n_clusters, alglib::real_2d_array xy, cluster &cluster_results){
  /*
  * @params n_clusters : then number of clusters
  * @params xy : matrix to do kmeans clustering on
  * @params cluster_results : object storing the output of k means clustering
  *
  * Performs K means ++ clustering.
  */
  alglib::clusterizerstate s;
  alglib::kmeansreport rep;
  cluster_results.n_clusters = n_clusters;
  
  int num_points = 6; //number of points in the data
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

void iterate_reads(bam1_t *r, IntervalTree &amplicons){
  /*
   * @param r : alignment object
   * @param track_haplotypes : vector containing haplotype objects per amplicon
   * 
   * In this function we encode the changes to the reference as follows:
   * A:0, C:1, G:2, T:4, Del:5, Ins:6
   */

  //get the cigar operation 
  uint32_t *cigar = bam_get_cigar(r);

  uint32_t i = 0;  

  //keep track of operation and length of operation
  uint32_t op;
  uint32_t op_len;

  //1 = A, 2 = C, 8 = G, 15 = N bits 
  //get a pointer to the sequence
  uint8_t *seq = bam_get_seq(r);

  //get pointer to MD:Z tag which tells us where substitutions are
  uint8_t *aux = bam_aux_get(r, "MD");
 
  //temp variable to count positions
  int start = 0;
  //these need to be tied to the amplicon object
  std::vector<int> haplotypes;
  std::vector<int> positions;
  std::cout << "\n";
  parse_md_tag(aux, haplotypes, positions); //total bases accounted in md tag
  std::cout << aux << std::endl;                                                                                                                
 
  //this refers to the start position relative to the reference
  //bool reverse = bam_is_rev(r);
  int abs_start_pos = r->core.pos;
  int abs_end_pos  = bam_endpos(r) - 1;

  //iterate through cigar ops for this read
  while(i < r->core.n_cigar){   
    //std::cout << "\n"; 
    //gives us the left mose mapping position comnsimes the reference, sc doesnt consume ref
    //assign this read to an amplicon
    //std::cout << "pos " << r->core.pos << std::endl;
    //std::cout << "end pos " << bam_endpos(r) << std::endl;

    op = bam_cigar_op(cigar[i]); //cigar operation
    op_len = bam_cigar_oplen(cigar[i]); //cigar length
    start += op_len; //total positions iterated
                     
    //these are the only operations we care about
    if(op == 1){
      std::cout << op << " " << op_len << std::endl;
      std::cout << aux << std::endl;
      std::cout << "bam seq i " << bam_seqi(seq, 22) << std::endl; 
    }
    i++;
  }
  //places haplotype on amplicon node
  amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions);
  //amplicons.print_amplicon_info();

}



//entry point for threshold determination
void determine_threshold(std::string bam, std::string bed, std::string pair_info, int32_t primer_offset = 0){
  /*
   * @param bam : path to the bam file
   * @param bed : path to the bed file
   * @param pair_info : path to the primer pair .tsv file
   * @param primer_offset : 
   */

  //initialize haplotype data structure
  //std::vector<haplotype> track_haplotypes;
  std::vector<primer> primers;
  IntervalTree amplicons;
  //populate primer, and primer pairs
  primers = populate_from_file(bed, primer_offset);
  amplicons = populate_amplicons(pair_info, primers);
 
  samFile *in = hts_open(bam.c_str(), "r");
  hts_idx_t *idx = sam_index_load(in, bam.c_str());
  bam_hdr_t *header = sam_hdr_read(in);
  bam1_t *aln = bam_init1();
  hts_itr_t *iter = NULL;
  std::string region_;
  //region refers to reference
  region_.assign(header->target_name[0]);
  
  iter  = sam_itr_querys(idx, header, region_.c_str());
  //test variable
  int j = 0;

  //this iterates over the reads and assigns them to an amplicon
  while(sam_itr_next(in, iter, aln) >= 0) {
    //pull out the relevant diff from reference
    if (j < 107000){
      j++;
      continue;
    }
    iterate_reads(aln, amplicons);
    j++;
    std::cout << j<< std::endl;
    if(j > 107020){
      break;
    }
  }
  
  //extract those reads into a format useable in the clustering
  //frequency_matrix = create_frequency_matrix(amplicons)
  
}

