#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdint.h>
#include <algorithm>
#include <string.h>
#include <sstream>
#include <vector>
#include "ap.h"
#include "dataanalysis.h"
#include "interval_tree.h"
#include "allele_functions.h"
#include "primer_bed.h"
#include "stdafx.h"
#include "clustering.h"
using namespace alglib;

std::vector<std::vector<uint32_t>> transpose(const std::vector<std::vector<uint32_t>> data) {
  // this assumes that all inner vectors have the same size and
  // allocates space for the complete result in advance
  std::vector<std::vector<uint32_t> > result(data[0].size(),
                                        std::vector<uint32_t>(data.size()));
  for (std::vector<uint32_t>::size_type i = 0; i < data[0].size(); i++) 
      for (std::vector<uint32_t>::size_type j = 0; j < data.size(); j++) {
          result[i][j] = data[j][i];
      }
  return result;
}


std::vector<std::vector<int>> transpose(const std::vector<std::vector<int>> data) {
  // this assumes that all inner vectors have the same size and
  // allocates space for the complete result in advance
  std::vector<std::vector<int> > result(data[0].size(),
                                        std::vector<int>(data.size()));
  for (std::vector<int>::size_type i = 0; i < data[0].size(); i++) 
      for (std::vector<int>::size_type j = 0; j < data.size(); j++) {
          result[i][j] = data[j][i];
      }
  return result;
}

std::vector<uint32_t> flatten(std::vector<std::vector<uint32_t>> const &vec)
{
    std::vector<uint32_t> flattened;
    for (auto const &v: vec) {
        flattened.insert(flattened.end(), v.begin(), v.end());
    }
    return flattened;
}

void zip(
    const std::vector<int> &haplotypes, 
    const std::vector<uint32_t> &positions, 
    std::vector<std::pair<uint32_t,int>> &zipped)
{
    for(size_t i=0; i< positions.size(); ++i)
    {
        zipped.push_back(std::make_pair(positions[i], haplotypes[i]));
    }
}

void unzip(const std::vector<std::pair<uint32_t, int>> &zipped, std::vector<int> &haplotypes, std::vector<uint32_t> &positions){
  for(size_t i=0; i < positions.size(); i++){
    positions[i] = zipped[i].first;
    haplotypes[i] = zipped[i].second;
  }
  positions.resize(zipped.size());
  haplotypes.resize(zipped.size());
}

void reorder_haplotypes(std::vector<int> &haplotypes, std::vector<uint32_t> &positions){
  /*
   * @param haplotypes : the hapltyope encoded nts
   * @param positions : the position relative to the reference
   *
   * Function reorders both vectors so that the positions are in relative order.
   */

    // Zip the vectors together
    std::vector<std::pair<std::uint32_t, int>> zipped;
    zip(haplotypes, positions, zipped);

    // Sort the vector of pairs
    std::sort(std::begin(zipped), std::end(zipped));

    // Write the sorted pairs back to the original vectors
    unzip(zipped, haplotypes, positions);

}

//print the cluster infor
void print_cluster_info(cluster cluster_results){
  /*
   * @param cluster_results : the cluster object containing clustering stats
   */
  std::cout << "N Clusters: " << cluster_results.n_clusters << std::endl;
  std::cout << "Sil Score: " << cluster_results.sil_score << std::endl;
  for(float x: cluster_results.centers){
    std::cout << "center @: " << x << std::endl;
  }
}

//get the average of a vector
float average(std::vector<float> x){
  float sumTotal = 0;
  for(float k=0; k < x.size(); ++k){
      sumTotal += x[k];            
  }
  return(sumTotal / x.size());
}

int encoded_nucs(std::string &tmp){
  /*
   * @param tmp : nucleotide to encode
   * @return encoded_nuc : nucleotide encoded as an int
   *
   * Encoded the char passed using the following system:
   * A:0, C:1, G:2, T:3, D:4, I:5, S:6
   */
  int encoded_nuc = -1;
  
  if (tmp == "A") {
    encoded_nuc = 0;
  }else if(tmp == "C"){
    encoded_nuc = 1;
  }else if (tmp == "G"){
    encoded_nuc = 2;
  }else if(tmp == "T"){
    encoded_nuc = 3;
  }else if(tmp == "D"){
    encoded_nuc = 4;
  }else if(tmp == "I"){
    encoded_nuc = 5;
  }else if(tmp == "S"){
    encoded_nuc = 6;
  }
  return(encoded_nuc);
}


void parse_md_tag(uint8_t *aux, std::vector<int> &haplotypes, std::vector<uint32_t> &positions,
    uint32_t abs_start_pos, std::vector<position> &all_positions,
    uint8_t *seq, uint32_t length, uint32_t correction_factor){
  /*
   * @param aux : the md tag
   * @param haplotypes : vector with encoded nuc haplotypes
   * @param positions : the positions of sub,ins, and dels that make the haplotype
   * @param abs_start_pos : the start position relative to the reference (is actually abs_end_pos) if reverse
   * @param reverse : whether this is a reverse read or not
   * @param ad : vector containing variant alleles
   * @param seq : pointer to the sequence
   *
   * Parse the MD tag, populating the haplotypes and positions of the amplicon with
   * substitutions.
   * Any discrepency between the cigar string and MD tag means insertions are in the read.
   * Cigar code handles insertions and deletions elsewhere.
   */
  if(correction_factor == 0){
    correction_factor = -1;
  }
  //std::cout << "aux " << aux << std::endl;
  bool deletion = false; //helping track the deletions
  std::string digits; //helping to track number operations 
  std::string nucs;
  std::string soft = "S";
  std::string del = "D";
  std::string nt;
  int i = 0;
  
  std::vector<std::string> nucleotides; //store the substituions & deletions

  do {
    char tmp = aux[i]; //this is the reference nuc   
    //std::cout << "tmp " << tmp << std::endl;
    if(isdigit(tmp)){ //on digit character
      if(deletion){
        deletion = false;
        abs_start_pos += std::stoi(digits) + nucs.length();
        positions.push_back(abs_start_pos);
        nt = seq_nt16_str[bam_seqi(seq, abs_start_pos+correction_factor - length)];
        haplotypes.push_back(encoded_nucs(nt));
        nucleotides.push_back(nt);
        //std::cout << "correction factor " << correction_factor << " " << length << std::endl;
        //std::cout << abs_start_pos+correction_factor - length << std::endl;
        //std::cout << abs_start_pos << " " << nt << std::endl;
        nucs.clear();
      }
      deletion = false;
      digits += tmp;
    }else if (isalpha(tmp)){ //on alpha character, not the very first one
      //std::cout << "here nucs " << nucs << " digits " << digits << " tmp " << tmp << std::endl;
      if(deletion){ //we already know this will be a deletion
        nucs += tmp;
        i++;
        continue;
      }
      if(digits != ""){
        //check if the position already got soft clipped
        abs_start_pos += std::stoi(digits) + nucs.length();
        positions.push_back(abs_start_pos);
        nt = seq_nt16_str[bam_seqi(seq, abs_start_pos+correction_factor - length)];
        nucleotides.push_back(nt);
        haplotypes.push_back(encoded_nucs(nt));
        //std::cout << "correction factor " << correction_factor << " " << length << std::endl;
        //std::cout << abs_start_pos+correction_factor - length << std::endl;
        //std::cout << abs_start_pos << " " << nt << std::endl;
        digits.clear();
        nucs.clear();
        nucs += tmp;

      }else{
        nucs += tmp;
      }
    }else if (tmp == '^'){
      //std::cout << "is del" << std::endl;
      deletion = true; 
    }
    i++;
  } while(aux[i] != '\0');
  update_allele_depth(all_positions, nucleotides, positions);
  std::cout << "\n";
}

//calculate the cluster centers
void calculate_cluster_centers(alglib::real_2d_array X, alglib::kmeansreport rep, int n_clusters, cluster &cluster_results){
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
  cluster_results.centers = centers;

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
void calculate_sil_score(alglib::real_2d_array X, alglib::kmeansreport rep, 
    int n_clusters, cluster &cluster_results){
  /*
  * @params X : array containing data points
  * @params rep : kmeans reported object
  * @params n_clusters : the number of clusters
  * @params cluster_results : the object for storing clustering metrics
  *
  * Calculate the silohuette score for each sample, and then average them
  * together to find the cumulative score.
  * (b - a) / max(a, b) sil score formula
  */
  //std::cout << "Calculating silohuette score." << std::endl;
  int rows = X.rows();
  float point = 0;
  float center = 0;
  std::vector<float> sil_scores;
  float tmp = 0;

  for (int i = 0; i < rows; i++) {
    point = X[i][0];
    center = rep.cidx[i];
    tmp = cluster_point_distances(X, rep, point, center, n_clusters);
    sil_scores.push_back(tmp);
  }

  //store the cumulative results
  cluster_results.sil_score = average(sil_scores);
  cluster_results.sil_scores = sil_scores;
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
  int num_points = xy.rows(); //number of points in the data
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
    std::cout << "Error in clustering haplotypes" << std::endl;
    exit(1);
  }

  /*
  int i = 0;
  while(i < n_clusters){
    std::cout << "center "<< i << " "  << rep.c[i][0] << std::endl;
    i++;
  }*/
 
  calculate_sil_score(xy, rep, n_clusters, cluster_results);
  calculate_cluster_centers(xy, rep, n_clusters, cluster_results);
}

void iterate_reads(bam1_t *r, IntervalTree &amplicons, std::vector<position> &all_positions){
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

  //we record the range of the read because if thye don't overlap perfectly could
  //accidentally falsely create multiple haplotypes
  std::vector<uint32_t> range;

  //1 = A, 2 = C, 8 = G, 15 = N bits 
  //get a pointer to the sequence
  uint8_t *seq = bam_get_seq(r);

  //get pointer to MD:Z tag which tells us where substitutions are
  uint8_t *aux = bam_aux_get(r, "MD");
  //temp variable to count positions
  uint32_t start = 0;
  //this refers to the start position relative to the reference
  bool reverse = bam_is_rev(r);
  uint32_t abs_start_pos = r->core.pos; //leftmost coordinate on ref
  uint32_t abs_end_pos  = bam_endpos(r); //rightmost coordinate on ref
  range.push_back(abs_start_pos);
  range.push_back(abs_end_pos);

  //these will later be place in the amplicon object
  std::vector<int> haplotypes;
  std::vector<uint32_t> positions;
  uint32_t insertion_pos = 0;
  std::string nucs;
  uint32_t correction_factor = 0;
  bool first_pass = true;

  //test lines
  std::cout << "abs start " << abs_start_pos << " abs end " << abs_end_pos << std::endl;
  //remember 0 is false in c++
  std::cout << "reverse " << reverse << std::endl;

  //iterate through cigar ops for this read
  while(i < r->core.n_cigar){   
    op = bam_cigar_op(cigar[i]); //cigar operation
    op_len = bam_cigar_oplen(cigar[i]); //cigar length
    
    //test line
    //std::cout << op << " " << op_len << std::endl;
    if(op == 4 && first_pass){
      correction_factor = op_len + 1;
      first_pass = false;
    }
    if(op != 4){ first_pass = false;}

    //these are the only operations we care about
    if(op == 1){
      if(reverse){
        insertion_pos = abs_end_pos - start;
      } else {
        insertion_pos = abs_start_pos + start;
      }
      //go get each nt in the insertion region
      for(uint32_t x = 0; x < op_len; x++){
        //the nucelotide at the insertion poi
        char nt = seq_nt16_str[bam_seqi(seq, start+x)];
        nucs += nt;
        //in this process, the positions are not longer in numeric order
        haplotypes.push_back(encoded_nucs(nucs));
        nucs.clear();
        positions.push_back(abs_start_pos+start+x);
      }
    }
    if (bam_cigar_type(op) & 2){
      start += op_len; //total positions iterated relation to ref
    }
                   
    i++;
  }
  
  parse_md_tag(aux, haplotypes, positions, abs_start_pos, all_positions, seq, abs_start_pos, correction_factor);
  //reoder the positions to be consistene
  reorder_haplotypes(haplotypes, positions);
  //places haplotype on amplicon node
  amplicons.find_amplicon_per_read(abs_start_pos, abs_end_pos, haplotypes, positions, reverse, range);
  //amplicons.print_amplicon_info();

}

std::vector<float> create_frequency_matrix(IntervalTree &amplicons){
  /*
   * @param amplicons : data strucuture containing read count, haplotype nt, and positions
   * @return frequencies : a flat vector containing all haplotype frequencies
   *
   * Function calculates the frequency of unique haplotypes on a per amplicon basis. In
   * the process, handles combining haplotypes that originate from reads that don't overlap.
   * Stores unique haplotypes, the associated positions, and frequency to amplicon object.
   */ 

  std::vector<float> frequencies;
  ITNode *node = amplicons.iterate_nodes();
  int read_count=0; //total reads in amplicon
  std::vector<std::vector<uint32_t>> positions;
  std::vector<std::vector<uint32_t>> ranges;
  std::vector<std::vector<int>> haplotypes;
  std::vector<uint32_t> position;
  std::vector<uint32_t> range;
  std::vector<int> haplotype;

  std::vector<float> unique_counts;
  std::vector<std::vector<uint32_t>> unique_positions;
  std::vector<std::vector<int>> unique_haplotypes;
  std::vector<std::vector<uint32_t>> unique_ranges;
  
  //loop through all the amplicons
  while(node != NULL){
    node = amplicons.iterate_nodes(node->right);
    if(node == NULL){
      break;
    }
    read_count = node->read_count;
    if(read_count == 0){
      continue;
    }
    positions = node->positions;
    haplotypes = node->haplotypes;
    ranges = node->ranges;
    int length = positions.size(); //total
    //std::cout << "\nLow: " << node->data->low << " High: " << node->data->high << std::endl;    
    //pool all the positions that have been modified in order to create a table
    std::vector<uint32_t> flattened = flatten(positions);
    std::sort(flattened.begin(), flattened.end());
    std::vector<uint32_t>::iterator ip = std::unique(flattened.begin(), flattened.end());
    flattened.resize(std::distance(flattened.begin(), ip));
    
    std::vector<std::vector<int>> haplo_1;
    std::vector<std::vector<uint32_t>> pos_1;

    //loop through all the haplotypes in the amplicon and find unqiue ones
    for(int i=0; i < length; i++){
      position = positions[i];
      haplotype = haplotypes[i];
      range = ranges[i];

      //fill out this haplotype with -1 for the things that aren't covered or soft clipped
      for(uint32_t p = 0; p < flattened.size(); p++){
        uint32_t tmp = flattened[p];
        //search the current haplotype for all positions seen on this amplicon
        if(std::find(position.begin(), position.end(), tmp) != position.end()){
          continue;
        }else{
          position.push_back(tmp);
          if(range[0] < tmp && tmp < range[1]){
            haplotype.push_back(-2);
          }else{
            haplotype.push_back(-1);
          }
        }
      }
      reorder_haplotypes(haplotype, position);
      haplo_1.push_back(haplotype);
      pos_1.push_back(position);
    }

    std::vector<std::vector<int>> final_haplotypes;
    std::vector<std::vector<uint32_t>> final_positions;
 
    //now we are going to remove positions that aren't relevant to due soft clipping/matching
    std::vector<std::vector<int>> transposed_vector = transpose(haplo_1);
    std::vector<std::vector<uint32_t>> transposed_positions = transpose(pos_1);

    for(uint32_t y = 0; y < transposed_vector.size(); y++){
      //check if all values are negative
      bool zeros = std::all_of(transposed_vector[y].begin(), transposed_vector[y].end(), [](int i) { return i<0; });
      //counts the depth of this allele in the haplotype
      //int count = std::count_if(transposed_vector[y].begin(), transposed_vector[y].end(), [] (int i) { return (i>=0);});
      //int read_min = 10;
      if((!zeros)){
        final_haplotypes.push_back(transposed_vector[y]);
        final_positions.push_back(transposed_positions[y]);
      }
    }
    if (final_positions.size() == 0){
      continue;
    }
    std::vector<std::vector<int>> transposed_vector_f = transpose(final_haplotypes);
    std::vector<std::vector<uint32_t>> transposed_positions_f = transpose(final_positions);
    
    bool found = false;
    //loop haplotypes, find the unique ones
    for(std::vector<int> haplo_2: transposed_vector_f){
      //check if we've seen this haplotype
      auto it = std::find(unique_haplotypes.begin(), unique_haplotypes.end(), haplo_2);
     
      //if we have seen this haplotype
      if(it != unique_haplotypes.end()) {
        //std::cout << "found in unique\n";
        /*for(std::vector<int>x:unique_haplotypes){
          for(int y:x){
            std::cout << y << ",";
          }
          std::cout << "\n";
        }*/
        int index = it - unique_haplotypes.begin();
        //make sure the nucleotides at these positions also match
        if(unique_haplotypes[index] == haplo_2){
          unique_counts[index] += 1;
        }else{
          //std::cout << "seen\n";
          //for(int u:haplo_2){std::cout << u << ", ";}
          unique_haplotypes.push_back(haplo_2);
          unique_counts.push_back(1);
        }
      }else{ //we haven't seen this haplotype before
        //std::cout << "this\n";
        //for(int u:haplo_2){std::cout  << u << ", ";}

        if(haplo_2.size() > 1){
          //lets see if we can find a neighbor, not including the soft clipped/not covered pos
          std::vector<int> results; //holds soft clipped indices
          auto it = std::find_if(std::begin(haplo_2), std::end(haplo_2), [](int i){return i == -1;});
      
          while (it != std::end(haplo_2)) {
            results.push_back(std::distance(std::begin(haplo_2), it));
            it = std::find_if(std::next(it), std::end(haplo_2), [](int i){return i == -1;});
          } 
          int j = 0;
          for(std::vector<int> uni_haplo : unique_haplotypes){
            //make a deep copy for removing SC/not covered bases
            std::vector<int> vect2;
            std::copy(uni_haplo.begin(), uni_haplo.end(), back_inserter(vect2));
            
           //replace problem indicies with 
            for(int remove:results){
              vect2[remove] = -1;
            }
            if(vect2 == haplo_2){
              found = true;
            }
            j++;
          }
        }else if((haplo_2.size() == 1) && (haplo_2[0] == -1)){
          read_count -= 1;
          continue;
        }else{
          found = false;
        }
        if(!found){
          //add it to the unique haplotypes
          unique_haplotypes.push_back(haplo_2);
          unique_counts.push_back(1);
          found=false;
        }
      }
    }
     /* for(uint32_t y:transposed_positions_f[0]){
        std::cout << y << ", ";
      }
      std::cout << "\n";
    
    
    for(std::vector<int> hap:unique_haplotypes){
      for(int h:hap){
        std::cout << h << ", ";
      }
      std::cout << "\n";
    }*/
    //save this info to the amplicon
    node->final_haplotypes = unique_haplotypes;
    node->final_positions = transposed_positions_f[0];
    for(float d: unique_counts){
      //std::cout << "uc " << d << " read count " << read_count <<std::endl;
      node->frequency.push_back(d / read_count);
      frequencies.push_back(d / read_count);
    }
    unique_haplotypes.clear();
    unique_positions.clear();
    unique_counts.clear();
    haplo_1.clear();
    pos_1.clear();
  }
  return(frequencies);
}

//entry point for threshold determination
void determine_threshold(std::string bam, std::string bed, std::string pair_info, int32_t primer_offset = 0){
  /*
   * @param bam : path to the bam file
   * @param bed : path to the bed file
   * @param pair_info : path to the primer pair .tsv file
   * @param primer_offset : 
   */

  //make sure the alleles get recorded
  std::vector<position> all_positions;

  std::string output_amplicon = "amplicon.txt";
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
  
  iter = sam_itr_querys(idx, header, region_.c_str());

  //fill md tag
  //bam_fillmd1_core(const char *ref_name, aln, char *ref, int flag, int max_nm)

  //this iterates over the reads and assigns them to an amplicon
  while(sam_itr_next(in, iter, aln) >= 0) {
    std::cout << bam_get_qname(aln) << std::endl;
    //pull out the relevant diff from reference
    iterate_reads(aln, amplicons, all_positions);
  }
  
  //extract those reads into a format useable in the clustering
  std::vector<float> all_frequencies = create_frequency_matrix(amplicons);
  for(float freq:all_frequencies){
    std::cout << freq << std::endl;
  }
  amplicons.print_amplicon_summary();  
  //amplicons.dump_amplicon_summary(output_amplicon);

  //reshape it into a real 2d array for alglib
  //alglib::real_2d_array xy;
  //xy.setlength(all_frequencies.size(), 1);
  //for(uint32_t i=0; i < all_frequencies.size(); i++){
  //  xy(i,0) = all_frequencies[i];
  //}
   
  //call kmeans clustering
  //cluster cluster_results;
  //for (int n =2; n <= 8; n++){
  //  k_means(n, xy, cluster_results);
  //  std::cout << "n " << n << " sil " << cluster_results.sil_score << std::endl;
  //}
  
  /*std::cout << cluster_results.sil_score << std::endl;
  for(float x: cluster_results.centers){
    std::cout << x << std::endl;
  }*/

}

