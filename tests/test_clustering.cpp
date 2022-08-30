#include<iostream>
#include "ap.h"
#include "../src/interval_tree.h"
#include "../src/primer_bed.h"
#include "../src/clustering.h"

//test the clustering function
void _unit_test_kmeans(){
  k_means(3);
}

int main(){
  //determine_threshold();
  //_unit_test_kmeans();
  determine_threshold("/Users/caceves/Desktop/file_227.final.bam", "/Users/caceves/Desktop/sarscov2_v2_primers.bed", "/Users/caceves/Desktop/primer_pairs.tsv", 0);
  return(0);
}
