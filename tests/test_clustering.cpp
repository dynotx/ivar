#include<iostream>
#include "ap.h"
#include "../src/clustering.h"

//test the clustering function
void _unit_test_kmeans(){
  k_means(3);
}

int main(){
  //determine_threshold();
  //_unit_test_kmeans();
  determine_threshold("/Users/caceves/Desktop/file_227.final.bam");
  return(0);
}
