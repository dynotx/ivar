#include<iostream>
#include "../src/clustering.h"

//test the clustering function
void _unit_test_kmeans(){
  k_means(3);
}

int main(){
  std::cout << "In clustering test code" << std::endl;
  _unit_test_kmeans();
  return(0);
}
