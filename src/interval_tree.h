#include <iostream>
#include "primer_bed.h"
using namespace std;

#ifndef interval_tree
#define interval_tree

// Structure to represent an interval
class Interval{   public:
  Interval(int val1, int val2, int val3, int val4): low(std::min(val1, val2)), high(std::max(val1, val2)), low_inner(val3), high_inner(val4) {}  // constructor
  int low, high, low_inner, high_inner;
};
// Structure to represent a node in Interval Search Tree
class ITNode{
  /*
    public:
    ITNode(Interval *values): data(value), left(nullptr), right(nullptr) {}  // constructor
    int max;
    // Getters - access member functions
    Interval getData()const;
    ITNode getLeft()const;
    ITNode getRight()const;
    // Setters - access member functions
    void setLeft(ITNode *node);
    void setRight(ITNode *node);
  */
public:
  ITNode(Interval value): data(new Interval(value)), left(nullptr), right(nullptr), max(value.high) {}  // constructor
  Interval *data;  // pointer to node's interval data object
  ITNode *left, *right; // pointer to node's left & right child node objects
  int max;
  std::vector<std::vector<int>> haplotypes;
  std::vector<std::vector<uint32_t>> positions;
  int read_count=0;

};

/////////////////////////////////////////////////////////////////////////////////////////
// IntervalTree class
class IntervalTree{
private:
  ITNode *_root;
  void insert(ITNode *root, Interval data);
  bool envelopSearch(ITNode *root, Interval data);
  void inOrder(ITNode * root); //prints amplicons starts and ends
  void print_amplicon_info(ITNode *root); //prints haplotypes and positions
  void find_amplicon_per_read(ITNode *root, int start, int end, std::vector<int> haplotypes, std::vector<uint32_t> positions, bool reverse); //places read info in amplicon

public:
  IntervalTree();  // constructor
  void insert(Interval data){ insert(_root, data);}
  bool envelopSearch(Interval data){ return envelopSearch(_root, data);}
  void inOrder() {inOrder(_root);}
  void print_amplicon_info() {print_amplicon_info(_root);}
  ITNode *iterate_nodes(ITNode *root); //used to returns nodes iteratively
  ITNode *iterate_nodes(){return iterate_nodes(_root);}
  void find_amplicon_per_read(int start, int end, std::vector<int> haplotypes, std::vector<uint32_t> positions, bool reverse){find_amplicon_per_read(_root, start, end, haplotypes, positions, reverse);}
};

IntervalTree populate_amplicons(std::string pair_info_file, std::vector<primer> &primers);

#endif
