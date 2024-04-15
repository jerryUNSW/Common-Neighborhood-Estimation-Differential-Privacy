#pragma once
#pragma once
#ifndef __UTILITY_H
#define __UTILITY_H

#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <chrono>
#include <atomic>
#include <mutex>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <queue>
#include <map>
// #include <string>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <iterator>
#include <cassert>
#include <queue> 
#include <random>

#include <omp.h>  // Add this line to include omp.h
#include <cmath>

// for solving equations:
#include <boost/math/tools/roots.hpp>
// #include <nlopt.hpp>

#include <gsl/gsl_multimin.h>

// for LDP: 
// #include "mt19937ar.h"
#include "MemoryOperation.h"
#include "include/stats.hpp"

typedef unsigned int vid_t;
typedef int num_t;

using namespace std;
/** Macros **/

#define MIN(a, b) (a <= b ? a : b)
#define MAX(a, b) (a >= b ? a : b)

extern long double power(long double base, int exponent) ;

extern double compute_epsilon(double exposure_risk, int m, int n1, int n2); 

extern long double add_geometric_noise(long double data, long double sensitivity, long double epsilon); 

extern double cosineSimilarity(const std::vector<double>& vectorA, const std::vector<double>& vectorB); 

extern double computeMSE( long double btf, const std::vector<double>& numbers); 

extern double computeStd(const std::vector<double>& numbers) ;

extern double computeCovariance(const std::vector<double>& x, const std::vector<double>& y);

extern double computeVariance(const std::vector<double>& numbers);

// extern double calculateMean(const std::vector<int>& vector);
// template <typename T>
// double calculateMean(const std::vector<T>& vector);

template <typename T>
long double calculateMean(const std::vector<T>& vector) {
    long double sum = 0;
    for (int i = 0; i < vector.size(); ++i) {
        sum += static_cast<double>(vector[i]);
    }
    return sum / vector.size();
}

template <typename T>
long double calculateRMSE(const std::vector<T>& predicted, const std::vector<T>& actual) {
    if (predicted.size() != actual.size() || predicted.empty()) {
        // Handle error: sizes don't match or vectors are empty
        return 0.0; // Or throw an exception
    }

    long double sumSquaredError = 0.0;
    for (size_t i = 0; i < predicted.size(); ++i) {
        sumSquaredError += static_cast<long double>((predicted[i] - actual[i]) * (predicted[i] - actual[i]));
    }

    return std::sqrt(sumSquaredError / predicted.size());
}

template <typename T>
long double calculateMAE(const std::vector<T>& predicted, const std::vector<T>& actual) {
    if (predicted.size() != actual.size() || predicted.empty()) {
        // Handle error: sizes don't match or vectors are empty
        return 0.0; // Or throw an exception
    }

    long double sumAbsoluteError = 0.0;
    for (size_t i = 0; i < predicted.size(); ++i) {
        sumAbsoluteError += std::abs(predicted[i] - actual[i]);
    }

    return sumAbsoluteError / predicted.size();
}

extern double calculateCorrelation(const std::vector<int>& vectorX, const std::vector<int>& vectorY);


extern double manhattanDistance(const std::vector<double>& vectorA, const std::vector<double>& vectorB); 

// linked list
struct Node
{
	int id;
	unordered_set<vid_t> vertices;
	Node *prev;
	Node *next;
	Node(){
		id=-1;
	}
};	
// doubly linked list, without a copy constructor or assignment operator
class List
{
	public:
	Node *head, *tail;
	vector<Node*> in_ptrs;
	int num_components;
	int num_nodes ;
	public:
	List(int n)
	{
		head=NULL;
		tail=NULL;
		num_components=0;
		num_nodes=0;
		in_ptrs.resize(n);
		for(int i=0;i<in_ptrs.size();i++){in_ptrs[i]=NULL;}
	}
	void insert_last(unordered_set<vid_t> data);	
	void display();
	void delete_first();
	void delete_last();
	void delete_node(Node* ptr); 
	void clear(){
		while(head){
			delete_first();
		}
	}
	~List(){
		clear();
	}
	int get_UB(int i ){
		return this->in_ptrs[i]->vertices.size();
	}
};
// linear heap 
class LinkNode {
	public:
	LinkNode *prev;
	LinkNode *next;
	int id; 
	int val;
	public:
	LinkNode();
};
class LinearHeap {
	public:
		int n, key_cap;
		int active;
		int max_key; 
		int min_key;
		LinkNode *nods;
		LinkNode *rank;
	public:
		LinearHeap(int n, int key_cap);
		// LinearHeap& operator=(const LinearHeap &obj);//assignment operator 
		LinearHeap(const LinearHeap &obj);//copy constructor, something wrong here
		void clear(){
			for(int i=0;i<n;i++){
				this->unlink(i); 
				this->nods[i].val=-1;
			}
		}
		void link(int v, int r);
		void unlink(int v);
		bool is_linked(int v);
		int get_rank(int v);
		bool has_rank(int r);
		int get_first_in_rank(int r);
		~LinearHeap();
		void print_heap();
		int get_min_key();
		int get_max_key();
		void update_key(int v, int new_rank);
};

struct node{
	friend bool operator< (node a, node b){
		return a.key > b.key;
	}
	node(int k, int id_){
		this->key=k;
		this->id=id_;
	}
	int key; 
	vid_t id;
};

class InputParser {
public:
	InputParser(int &argc, char **argv) {
		for (int i = 1; i < argc; ++i)
			this->tokens.push_back(std::string(argv[i]));
	}
	const std::string& getCmdOption(const std::string &option) const {
		std::vector<std::string>::const_iterator itr;
		itr = std::find(this->tokens.begin(), this->tokens.end(), option);
		if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
			return *itr;
		}
		static const std::string empty_string("");
		return empty_string;
	}
	bool cmdOptionExists(const std::string &option) const {
		return std::find(this->tokens.begin(), this->tokens.end(), option)
			!= this->tokens.end();
	}
private:
	std::vector <std::string> tokens;
};

// Computes n choose k
extern double n_choose_k(int n, int k);

enum class Index_update { withlimit, withlimit_base_opt, withlimit_parallel, withlimit_dfs, withlimit_dfs_parallel, withoutlimit };

extern std::ostream& operator<<(std::ostream& out, const Index_update value);

extern void get_anchor_groups(std::vector<vid_t> candidates, int budget, std::vector<std::vector<vid_t>>& result);

extern void zip( const std::vector<vid_t> &a, const std::vector<int> &b, std::vector<std::pair<vid_t,int>> &zipped);

extern void unzip( const std::vector<std::pair<vid_t, int>> &zipped, std::vector<vid_t> &a, std::vector<int> &b);

// make sure the vectors are sorted before applying these functions
extern vector<vid_t> vector_intersection(vector<vid_t> A, vector<vid_t> B);

extern vector<vid_t>  vector_diff(vector<vid_t> A, vector<vid_t> B);

extern vector<vid_t>  vector_union(vector<vid_t> A, vector<vid_t> B);

extern vector<vid_t> sort_by_left(vector<vid_t> old_ids, vector<int>& UB);


extern void print_dash(int n);

// template functions
template <typename T> void vector_show(vector<T> A)
{
	for(auto kkk:A){cout<<kkk<<" ";}cout<<endl;
}

#endif  /* __UTILITY_H */