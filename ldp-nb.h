#pragma once
#ifndef ABCORE_H
#define ABCORE_H
#include "bigraph.h"
#include <sstream>

void randomizedResponses(BiGraph& g, BiGraph& g2, 
    int query_vertex, int flipping_from, int flipping_to, double p);


// Define function fff(x, alpha)
double fff(double x, double alpha) ; 
// Define objective function to minimize
double objective(const gsl_vector *v, void *params) ;

std::vector<double> minimize(double d1__, double d2__, double epsilon__) ;

long double solve_equation(long double d, long double epsilon) ;
// Function to uniformly sample K vertices from a range
std::vector<int> uniformSampleKVertices(int K, bool useFirstRange, BiGraph& g); 

long double check_two_hop_path(int s, int v, int t, BiGraph& g2);

double locally_compute_f_given_q_and_x(int q, int x, BiGraph& g, BiGraph& g2) ; 

double compute_f_given_q_x(int q, int x, BiGraph& g2) ;

vector<double> locally_compute_num_two_paths_given_q(int q, BiGraph& g, BiGraph& g2);

vector<double> locally_compute_twoHop_given_q(int q, BiGraph& g, BiGraph& g2); 

long double locally_check_distance_two(int s, int t, BiGraph& g, BiGraph& g2) ; 


vector<double> compute_twoHop_given_q(int q, BiGraph& g2); 

long double one_round_even_partition(BiGraph& g, int k1, int k2, unsigned long seed);

long double one_round_biclique(BiGraph& g, unsigned long seed, int p__, int q__, int _switch);

std::vector<std::vector<int>> generateCombinations(const std::vector<int>& up_vs, int p) ;

void private_estimate_of_degrees(BiGraph& g); 

long double BFC_EVP(BiGraph& g);

bool satisfy_bound(long double upper, long double lower, int u, int x, int v, int w, BiGraph& g); 

double my_genrand_real2();

void my_init_genrand(unsigned long seed);

// layer based approach:
long double layer_based_two_round(BiGraph& g, unsigned long seed); 

// u>x && v>w && u >v , each btf is counted once.
long double reduce_redundancy_two_round(BiGraph& g, unsigned long seed); 

// with edge dependency.
long double reduce_redundancy_two_round_dependency(BiGraph& g, unsigned long seed);

// new method two round
long double new_method(BiGraph& g, unsigned long seed); 

// long double BFC_EVP_sample(BiGraph& g, long double sampling_prob); 

void compute_m3_m2(long double &m4, long double &m3, long double &m2, long double &m1, long double &m0, BiGraph& g2);

long double one_round_btf(BiGraph& g, unsigned long seed);

long double two_round_btf(BiGraph& g, unsigned long seed); 

// naive noisy motif counts
void get_noisy_naive(BiGraph& g, BiGraph& g2, long double& local_btfs, long double& local_cate, long double& res ); 
// process butterflies in batch 
void BFC_EVP_noisy(BiGraph& g, BiGraph& g2, long double& BTF, long double& cate, long double& res );

void approx_abcore(BiGraph& g, long double Eps0); 

void construct_noisy_graph(BiGraph& g, BiGraph& g2, unsigned long seed);

void construct_noisy_graph_2(int upper, int lower, BiGraph& g, BiGraph& g2, unsigned long seed);

void construct_noisy_graph_3(int upper, int lower, BiGraph& g, BiGraph& g2, unsigned long seed); 

long double get_wedges(BiGraph& g);

long double get_cate(BiGraph& g);

int findPercentile(std::vector<int>& data); 

void test_random_number(unsigned long seed); 

long double get_laplace(long double parameter);

// void BFS(vid_t v, vector<bool>& visted, BiGraph& g, unordered_set<vid_t> &result);
// void get_connected_components(BiGraph& g, List &list_);
void compute_m3_m2_2(long double &m4, long double &m3, long double &m2, long double &m1, long double &m0, BiGraph& g2);


// Custom hash function for std::vector<int>
struct VectorHasher {
    std::size_t operator()(const std::vector<int>& v) const {
        std::size_t seed = v.size();
        for (const auto& i : v) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

std::vector<std::vector<int>> generateDistinctMotifs(int p, int q, BiGraph& g, int num_motifs) ;
#endif