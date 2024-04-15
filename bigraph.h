#pragma once
#ifndef __BIGRAPH_H
#define __BIGRAPH_H
#include "utility.h"

using namespace std;
using namespace std::chrono; 

// std::mt19937 rng(std::random_device{}()); // for seeding

// for each root, store the other vertices in the CC of root. 


class BiGraph{
public:
	// BiGraph(std::string dir);
	BiGraph(string dir, bool is_bigraph );
	BiGraph(const BiGraph& other); 
	BiGraph();
	~BiGraph() {
		// destroy();
	}
	// note that addEdgeRaw is processing raw ids.
	void addEdgeRaw(vid_t u, vid_t v);

	// addEdge is processing real ids. 
	void addEdge(vid_t u, vid_t v);
	// deleteEdge is processing real ids. 
	void deleteEdge(vid_t u, vid_t v);

public:

    // Function to sample K pairs of vertices from the same layer in a bipartite graph
    vector<pair<vid_t, vid_t>> sample_bipartite_pairs(int K, int balancing_factor) {
        vector<pair<vid_t, vid_t>> pairs;
        vector<vid_t> upper_layer_vertices, lower_layer_vertices;

        // Collect vertices from upper and lower layers
        for (vid_t u = 0; u < num_v1; ++u) {
            upper_layer_vertices.push_back(u);
        }
        for (vid_t v = num_v1; v < num_nodes(); ++v) {
            lower_layer_vertices.push_back(v);
        }

        // Randomly sample K pairs of vertices, ensuring both vertices in each pair are from the same layer
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> upper_dist(0, upper_layer_vertices.size() - 1);
        uniform_int_distribution<> lower_dist(0, lower_layer_vertices.size() - 1);

        // for (int i = 0; i < K; ++i) {
		while(pairs.size()<K){
            vid_t vertex1, vertex2;
            if (uniform_int_distribution<>(0, 1)(gen) == 0) {
                // Sample from upper layer
                int idx1 = upper_dist(gen);
                vertex1 = upper_layer_vertices[idx1];
                int idx2 = upper_dist(gen);
                vertex2 = upper_layer_vertices[idx2];
            } 
			else {
                // Sample from lower layer
                int idx1 = lower_dist(gen);
                vertex1 = lower_layer_vertices[idx1];
                int idx2 = lower_dist(gen);
                vertex2 = lower_layer_vertices[idx2];
            }

			if(balancing_factor==0){
				pairs.push_back(make_pair(vertex1, vertex2));
			}else{
				int multiples = pow(10, balancing_factor); 
				// int multiples = balancing_factor;
				if (degree[vertex1] > 0 && degree[vertex2] > 0 &&
					(degree[vertex1] > degree[vertex2] * multiples ||
					degree[vertex2] > degree[vertex1] * multiples)) {
					pairs.push_back(std::make_pair(vertex1, vertex2));
				}
			}
        }
        return pairs;
    }


	void computePriority();
	void init(unsigned int num_v1, unsigned int num_v2);

	// for bipartite graphs:
	void loadGraph(std::string dir);

	// for general graphs
	void load_general_graph(std::string);

	void print_graph();
	void show_nodes();
	// vector<LinearHeap> make_heaps();
	bool same_layer(vid_t u, vid_t v ){
		if(is_upper(u) && is_upper(v)) return true; 
		if(is_lower(u) && is_lower(v)) return true; 
		return false ; 
	}

	bool diff_layer(vid_t u, vid_t v ){
		return !same_layer(u,v);
	}

	bool is_upper(vid_t u){ 
		return (u < this->num_v1);
	}
	bool is_lower(vid_t u){
		return (u>= this->num_v1);
	}
	// void generalGraph(string dir);

	bool com_p(vid_t u, vid_t v){
		// return prio[u] > prio[v]; 
		return prio[u] < prio[v]; 
	}

	bool is_active(int vertex){ 
		return degree[vertex]>0;
	}
	int num_nodes(){ 
		if(is_bipartite){
			return this->num_v1+this->num_v2;
		}else{
			return num_vertices;
		}
		// return this->num_v1+this->num_v2;
	}

	bool has(int i, int j){
		return std::find(neighbor[i].begin(), neighbor[i].end(), j) != neighbor[i].end(); 
	}
	bool hasEdge(int u, int v) {
		if (u < 0 || u >= num_nodes() || v < 0 || v >= num_nodes()) {
			cout<<"Invalid vertex index\n";
			// return false;
			exit(1);
		}
		// Check if u and v belong to different layers
		if (same_layer(u,v)) {
			// Invalid edge between upper and lower vertices
			cout<<"Invalid edge between upper and lower vertices"<<endl;
			// return false;
			exit(1);
		}
		// Get the degrees of vertices u and v
		int degreeU = degree[u];
		int degreeV = degree[v];
		
		// Determine which vertex to scan first based on their degrees
		int vertexToScanFirst = (degreeU < degreeV) ? u : v;
		int vertexToScanSecond = (vertexToScanFirst == u) ? v : u;
	
		// Check if vertexToScanFirst is connected to vertexToScanSecond
		for (const auto& neighbor : neighbor[vertexToScanFirst]) {
			if (neighbor == vertexToScanSecond) {
				// Edge found
				return true;
			}
		}
		// Edge not found
		return false;
	}

	// void show_amat(); 
	
	// we also need a function that returns all vertices that are two-hop reachable from q. 
    vector<int> twoHopNeighbors(int q, vector<int>& two_h_paths) {
        vector<int> twoHop;
        unordered_set<int> visited;

        for (int u : neighbor[q]) {
            for (int x : neighbor[u]) { 
				// if x is not a one-hop neighbor of q AND x != q. 
				// I feel like there is something wrong with this code
                // if (std::find(neighbor[q].begin(), neighbor[q].end(), x) == neighbor[q].end() && (x!= q) ) {
				if ( x!= q ) {
					// cout<<"visiting "<<q<<" "<<u<<" "<<x<<endl;
					// if we want to returns exactly the vertices such that distance(q, x) = 2.
					// we need std::find(neighbor[q].begin(), neighbor[q].end(), x) == neighbor[q].end()
					if(visited.find(x) == visited.end()){
						// if not visited. 
						twoHop.push_back(x);
						visited.insert(x);
					}
					// in here, we still visit all 2-hop paths. 
					two_h_paths[x]++;
                }
            }
        }
        return twoHop;
    }

	vector<int> threeHopNeighbors(int q, vector<int>& three_h_paths) {
		vector<int> threeHop;
		unordered_set<int> visited;

		// Iterate over the neighbors of vertex q
		for (int u : neighbor[q]) {
			for (int v : neighbor[u]) {
				if(v==q) continue;
				for (int x : neighbor[v]) {
					if(x==u) continue;
					if(x==q) continue;

					if(visited.find(x) == visited.end()){
						// if not visited. 
						threeHop.push_back(x);
						visited.insert(x);
					}

					three_h_paths[x]++;
				}
			}
		}
		return threeHop;
	}



	std::string dir;

	num_t num_v1;
	num_t num_v2;

	long num_edges;

	long num_vertices;

	vector<vector<vid_t>> neighbor; 
	vector<int> degree; 

	vector<int> prio;

	vector<int> U, V;

	vector<int> nodes; 

public:
	// max and min degrees
	int v1_max_degree;
	int v2_max_degree;

	// the switch determing if the incoming graph is bipartite
	bool is_bipartite = true; 
};

extern void computeCore(int ALPHA, int BETA, BiGraph& g);

extern void peeling(LinearHeap& degree_heap, int ENG, unordered_set<vid_t>& affected_vertices, BiGraph& g);

extern void update_engagement(LinearHeap& upper_heap, LinearHeap& lower_heap, BiGraph& g, unordered_set<vid_t>& affected_nodes); 


#endif  /* __BIGRAPH_H */