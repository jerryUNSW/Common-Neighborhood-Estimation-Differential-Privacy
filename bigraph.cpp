#include "bigraph.h"

using namespace std;

extern bool one_round; 

// a copy constructor: this only copies the skeleton of the graph, not the meaning members like neighbors, et.c 
BiGraph::BiGraph(const BiGraph& other) {

	// cout<<"cp constructor called "<<endl;

	this->is_bipartite = other.is_bipartite;
	if(this->is_bipartite){
		this->dir = other.dir;

		this->init(other.num_v1, other.num_v2);

		// this->num_edges = 0;
		this->v1_max_degree = 0;
		this->v2_max_degree = 0;
	}else{

		// cout<<"copying a general graph"<<endl;

		this->dir = other.dir;

		this-> num_vertices = other.num_vertices ; 

		neighbor.resize(this-> num_vertices);
		degree.resize(this-> num_vertices);
		fill_n(degree.begin(), this-> num_vertices, 0);

		// this->num_edges = 0;
		this->v1_max_degree = 0;

	}
}

// add a switch to this to distinguish between bigraph and general graph
BiGraph::BiGraph(string dir, bool is_bigraph )
{
	num_v1 = 0;
	num_v2 = 0;
	num_edges = 0;
	v1_max_degree=0; 
	v2_max_degree=0;
	neighbor.clear();
	degree.clear();
	this->dir = dir;

	// do this for bipartite graphs: 
	if(is_bigraph){
		cout<<"processing bipartite graphs"<<endl;
		loadGraph(dir);
	}else{
		cout<<"processing general graphs"<<endl;
		this->is_bipartite =false;
		load_general_graph(dir);
	}
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
        return std::hash<T1>()(p.first) ^ std::hash<T2>()(p.second);
    }
};

// default constructor 
BiGraph::BiGraph() {
	dir = "";
	num_v1 = 0;
	num_v2 = 0;
	num_edges = 0;
	v1_max_degree=0; 
	v2_max_degree=0;
	neighbor.clear();
	degree.clear();
}

void BiGraph::print_graph()
{
	print_dash(50);
	cout<<"print_graph() called\n";
	for(int i=0;i<degree.size();i++){
		vector<vid_t> NB = neighbor[i];
		cout<<i<<": ";
		vector_show<vid_t>(NB); 
	}
}
void BiGraph::show_nodes()
{
	cout<<"show_nodes() called\n";
	// cout<<"upper nodes: "; 
	for(int i=0;i<num_v1+num_v2;i++){
		if(degree[i]>0){
			cout<<i<<" ";
		}
	}
	cout<<endl;
}
// initialize the BiGraph with size requirements  
void BiGraph::init(unsigned int num1, unsigned int num2)
{
	num_v1 = num1;
	num_v2 = num2;
	num_edges = 0;
	neighbor.resize(num_v1+num_v2);
	degree.resize(num_v1+num_v2);

	prio.resize(num_v1+num_v2);

	fill_n(degree.begin(), num_v1+num_v2, 0);
	// neighborHash.resize(num_v1+num_v2);
}

void BiGraph::computePriority(){
	std::vector<std::pair<int, int>> vertexDegrees(num_nodes());
	for (int i = 0; i < num_nodes(); i++) vertexDegrees[i] = std::make_pair(i, degree[i]);
	// Sort the vertexDegrees based on degrees and IDs
	std::sort(vertexDegrees.begin(), vertexDegrees.end(), [](const auto& a, const auto& b) {
		// the vertex with higher priority has lower degree.
		if (a.second != b.second) {
			return a.second > b.second; 
		} else {
			return a.first > b.first; 
		}
	});
	for(int i=0; i<num_nodes() ;i++)
	{
		neighbor[i].shrink_to_fit();
		sort(neighbor[i].begin(), neighbor[i].end());
		// cout<<"rank = "<<i<<",  id = "<<vertexDegrees[i].first<<",  deg = "<<vertexDegrees[i].second<<endl;
		prio[vertexDegrees[i].first] = i; 
	}
}

void BiGraph::loadGraph(string dir)
{

	// long double upper__ = 60 ;
	// long double lower__ = 50 ;
	unsigned int n1, n2;
	unsigned int edges = 0;
	int u, v, r;
	string metaFile = dir + ".meta";
	string edgeFile = dir + ".e";
	FILE * metaGraph = fopen(metaFile.c_str(), "r");
	FILE * edgeGraph = fopen(edgeFile.c_str(), "r");

	// bool include_amat = true ;
	// scan the meta file and read number of nodes 
	if (fscanf(metaGraph, "%d\n%d", &n1, &n2) != 2)
	{
		fprintf(stderr, "Bad file format: n1 n2 incorrect\n");
		exit(1);
	}
	fprintf(stdout, "n1: %d, n2: %d\n", n1, n2);
	init(n1, n2);

	// if(include_amat){a_mat = new map<int, int>[n1+n2+1];}

	while ((r = fscanf(edgeGraph, "%d %d", &u, &v)) != EOF)
	{
		//fprintf(stderr, "%d, %d\n", u, v);
		if (r != 2)
		{
			fprintf(stderr, "Bad file format: u v incorrect\n");
			exit(1);
		}
		addEdgeRaw(u,v);
		// neighbor[u].push_back(v+num_v1);
		// neighbor[v+num_v1].push_back(u);

		// if(include_amat){a_mat[u][v+num_v1] = a_mat[v+num_v1][u] = 1;}
	}
	// for(int u=0;u<num_nodes();u++){
	// 	degree[u] = neighbor[u].size();
	// 	if(u<num_v1){
	// 		v1_max_degree = v1_max_degree > degree[u] ? v1_max_degree : degree[u]; 
	// 	}else{
	// 		v2_max_degree = v2_max_degree > degree[u] ? v2_max_degree : degree[u]; 
	// 	}
	// 	num_edges+= degree[u] ;
	// }
	// num_edges/=2;
	
	cout<<"|E| = "<<num_edges<<endl;

	cout<<"fill_rate = "<<num_edges*1.0/(num_v1 * num_v2)<<endl;

	fclose(metaGraph);
	fclose(edgeGraph);

	computePriority();	
}

void BiGraph::load_general_graph(string dir)
{

	unsigned int n__, m__;

	int u, v, r;

	string edgeFile = dir; 

	FILE * edgeGraph = fopen(edgeFile.c_str(), "r");


	if (fscanf(edgeGraph, "%d\n%d", &n__, &m__) != 2)
	{
		fprintf(stderr, "Bad file format: n m incorrect\n");
		exit(1);
	}

	this-> num_vertices = n__ ; 
	this-> num_edges    = m__ ; 

	// cout<<"this-> num_vertices = "<< this-> num_vertices <<endl;
	// cout<<"this-> num_edges = "<< this-> num_edges <<endl;

	v1_max_degree = 0; // on general graphs, we use this to represent the maximum degree

	neighbor.resize(this-> num_vertices);
	degree.resize(this-> num_vertices);
	fill_n(degree.begin(), this-> num_vertices, 0);

	int edge_counter = 0;
	while ((r = fscanf(edgeGraph, "%d %d", &u, &v)) != EOF)
	{
		if (r != 2)
		{
			fprintf(stderr, "Bad file format: u v incorrect\n");
			exit(1);
		}
		// keep record of neighbors and degrees
		neighbor[u].push_back(v);
		neighbor[v].push_back(u);		
		// counting the visited edges
		edge_counter++;
	}

	// computing the maximum degree.
	for(int u=0;u<num_nodes();u++){
		degree[u] = neighbor[u].size();
		v1_max_degree = v1_max_degree > degree[u] ? v1_max_degree : degree[u]; 
	}

	// Sort the neighbor lists by vertex degrees
	for (int i = 0; i < num_nodes(); ++i) {
		std::sort(this->neighbor[i].begin(), this->neighbor[i].end(), [&](const int& a, const int& b) {
			return this->degree[a] > this->degree[b]; // Sort by degree in non-increasing order
		});
	}

	// verify the number of edges.
	assert(edge_counter == this-> num_edges ); 

	cout<<"|E| = "<<num_edges<<endl;

	fclose(edgeGraph);
}

bool satisfy_bound(long double upper, long double lower, int u, int x, int v, int w, BiGraph& g){
	// return (g.degree[u] <= upper) & (g.degree[x] <= upper) & (g.degree[v] <= lower) & (g.degree[w] <= lower); 

	return (g.degree[u] >= upper) & (g.degree[x] >= upper) & (g.degree[v] >= lower) & (g.degree[w] >= lower); 
	// all heavy vertices. 
}

// u,v are raw ids, convert them into real ids, and then update neighbor and degree
void BiGraph::addEdgeRaw(vid_t u, vid_t v)
{
	neighbor[u].push_back(v+num_v1);

	degree[u]++; 

	neighbor[v+num_v1].push_back(u);

	degree[v+num_v1]++;

	num_edges++;
	v1_max_degree = v1_max_degree > degree[u] ? v1_max_degree : degree[u]; 
	v2_max_degree = v2_max_degree > degree[v+num_v1] ? v2_max_degree : degree[v+num_v1]; 
}

// maybe this function takes too long? 
void BiGraph::addEdge(vid_t u, vid_t v)
{	
	// it looks like only these two lines are necessary 
	neighbor[u].push_back(v);
	neighbor[v].push_back(u);

	// the degrees, num_edges and max degrees can be computed on O(n) time. 
	// currently we are 

	degree[u]++; 
	degree[v]++;
	num_edges++;
	v1_max_degree = v1_max_degree > degree[u] ? v1_max_degree : degree[u]; 
	v2_max_degree = v2_max_degree > degree[v] ? v2_max_degree : degree[v]; 
}

// u1, u2 are real ids
void BiGraph::deleteEdge(vid_t u1, vid_t u2)
{
	vid_t upper, lower; 
	if(same_layer(u1,u2)){
		fprintf(stderr, "Bad input for delete edge: u v on the same layer\n");
		exit(1);		
	}
	for (int i = 0; i < degree[u1]; ++i)
	{
		int vv = neighbor[u1][i];
		if (vv == u2)
		{	
			swap(neighbor[u1][i], neighbor[u1][degree[u1] - 1]);
			--degree[u1];
			neighbor[u1].pop_back();
			num_edges--;
			break;
		}
	}
	for (int i = 0; i < degree[u2]; ++i)
	{
		int uu = neighbor[u2][i];
		if (uu == u1)
		{
			swap(neighbor[u2][i], neighbor[u2][degree[u2] - 1]);
			--degree[u2];
			neighbor[u2].pop_back();
			break;
		}
	}
}
// this function tests the validity of the deleteEdge function.

void computeCore(int ALPHA, int BETA, BiGraph& g)
{
	// make heap for calculations. this may not be necessary
	// printf("make heaps \n");
	LinearHeap upper_heap(g.num_v1, g.v1_max_degree);
	LinearHeap lower_heap(g.num_nodes(), g.v2_max_degree);

	// printf("lower cap = %d \n", lower_heap.key_cap);
	// printf("link \n");
	for(int i=0;i<g.num_nodes();i++){
		if(g.degree[i]==0) continue;
		if(g.is_upper(i)) {
			// cout<<"up "<< i <<" "<<g.degree[i] <<endl; 
			upper_heap.link(i, g.degree[i]);
		}
		if(g.is_lower(i)){
			// cout<<"lo "<< i <<" "<<g.degree[i] <<endl; 
			lower_heap.link(i, g.degree[i]);
		}
	}
	// printf("peeling \n");
	printf("upper_heap.get_min_key()  = %d  alpha = %d \n", upper_heap.get_min_key(), ALPHA);
	printf("lower_heap.get_min_key()  = %d  beta = %d \n", lower_heap.get_min_key(), BETA);
	// peeling
	while((upper_heap.get_min_key()<ALPHA)||(lower_heap.get_min_key()<BETA))
	{
		unordered_set<vid_t> affected_vertices;

		// printf("peeling upper \n"); 
		peeling(upper_heap,ALPHA,affected_vertices,g);
		// printf("peeling lower \n"); 
		peeling(lower_heap,BETA, affected_vertices,g);

		update_engagement(upper_heap, lower_heap, g, affected_vertices);
		if(g.num_edges==0){
			// cout<<"Resulting in empty graph.\n";
			break;
		}
	}
}
void peeling(LinearHeap& degree_heap, int ENG, unordered_set<vid_t>& affected_vertices, BiGraph& g)
{
	for(int key=degree_heap.min_key;key<ENG;key++){
		if(!degree_heap.has_rank(key)){
			// printf("\t no rank : %d \n", key);
			continue;
		}

		// printf("processing rank : %d \n", key);
		LinkNode *ptr = degree_heap.rank+key;
		ptr=ptr->next;
		while(ptr){
			vid_t u = ptr->id;
			// printf("deleting ptr->id = %d , deg = %d \n", u, g.degree[u]);
			vector<vid_t> NB = g.neighbor[u]; 
			for(vid_t v:NB){
				g.deleteEdge(u,v);
				affected_vertices.insert(v) ;
			}
			ptr=ptr->next;
			affected_vertices.insert(u) ;
		}
	}
	// if(g.num_edges<=724) cout<<" ******* iteration  *******"<<endl;
}
// make sure do not include anchors as affected vertices
void update_engagement(LinearHeap& upper_heap, LinearHeap& lower_heap, BiGraph& g, unordered_set<vid_t>& affected_nodes)
{
	for(auto i:affected_nodes)		
	{
		if(i<g.num_v1){
			if(g.degree[i]>0){
				upper_heap.update_key(i,g.degree[i]);
			}else{
				upper_heap.unlink(i); 
			}	
		}else
		{
			if(g.degree[i]>0){
				lower_heap.update_key(i,g.degree[i]);
			}else{
				lower_heap.unlink(i); 
			}	
		}
	}	
}
