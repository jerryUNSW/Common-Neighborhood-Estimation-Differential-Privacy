#include "ldp-nb.h"
#include "mt19937ar.h"

using namespace std;

vector<int> upper_sample, lower_sample;

double sample_ratio = 1.0; 

bool samling_one_round = false;

// long double verified = 0 , not_verified = 0;

extern long double Eps, Eps0, Eps1, Eps2, p; 

// m3__, m2__, m1__, m0__;
// priv deg related. 
extern vector<int> priv_deg; 
extern vector<long double> naive_estis;
extern int priv_dmax_1, priv_dmax_2; 

extern long double communication_cost; 

extern double gamma__ ;

extern int iteration;

extern double d1, d2, epsilon;

stats::rand_engine_t engine(std::time(0)); // used to be 1776

bool one_round = false, edge_clipping = true;

bool naive_switch;

bool sampling_noisy_graph = false ; 

double p____ = 0.5; // this is the sampling ratio.

// int alpha = 10; // this is used to cap two-round two-hop neighbors as well
int alpha = 3;

bool clipping = true; // this is our new switch 

bool eva_comm = false;

extern int q_vertex , x_vertex; 

extern long double RR_time, server_side_time, naive_server_side;

void private_estimate_of_degrees(BiGraph& g){
	// Eps0 = 0.05;

	// private estimate degrees. 
	priv_dmax_1 = 0 ;       
	priv_dmax_2 = 0 ; 
	priv_deg.resize(g.num_nodes());

	for(int i=0;i<g.num_nodes();i++){

		// priv_deg[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
		// need to re-run two-round algorithms with 1 --> 2 
		priv_deg[i] = g.degree[i]+stats::rlaplace(0.0, 2/(Eps0), engine); 
		
		if(edge_clipping){
			priv_deg[i]+=alpha;
		}
		// long double tmp = add_geometric_noise(g.degree[i], 1, Eps0);  
		// we found that geometric 
		// cout<<"deg = "<<g.degree[i]<<",\t"<<priv_deg[i]<<",\t"; cout<<tmp<<endl; 
		if(g.is_upper(i)){
			priv_dmax_1 = priv_dmax_1 > priv_deg[i] ? priv_dmax_1 : priv_deg[i] ; 
		}else{
			priv_dmax_2 = priv_dmax_2 > priv_deg[i] ? priv_dmax_2 : priv_deg[i] ; 
		}
		// if(choose_upper && g.is_upper(i))x1+=priv_deg[i];
		// if(!choose_upper && g.is_lower(i))x1+=priv_deg[i];
	}

	// cout<<"esitmated  = "<<x1<<endl;
	// g.num_edges = x1;
}

double my_genrand_real2(){return genrand_real2(); }

void my_init_genrand(unsigned long seed){init_genrand(seed);}

// Function to uniformly sample K vertices from a range
std::vector<int> uniformSampleKVertices(int K, bool useFirstRange, BiGraph& g) {
    int range_from, range_to;
    if (useFirstRange) {
        range_from = 0;
        range_to = g.num_v1 - 1;
    } else {
        range_from = g.num_v1;
        range_to = g.num_nodes() - 1;
    }

    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<int> distr(range_from, range_to);

    std::vector<int> sampledVertices;
    for (int i = 0; i < K; ++i) {
        int vertex = distr(generator);
        sampledVertices.push_back(vertex);
    }

    return sampledVertices;
}


// Define the function
void randomizedResponses(BiGraph& g, BiGraph& g2, int query_vertex, int flipping_from, int flipping_to, double p) {
    // Loop through vertices within the specified range
    // for(int vvv = flipping_from; vvv < flipping_to; vvv++) {
	for(int vvv = flipping_from; vvv <= flipping_to; vvv++) {
        // Check if vvv is a neighbor of query_vertex
        if(std::find(g.neighbor[query_vertex].begin(), g.neighbor[query_vertex].end(), vvv) != g.neighbor[query_vertex].end()) {
            // Case: vvv is a neighbor
            if(my_genrand_real2() >= p) { // Probability check
                g2.addEdge(query_vertex, vvv); // Add edge
            }
        } else {
            // Case: vvv is not a neighbor
            if(my_genrand_real2() < p) { // Probability check
                g2.addEdge(query_vertex, vvv); // Add edge
            }
        }
    }
}

void construct_noisy_graph(BiGraph& g, BiGraph& g2, unsigned long seed){

	if(g.is_bipartite){
		// cout<<"is bigraph"<<endl;
		const int range_from  = g2.num_v1;
		const int range_to    = g2.num_nodes()-1;
		std::random_device                  rand_dev;
		std::mt19937                        generator(rand_dev());
		std::uniform_int_distribution<int>  distr(range_from, range_to);

		init_genrand(seed);
		int flip1 = 0; 
		int visited_vertices = 0;
		int total_vertices = g2.num_v1;
		int ten_percent = total_vertices / 5;
		long double max_time_per_user = -1;

		for(int i=0;i<g2.num_v1;i++){
			visited_vertices++; 
			// if (visited_vertices % ten_percent == 0) {
			// 	int progress = visited_vertices * 100 / total_vertices;
			// 	cout << "Processed " << progress << "% of vertices" << endl;
			// }
			double tx = omp_get_wtime();
			for(int j=g2.num_v1;j<g2.num_nodes();j++){
				if(std::find(g.neighbor[i].begin(), g.neighbor[i].end(), j) != g.neighbor[i].end() ){
					if(genrand_real2() >= p ){ // 1  --> 1 
						g2.addEdge(i,j);
						flip1++;
					}	
				}else{
					if(genrand_real2() < p){   // 0 --> 1
						g2.addEdge(i,j);
						flip1++;	
					}
				}
			}
			double ty = omp_get_wtime();
			max_time_per_user = max_time_per_user > (ty-tx) ? max_time_per_user : (ty-tx);
		}
		cout<<"noisy edges = "<<flip1<<endl;

		communication_cost+= flip1*sizeof(int); 

		long double expected_E = g.num_edges*(1-p) +(g.num_v1*g.num_v2-g.num_edges)*p; 

		cout<<"expected E = "<< expected_E<<endl;

		g2.computePriority();

	}else{
		cout<<"constructing a noisy graph for a general graph"<<endl;
		// constructing a noisy graph for a general graph
		const int range_from  = 0;
		const int range_to    = g2.num_nodes()-1;
		std::random_device                  rand_dev;
		std::mt19937                        generator(rand_dev());
		std::uniform_int_distribution<int>  distr(range_from, range_to);
		init_genrand(seed);
		//
		int flip1 = 0; 
		int visited_vertices = 0;
		int total_vertices = g2.num_nodes();
		int ten_percent = total_vertices / 2;
		// long double max_time_per_user = -1;
		for(int i=0;i<g2.num_nodes();i++){
			visited_vertices++; 
			if (visited_vertices % ten_percent == 0) {
				int progress = visited_vertices * 100 / total_vertices;
				cout << "Processed " << progress << "% of vertices" << endl;
			}
			// double tx = omp_get_wtime();

			// for(int j=g2.num_v1;j<g2.num_nodes();j++){
			for(int j=0; j<i;j++){
				// the lower left triangle
				if(std::find(g.neighbor[i].begin(), g.neighbor[i].end(), j) != g.neighbor[i].end() ){
					// if edge (i,j) ==1
					// cout<<"edge "<<i<<" "<<j<<" = 1"<<endl;
					if(genrand_real2() >= p ){ // 1  --> 1, keep with probability 1-p
						g2.addEdge(i,j);
						flip1++;
					}	
				}else{
					// if edge (i,j) ==0
					// cout<<"edge "<<i<<" "<<j<<" = 0"<<endl;
					if(genrand_real2() < p){   // 0 --> 1, flip with probability p.
						g2.addEdge(i,j);
						flip1++;	
					}
				}
			}
			double ty = omp_get_wtime();
			// max_time_per_user = max_time_per_user > (ty-tx) ? max_time_per_user : (ty-tx);
		}
		long double possible_edges = g.num_vertices*(g.num_vertices-1)/2; 

		cout<<"# possible edges = "<<possible_edges <<endl;

		long double expected_E = g.num_edges*(1-p) +(possible_edges -g.num_edges) *p; 

		cout<<"expected E = "<< expected_E<<endl;

		cout<<"# noisy edges = "<<flip1<<endl;

	}
}


// is there a path s -- v --- t in G ? 
// these functions are returning two-hop reachable vertices. 
// two hop neighbors = immediate neighbors + two-hop reachable vertices.
// works for all graphs
// this function takes constant time complexity.
long double check_two_hop_path(int s, int v, int t, BiGraph& g2) {

	// cout<<"checking "<<s<<" "<<v<<" "<<t<<endl;

    long double result1, result2; 
    result1 = g2.has(s, v) ? 1 : 0 ; 
    result2 = g2.has(v, t) ? 1 : 0 ; 

	// do we need to minus p from result1 and result2 ? 
	result1-= p;
	result2-= p;
	// these two lines are very important

    result1 /= (1 - 2 * p);
    result2 /= (1 - 2 * p); 
    return result1 * result2;
}

// does there exist a path of length two between s and t in G? // this works
long double check_distance_two(int s, int t, BiGraph& g2) {
    long double result = 1.0; // Initialize result to 1.0

	int start, end; 

	if(g2.is_bipartite){
		// Determine the range of vertices to iterate based on the position of s
		start = (g2.is_upper(s)) ? g2.num_v1 : 0;
		end = (g2.is_upper(s)) ? g2.num_nodes() : g2.num_v1;
	}else{
		start = 0;
		end = g2.num_nodes();
	}

	// Iterate over the vertices to compute the result
	// this is very slow, maybe we could parallelize this part
	#pragma omp parallel for reduction(*:result)
	for (int xxx = start; xxx < end; xxx++) {

		// Ensure that the vertex xxx is in the correct partition
		if(g2.is_bipartite){
			assert((g2.is_upper(s) && !g2.is_upper(xxx)) || (!g2.is_upper(s) && g2.is_upper(xxx)));
		}else{
			if(xxx==s) 
				continue;
			if(xxx==t) 
				continue;
		}
		// Calculate the probability of a two-hop path
		long double tmp = check_two_hop_path(s, xxx, t, g2);

		// Update the result using the probability
		result *= (1 - tmp);
	}
	return 1.0 - result;
}


double compute_f_given_q_x(int q, int x, BiGraph& g2) 
{
	int start, end;
	if(g2.is_bipartite){

		start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
		end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1 ;

		// cout<<"start = "<<start <<endl;
		// cout<<"end = "<<end <<endl;
	}else{
		start = 0;
		end = g2.num_nodes();
	}

	double res = 0;
	// Iterate over the vertices to compute the result
	// cout<<"Enumerating all 2 hop paths..."<<endl;

	// this can be made a lot faster.
	for (int vvv = start; vvv < end; vvv++) {
		if(vvv == q) continue;
		if(vvv == x) continue;
		// accumulate res[xxx]
		res += check_two_hop_path( q, vvv, x,  g2) ; 
	}

	


	return res;
}

// given a vertex q, we iterature all vertices that could be its two hop neighbors and 
vector<double> compute_twoHop_given_q(int q, BiGraph& g2) {

	vector<double> res(g2.num_nodes() , -1); 

	int start, end;
	if(g2.is_bipartite){
		start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
		end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;
	}else{
		start = 0;
		end = g2.num_nodes();
	}

	// Iterate over the vertices to compute the result
	// cout<<"Enumerating all 2 hop paths..."<<endl;
	// for each two hop candidate: O(N) 
	// this loop cannot be sampled because each vertex needs a score.
	for (int xxx = start; xxx < end; xxx++) {

		// do not consider itself as a two hop neighbor
		if(xxx == q) 
			continue; 


		if(xxx==x_vertex){
			// do somnthing here. 

			for (int vvv = start; vvv < end; vvv++) {
				if(vvv == q) 
					continue;
				if(vvv == xxx) 
					continue;
				// accumulate res[xxx]
				res[xxx] += check_two_hop_path( q, vvv, xxx,  g2) ; 
			}

		}else{
			continue;
		}

		// res[xxx] = check_distance_two(q, xxx, g2); // this takes O(N)

		// cout<<"ver: "<<xxx<<" ";
		// cout<<res[xxx] <<endl;
	}
	// the overall complexity: O(N^2). 
	// how can we reduce this complexity? 

	return res;
}


bool selectively_visit_neighbors = true; 


long double locally_check_distance_two(int s, int t, BiGraph& g, BiGraph& g2) {
	// cout<<"s = "<<s <<endl;
	// cout<<"t = "<<t <<endl; 

    long double result = 1.0; // Initialize result to 1.0

	// only visit the middle vertices that are neighbors of s in G
	// pruning: onsly visit top K vertices in here.
	int cnt = 0; 
	for(auto xxx : g.neighbor[s]){
		if(clipping && cnt == alpha){
			break;
		}
		// Calculate the probability of a two-hop path
		long double tmp = g2.has(xxx, t) ? 1 : 0 ;   
		tmp -= p; 
		tmp /= (1 - 2 * p);
		// Update the result using the probability
		result *= (1 - tmp);
		cnt++;
	}
	return 1.0 - result;
}

// given a vertex q, we iterature all vertices that could be its two hop neighbors and 
vector<double> locally_compute_twoHop_given_q(int q, BiGraph& g, BiGraph& g2) {

	vector<double> res; 
	res.resize(g2.num_nodes()); 
	fill(res.begin(), res.end(), -1); 

	int start, end;
	if(g2.is_bipartite){
		start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
		end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;
	}else{
		start = 0;
		end = g2.num_nodes();
	}
	// int start = ( g2.is_upper(q) ) ?   0 : g2.num_v1 ;
	// int end =   ( g2.is_upper(q) ) ?   g2.num_v1 : g2.num_nodes();

	for (int xxx = start; xxx < end; xxx++) {

		if(xxx == q) 
			continue; // we compute res for all vertices that are not q.

		// we assume the middle vertex is the neighbor of q. 
		res[xxx] = locally_check_distance_two(q, xxx, g, g2);

		long double global_sensitivity; 

		if(clipping){
			global_sensitivity = pow(gamma__, 1.0*(alpha+1));
		}else{
			global_sensitivity = pow(gamma__, 1.0*(g.degree[q]+1));
		}

		// Add Laplacian noise here:

		long double noise = stats::rlaplace(0.0, (global_sensitivity/Eps2), engine); 

		// cout<<"deg = "<<g.degree[q] <<endl;
		// cout<<"gs = "<<global_sensitivity<<" ";
		// cout<<" noise = "<<noise <<endl;

		res[xxx] += noise; 
		// this noise is too much because the global sensitivity is exponential
		// idea: only visit top-K vertices in the neighbors of q.
		// these vertices contribute more 
	}

	return res;
}

double locally_compute_f_given_q_and_x(int q, int x, BiGraph& g, BiGraph& g2) {

	double res =-1;

	int start, end;
	if(g2.is_bipartite){
		start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
		end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;
	}else{
		start = 0;
		end = g2.num_nodes();
	}

	double Nx_cap_Nq_minus_x =0, Nq_minus_Nx_minus_x =0;


	// it looks like edge clipping makes it worse
	bool x_is_a_nb_of_q = false;

	for(auto nb: g.neighbor[q]){
		if(nb == x){
			x_is_a_nb_of_q = true;
		}

		if(g2.has(nb, x)){
			Nx_cap_Nq_minus_x++;
		}else{
			Nq_minus_Nx_minus_x++;
		}
	}

	// if x is not a neighbor of q, then 
	if(g2.is_bipartite){
		// in the case of bipartite graphs
		if(x_is_a_nb_of_q){
			cout<<"x cannot be a neighbor of q. They are from the same layer."<<endl;
			exit(1);
		}
	}else{
		// in the case of general graphs 
		if(x_is_a_nb_of_q){
			assert(Nq_minus_Nx_minus_x == g.degree[q] - Nx_cap_Nq_minus_x-1);
		}else{
			assert(Nq_minus_Nx_minus_x == g.degree[q] - Nx_cap_Nq_minus_x); 
		}
	}

	// locally the degree of q is known. 
	res = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 

	long double noise = stats::rlaplace(0.0, (gamma__/Eps2), engine); 

	res += noise; 

	return res;
}


// given a vertex q, we iterature all vertices that could be its two hop neighbors and 
vector<double> locally_compute_num_two_paths_given_q(int q, BiGraph& g, BiGraph& g2) {

	vector<double> res; 
	res.resize(g2.num_nodes()); 
	fill(res.begin(), res.end(), -1); 

	int start, end;
	if(g2.is_bipartite){
		start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
		end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;
	}else{
		start = 0;
		end = g2.num_nodes();
	}

	for (int xxx = start; xxx < end; xxx++) {

		if(q == q_vertex  && xxx != x_vertex){
			continue;
		}
		if(q ==x_vertex && xxx != q_vertex){
			continue;
		}

		if(xxx == q) 
			continue; // we compute res for all vertices that are not q.

		double Nx_cap_Nq_minus_x =0;

		double Nq_minus_Nx_minus_x =0;


		// it looks like edge clipping makes it worse
		int cnt = 0; 
		bool x_is_a_nb_of_q = false;
		for(auto nb: g.neighbor[q]){
			
			if(nb == xxx){
				x_is_a_nb_of_q = true;
				continue;
			}
			// if(clipping && cnt == alpha){
			// 	break;
			// }

			if(g2.has(nb, xxx)){
				Nx_cap_Nq_minus_x++;
			}else{
				Nq_minus_Nx_minus_x++;
			}
			cnt++;
		}

		// if x is not a neighbor of q, then 
		if(x_is_a_nb_of_q){
			assert(Nq_minus_Nx_minus_x == g.degree[q] - Nx_cap_Nq_minus_x-1);
		}else{
			assert(Nq_minus_Nx_minus_x == g.degree[q] - Nx_cap_Nq_minus_x); 
		}

		// it looks like the laplacian noise is not even the issue here.

		// we assume the middle vertex is the neighbor of q. 
		// locally the degree of q is known. 

		if(clipping){
			res[xxx] = Nx_cap_Nq_minus_x * gamma__ + (alpha - Nx_cap_Nq_minus_x)*(-p)/(1-2*p); 
		}else{
			res[xxx] = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x*(-p)/(1-2*p); 
		}
		// if we rank vertices based on this score, it is effectively ranking them by the number of noisy edges

		if(q == q_vertex && xxx == x_vertex){
			cout<<"q = "<<q <<endl;
			cout<<"Nx_cap_Nq_minus_x = "<<Nx_cap_Nq_minus_x <<endl;
			cout<<"Nq_minus_Nx_minus_x = " <<Nq_minus_Nx_minus_x <<endl;
		}

		if(q == x_vertex && xxx == q_vertex){
			cout<<"q = "<<q <<endl;
			cout<<"Nx_cap_Nq_minus_x = "<<Nx_cap_Nq_minus_x <<endl;
			cout<<"Nq_minus_Nx_minus_x = " <<Nq_minus_Nx_minus_x <<endl;
		}


		long double noise = stats::rlaplace(0.0, (gamma__/Eps2), engine); 

		res[xxx] += noise; 

		// even withou the noise, estimate of f is not very accurate
	}

	return res;
}


long double get_laplace(long double parameter){
	return stats::rlaplace(0.0, parameter, engine); 
}


std::vector<std::vector<int>> generateCombinations(const std::vector<int>& up_vs, int p) {
    std::vector<std::vector<int>> combinations;
    std::vector<bool> selected(up_vs.size());

    std::fill(selected.begin(), selected.begin() + p, true);

    do {
        std::vector<int> combination;
        for (size_t i = 0; i < up_vs.size(); i++) {
            if (selected[i]) {
                combination.push_back(up_vs[i]);
            }
        }
        combinations.push_back(combination);
    } while (std::prev_permutation(selected.begin(), selected.end()));

    return combinations;
}

// this uses Newton-Ralph's method to find the global minimum.
long double solve_equation(long double d, long double epsilon) {
	double target = 0;
	auto target_func = [d, epsilon, &target](double x) {
		// f0 is the equation to be solved
		double f0 = d*pow(epsilon, 3)*exp(x) + d*pow(epsilon, 3) - 3*d*pow(epsilon, 2)*x*exp(x) - 
			3*d*pow(epsilon, 2)*x + 3*d*epsilon*pow(x, 2)*exp(x) + 3*d*epsilon*pow(x, 2) - 
			d*pow(x, 3)*exp(x) - d*pow(x, 3) + 8*epsilon*exp(x) - 8*x*exp(x) - 8*exp(2*x) 
			+ 8*exp(x);
		// f1 is the derivative of f0
		double f1 = d*pow(epsilon, 3)*exp(x) - 3*d*pow(epsilon, 2)*x*exp(x) - 3*d*pow(epsilon, 2)*exp(x) - 
			3*d*pow(epsilon, 2) + 3*d*epsilon*pow(x, 2)*exp(x) + 6*d*epsilon*x*exp(x) + 6*d*epsilon*x - 
			d*pow(x, 3)*exp(x) - 3*d*pow(x, 2)*exp(x) - 3*d*pow(x, 2) + 8*epsilon*exp(x) - 8*x*exp(x) 
			- 16*exp(2*x);
		return std::pair<double, double>(f0 - target, f1);
	};
	double guess = 0.5*epsilon;
	double min = 0.01 * epsilon;
	double max = 0.99 * epsilon;
	int sf = 10; // number of significant figures

	double root = boost::math::tools::newton_raphson_iterate(target_func, guess, min, max, sf);
	// double true_value = std::atanh(target);
    return root;
}
// Define function fff(x, alpha)
double fff(double x, double alpha) {
    double p = 1 / (1 + exp(x));
    double s1 = p * (1 - p) / pow(1 - 2 * p, 2);
    double s2 = pow(alpha, 2) * d1 + pow(1 - alpha, 2) * d2;
    double s3 = (2 * pow(1 - p, 2)) / pow(1 - 2 * p, 2) / pow(epsilon - x, 2);
    double s4 = pow(alpha, 2) + pow(1 - alpha, 2);
    return s1 * s2 + s3 * s4;
}

// Define objective function to minimize
double objective(const gsl_vector *v, void *params) {
    double x = gsl_vector_get(v, 0);
    double alpha = gsl_vector_get(v, 1);

    // Check constraints
    if (x <= 0 || x >= epsilon || alpha < 0 || alpha > 1) {
        // Penalize the objective function if constraints are violated
        return std::numeric_limits<double>::max();
    }

    return fff(x, alpha);
}
std::vector<double> minimize(double d1__, double d2__, double epsilon__) {
    const size_t dim = 2;

	// pass these to global variable
	d1 = d1__; 
	d2 = d2__; 
	epsilon = epsilon__;

    gsl_vector *x = gsl_vector_alloc(dim);
    gsl_vector_set(x, 0, 0.5 * epsilon); // Set x
    gsl_vector_set(x, 1, 0.5);           // Set alpha

    gsl_multimin_function obj_func;
    obj_func.n = 2; // Number of dimensions
    obj_func.f = &objective; // Objective function
    obj_func.params = nullptr; // No additional parameters needed

    gsl_multimin_fminimizer *minimizer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, dim);

    gsl_vector *step_size = gsl_vector_alloc(dim);
    gsl_vector_set_all(step_size, 0.01); // Set the step size for all dimensions

    gsl_multimin_fminimizer_set(minimizer, &obj_func, x, step_size);

    int status;
    size_t iter = 0;
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(minimizer);
        if (status) break;
        status = gsl_multimin_test_size(minimizer->size, 1e-6);
    } while (status == GSL_CONTINUE && iter < 1000);

    std::vector<double> result(2);
    if (status == GSL_SUCCESS) {
        result[0] = gsl_vector_get(minimizer->x, 0); // argmin x
        result[1] = gsl_vector_get(minimizer->x, 1); // argmin alpha
    } else {
		cout<<"Optimization failed!" << std::endl;
    }
	
	assert(result[0]>0 && result[0] < epsilon);
	assert(result[1]>=0 && result[1] <= 1);

    gsl_multimin_fminimizer_free(minimizer);
    gsl_vector_free(x);
    gsl_vector_free(step_size);

    return result;
}