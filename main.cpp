#include "ldp-nb.h"
using namespace std;
using namespace std::chrono; 

int num_rounds;
vector<long double> estis; // we do not use this.
vector<long double> naive_estis, bs_estis, adv_estis;

vector<long double> bs_estis2;

extern bool one_round, count_cate, naive_switch ; 

long double p; // this is the flipping probability, corresponding to Eps1.

double gamma__ ;

long double Eps, Eps0, Eps1, Eps2;  
// private degrees: 
vector<int> priv_deg; 
int priv_dmax_1, priv_dmax_2; 

int iteration;

extern int alpha;

long double real = 0;

extern bool eva_comm ;


extern stats::rand_engine_t engine;  // Declare the engine as extern

// efficiency evaluation
long double RR_time, server_side_time, naive_server_side;

long double communication_cost=0;

int q_vertex = -1, x_vertex = -1; 

double d1, d2, epsilon; // these are global variables

int K2 = 1;

int main(int argc, char *argv[]) {

	// input parser:
	long double Eps =stold(argv[1]);

	string dataset = argv[2]; 

	num_rounds = atoi(argv[3]); // by default = 4

	int balancing_factor = atoi(argv[4]); // kappa, by default set to zero. 

	int K = atoi(argv[5]);// number of q vertices to sample 
	
	// initialize time
	RR_time=0, server_side_time=0, naive_server_side=0;

	std::mt19937 rng(std::random_device{}()); // for seeding

	bool is_bipartite  = true ;

	BiGraph g(dataset, is_bipartite);

	vector<long double> actuals, all_rmse,  allmae_bs, allmae2, allmae3; 

	// by default the balancing factor should be zero.
	cout<<"balancing_factor = "<<balancing_factor <<endl;

	// sample vertex pairs randomly
	vector<pair<vid_t, vid_t>> sampled_pairs = g.sample_bipartite_pairs(K, balancing_factor); 

    // Print the sampled pairs
    cout << "Sampled pairs of vertices from the same layer:" << endl;
    for (const auto& pair : sampled_pairs) {
		int uuu = pair.first, vvv = pair.second ; 
		vector<int> tmp; 
		set_intersection(
			g.neighbor[uuu].begin(), g.neighbor[uuu].end(), 
			g.neighbor[vvv].begin(), g.neighbor[vvv].end(), 
			back_inserter(tmp)
		);
		long double Nu_cap_Nv = tmp.size();
		actuals.push_back(Nu_cap_Nv);

    }



		// if we investigate bipartite graphs: 
		cout<<"Eps = "<<Eps<<endl;
		double t0 = omp_get_wtime();

		// unsigned int seed = rng(); 
		// cout<<"seed = "<<seed<<endl;

		bool naive_switch = false ;
		bool one_round_switch = false ;
		bool multi_ss_switch = false ;
		bool multi_ds_basic_switch = false ;
		bool multi_ds_switch = false ;
		bool multi_ds_im_switch = false ;

		long double alpha___ = 0;
		
		vector<long double> naive_time(num_rounds, 0), bs_time(num_rounds,0), adv_time(num_rounds,0);
		vector<long double> naive_com(num_rounds, 0), adv_com(num_rounds,0);
		server_side_time=0;

		for (int iteration = 0; iteration < num_rounds; iteration++) {
			cout<<"running default experiments"<<endl;
			multi_ds_basic_switch = false;
			switch(iteration) {
				case 0:
					naive_switch = true;
					one_round_switch = true;
					break;
				case 1:
					naive_switch = false ; 
					one_round_switch = false ;
					multi_ss_switch = true;
					break;
				// we do not need this version of MultiDS
				case 2:
					multi_ss_switch = false ; 
					break;
				case 3:
					multi_ds_im_switch = true;
					break;
				default:
					break;
			}
				
			naive_estis.clear();
			bs_estis.clear();
			adv_estis.clear();



			// for each sampled vertex pair
			for (int itr_sam = 0; itr_sam < sampled_pairs.size(); ++itr_sam) {
				const auto& pair = sampled_pairs[itr_sam];

				// each vertex deserves a unique seed 
				unsigned int seed = rng(); 
				cout<<"seed = "<<seed<<endl;

				cout << "(" << pair.first << ", " << pair.second << ")" << endl;

				int query_vertex = pair.first, x_vertex = pair.second ; 

				cout<<"deg 1 = "<<g.degree[query_vertex] <<endl;
				cout<<"deg 2 = "<<g.degree[x_vertex] <<endl;

				my_init_genrand(seed);

				int flipping_from, flipping_to; 
				if (g.is_lower(query_vertex)) {
					flipping_from = 0;
					flipping_to = g.num_v1 - 1;
				} else {
					flipping_from = g.num_v1;
					flipping_to = g.num_nodes() - 1;
				}
				double tt0 = omp_get_wtime();

				if(naive_switch || one_round_switch){
					if(naive_switch){
						cout<<"Running Naive for rounds"<<endl;
					}else{
						cout<<"Running Baseline for rounds"<<endl;
					}

					// in each round, compute a mae, 
					// compute the average of the mae and report

					vector<long double> naive_tmp, bs_tmp, reals;
					for(int xxx =0;xxx <K2;xxx++){
						reals.push_back(actuals[itr_sam]);
						// refresh the seed to get different noisy graph each time 
						my_init_genrand(rng());
						// 
						p = 1.0 / (exp(Eps) + 1.0);
						// cout<<"p = "<<p <<endl;
						BiGraph g2(g); // make a noisy graph
						// randomized responses from q.
						randomizedResponses(g, g2, query_vertex, flipping_from, flipping_to, p) ;
						randomizedResponses(g, g2, x_vertex, flipping_from, flipping_to, p) ;
						double tt1 = omp_get_wtime();
						// naive_time[iteration] += tt1-tt0; 
						// bs_time[iteration] += tt1-tt0; 
						// naive: server side:
						vector<int> Nq_cap_Nxx; 
						set_intersection(g2.neighbor[query_vertex].begin(), 
							g2.neighbor[query_vertex].end(), g2.neighbor[x_vertex].begin(), 
							g2.neighbor[x_vertex].end(), back_inserter(Nq_cap_Nxx));

						naive_tmp.push_back(Nq_cap_Nxx.size()); 
						// cout<<"naive estimate = "<<Nq_cap_Nxx.size() <<endl;
						double tt2 = omp_get_wtime();
						// naive_time[iteration] += tt2-tt1; 
						if(one_round_switch){
							// obtain the baseline estimate
							vector<int> Nq_union_Nxx; 
							set_union(g2.neighbor[query_vertex].begin(), 
								g2.neighbor[query_vertex].end(), 
								g2.neighbor[x_vertex].begin(), 
								g2.neighbor[x_vertex].end(), back_inserter(Nq_union_Nxx));
							long double bs_estimate = 0;
							long double param1 = std::pow((1 - p), 2) / std::pow((1 - 2 * p), 2);
							long double param2 = (1 - p) * (-p) / std::pow((1 - 2 * p), 2);
							long double param3 = std::pow(p, 2) / std::pow((1 - 2 * p), 2);
							bs_estimate += Nq_cap_Nxx.size() * param1; 
							bs_estimate += (Nq_union_Nxx.size() - Nq_cap_Nxx.size() ) * param2; 
							bs_estimate += (flipping_to - flipping_from +1 - Nq_union_Nxx.size()) * param3; 
							// bs_estis.push_back(bs_estimate);
							bs_tmp.push_back(bs_estimate);
							cout<<"baseline estimate = "<<bs_estimate <<endl;
							double tt3 = omp_get_wtime();
							// bs_time[iteration] += tt3-tt2; 	
						}
					}
					cout<<"naive estis : "<<naive_tmp[0] <<endl;
					cout<<"bs estis : "<<bs_tmp[0] <<endl;

					long double naive_mae = calculateMAE(naive_tmp, reals ); 
					long double bs_mae =    calculateMAE(bs_tmp, reals ); 
					naive_estis.push_back(naive_mae);
					bs_estis.push_back(bs_mae);
				}


				// need to compute adv_estis:
				double adv_time0 = omp_get_wtime();
				if(multi_ss_switch){
					// single source algorithm
					cout<<"MultiR-SS"<<endl;
					Eps0 = 0;
					Eps1 = Eps * 0.5;
				}
				// if running AVG basic algorithm, then run it for 8 iterations.
				if(multi_ds_basic_switch ){
					// AVG basic
					cout<<"MultiR-DS-Basic"<<endl;
					Eps0 = 0; // for double-source-basic, we do not need the first round.
					Eps1 = Eps * 0.5;
				}
				
				if(multi_ds_switch){
					// AVG smart
					cout<<"MultiR-DS"<<endl;
					Eps0 = Eps * 0.05; 
					long double deg1 = g.degree[query_vertex] + stats::rlaplace(0.0, 1/(Eps0), engine); 
					long double deg2 = g.degree[x_vertex] + stats::rlaplace(0.0, 1/(Eps0), engine); 
					//if(eva_comm) adv_com[iteration] += 2* sizeof(int);

					long double D = deg1 + deg1; 
					cout<<"d1 = "<<deg1 <<endl;
					cout<<"d2 = "<<deg1 <<endl;
					cout<<"D = "<<D<<endl;
					long double res_eps = Eps-Eps0;
					if(D>0){					
						// double solve_time1 = omp_get_wtime();
						// server_side_time+=solve_time1-solve_time0;
						Eps1 = solve_equation(D, res_eps); 
						std::cout << "Optimal Epsilon 1 computed: " << Eps1 << std::endl;
					}else{
						cout<<"regress to basic"<<endl;
						Eps1 = 0.5*res_eps ;
					}
				}
				// if( (!single_source_switch) && (!avg_double_source_smart) && smart_double_source){
				if(multi_ds_im_switch){
					cout<<"MultiR-DS-Im"<<endl;
					Eps0 = Eps * 0.05;
					// using the avg degree approach benefits the imbalanced cases 
					long double avgdeg = 0;
					long double deg1 =0, deg2 = 0;
					int startIdx = (g.is_upper(query_vertex)) ? 0 : g.num_v1;
					int endIdx = (g.is_upper(query_vertex)) ? g.num_v1 : g.num_nodes();
					for (int i = startIdx; i < endIdx; i++) {
						long double noisy_deg  = g.degree[i] + stats::rlaplace(0.0, 1 / Eps0, engine);
						if(i==query_vertex){
							deg1 = noisy_deg; 
						}
						if(i==x_vertex){
							deg2 = noisy_deg; 
						}						
						avgdeg += noisy_deg; 
					}
					avgdeg /= (g.is_upper(query_vertex)) ? g.num_v1 : g.num_v2;
					cout<<"avg = "<<avgdeg<<endl;
					long double res_eps = Eps-Eps0;
					cout<<"d1 = "<<deg1<<endl;
					cout<<"d2 = "<<deg2<<endl;
					if(deg1<=0){
						deg1 = avgdeg; 
					}
					if(deg2<=0){
						deg2 = avgdeg; 
					}
					std::vector<double> result = minimize(deg1, deg2, res_eps);
					Eps1 = result[0]; 
					alpha___ = result[1]; 
					std::cout << "argmin Eps1: " << Eps1 << ", argmin alpha: " << alpha___ << std::endl;
				}

				cout<<"Epsilon 1 = "<<Eps1 <<endl;

				p = 1.0 / (exp(Eps1) + 1.0); // computing the flipping probability

				BiGraph g2__(g); // make a different noisy graph

				// randomized responses from q and x
				randomizedResponses(g, g2__, query_vertex, flipping_from, flipping_to, p) ;
				randomizedResponses(g, g2__, x_vertex, flipping_from, flipping_to, p) ;
				// round 2: perform local graph analysis:
				Eps2 = Eps - Eps1 - Eps0; 
				cout<<"Epsilon 2 = "<<Eps2 <<endl;
				gamma__ = (1-p) / (1-2*p);
				long double adv_estimate=0; 
				// single source algorithms
				if(multi_ss_switch){// just use q.
					adv_estimate = locally_compute_f_given_q_and_x(query_vertex, x_vertex, g, g2__);
					//if(eva_comm) adv_com[iteration] += sizeof(long double); 
				}
				// double source algorithms
				if(multi_ds_im_switch){
					// need to solve for alpha___
					long double res_1 = locally_compute_f_given_q_and_x(query_vertex, x_vertex, g, g2__);
					long double res_2 = locally_compute_f_given_q_and_x(x_vertex, query_vertex, g, g2__);
					adv_estimate = alpha___ * res_1 + (1- alpha___) *res_2 ;
				}

				adv_estis.push_back(adv_estimate) ; 
				// cout<<"ADV estimate = "<<adv_estimate <<endl<<endl;
				double adv_time1 = omp_get_wtime();

				adv_time[iteration] += adv_time1- adv_time0; 
			}
			if(iteration ==0){

				long double mae1 = calculateMean(naive_estis ); 
				cout<<"naive MAE : "<<mae1 <<endl;

				long double mae2 = calculateMean(bs_estis ); 
				cout<<"oneR MAE : "<<mae2 <<endl;
			
			}else{
				long double mae3 = calculateMAE(adv_estis,   actuals ); 
				if(iteration ==1){
					cout<<"MultuiR single-source MAE: "<<mae3 <<endl;
				}
				if(iteration ==2){
					cout<<"MultuiR double-source MAE: "<<mae3 <<endl;
				}
				
			}

		}
		// efficiency evaluations:
		double t1 = omp_get_wtime();
		double seconds = t1 - t0;
		printf("time:%f\n", seconds);

	// }

	return 0;
}