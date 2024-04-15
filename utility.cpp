// experiment results:
#include "utility.h"
// #include "mt19937ar.h"

long double add_geometric_noise(long double data, long double sensitivity, long double epsilon) {
    long double p__ = 1 - std::exp(-epsilon / sensitivity);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::geometric_distribution<int> geometric_dist(p__);

    int noise = geometric_dist(gen) - geometric_dist(gen);
    long double noisy_data = data + noise;

    return noisy_data;
}
long double power(long double base, int exponent) {
    long double result = 1.0;

    if (exponent >= 0) {
        for (int i = 0; i < exponent; i++) {
            result *= base;
        }
    } else {
        for (int i = 0; i < -exponent; i++) {
            result /= base;
        }
    }
    return result;
}
double compute_epsilon(double exposure_risk, int m, int n1, int n2) {

	long double m_ = m, n1_ = n1, n2_ = n2; 

	long double num = exposure_risk * (m_ - n1_ * n2_) ; 

	long double demon = (m_ * (exposure_risk - 1)) ; 

	// cout<<"num = "<<num<<endl;

	// cout<<"demon = "<<demon<<endl;
	assert(num<0); 
	assert(demon<0);

    return log(num/demon );
}

double computeCovariance(const std::vector<double>& x, const std::vector<double>& y) {
    double covariance = 0.0;
    int size = x.size();

    // Calculate the means
    double meanX = 0.0;
    double meanY = 0.0;
    for (int i = 0; i < size; ++i) {
        meanX += x[i];
        meanY += y[i];
    }
    meanX /= size;
    meanY /= size;

    // Calculate the covariance
    for (int i = 0; i < size; ++i) {
        covariance += (x[i] - meanX) * (y[i] - meanY);
    }
    covariance /= size-1;

    return covariance;
}

double computeVariance(const std::vector<double>& numbers) {

    double mean = 0.0;
    double variance = 0.0;
    int size = numbers.size();

    // Calculate the mean
    for (const double& number : numbers)  mean += number;
    mean /= size;

    // Calculate the variance
    for (const double& number : numbers) {
        variance += (number - mean) * (number - mean);
    }
    variance /= size-1;
	
    return variance;
}

double computeStd(const std::vector<double>& numbers) {
    double mean = 0.0;
    double variance = 0.0;
    int size = numbers.size();

    // Calculate the mean
    for (const double& number : numbers)
        mean += number;
    mean /= size;

    // Calculate the variance
    for (const double& number : numbers)
        variance += (number - mean) * (number - mean);
    variance /= size;
    // Calculate the standard deviation
    double stdDeviation = std::sqrt(variance);
    return stdDeviation;
}

// Calculate L2 loss and MSE cost function
double computeMSE( long double btf, const std::vector<double>& numbers) {

	double mse = 0; 
    int size = numbers.size();
    // Calculate the mean
    for (const double& number : numbers){
		mse += (number - btf)*(number - btf);
	}
    mse /= size;
    return mse;
}

// Function to calculate the Manhattan distance between two vectors
double manhattanDistance(const std::vector<double>& vectorA, const std::vector<double>& vectorB) {
    if (vectorA.size() != vectorB.size()) {
        std::cerr << "Error: Vectors must have the same dimensions." << std::endl;
        return 0.0;
    }

    double distance = 0.0;
    double sumA = 0.0;
    for (size_t i = 0; i < vectorA.size(); ++i) {
        distance += std::abs(vectorA[i] - vectorB[i]);
        sumA += vectorA[i];
    }
    return distance / sumA;
}

double calculateCorrelation(const std::vector<int>& vectorX, const std::vector<int>& vectorY) {
	
    int n = vectorX.size();
    double meanX = calculateMean(vectorX);
    double meanY = calculateMean(vectorY);

    double numerator = 0.0;
    double sumXDiffSquared = 0.0;
    double sumYDiffSquared = 0.0;

    for (int i = 0; i < n; ++i) {
        double diffX = vectorX[i] - meanX;
        double diffY = vectorY[i] - meanY;
        numerator += diffX * diffY;
        sumXDiffSquared += diffX * diffX;
        sumYDiffSquared += diffY * diffY;
    }
    double denominator = sqrt(sumXDiffSquared * sumYDiffSquared);

    if (denominator == 0.0) {
        // Handle division by zero error or undefined correlation
        return 0.0;
    }

    double correlation = numerator / denominator;

    return correlation;
	
}

double n_choose_k(int n, int k) 
{
    if (k > n) {
        return 0; // Cannot choose more items than available
    }
    int num = 1, den = 1;
    
    // Compute numerator: n * (n-1) * ... * (n-k+1)
    for (int i = n; i >= n - k + 1; i--) {
        num *= i;
    }
    
    // Compute denominator: k!
    for (int i = 2; i <= k; i++) {
        den *= i;
    }
    
    double c = num / den;
    return c;
}


// Function to calculate the cosine similarity between two vectors
double cosineSimilarity(const std::vector<double>& vectorA, const std::vector<double>& vectorB) {
    if (vectorA.size() != vectorB.size()) {
        std::cerr << "Error: Vectors must have the same dimensions." << std::endl;
        return 0.0;
    }

    double dotProduct = 0.0;
    double magnitudeA = 0.0;
    double magnitudeB = 0.0;

    for (size_t i = 0; i < vectorA.size(); ++i) {
        dotProduct += vectorA[i] * vectorB[i];
        magnitudeA += vectorA[i] * vectorA[i];
        magnitudeB += vectorB[i] * vectorB[i];
    }

    magnitudeA = std::sqrt(magnitudeA);
    magnitudeB = std::sqrt(magnitudeB);

    if (magnitudeA == 0.0 || magnitudeB == 0.0) {
        std::cerr << "Error: One or both vectors have zero magnitude." << std::endl;
        return 0.0;
    }

    return dotProduct / (magnitudeA * magnitudeB);
}

// linked list functions
void List::insert_last(unordered_set<vid_t> data) 
{	
	// cout<<"inserting last\n";
	Node *temp=new Node;
	temp->vertices=data;

	// correctly assign id and in_ptrs
	int min=in_ptrs.size();
	for(auto v:data){
		min = min<v ? min:v;
		in_ptrs[v]=temp;
		num_nodes++;
	}
	temp->id=min;

	temp->prev=NULL;
	temp->next=NULL;
	if(head==NULL){ // the list is currently empty
		head=temp;
		tail=temp;
		temp=NULL;
	}
	else{	// not empty
		tail->next=temp;
		temp->prev=tail;
		tail=temp;
	}
	num_components++;
}
void List::display() // this function does not change
{
	Node *temp=head;
	if(head==NULL){
		cout<<"the list is empty\n";
		return;
	}
	while(temp!=NULL)
	{
		cout<<temp->id<<":";
		// for(auto xxx:temp->vertices){
		// 	cout<<xxx<<" ";
		// }
		cout<<temp->vertices.size();
		cout<<endl;
		temp=temp->next;
	}
} 
void List::delete_first()
{
	Node *temp=head;
	head=head->next;
	if(head!=NULL){
		head->prev=NULL;
	}else{
		tail=NULL;
	}
	for(auto v:temp->vertices){
		in_ptrs[v]=NULL;
	}
	delete temp;
	num_nodes--;
}
void List::delete_last()
{
	Node *temp=tail;
	tail=tail->prev;
	if(tail!=NULL){
		tail->next=NULL;
	}else{
		head=NULL;
	}
	for(auto v:temp->vertices){
		in_ptrs[v]=NULL;
	}
	delete temp;
	num_nodes--;
}
// it is assumed that ptr must be a node in the list
void List::delete_node(Node* ptr) 
{
	if((ptr==NULL)||(head==NULL)){
		return;
	}
	if(ptr==head){
		delete_first();
		return;
	}
	if(ptr==tail){
		delete_last();
		return;
	}
	// at this point ptr is in the middle
	Node* before = ptr->prev;
	Node* after = ptr->next;
	before->next=after;
	after->prev=before;
	for(auto v:ptr->vertices){
		in_ptrs[v]=NULL;
	}
	delete ptr;
	num_nodes--;
}

// linkedHeap functions
LinkNode::LinkNode() {
	prev = NULL;
	next = NULL;
	id = -1;
	val = -1;
}
// n is the number of nodes
LinearHeap::LinearHeap(int n, int key_cap) {
	this->max_key=0;
	this->active=0;
	this->min_key=key_cap;
	this->n = n;
	this->key_cap=key_cap;
	nods = new LinkNode[n];
	for( int i = 0; i < n; ++i) nods[i].id = i;
	rank = new LinkNode[key_cap+1];
}
// copy constructor 
LinearHeap::LinearHeap(const LinearHeap &obj){
	this->n = obj.n;
	this->key_cap=obj.key_cap;
	this->max_key=obj.max_key; 
	this->min_key=obj.min_key;
	this->nods = new LinkNode[n];
	this->rank = new LinkNode[key_cap+1];
	for( int i = 0; i < n; ++i){
		this->nods[i].id = i; 
		if(obj.nods[i].val!=-1){
			this->link(i, obj.nods[i].val);
		}
	}
}
bool LinearHeap::is_linked(int v) {return nods[v].prev != NULL;}
LinearHeap::~LinearHeap() {
	delete[] nods;
	delete[] rank;
}
// when unlink a node, reset its value
void LinearHeap::unlink(int v) {
	if(!is_linked(v)) return;
	LinkNode *nod = nods + v;
	if(nod->prev) nod->prev->next = nod->next;
	if(nod->next) nod->next->prev = nod->prev;
	nod->prev = NULL;
	nod->next = NULL;
	this->nods[v].val=-1;// important
	active--;
}
// v th nod with rank r
void LinearHeap::link(int v, int r) {
	if (r >key_cap){
		printf("ver = %d, rank = %d , cap = %d \n", v, r,  key_cap); 
		fprintf(stderr, "Invalid rank: rank exceeded kay cap\n"); exit(1);
	}	
	LinkNode *nod = nods + v;
	nod->val = r; // important
	max_key = max_key > r ? max_key : r; 
	min_key = min_key < r ? min_key : r;
	if(nod->prev) nod->prev->next = nod->next;
	if(nod->next) nod->next->prev = nod->prev;
	nod->prev = rank+r;
	nod->next = rank[r].next;
	if(rank[r].next) rank[r].next->prev = nod;
	rank[r].next = nod;
	active++;
}
int LinearHeap::get_rank(int v) {
	return nods[v].val;
}
bool LinearHeap::has_rank(int r) {
	if(r>key_cap){
		return false;cout<<"exceeding max rank\n";
	}
	return rank[r].next != NULL;
}
int LinearHeap::get_first_in_rank(int r) {
	return rank[r].next == NULL ? -1 : rank[r].next->id;
}
void LinearHeap::print_heap()
{
	cout<<"printing heap\n";
	for(int key=this->min_key; key<=this->max_key; key++){
		if(!this->has_rank(key)){
			continue;
		}
		LinkNode *ptr = this->rank+key;
		ptr=ptr->next;
		cout<<key<<":";
		while(ptr){
			cout<<ptr->id<<" ";
			ptr=ptr->next;
		}
		cout<<endl;
	}
}
int LinearHeap::get_min_key(){	
	// int old_key = this->min_key;
	for(int key=0;key<=this->key_cap; key++){
		if(has_rank(key)){
			min_key = key;
			break;
		}
	}
	return min_key;
}
int LinearHeap::get_max_key()
{
	// int old_key = this->max_key;
	for(int key=this->key_cap;key>=0; key--){
		if(has_rank(key)){
			max_key=key;
			break;
		}
	}
	// cout<<"the max key is "<<max_key<<endl;
	return max_key;
}

void LinearHeap::update_key(int v, int new_rank){
	unlink(v);
	link(v,new_rank);
}


std::ostream& operator<<(std::ostream& out, const Index_update value) {
	switch (value)
	{
	case Index_update::withlimit: return out << "withlimit";
	case Index_update::withlimit_base_opt: return out << "withlimit_base_opt";
	case Index_update::withlimit_dfs: return out << "withlimit_dfs";
	case Index_update::withlimit_dfs_parallel: return out << "withlimit_dfs_parallel";
	case Index_update::withlimit_parallel: return out << "withlimit_parallel";
	case Index_update::withoutlimit: return out << "withoutlimit";
	};
	return out << static_cast<std::uint16_t>(value);
}

void get_anchor_groups(std::vector<vid_t> candidates, int budget, std::vector<std::vector<vid_t>>& result)
{
	int K=budget,N=candidates.size();
	std::string bitmask(K, 1); 
	bitmask.resize(N, 0); 
	do {
		std::vector<vid_t> tmp;
		for (int i = 0; i < N; ++i) // [0..N-1] integers
		{
			if (bitmask[i]){
				tmp.push_back(candidates[i]);
			}
		}
		sort(tmp.begin(),tmp.end());
		result.push_back(tmp);
	} while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}

void zip(
    const std::vector<vid_t> &a, 
    const std::vector<int> &b, 
    std::vector<std::pair<vid_t,int>> &zipped)
{
    for(size_t i=0; i<a.size(); ++i)
    {
        zipped.push_back(std::make_pair(a[i], b[i]));
    }
}
void unzip(
    const std::vector<std::pair<vid_t, int>> &zipped, 
    std::vector<vid_t> &a, 
    std::vector<int> &b)
{
    for(size_t i=0; i<a.size(); i++)
    {
        a[i] = zipped[i].first;
        b[i] = zipped[i].second;
    }
}

vector<vid_t> vector_intersection(vector<vid_t> A, vector<vid_t> B)
{
	vector<vid_t> tmp; 
	set_intersection(A.begin(), A.end(), B.begin(), B.end(), back_inserter(tmp));
	return tmp;
}
vector<vid_t>  vector_diff(vector<vid_t> A, vector<vid_t> B){
	vector<vid_t> diff; 
	set_difference(A.begin(),A.end(),B.begin(),B.end(),inserter(diff, diff.begin()));
	return diff;
}
vector<vid_t>  vector_union(vector<vid_t> A, vector<vid_t> B){
	vector<vid_t> U_; 
	set_union(A.begin(),A.end(),B.begin(),B.end(),inserter(U_, U_.begin()));
	return U_;
}
// this function pass by value
vector<vid_t> sort_by_left(vector<vid_t> old_ids, vector<int>& UB)
{
	vector<pair<vid_t,int>> _zipped;
	vector<vid_t> old_ids_ = old_ids;
	zip(old_ids_, UB, _zipped);
	sort(begin(_zipped), end(_zipped), 
	[&](const auto& a, const auto& b)
	{	
		if(a.second == b.second){ return a.first > b.first;}
		return a.second > b.second;
	});
	unzip(_zipped, old_ids_, UB);	
	return old_ids_;
}

void print_dash(int n){
	for(int i=0;i<n;i++){cout<<"-";}
	cout<<endl;
}
