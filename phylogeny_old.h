#ifndef phylogenyOld_H
#define phylogenyOld_H
//#include <iomanip>
#include <algorithm>
#include <random> 
#include <iostream>
#include "rootedPhylogeny.h"
class phylogeny_vertex{
public:
	string name;
	string newickLabel = "";
	int timesVisited = 0;
	double sumOfEdgeLogLikelihoods = 0;
	double logLikelihoodForObservingSequence = 0;
	int degree = 0;
	int id;
	vector <phylogeny_vertex*> neighbors;
	phylogeny_vertex(int idToAdd, string nameToAdd){
		id = idToAdd;
		name = nameToAdd;
	}
	~phylogeny_vertex(){		
	}
};

class phylogeny_tree{
private:
	pair <phylogeny_vertex*,phylogeny_vertex*> edgeForRooting;	
	string DNA_alphabet = "ACGT";
public:
	int sequenceLength;
	double logLikelihoodForRooting = -1000000000000;
	double maxLogLikelihood = logLikelihoodForRooting - 10;
	default_random_engine generator;
	phylogeny_vertex * vertexForRooting;
	map <int,phylogeny_vertex*>* vertexMap;
	vector <pair<int,int>> * edgeList;
	map <pair<int,int>,double>* edgeLogLikelihoodMap; // bidirected edges
//	map <pair<int,int>,vector<unsigned char>>* rootSequencePerEdge; // undirected edges
//	map <pair<int,int>,double>* dyadLogLikelihoodMap; // undirected edges
//	map <pair<int,int>,pair<double,double>>* dyadSeperateLogLikelihoodMap; // undirected edges
//	map <pair<int,int>,array<array<float,4>,4>>* transitionMatricesForEdgesMap; // bidirected edges
//	map <pair<int,int>,pair<array<array<float,4>,4>,array<array<float,4>,4>>>* transitionMatricesForDyadMap; // undirected edges
//	map <pair<int,int>,array<float,4>>* rootProbabilityForDyadMap; // undirected edges
//	map <pair<int,int>,float>* edgeLengthsMap; // undirected edges
//	map <pair<int,int>,pair<float,float>>* edgeLengthsForDyadMap; // for each undirected edge {u,v} compute length of edges (r,u) and (r,v), 
	vector <pair<int,int>>* directedEdgeList;
//	void ComputeDirectedEdgeList();
	void AddEdge(int u_id, vector <unsigned char> seq_u, int v_id, vector <unsigned char> seq_v);
//	void AddEdge(int u_id, vector <unsigned char> seq_u, int v_id, vector <unsigned char> seq_v, vector <unsigned char> seq_r);
	void AddVertex(int v_id, string v_name, vector <unsigned char> seq_v);
//	void AddVertex(int v_id, string v_name);
	bool ContainsVertex(int v_id);
//	void ComputeDyadLogLikelihood(int u_id, vector<unsigned char> seq_u, int v_id, vector<unsigned char> seq_v);
	void ComputeEdgeLengthAndBidirectionalTransitionMatricesAndBidirectionalLogLikelihoods(int u_id, vector<unsigned char> seq_u, int v_id, vector<unsigned char> seq_v);
//	void ComputeDyadEdgeLengthsAndDyadTransitionMatricesAndDyadRootProbabilityAndDyadLogLikelihoods(int u_id, vector<unsigned char> seq_u, int v_id, vector<unsigned char> seq_v);
//	void SelectEdgeForRooting();
	void SelectVertexForRooting();
	void RootTreeAtEdge(int u_id, int v_id);
	void RootTreeAtVertex(int v_id);
//	void WriteUnrootedTreeToFile(string sequenceFileName);
//	void WriteOutputFiles(string sequenceFileName);
	float GetEdgeLength(int v_id, int n_id);
	
	phylogeny_tree(){
		this->vertexMap = new map<int,phylogeny_vertex*>;
//		rootSequencePerEdge = new map <pair<int,int>,vector<unsigned char>>;
		this->edgeLogLikelihoodMap = new map <pair<int,int>, double>;
		this->edgeList = new vector <pair<int,int>>;
//		dyadLogLikelihoodMap = new map <pair<int,int>, double>;
//		dyadSeperateLogLikelihoodMap = new map <pair<int,int>, pair<double,double>>;
//		edgeLengthsMap = new map <pair<int,int>, float>;
//		edgeLengthsForDyadMap = new map <pair<int,int>, pair<float,float>>;
//		transitionMatricesForDyadMap = new map <pair<int,int>,pair<array<array<float,4>,4>,array<array<float,4>,4>>>;
//		transitionMatricesForEdgesMap = new map <pair<int,int>,array<array<float,4>,4>>;
//		rootProbabilityForDyadMap = new map <pair<int,int>,array<float,4>>;
		this->directedEdgeList = new vector <pair<int,int>>;
	}
	~phylogeny_tree(){
		this->edgeLogLikelihoodMap->clear();
//		dyadLogLikelihoodMap->clear();
//		dyadSeperateLogLikelihoodMap->clear();
//		edgeLengthsMap->clear();
//		edgeLengthsForDyadMap->clear();
		for (pair<int,phylogeny_vertex*> vIdAndVertexPtr: *this->vertexMap){
			delete vIdAndVertexPtr.second;
		}
		this->vertexMap->clear();
//		transitionMatricesForDyadMap->clear();
//		transitionMatricesForEdgesMap->clear();
		delete this->vertexMap;
		delete this->edgeList;
		delete this->edgeLogLikelihoodMap;
//		delete dyadLogLikelihoodMap;
//		delete dyadSeperateLogLikelihoodMap;
//		delete edgeLengthsMap;
//		delete edgeLengthsForDyadMap;
//		delete transitionMatricesForDyadMap;
//		delete transitionMatricesForEdgesMap;
//		delete rootProbabilityForDyadMap;
//		delete rootSequencePerEdge;
		delete this->directedEdgeList;
	}
};

//float phylogeny_tree::GetEdgeLength(int v_id, int n_id){
//	if (v_id < n_id){
//		return ((*this->edgeLengthsMap)[pair<int,int>(v_id,n_id)]);
//	} else {
//		return ((*this->edgeLengthsMap)[pair<int,int>(n_id,v_id)]);
//	}
//}

//void phylogeny_tree::WriteUnrootedTreeToFile(string sequenceFileName){
//	ofstream phylogenyFile;
//	phylogenyFile.open(sequenceFileName + ".unrootedPhylogeny_" + to_string((*this->vertexMap).size()));
//	phylogeny_vertex* v;
//	for (pair<int,phylogeny_vertex*> vIdPtr : (*this->vertexMap)){
//		v = vIdPtr.second;
//		for (phylogeny_vertex* n : v->neighbors){
//			if (v->id < n->id){
//				phylogenyFile << v->name << "\t" << n->name << "\t" << this->GetEdgeLength(v->id, n->id) << endl;
//			}
//		}
//	}
//	phylogenyFile.close();
//}
	

void phylogeny_tree::ComputeEdgeLengthAndBidirectionalTransitionMatricesAndBidirectionalLogLikelihoods(int u_id, vector<unsigned char> seq_u, int v_id, vector<unsigned char> seq_v){
//	float edgeLength = 0;
	int sequenceLength = int(seq_u.size());
//	for (int sitePos = 0; sitePos < sequenceLength; sitePos++){
//		if (seq_u[sitePos] != seq_v[sitePos]){
//			edgeLength += 1;
//		}
//	}
//	edgeLength /= float(sequenceLength);
//	if (u_id < v_id){
//		edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(u_id,v_id),edgeLength));
//	} else {
//		edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(v_id,u_id),edgeLength));
//	}
	double edgeLogLikelihoodUToV = 0;
	double edgeLogLikelihoodVToU = 0;
	array<array<float,4>,4> transitionMatrixUToV;
	array<array<float,4>,4> transitionMatrixVToU;
	for (unsigned char char_u = 0; char_u < 4; char_u++){
		for (unsigned char char_v = 0; char_v < 4; char_v++){
			transitionMatrixUToV[char_u][char_v] = 0.0;
			transitionMatrixVToU[char_v][char_u] = 0.0;
		}
	}
	unsigned char char_u; unsigned char char_v;
	for (int sitePos =0; sitePos < sequenceLength; sitePos++){
		char_u = seq_u[sitePos];
		char_v = seq_v[sitePos];
		transitionMatrixUToV[char_u][char_v] += 1;
		transitionMatrixVToU[char_v][char_u] += 1;
	}
	float rowSum_1; float rowSum_2;
	unsigned char char_p; unsigned char char_c;
	for (char_p = 0; char_p < 4; char_p++){
		rowSum_1 = 0; rowSum_2 = 0;
		for (char_c = 0; char_c < 4; char_c++){
			rowSum_1 += transitionMatrixUToV[char_p][char_c];
			rowSum_2 += transitionMatrixVToU[char_p][char_c];
		}
		for (char_c = 0; char_c < 4; char_c++){
			transitionMatrixUToV[char_p][char_c] /= float(rowSum_1);
			transitionMatrixVToU[char_p][char_c] /= float(rowSum_2);
		}
	}
	for (int sitePos = 0; sitePos < sequenceLength; sitePos++){
		char_u = seq_u[sitePos];
		char_v = seq_v[sitePos];
		edgeLogLikelihoodUToV += log(transitionMatrixUToV[char_u][char_v]);
		edgeLogLikelihoodVToU += log(transitionMatrixVToU[char_v][char_u]);
	}
//	cout << edgeLogLikelihoodUToV << "\tEdge log-likelihood_u_to_v" << endl;
//	cout << edgeLogLikelihoodVToU << "\tEdge log-likelihood_v_to_u" << endl;
//	cout << edgeLength << "\tedge length" << endl;
	edgeLogLikelihoodMap->insert(pair<pair<int,int>,double>(pair<int,int>(u_id,v_id),edgeLogLikelihoodUToV));
	edgeLogLikelihoodMap->insert(pair<pair<int,int>,double>(pair<int,int>(v_id,u_id),edgeLogLikelihoodVToU));
//	transitionMatricesForEdgesMap->insert(pair<pair<int,int>,array<array<float,4>,4>>(pair<int,int>(u_id,v_id),transitionMatrixUToV));
//	transitionMatricesForEdgesMap->insert(pair<pair<int,int>,array<array<float,4>,4>>(pair<int,int>(v_id,u_id),transitionMatrixVToU));
}

bool phylogeny_tree::ContainsVertex(int v_id){
	return (vertexMap->find(v_id) != vertexMap->end());
}

void phylogeny_tree::AddVertex(int v_id, string v_name, vector <unsigned char> seq_v){
	phylogeny_vertex * v = new phylogeny_vertex(v_id,v_name);
	vertexMap->insert(pair<int,phylogeny_vertex*>(v_id,v));
	array <float, 4> rootProbability;
	for (int i = 0; i < 4; i++){
		rootProbability[i] = 0;
	}
	int sequenceLength = seq_v.size();
	for (int sitePos = 0; sitePos < sequenceLength; sitePos++){
		rootProbability[seq_v[sitePos]] += 1;
	}
	for (unsigned char dna = 0; dna < 4; dna ++){
		rootProbability[dna] /= float(sequenceLength);
	}
	for (int sitePos = 0; sitePos < sequenceLength; sitePos++){
		v->logLikelihoodForObservingSequence += log(rootProbability[seq_v[sitePos]]);
	}
//	for (int i = 0; i < 4; i++){
//		cout << rootProbability[i] << "\t";
//	}
//	cout << endl;
}

//void phylogeny_tree::AddVertex(int v_id, string v_name){
//	phylogeny_vertex * v = new phylogeny_vertex(v_id,v_name);
//	vertexMap->insert(pair<int,phylogeny_vertex*>(v_id,v));
//}

void phylogeny_tree::AddEdge(int u_id, vector <unsigned char> seq_u, int v_id, vector <unsigned char> seq_v){
	phylogeny_vertex * u = (*this->vertexMap)[u_id];
	phylogeny_vertex * v = (*this->vertexMap)[v_id];
	u->neighbors.push_back(v);
	u->degree += 1;
	v->neighbors.push_back(u);
	v->degree +=1;
//	if (u_id == 163){
//		cout << "{" << u_id << "," << v_id << "}" << endl;
//	}
//	if (v_id == 163){
//		cout << "{" << u_id << "," << v_id << "}" << endl;
//	}
	this->edgeList->push_back(pair<int,int>(u_id, v_id));
//	ComputeDyadEdgeLengthsAndDyadTransitionMatricesAndDyadRootProbabilityAndDyadLogLikelihoods(u_id, seq_u, v_id, seq_v);
	ComputeEdgeLengthAndBidirectionalTransitionMatricesAndBidirectionalLogLikelihoods(u_id, seq_u, v_id, seq_v);
}

//void phylogeny_tree::AddEdge(int u_id, vector <unsigned char> seq_u, int v_id, vector <unsigned char> seq_v, vector <unsigned char> seq_r){
//	phylogeny_vertex * u = (*this->vertexMap)[u_id];
//	phylogeny_vertex * v = (*this->vertexMap)[v_id];
//	u->neighbors.push_back(v);
//	u->degree += 1;
//	v->neighbors.push_back(u);
//	v->degree +=1;
//	if (u_id < v_id){
//		(*this->rootSequencePerEdge)[pair<int,int>(u_id,v_id)] = seq_r;		
//	} else {
//		(*this->rootSequencePerEdge)[pair<int,int>(v_id,u_id)] = seq_r;		
//	}	
//	ComputeDyadEdgeLengthsAndDyadTransitionMatricesAndDyadRootProbabilityAndDyadLogLikelihoods(u_id, seq_u, v_id, seq_v);
//	ComputeDyadLogLikelihood(u_id, seq_u, v_id, seq_v);
//	ComputeEdgeLengthAndBidirectionalTransitionMatricesAndBidirectionalLogLikelihoods(u_id, seq_u, v_id, seq_v);
//}
//
//void phylogeny_tree::ComputeDyadEdgeLengthsAndDyadTransitionMatricesAndDyadRootProbabilityAndDyadLogLikelihoods(int u_id, vector<unsigned char> seq_u, int v_id, vector<unsigned char> seq_v){
////	 E-M for root state, root probability, and edge transition matrices
////	 Initialize root state
//	double dyadLogLikelihood;
//	double logLikelihoodRootToV;
//	double logLikelihoodRootToU;
//	double dyadLogLikelihood_prev;	
//	float edgeLengthFromRootToU = 0;
//	float edgeLengthFromRootToV = 0;
//	array<array<float,4>,4> transitionMatrixFromRootToU;
//	array<array<float,4>,4> transitionMatrixFromRootToV;	
//	array <float, 4> rootProbability;
//	int sequenceLength = int(seq_u.size());
//	vector<unsigned char> compressed_seq_u;
//	vector<unsigned char> compressed_seq_v;
//	vector<unsigned char> compressed_seq_r;
//	unsigned char char_u; unsigned char char_v;
//	float rowSum_U; float rowSum_V;
//	map <pair<unsigned char,unsigned char>, int> edgeSitePattern;
//	vector <int> sitePatternWeight;
//	double siteLikelihood;
//	int choice;
//	for (unsigned char char_u = 0; char_u < 4; char_u++){
//		for (unsigned char char_v = 0; char_v < 4; char_v++){
//			edgeSitePattern[pair<unsigned char,unsigned char>(char_u,char_v)] = 0;
//		}
//	}
//	for (int sitePos = 0; sitePos < int(seq_u.size()); sitePos++){
//		char_u = seq_u[sitePos];
//		char_v = seq_v[sitePos];
//		if (edgeSitePattern[pair<unsigned char,unsigned char>(char_u,char_v)] == 0){
//			compressed_seq_u.push_back(char_u);
//			compressed_seq_v.push_back(char_v);
//		}
//		edgeSitePattern[pair<unsigned char,unsigned char>(char_u,char_v)] += 1;
//	}	
//	for (int sitePos = 0; sitePos < int(compressed_seq_u.size()); sitePos++){
//		char_u = compressed_seq_u[sitePos];
//		char_v = compressed_seq_v[sitePos];
//		sitePatternWeight.push_back(edgeSitePattern[pair<unsigned char,unsigned char>(char_u,char_v)]);
//	}
//	
//	for (int sitePos = 0; sitePos < int(compressed_seq_u.size()); sitePos++){
//		char_u = compressed_seq_u[sitePos];
//		char_v = compressed_seq_v[sitePos];
//		if (char_u == char_v){
//			compressed_seq_r.push_back(char_u);
//		} else {
//			uniform_int_distribution <int> distribution(0,1);
//			choice = distribution(generator);				
//			if (choice == 0){
//				compressed_seq_r.push_back(char_u);
//			} else {
//				compressed_seq_r.push_back(char_v);
//			}
//		}
//	}	
//	dyadLogLikelihood = -100000000;
//	dyadLogLikelihood_prev = dyadLogLikelihood - 10;
//	int iteration = 0;
//	float currentProbability; unsigned char MAPState; 
//	float maxProbability = 0.0;
//	while (dyadLogLikelihood > dyadLogLikelihood_prev and iteration < 100){
//		iteration += 1;
//		dyadLogLikelihood_prev = dyadLogLikelihood;
//		dyadLogLikelihood = 0.0;
////		M-step root probability
//		unsigned char char_r;
//		for (int i=0; i< 4; i++){
//			rootProbability[i] = 0.0;
//		}
//		for (int sitePos = 0; sitePos < int(compressed_seq_u.size()); sitePos++){
//			char_r = compressed_seq_r[sitePos];
//			rootProbability[char_r] += sitePatternWeight[sitePos];
//		}
//		for (unsigned char dna = 0; dna < 4; dna ++){
//			rootProbability[dna] /= float(sequenceLength);
//		}
////		M-step transition matrices
////		MLE of transition matrix from root to u, and ,MLE of transition matrix from root to u.
//		for (unsigned char char_u = 0; char_u < 4; char_u++){
//			for (unsigned char char_v = 0; char_v < 4; char_v++){
//				transitionMatrixFromRootToU[char_u][char_v] = float(0);
//				transitionMatrixFromRootToV[char_u][char_v] = float(0);
//			}
//		}
//		for (int sitePos = 0; sitePos < int(compressed_seq_u.size()); sitePos++){
//			char_r = compressed_seq_r[sitePos];
//			char_u = compressed_seq_u[sitePos];
//			char_v = compressed_seq_v[sitePos];
//			transitionMatrixFromRootToU[char_r][char_u] += sitePatternWeight[sitePos];
//			transitionMatrixFromRootToV[char_r][char_v] += sitePatternWeight[sitePos];
//		}
//		for (unsigned char dna_ind_p = 0; dna_ind_p < 4; dna_ind_p ++){
//			rowSum_U = 0.0;
//			rowSum_V = 0.0;
//			for (unsigned char dna_ind_c = 0; dna_ind_c < 4; dna_ind_c ++){
//				rowSum_U += transitionMatrixFromRootToU[dna_ind_p][dna_ind_c];
//				rowSum_V += transitionMatrixFromRootToV[dna_ind_p][dna_ind_c];
//			}
//			for (unsigned char dna_ind_c = 0; dna_ind_c < 4; dna_ind_c ++){
//				transitionMatrixFromRootToU[dna_ind_p][dna_ind_c]/=float(rowSum_U);
//				transitionMatrixFromRootToV[dna_ind_p][dna_ind_c]/=float(rowSum_V);				
//			}
//		}
////		E-step root state
//		for (int sitePos = 0; sitePos < int(compressed_seq_u.size()); sitePos++){
//			siteLikelihood = 0.0;
//			maxProbability = 0.0;
//			char_u = compressed_seq_u[sitePos];
//			char_v = compressed_seq_v[sitePos];
//			MAPState = -1;
//			for (unsigned char dna = 0; dna < 4; dna++){
//				currentProbability = rootProbability[dna] * transitionMatrixFromRootToU[dna][char_u] * transitionMatrixFromRootToV[dna][char_v];
//				siteLikelihood += currentProbability;
//				if (currentProbability > maxProbability){
//					maxProbability = currentProbability;
//					MAPState = dna;
//				}
//			}
//			if (MAPState == -1){
//				cout << "error computing MAP state for root sequence";
//			} else {
//				compressed_seq_r[sitePos] = MAPState;
//			}			
//			dyadLogLikelihood += log(siteLikelihood) * sitePatternWeight[sitePos];
//			if (compressed_seq_r[sitePos] != compressed_seq_u[sitePos]){
//				edgeLengthFromRootToU += sitePatternWeight[sitePos];
//			}
//			if (compressed_seq_r[sitePos] != compressed_seq_v[sitePos]){
//				edgeLengthFromRootToV += sitePatternWeight[sitePos];
//			}
//		}
//	}
//	edgeLengthFromRootToU /= float(sequenceLength);
//	edgeLengthFromRootToV /= float(sequenceLength);
//	if (u_id < v_id){
//		edgeLengthsForDyadMap->insert(pair<pair<int,int>,pair<float,float>>(pair<int,int>(u_id,v_id),pair<float,float>(edgeLengthFromRootToU,edgeLengthFromRootToV)));
//		rootProbabilityForDyadMap->insert(pair<pair<int,int>,array<float,4>>(pair<int,int>(u_id,v_id),rootProbability));
////		dyadLogLikelihoodMap->insert(pair<pair<int,int>,double>(pair<int,int>(u_id,v_id),dyadLogLikelihood));
//		dyadSeperateLogLikelihoodMap->insert(pair<pair<int,int>,pair<double,double>>(pair<int,int>(u_id,v_id),pair<double,double>(logLikelihoodRootToU,logLikelihoodRootToV)));
//		transitionMatricesForDyadMap->insert(pair<pair<int,int>,pair<array<array<float,4>,4>,array<array<float,4>,4>>>(pair<int,int>(u_id,v_id),pair<array<array<float,4>,4>,array<array<float,4>,4>>(transitionMatrixFromRootToU,transitionMatrixFromRootToV)));
//	} else {
//		edgeLengthsForDyadMap->insert(pair<pair<int,int>,pair<float,float>>(pair<int,int>(v_id,u_id),pair<float,float>(edgeLengthFromRootToV,edgeLengthFromRootToU)));
//		rootProbabilityForDyadMap->insert(pair<pair<int,int>,array<float,4>>(pair<int,int>(v_id,u_id),rootProbability));
////		dyadLogLikelihoodMap->insert(pair<pair<int,int>,double>(pair<int,int>(v_id,u_id),dyadLogLikelihood));
//		dyadSeperateLogLikelihoodMap->insert(pair<pair<int,int>,pair<double,double>>(pair<int,int>(v_id,u_id),pair<double,double>(logLikelihoodRootToU,logLikelihoodRootToV)));
//		transitionMatricesForDyadMap->insert(pair<pair<int,int>,pair<array<array<float,4>,4>,array<array<float,4>,4>>>(pair<int,int>(v_id,u_id),pair<array<array<float,4>,4>,array<array<float,4>,4>>(transitionMatrixFromRootToV,transitionMatrixFromRootToU)));
//	}
//}
//
//void phylogeny_tree::ComputeDyadLogLikelihood(int u_id, vector<unsigned char> seq_u, int v_id, vector<unsigned char> seq_v){
//    vector <unsigned char> seq_r;
//	if (u_id < v_id){
//		seq_r = (*this->rootSequencePerEdge)[pair<int,int>(u_id,v_id)];
//	} else {
//		seq_r = (*this->rootSequencePerEdge)[pair<int,int>(v_id,u_id)];
//	}
//	int sequenceLength = seq_r.size();
//	// Compute rootProbability
//	array <float, 4> rootProbability;
//	for (int i=0; i< 4; i++){
//		rootProbability[i] = 0;
//	}
//	for (int sitePos = 0; sitePos < sequenceLength; sitePos++){
//		rootProbability[seq_r[sitePos]] += 1;
//	}
//	for (unsigned char dna = 0; dna < 4; dna ++){
//		rootProbability[dna] /= float(sequenceLength);
//	}
//	// Compute MLE of transition matrices
//	array<array<float,4>,4> transitionMatrixFromRootToU;
//	array<array<float,4>,4> transitionMatrixFromRootToV;	
//	float rowSum_U; float rowSum_V;
//	for (unsigned char char_u = 0; char_u < 4; char_u++){
//		for (unsigned char char_v = 0; char_v < 4; char_v++){
//			transitionMatrixFromRootToU[char_u][char_v] = 0;
//			transitionMatrixFromRootToV[char_u][char_v] = 0;
//		}
//	}
//	unsigned char char_r; unsigned char char_u; unsigned char char_v;
//	for (int sitePos = 0; sitePos < sequenceLength; sitePos++){
//		char_r = seq_r[sitePos];
//		char_u = seq_u[sitePos];
//		char_v = seq_v[sitePos];
//		transitionMatrixFromRootToU[char_r][char_u] += 1;
//		transitionMatrixFromRootToV[char_r][char_v] += 1;
//	}
//	for (unsigned char dna_ind_p = 0; dna_ind_p < 4; dna_ind_p ++){
//		rowSum_U = 0;
//		rowSum_V = 0;
//		for (unsigned char dna_ind_c = 0; dna_ind_c < 4; dna_ind_c ++){
//			rowSum_U += transitionMatrixFromRootToU[dna_ind_p][dna_ind_c];
//			rowSum_V += transitionMatrixFromRootToV[dna_ind_p][dna_ind_c];
//		}
//		for (unsigned char dna_ind_c = 0; dna_ind_c < 4; dna_ind_c ++){
//			transitionMatrixFromRootToU[dna_ind_p][dna_ind_c]/=float(rowSum_U);
//			transitionMatrixFromRootToV[dna_ind_p][dna_ind_c]/=float(rowSum_V);				
//		}
//	}	
//	double dyadLogLikelihood = 0;
//	for (int sitePos = 0; sitePos < sequenceLength; sitePos++){
//		char_r = seq_r[sitePos];
//		char_u = seq_u[sitePos];
//		char_v = seq_v[sitePos];
////		dyadLogLikelihood += log(rootProbability[seq_r[sitePos]]);
////		dyadLogLikelihood += log(transitionMatrixFromRootToU[char_r][char_u]);
////		dyadLogLikelihood += log(transitionMatrixFromRootToV[char_r][char_v]);		
//	}
////	double logLikelihoodRootToV;
////	double logLikelihoodRootToU;
////	double dyadLogLikelihood_prev;	
////	float edgeLengthFromRootToU = 0;
////	float edgeLengthFromRootToV = 0;
////	
////	vector<unsigned char> compressed_seq_u;
////	vector<unsigned char> compressed_seq_v;
////	vector<unsigned char> compressed_seq_r;
////	unsigned char char_u; unsigned char char_v;
////	float rowSum_U; float rowSum_V;
////	map <pair<unsigned char,unsigned char>, int> edgeSitePattern;
////	vector <int> sitePatternWeight;
////	double siteLikelihood;
////	int choice;
////	for (unsigned char char_u = 0; char_u < 4; char_u++){
////		for (unsigned char char_v = 0; char_v < 4; char_v++){
////			edgeSitePattern[pair<unsigned char,unsigned char>(char_u,char_v)] = 0;
////		}
////	}
////	for (int sitePos = 0; sitePos < int(seq_u.size()); sitePos++){
////		char_u = seq_u[sitePos];
////		char_v = seq_v[sitePos];
////		if (edgeSitePattern[pair<unsigned char,unsigned char>(char_u,char_v)] == 0){
////			compressed_seq_u.push_back(char_u);
////			compressed_seq_v.push_back(char_v);
////		}
////		edgeSitePattern[pair<unsigned char,unsigned char>(char_u,char_v)] += 1;
////	}	
////	for (int sitePos = 0; sitePos < int(compressed_seq_u.size()); sitePos++){
////		char_u = compressed_seq_u[sitePos];
////		char_v = compressed_seq_v[sitePos];
////		sitePatternWeight.push_back(edgeSitePattern[pair<unsigned char,unsigned char>(char_u,char_v)]);
////	}
////	
////	for (int sitePos = 0; sitePos < int(compressed_seq_u.size()); sitePos++){
////		char_u = compressed_seq_u[sitePos];
////		char_v = compressed_seq_v[sitePos];
////		if (char_u == char_v){
////			compressed_seq_r.push_back(char_u);
////		} else {
////			uniform_int_distribution <int> distribution(0,1);
////			choice = distribution(generator);				
////			if (choice == 0){
////				compressed_seq_r.push_back(char_u);
////			} else {
////				compressed_seq_r.push_back(char_v);
////			}
////		}
////	}
////	
////	dyadLogLikelihood = -100000000;
////	dyadLogLikelihood_prev = dyadLogLikelihood - 10;
////	int iteration = 0;
////	float currentProbability; unsigned char MAPState; 
////	float maxProbability = 0.0;
////	while (dyadLogLikelihood > dyadLogLikelihood_prev and iteration < 100){
////		iteration += 1;
////		dyadLogLikelihood_prev = dyadLogLikelihood;
////		dyadLogLikelihood = 0.0;
//////		M-step root probability
////		unsigned char char_r;
////		for (int i=0; i< 4; i++){
////			rootProbability[i] = 0.0;
////		}
////		for (int sitePos = 0; sitePos < int(compressed_seq_u.size()); sitePos++){
////			char_r = compressed_seq_r[sitePos];
////			rootProbability[char_r] += sitePatternWeight[sitePos];
////		}
////		for (unsigned char dna = 0; dna < 4; dna ++){
////			rootProbability[dna] /= float(sequenceLength);
////		}
//////		M-step transition matrices
//////		MLE of transition matrix from root to u, and ,MLE of transition matrix from root to u.
////		for (unsigned char char_u = 0; char_u < 4; char_u++){
////			for (unsigned char char_v = 0; char_v < 4; char_v++){
////				transitionMatrixFromRootToU[char_u][char_v] = float(0);
////				transitionMatrixFromRootToV[char_u][char_v] = float(0);
////			}
////		}
////		for (int sitePos = 0; sitePos < int(compressed_seq_u.size()); sitePos++){
////			char_r = compressed_seq_r[sitePos];
////			char_u = compressed_seq_u[sitePos];
////			char_v = compressed_seq_v[sitePos];
////			transitionMatrixFromRootToU[char_r][char_u] += sitePatternWeight[sitePos];
////			transitionMatrixFromRootToV[char_r][char_v] += sitePatternWeight[sitePos];
////		}
////		for (unsigned char dna_ind_p = 0; dna_ind_p < 4; dna_ind_p ++){
////			rowSum_U = 0.0;
////			rowSum_V = 0.0;
////			for (unsigned char dna_ind_c = 0; dna_ind_c < 4; dna_ind_c ++){
////				rowSum_U += transitionMatrixFromRootToU[dna_ind_p][dna_ind_c];
////				rowSum_V += transitionMatrixFromRootToV[dna_ind_p][dna_ind_c];
////			}
////			for (unsigned char dna_ind_c = 0; dna_ind_c < 4; dna_ind_c ++){
////				transitionMatrixFromRootToU[dna_ind_p][dna_ind_c]/=float(rowSum_U);
////				transitionMatrixFromRootToV[dna_ind_p][dna_ind_c]/=float(rowSum_V);				
////			}
////		}
//////		E-step root state
////		for (int sitePos = 0; sitePos < int(compressed_seq_u.size()); sitePos++){
////			siteLikelihood = 0.0;
////			maxProbability = 0.0;
////			char_u = compressed_seq_u[sitePos];
////			char_v = compressed_seq_v[sitePos];
////			MAPState = -1;
////			for (unsigned char dna = 0; dna < 4; dna++){
////				currentProbability = rootProbability[dna] * transitionMatrixFromRootToU[dna][char_u] * transitionMatrixFromRootToV[dna][char_v];
////				siteLikelihood += currentProbability;
////				if (currentProbability > maxProbability){
////					maxProbability = currentProbability;
////					MAPState = dna;
////				}
////			}
////			if (MAPState == -1){
////				cout << "error computing MAP state for root sequence";
////			} else {
////				compressed_seq_r[sitePos] = MAPState;
////			}			
////			dyadLogLikelihood += log(siteLikelihood) * sitePatternWeight[sitePos];
////			if (compressed_seq_r[sitePos] != compressed_seq_u[sitePos]){
////				edgeLengthFromRootToU += sitePatternWeight[sitePos];
////			}
////			if (compressed_seq_r[sitePos] != compressed_seq_v[sitePos]){
////				edgeLengthFromRootToV += sitePatternWeight[sitePos];
////			}
////		}
////	}
////	edgeLengthFromRootToU /= float(sequenceLength);
////	edgeLengthFromRootToV /= float(sequenceLength);
//	if (u_id < v_id){
////		edgeLengthsForDyadMap->insert(pair<pair<int,int>,pair<float,float>>(pair<int,int>(u_id,v_id),pair<float,float>(edgeLengthFromRootToU,edgeLengthFromRootToV)));
////		rootProbabilityForDyadMap->insert(pair<pair<int,int>,array<float,4>>(pair<int,int>(u_id,v_id),rootProbability));
//		dyadLogLikelihoodMap->insert(pair<pair<int,int>,double>(pair<int,int>(u_id,v_id),dyadLogLikelihood));
////		dyadSeperateLogLikelihoodMap->insert(pair<pair<int,int>,pair<double,double>>(pair<int,int>(u_id,v_id),pair<double,double>(logLikelihoodRootToU,logLikelihoodRootToV)));
////		transitionMatricesForDyadMap->insert(pair<pair<int,int>,pair<array<array<float,4>,4>,array<array<float,4>,4>>>(pair<int,int>(u_id,v_id),pair<array<array<float,4>,4>,array<array<float,4>,4>>(transitionMatrixFromRootToU,transitionMatrixFromRootToV)));
//	} else {
////		edgeLengthsForDyadMap->insert(pair<pair<int,int>,pair<float,float>>(pair<int,int>(v_id,u_id),pair<float,float>(edgeLengthFromRootToV,edgeLengthFromRootToU)));
////		rootProbabilityForDyadMap->insert(pair<pair<int,int>,array<float,4>>(pair<int,int>(v_id,u_id),rootProbability));
//		dyadLogLikelihoodMap->insert(pair<pair<int,int>,double>(pair<int,int>(v_id,u_id),dyadLogLikelihood));
////		dyadSeperateLogLikelihoodMap->insert(pair<pair<int,int>,pair<double,double>>(pair<int,int>(v_id,u_id),pair<double,double>(logLikelihoodRootToU,logLikelihoodRootToV)));
////		transitionMatricesForDyadMap->insert(pair<pair<int,int>,pair<array<array<float,4>,4>,array<array<float,4>,4>>>(pair<int,int>(v_id,u_id),pair<array<array<float,4>,4>,array<array<float,4>,4>>(transitionMatrixFromRootToV,transitionMatrixFromRootToU)));
//	}
//}

//
//void phylogeny_tree::SelectEdgeForRooting(){	
//	vector <phylogeny_vertex*> verticesToVisit;
//	for (pair<int,phylogeny_vertex*> idAndVertex: *this->vertexMap){
//		if (idAndVertex.second->degree == 1){
//			verticesToVisit.push_back(idAndVertex.second);			
//		}
//	}
//	phylogeny_vertex* v;
//	vector <phylogeny_vertex*> verticesVisited;
//	unsigned int numberOfVerticesToVisit = verticesToVisit.size();
//	while (numberOfVerticesToVisit > 0){
//		v = verticesToVisit[numberOfVerticesToVisit-1];
//		verticesToVisit.pop_back();
//		verticesVisited.push_back(v);
//		numberOfVerticesToVisit -= 1;
//		for (phylogeny_vertex* n: v->neighbors){
//			if (find(verticesVisited.begin(),verticesVisited.end(),n)==verticesVisited.end()){
//				n->timesVisited += 1;
//				n->sumOfEdgeLogLikelihoods += (*this->edgeLogLikelihoodMap)[pair<int,int>(n->id,v->id)] + v->sumOfEdgeLogLikelihoods;
//				if ((n->degree - n->timesVisited) == 1){
//					verticesToVisit.push_back(n);
//					numberOfVerticesToVisit +=1;
//				}
//			}
//		}		
//	}
//	verticesVisited.clear();
//	verticesToVisit.clear();
//	verticesToVisit.push_back(v);
//	numberOfVerticesToVisit = verticesToVisit.size();	
//	while (numberOfVerticesToVisit > 0){
//		v = verticesToVisit[numberOfVerticesToVisit-1];
//		verticesToVisit.pop_back();
//		verticesVisited.push_back(v);
//		numberOfVerticesToVisit -= 1;
//		for (phylogeny_vertex* n: v->neighbors){
//			if (find(verticesVisited.begin(),verticesVisited.end(),n) == verticesVisited.end()){
//				this->logLikelihoodForRooting = v->sumOfEdgeLogLikelihoods - (*this->edgeLogLikelihoodMap)[pair<int,int>(v->id,n->id)];
//				if (v->id < n->id){
//					this->logLikelihoodForRooting += (*this->dyadLogLikelihoodMap)[pair<int,int>(v->id,n->id)];
//				} else {
//					this->logLikelihoodForRooting += (*this->dyadLogLikelihoodMap)[pair<int,int>(n->id,v->id)];
//				}
//				if (this->logLikelihoodForRooting > this->maxLogLikelihood){
//					this->maxLogLikelihood = this->logLikelihoodForRooting;
//					if (v->id < n->id){
//						edgeForRooting.first = v;
//						edgeForRooting.second = n;
//					} else {
//						edgeForRooting.first = n;
//						edgeForRooting.second = v;
//					}
//				}
//				n->sumOfEdgeLogLikelihoods += (v->sumOfEdgeLogLikelihoods - n->sumOfEdgeLogLikelihoods - (*this->edgeLogLikelihoodMap)[pair<int,int>(v->id,n->id)]);
//				n->sumOfEdgeLogLikelihoods += (*this->edgeLogLikelihoodMap)[pair<int,int>(n->id,v->id)];
//				verticesToVisit.push_back(n);
//				numberOfVerticesToVisit += 1;
//			}
//		}
//	}
//}

void phylogeny_tree::SelectVertexForRooting(){
	vector <phylogeny_vertex*> verticesToVisit;
	for (pair<int,phylogeny_vertex*> idAndVertex: *this->vertexMap){
		if (idAndVertex.second->degree == 1){
			verticesToVisit.push_back(idAndVertex.second);			
		}
	}
	phylogeny_vertex* v;
	vector <phylogeny_vertex*> verticesVisited;
	// Traverse from leaves to central vertex
	unsigned int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0){
		v = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		verticesVisited.push_back(v);
		numberOfVerticesToVisit -= 1;
		for (phylogeny_vertex* n: v->neighbors){
			if (find(verticesVisited.begin(),verticesVisited.end(),n)==verticesVisited.end()){
				n->timesVisited += 1;
				n->sumOfEdgeLogLikelihoods += (*this->edgeLogLikelihoodMap)[pair<int,int>(n->id,v->id)] + v->sumOfEdgeLogLikelihoods;
				if ((n->degree - n->timesVisited) == 1){
					verticesToVisit.push_back(n);
					numberOfVerticesToVisit +=1;
				}
			}
		}		
	}
	verticesVisited.clear();
	verticesToVisit.clear();
	verticesToVisit.push_back(v);
	// Traverse from central vertex to leaves
	numberOfVerticesToVisit = verticesToVisit.size();	
	this->logLikelihoodForRooting = v->sumOfEdgeLogLikelihoods + v->logLikelihoodForObservingSequence;	
	this->maxLogLikelihood = this->logLikelihoodForRooting;
	this->vertexForRooting = v;	
	while (numberOfVerticesToVisit > 0){
		v = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		verticesVisited.push_back(v);
		numberOfVerticesToVisit -= 1;
		for (phylogeny_vertex* n: v->neighbors){
			if (find(verticesVisited.begin(),verticesVisited.end(),n) == verticesVisited.end()){
				n->sumOfEdgeLogLikelihoods += v->sumOfEdgeLogLikelihoods - n->sumOfEdgeLogLikelihoods - (*this->edgeLogLikelihoodMap)[pair<int,int>(v->id,n->id)];
				n->sumOfEdgeLogLikelihoods += (*this->edgeLogLikelihoodMap)[pair<int,int>(n->id,v->id)];
				this->logLikelihoodForRooting = n->sumOfEdgeLogLikelihoods + n->logLikelihoodForObservingSequence;
				if (this->logLikelihoodForRooting > this->maxLogLikelihood){
					this->maxLogLikelihood = this->logLikelihoodForRooting;
					this->vertexForRooting = n;
				}
				verticesToVisit.push_back(n);
				numberOfVerticesToVisit += 1;
			}
		}
	}
	cout << "logLikelihoodForObservingSequenceIs" << endl;
	cout << this->vertexForRooting->logLikelihoodForObservingSequence << endl;
	// Set directed edge list
	phylogeny_vertex * root = this->vertexForRooting;
	verticesVisited.clear();
	verticesToVisit.clear();
	verticesToVisit.push_back(root);
	numberOfVerticesToVisit = 1;
	while (numberOfVerticesToVisit > 0){
		v = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		verticesVisited.push_back(v);
		numberOfVerticesToVisit -= 1;
		for (phylogeny_vertex * n: v->neighbors){
			if (find(verticesVisited.begin(),verticesVisited.end(),n)==verticesVisited.end()){
				this->directedEdgeList->push_back(pair<int,int>(v->id,n->id));
				verticesToVisit.push_back(n);
				numberOfVerticesToVisit += 1;
			}
		}
	}
}

void phylogeny_tree::RootTreeAtEdge(int u_id, int v_id){	
	phylogeny_vertex * v;
	this->directedEdgeList->clear();
	this->directedEdgeList->push_back(pair<int,int>(-1,u_id));
	this->directedEdgeList->push_back(pair<int,int>(-1,v_id));
	vector<phylogeny_vertex *> verticesVisited;
	vector<phylogeny_vertex *> verticesToVisit;
	verticesToVisit.push_back((*this->vertexMap)[u_id]);
	verticesToVisit.push_back((*this->vertexMap)[v_id]);
	verticesVisited.push_back((*this->vertexMap)[u_id]);
	verticesVisited.push_back((*this->vertexMap)[v_id]);
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0){
		v = verticesToVisit[numberOfVerticesToVisit-1];
		numberOfVerticesToVisit -= 1;
		verticesToVisit.pop_back();
		for (phylogeny_vertex * n : v->neighbors){
			if (find(verticesVisited.begin(),verticesVisited.end(),n)==verticesVisited.end()){
				verticesToVisit.push_back(n);
				numberOfVerticesToVisit += 1;
				verticesVisited.push_back(n);
				this->directedEdgeList->push_back(pair<int,int>(v->id,n->id));
			}
		}
	}
}

void phylogeny_tree::RootTreeAtVertex(int v_id){
	this->directedEdgeList->clear();
	phylogeny_vertex * root = (*this->vertexMap)[v_id];
	phylogeny_vertex * v;
	vector <phylogeny_vertex *> verticesVisited;
	vector <phylogeny_vertex *> verticesToVisit;
	verticesToVisit.push_back(root);
	int numberOfVerticesToVisit = 1;
	while (numberOfVerticesToVisit > 0){
		v = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		verticesVisited.push_back(v);
		numberOfVerticesToVisit -= 1;
		for (phylogeny_vertex * n: v->neighbors){
			if (find(verticesVisited.begin(),verticesVisited.end(),n)==verticesVisited.end()){
				this->directedEdgeList->push_back(pair<int,int>(v->id,n->id));
				verticesToVisit.push_back(n);
				numberOfVerticesToVisit += 1;
			}
		}
	}	
}
//
//void phylogeny_tree::WriteOutputFiles(string sequenceFileName){		
//	string newickLabel = "";
//	vector <phylogeny_vertex*> verticesToVisit;
//	for (pair<int,phylogeny_vertex*> idAndVertex: *this->vertexMap){
//		if (idAndVertex.second->degree == 1){
//			idAndVertex.second->newickLabel = idAndVertex.second->name;
//			if (idAndVertex.second->id != this->edgeForRooting.first->id and idAndVertex.second->id != this->edgeForRooting.second->id){
//				verticesToVisit.push_back(idAndVertex.second);
//			}
//		} else {
//			idAndVertex.second->timesVisited = 0;
//		}
//	}
//	vector <phylogeny_vertex*> verticesVisited;
//	phylogeny_vertex* v;
//	unsigned int numberOfVerticesToVisit = verticesToVisit.size();
//	int u_id = edgeForRooting.first->id;
//	int v_id = edgeForRooting.second->id;
//	while (numberOfVerticesToVisit > 0){
//		v = verticesToVisit[numberOfVerticesToVisit-1];
//		verticesToVisit.pop_back();
//		verticesVisited.push_back(v);
//		numberOfVerticesToVisit -= 1;
//		float edgeLength; int n_ind;
//		for (phylogeny_vertex* n: v->neighbors){
//			if (find(verticesVisited.begin(),verticesVisited.end(),n)==verticesVisited.end()){
//				n->timesVisited += 1;
//				if (n->id < v->id){
//					edgeLength = (*this->edgeLengthsMap)[pair<int,int>(n->id,v->id)];
//				} else {
//					edgeLength = (*this->edgeLengthsMap)[pair<int,int>(v->id,n->id)];
//				}
//				if ((n->degree - n->timesVisited) == 1){
//					n->newickLabel += "," + v->newickLabel + ":" + to_string(edgeLength) + ")";
//					if (n->id != u_id and n->id != v_id){
//						verticesToVisit.push_back(n);
//						numberOfVerticesToVisit +=1;
//					}
//				} else {
//					n->newickLabel += "(" + v->newickLabel + ":" + to_string(edgeLength);
//				}
//			} else {
//				n_ind = find(verticesVisited.begin(),verticesVisited.end(),n) - verticesVisited.begin();
//				verticesVisited.erase(verticesVisited.begin()+n_ind);
//			}
//		}
//	}
//	float edgeLength_first; float edgeLength_second;	
//	tie(edgeLength_first, edgeLength_second) = (*this->edgeLengthsForDyadMap)[pair<int,int>(u_id,v_id)];
//	newickLabel = "(" + edgeForRooting.first->newickLabel + ":" + to_string(edgeLength_first);
//	newickLabel += "," + edgeForRooting.second->newickLabel + ":" + to_string(edgeLength_second)+");";
//	ofstream newickFile;
//	newickFile.open(sequenceFileName+".newick");
//	newickFile << newickLabel << endl;
//	newickFile.close();

//	bool writeGMMToFile = 0;
//	
//	if (writeGMMToFile){
//		ofstream GMMFile;
//		GMMFile.open(sequenceFileName+".gmm");
//		GMMFile << "Root probability" << endl;	
//		array <float,4> rootProbability;
//
//		rootProbability = (*this->rootProbabilityForDyadMap)[pair<int,int>(u_id,v_id)];	
//		
//		for (unsigned char dna = 0; dna < 4; dna ++){
//			GMMFile << "P(" << this->DNA_alphabet[dna] << ")" << "\t";
//		}
//		GMMFile << endl;
//		for (unsigned char dna = 0; dna < 4; dna ++){
//			GMMFile << rootProbability[dna] << "\t";
//		}
//		GMMFile << endl;
//		GMMFile << "Transition matrices" << endl;
//		GMMFile << "Vertex_from" << "\t" << "Vertex_to" << "\t";
//		for (unsigned char dna_from = 0; dna_from < 4; dna_from ++){
//			for (unsigned char dna_to = 0; dna_to < 4; dna_to ++){
//				GMMFile << "P(" << this->DNA_alphabet[dna_from] << "->" << this->DNA_alphabet[dna_to] << ")" << "\t";
//			}		
//		}
//		GMMFile << endl;
//		array <array<float,4>,4> transitionMatrixFromRootToU;
//		array <array<float,4>,4> transitionMatrixFromRootToV;
//		
//		tie (transitionMatrixFromRootToU,transitionMatrixFromRootToV) = (*this->transitionMatricesForDyadMap)[pair<int,int>(u_id,v_id)];
//		
//		GMMFile << "h_root" << "\t" << this->edgeForRooting.first->name << "\t";
//		for (unsigned char dna_from = 0; dna_from < 4; dna_from ++){
//			for (unsigned char dna_to = 0; dna_to < 4; dna_to ++){
//				GMMFile << transitionMatrixFromRootToU[dna_from][dna_to] << "\t";
//			}		
//		}
//		GMMFile << endl;
//		GMMFile << "h_root" << "\t" << this->edgeForRooting.second->name << "\t";
//		for (unsigned char dna_from = 0; dna_from < 4; dna_from ++){
//			for (unsigned char dna_to = 0; dna_to < 4; dna_to ++){
//				GMMFile << transitionMatrixFromRootToV[dna_from][dna_to] << "\t";
//			}		
//		}
//		GMMFile << endl;
//		verticesToVisit.clear();
//		verticesVisited.clear();
//		verticesToVisit.push_back(this->edgeForRooting.first);
//		verticesToVisit.push_back(this->edgeForRooting.second);
//		verticesVisited.push_back(this->edgeForRooting.first);
//		verticesVisited.push_back(this->edgeForRooting.second);
//		numberOfVerticesToVisit = int(verticesToVisit.size());
//		array <array <float, 4>,4> transitionMatrixFromVToN;
//		while (numberOfVerticesToVisit > 0){
//			v = verticesToVisit[numberOfVerticesToVisit-1];
//			verticesToVisit.pop_back();
//			verticesVisited.push_back(v);
//			numberOfVerticesToVisit -= 1;
//			for (phylogeny_vertex* n: v->neighbors){
//				if (find(verticesVisited.begin(),verticesVisited.end(),n) == verticesVisited.end()){
//					verticesToVisit.push_back(n);
//					numberOfVerticesToVisit += 1;
//					GMMFile << v->name << "\t" << n->name << "\t";				
//					transitionMatrixFromVToN = (*this->transitionMatricesForEdgesMap)[pair<int,int>(v->id,n->id)];
//					for (unsigned char dna_from = 0; dna_from < 4; dna_from ++){
//						for (unsigned char dna_to = 0; dna_to < 4; dna_to ++){
//							GMMFile << transitionMatrixFromVToN[dna_from][dna_to] << "\t";
//						}
//					}
//					GMMFile << endl;
//				}
//			}
//		}
//		GMMFile.close();
//	}
//}

#endif