#ifndef rootedPhylogeny_H
#define rootedPhylogeny_H
#include <iostream>
#include "utilities.h"
#include <boost/algorithm/string.hpp>
#include <boost/math/tools/minima.hpp>
#include <random>
#include <chrono>
#include <iomanip>
#include <math.h>
#include "brent.h"
#include "asa047.hpp"
using namespace Eigen;

//using namespace boost;
//class LogLikelihood : public Problem<double> {
//  public:
//    using typename cppoptlib::Problem<double>::Scalar;
//    using typename cppoptlib::Problem<double>::TVector;
//
//    double value(vector<float> parameters) {
//        
//        return   t1 * t1 + 100 * t2 * t2;
//    }
//    void gradient(const TVector &x, TVector &grad) {
//        grad[0]  = -2 * (1 - x[0]) + 200 * (x[1] - x[0] * x[0]) * (-2 * x[0]);
//        grad[1]  = 200 * (x[1] - x[0] * x[0]);
//    }
//};

class rootedPhylogeny_vertex{	
public:
	int parent_id = -1;	
	int numberOfDescendants = 0;
	int timesVisited = 0;
	int rateCategory = -1;
	vector <int> children_id;	
	int id;
	string name;	
	string newickLabel;
	vector <unsigned char> sequence;
	vector <unsigned char> compressedSequence;
	array <float, 4> baseFreq;
//	vector <array<int,4>> conditionalLikelihood_int;
	vector <array<double,4>> conditionalLikelihood;
	void AddParent(int p_id);
	void AddChild(int c_id);
	void AddNumberOfObservedSequences(int numberOfObservedSequencesToAdd);
	rootedPhylogeny_vertex(int idToAdd, string nameToAdd, vector <unsigned char> sequenceToAdd){
		id = idToAdd;
		name = nameToAdd;
		sequence = sequenceToAdd;
		for (int dna = 0; dna < 4; dna++){
			baseFreq[dna] = 0;
		}
		newickLabel = "";
	}
	~rootedPhylogeny_vertex(){
		
	}
};


void rootedPhylogeny_vertex::AddParent(int p_id){
	this->parent_id = p_id;
}

void rootedPhylogeny_vertex::AddChild(int c_id){
	this->children_id.push_back(c_id);
}

class rootedPhylogeny_tree {
private:
	int numberOfObservedSequences;
	default_random_engine generator;		
	map <int, float> * scalingFactorForRateCategory;
	map <int, MatrixXf> * stationaryDistributionForCategory;
	map <int, int> * changePointForRateCategory;
	rootedPhylogeny_vertex * root;
	map <rootedPhylogeny_vertex*, float> * baseFreqChangeForEachDescendant;
	vector <rootedPhylogeny_vertex*> * changePoints;	
	int rateCategoryForOptimization;
	vector <float> * thresholdList;
	MatrixXd searchDirection;
public:
	map <int, MatrixXd> * parametersPerRateCategory;
	map <int, Matrix4f> * rateMatrixPerRateCategory;
	float optimalThreshold;
	float optimal_BIC;
	double logLikelihood;
	array <double, 11> jacobian_dep;
	float scalingFactor;
	MatrixXd JacobianDep;
	MatrixXd HessianDep;
	array <float, 4> rootProbability;
	vector <int> siteWeights;
	Matrix4f rateMatrix;	
	VectorXd initialEstimateForFreeParametersOfQ;
	VectorXd initialEstimateForFreeParameters;
	VectorXd freeParametersExcBaseFreq;
	MatrixXd freeParametersIncBaseFreq;
	int numberOfRateCategories;	
	vector <rootedPhylogeny_vertex*> listOfChangePoints;
//	vector <int> * leaf_ids;
	map <pair<int, int>, float> * edgeLengthsMap;
	vector <pair<int, int>> * edgesForPostOrderTreeTraversal;
	vector <rootedPhylogeny_vertex*> * verticesForPreOrderTreeTraversal;
	map <int, rootedPhylogeny_vertex*>* vertexMap;
	vector <rootedPhylogeny_vertex *> * leaves;
//	vector <pair<int,int>> * directedEdgeList;
	float ComputeNTDiff(rootedPhylogeny_vertex* p, rootedPhylogeny_vertex* c);
	void AddDirectedEdges(vector <pair<int, int>> * directedEdgeList_ptr);
	void AddEdge(int p_id, int c_id);
	void AddEdgeLength(int p_id, int c_id, float edgeLength);
	void SetEdgeLength(int p_id, int c_id, float edgeLength);
	void RemoveEdges();
	float GetEdgeLength(int u_id, int v_id);
	void ComputeAndSetEdgeLength(int u_id, int v_id);
	void ComputeEdgeLengths();
	void AddVertex(int idToAdd, string nameToAdd, vector <unsigned char> sequenceToAdd);
	void AddVertex(int idToAdd, string nameToAdd, vector <unsigned char> sequenceToAdd, vector <unsigned char> compressdSequenceToAdd);
	void ContractEdge(int idOfVertexToKeep, int idOfVertexToRemove);
	bool IsEdgeContractionFeasbile(int vertex_id_1, int vertex_id_2);
	void ContractOutdegreeZeroLeafIncidentZeroLengthEdges();
	void ContractZeroLengthEdges();
	void ContractEdgesBasedOnAIC();
	void ContractEdgesBasedOnBIC();
	float ComputeScalingFactor(Matrix4f Q);
	MatrixXf ComputeStationaryDistribution(Matrix4f Q);
	void AddNumberOfObservedSequences(int numberOfObservedSequencesToAdd);
	void WriteEdgeList(string sequenceFileName);
	void WriteAncestralSequences(string sequenceFileName);
	void ComputeBaseFreq();
	void SetThresholds();
	double local_min_brent();
	static double powell (double x[4]);
	void testNelderMead_1();
	void testNelderMead_2();
	void WriteNewickFile(string sequenceFileName);
//	void SetLeafIds();
	void ComputeNumberOfDescendants();
	void ReadTreeFile(string treeFileName);
	void SetParameters(VectorXd x);
	void SetRateMatrixForRateCategory(Matrix4f Q, int rateCategory);
	void AddStationaryDistributionForCategory(MatrixXf stationaryDistribution ,int rateCategory);
	void AddRateCategoryForVertex(int v_id, int rateCategory);	
	void SetChangePointForRateCategory(int v_id, int rateCategory);
//	void ComputeLogLikelihoodUsingAncestralStates();
	MatrixXd GetFreeParametersDep(Matrix4f Q);
	Matrix4f GetRateMatrixForFreeParameters(MatrixXd B);
	void SetRateMatrixUsingFreeParametersDep(MatrixXd B);
	void InitializeModelParametersForNelderMead();
	void InitializeModelParametersForBFGS();
	void OptimizeModelParametersDep();
	void OptimizeModelParametersForRateCatDep(int rateCat);
	void ComputeMLEOfRootProbability();
	void OptimizeModelParametersForAGivenThresholdUsingNelderMead(float threshold);
	void EstimateAncestralSequencesByFittingTheGMM();
	void ComputeMPEstimateOfAncestralSequences();
	void ComputeMAPEstimateOfAncestralSequences();
	void ComputeInitialEstimateForRateMatrix();
	void ScaleEdgeLengths();
	int GetVertexId(string v_name);
	void SetMinLengthOfEdges();
	void ComputeMLEOfRateMatrices();
	void ComputeMLEOfRateMatricesForLeafLabeledTrees();
	void ComputeMLEOfEdgeLengths();
	void NelderMeadForPowell(int n, double start[], double xmin[], double *ynewlo,
		 double reqmin, double step[], int konvge, int kcount, 
		 int *icount, int *numres, int *ifault );
	void NelderMeadForOptimizingRateParametersForRateCat(int rateCat, int n, double start[], double xmin[], 
		 double *ynewlo, double reqmin, double step[], int konvge,
		 int kcount, int *icount, int *numres, int *ifault);
	void NelderMeadForOptimizingParametersOfRootProb(int n, double start[], double xmin[], 
		 double *ynewlo, double reqmin, double step[], int konvge,
		 int kcount, int *icount, int *numres, int *ifault);
	void ComputeJacobianForFullyLabeledTree();
	void ComputeHessianForRateMatrixParametersDep();
	void ComputeInitialEstimateForFreeParametersIncBaseFreq();
	void ComputeConditionalLikelihoods();
	MatrixXd GetJacobianForRateCategory(int rateCategoryForOptimization);
	void ComputeJacobianOfLogLikelihoodDep();
	MatrixXd GetJacobianOfLogLikelihoodForRootProbParameters();
	MatrixXd GetHessianOfLogLikelihoodForRootProbParameters();
	void ComputeLogLikelihoodUsingStoredConditionals();
	void ComputeLogLikelihood();
	double GetBIC();
	void PerformModelSelectionUsingNelderMead();
	void PerformModelSelectionUsingBFGS();
	double GetNegLogLikelihoodForStepSize (double stepSize);	
	double GetNegLogLikelihoodForRateParametersForRateCat (double x[], int rateCat);
	double GetNegLogLikelihoodForParametersOfRootProb (double x[]);
	float ComputeAbsDifferenceInBaseFreq(array <float, 4> baseFreq_1, array <float, 4> baseFreq_2);
	double SampleFunctionForMinimization(double x);
	double BrentLineSearchForSelectingStepSize(double a, double b, double t, double &x);
	double GetOptimalStepSizeUsingLineSearch(MatrixXd B_current);
	double ComputeNormOfJacobian(MatrixXd Jacobian);
	void ComputeLogLikelihoodForFullyLabeledTree();
	void SetParametersForRateMatrixForNelderMead(double x[], int rateCat);
	void SetEdgesForPostOrderTreeTraversal();
	void SetEdgeLengths();
	void SetLeaves();
	void SetVerticesForPreOrderTreeTraversal();
	void SetParametersForRateMatrixForBFGS(MatrixXd parameters, int rateCat);
	void SetRateCategories(float baseFreqThreshold);	
	void WriteRateCategoryPerVertexAndModelParameters(string sequenceFileName);
	void OptimizeModelParametersUsingNelderMead();
	void OptimizeModelParametersUsingBFGS();
	void OptimizeParametersForRateMatrixUsingNelderMead(int rateCat);
	void OptimizeParametersForRateMatrixUsingBFGS(int rateCat);
	void OptimizeEdgeLengthsForRateCat(int rateCat);
	void OptimizeModelParametersForRootProbUsingNelderMead();
	void OptimizeModelParametersForRootProbUsingNewtonRaphson();
	double ComputeEdgeLogLikelihood(rootedPhylogeny_vertex * p, rootedPhylogeny_vertex * c, Matrix4f P);
	void SetSiteWeights(vector <int> siteWeightsToSet);
	array <double, 4> GetLikelihoodArray(double elem_0, double elem_1, double elem_2, double elem_3);
	void InitializeConditionalLikelihoods();
	void ResetConditionalLikelihoodsForAncestors();	
	Matrix4f ComputeFirstDerivativeOfRateMatrix(int rateCat, int par);
	Matrix4f ComputeFirstDerivativeOfMatrixExponential(float t, int rateCat, int par); 	
	Matrix4f ComputeFirstDerivativeOfRateMatrixDep(int par);
	Matrix4f ComputeDerivativeOfMatrixExponentialDep(float t, int par);
	Matrix4f ComputeSecondDerivativeOfRateMatrix(int par_1, int par_2);
	Matrix4f ComputeSecondDerivativeOfMatrixExponential(float t, int par_1, int par_2);	
	rootedPhylogeny_tree(){
		this->optimal_BIC = pow(10,10);
		this->leaves = new vector <rootedPhylogeny_vertex *>;
		this->thresholdList = new vector <float>;
		this->initialEstimateForFreeParameters = VectorXd(11);
		this->freeParametersExcBaseFreq = VectorXd(11);
		this->rateMatrix = ArrayXXf::Zero(4,4);		
		this->JacobianDep = ArrayXXd::Zero(11,1);
		this->HessianDep = ArrayXXd::Zero(11,11);
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		this->generator = default_random_engine(seed);
		vector <unsigned char> root_sequence;
		this->verticesForPreOrderTreeTraversal = new vector <rootedPhylogeny_vertex*>;
		this->baseFreqChangeForEachDescendant = new map <rootedPhylogeny_vertex*, float>;
		this->root = new rootedPhylogeny_vertex(-1, "h_root", root_sequence);
		this->vertexMap = new map <int, rootedPhylogeny_vertex*>;
		this->changePoints= new vector <rootedPhylogeny_vertex*>;
		this->vertexMap->insert(pair<int, rootedPhylogeny_vertex*>(-1, this->root));
		this->edgesForPostOrderTreeTraversal = new vector<pair<int, int>>;
		this->edgeLengthsMap = new map <pair<int, int>, float>;
		this->parametersPerRateCategory = new map <int, MatrixXd>;
		this->rateMatrixPerRateCategory = new map <int, Matrix4f>;
		this->stationaryDistributionForCategory = new map <int, MatrixXf>;
		this->changePointForRateCategory = new map <int, int>;
		this->scalingFactorForRateCategory = new map <int, float>;
	}
	~rootedPhylogeny_tree(){
		for (pair<int,rootedPhylogeny_vertex*> idAndVertexPtrPair : (*this->vertexMap)){
			delete idAndVertexPtrPair.second;
		}
		delete this->vertexMap;
		delete this->thresholdList;
		delete this->leaves;
		delete this->edgesForPostOrderTreeTraversal;
		delete this->edgeLengthsMap;
		delete this->parametersPerRateCategory;
		delete this->rateMatrixPerRateCategory;
		delete this->stationaryDistributionForCategory;
		delete this->verticesForPreOrderTreeTraversal;
		delete this->changePointForRateCategory;
		delete this->scalingFactorForRateCategory;
		delete this->baseFreqChangeForEachDescendant;
		delete this->changePoints;
	}
};

void rootedPhylogeny_tree::SetLeaves(){
	this->leaves->clear();
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
		if (idPtrPair.second->children_id.size() == 0){
			this->leaves->push_back(idPtrPair.second);
		}
	}
}

void rootedPhylogeny_tree::WriteRateCategoryPerVertexAndModelParameters(string sequenceFileName){
	ofstream modelParametersFile;
	modelParametersFile.open(sequenceFileName+".modelParameters");
	modelParametersFile << "Vertex name" << "\t" << "Rate category" << endl;
	rootedPhylogeny_vertex * v;
	for (pair<int,rootedPhylogeny_vertex *> idPtrPair : (*this->vertexMap)){
		v = idPtrPair.second;
		modelParametersFile << v->name << "\t" << v->rateCategory << endl;
	}
	modelParametersFile << "==============================================" << endl;
	modelParametersFile << "Root probability:" << endl;
	// A : 0, C : 1, G : 2, T : 3
	string dnaString = "ACGT";
	for (int i = 0; i < 4; i++){
		modelParametersFile << "P(" << dnaString[i] << ") = " ;
		modelParametersFile << this->rootProbability[i] << "\t";
	}
	modelParametersFile << endl;
	modelParametersFile << "==============================================" << endl;
	modelParametersFile << "Rate matrix parameters:" << endl;
	int rateCat; Matrix4f Q;
	for (pair<int,Matrix4f> catAndRateMatrixPair : *this->rateMatrixPerRateCategory){
		rateCat = catAndRateMatrixPair.first;
		Q = catAndRateMatrixPair.second;
		Q /= this->ComputeScalingFactor(Q);
		modelParametersFile << "Rate matrix for category " << rateCat << " is " << endl;
		modelParametersFile << Q << endl;
	}
	modelParametersFile << "==============================================" << endl;
	modelParametersFile.close();	
}

float rootedPhylogeny_tree::ComputeNTDiff(rootedPhylogeny_vertex* p, rootedPhylogeny_vertex* c){
	float ntDiff = 0;
	float seqLength = 0;
	unsigned char dna_p; unsigned char dna_c;
	for (unsigned int site = 0; site < this->siteWeights.size(); site++){
		dna_p = p->compressedSequence[site];
		dna_c = c->compressedSequence[site];
		if (dna_p != dna_c){
			ntDiff += this->siteWeights[site];
		}
		seqLength += this->siteWeights[site];
	}
	ntDiff /= seqLength;
	return (ntDiff);
}

void rootedPhylogeny_tree::SetEdgeLengths(){
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	float edgeLength;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			edgeLength = this->ComputeNTDiff(p,c);
			if (p->id < c->id){
				(*this->edgeLengthsMap)[pair<int,int>(p->id,c->id)] = edgeLength;
			} else {
				(*this->edgeLengthsMap)[pair<int,int>(c->id,p->id)] = edgeLength;
			}
		}	
	}
}

void rootedPhylogeny_tree::PerformModelSelectionUsingBFGS(){
	// Compute ancestral states by fitting the GMM
	this->EstimateAncestralSequencesByFittingTheGMM();
	this->SetLeaves();
	// Set edge lengths
	this->SetEdgeLengths();	
//	cout << this->GetEdgeLength(this->GetVertexId("h_root"),this->GetVertexId("h_1")) << endl;
	this->SetMinLengthOfEdges();
	this->SetVerticesForPreOrderTreeTraversal();
	this->SetEdgesForPostOrderTreeTraversal();
	// Compute base freq changes
	this->ComputeBaseFreq();					
	// Compute list of thresholds
	this->SetThresholds();
	// Iterate over thresholds in order of decreasing value	
	int noOfThresholdsTried = 0;
	double current_BIC = 0;
	this->optimal_BIC = 0;
	for (float threshold : (*this->thresholdList)){
		noOfThresholdsTried += 1;
		this->SetRateCategories(threshold);
		this->InitializeModelParametersForBFGS();
		this->OptimizeModelParametersUsingBFGS();
		current_BIC = this->GetBIC();		
		if (this->optimal_BIC > current_BIC or noOfThresholdsTried == 1){			
			this->optimal_BIC = current_BIC;
			this->optimalThreshold = threshold;
			cout << "number of rate categories is " << this->numberOfRateCategories << endl;
			cout << "optimal BIC is " << this->optimal_BIC << endl;
		} else {
			break;
		}
	}
}

void rootedPhylogeny_tree::PerformModelSelectionUsingNelderMead(){
	// Compute ancestral states by fitting the GMM
	this->EstimateAncestralSequencesByFittingTheGMM();
	this->SetLeaves();
	// Set edge lengths
	this->SetEdgeLengths();	
//	cout << this->GetEdgeLength(this->GetVertexId("h_root"),this->GetVertexId("h_1")) << endl;
	this->SetMinLengthOfEdges();
	this->SetVerticesForPreOrderTreeTraversal();
	this->SetEdgesForPostOrderTreeTraversal();
	// Compute base freq changes
	this->ComputeBaseFreq();					
	// Compute list of thresholds
	this->SetThresholds();
	// Iterate over thresholds in order of decreasing value	
	int noOfThresholdsTried = 0;
	double current_BIC = 0;
	this->optimal_BIC = 0;
	for (float threshold : (*this->thresholdList)){
		noOfThresholdsTried += 1;
		this->SetRateCategories(threshold);
//		cout << "number of rate categories is " << this->numberOfRateCategories << endl;
		this->InitializeModelParametersForNelderMead();
		this->OptimizeModelParametersUsingNelderMead();
		current_BIC = this->GetBIC();				
		if (this->optimal_BIC > current_BIC or noOfThresholdsTried == 1){			
			this->optimal_BIC = current_BIC;
			this->optimalThreshold = threshold;
//			cout << "number of rate categories for optimal threshold is " << this->numberOfRateCategories << endl;
//			cout << "optimal BIC is " << this->optimal_BIC << endl;
		} else {
			break;
		}
	}
}


void rootedPhylogeny_tree::OptimizeModelParametersForAGivenThresholdUsingNelderMead(float threshold){
	this->EstimateAncestralSequencesByFittingTheGMM();
	// Set edge lengths
	this->SetEdgeLengths();
	this->SetMinLengthOfEdges();
	// Compute base freq changes
	this->ComputeBaseFreq();
	this->SetRateCategories(threshold);
	this->InitializeModelParametersForNelderMead();
	this->OptimizeModelParametersUsingNelderMead();
	
}
double rootedPhylogeny_tree::GetBIC(){
	float sequenceLength = 0;
	for (unsigned int site = 0; site < this->siteWeights.size(); site++){
		sequenceLength += this->siteWeights[site];
	}
	double numberOfParameters;
//	cout << "current logLikelihood is " << this->logLikelihood << endl;
//	cout << "number of rate categories is " << this->numberOfRateCategories << endl;
	// If no of vertices in the same category as the root is one then 
	rootedPhylogeny_vertex * c_l;
	rootedPhylogeny_vertex * c_r;
	c_l = (*this->vertexMap)[this->root->children_id[0]];
	c_r = (*this->vertexMap)[this->root->children_id[1]];
	if (this->root->rateCategory != c_l->rateCategory and this->root->rateCategory != c_r->rateCategory){
		numberOfParameters = (11.0 * (this->numberOfRateCategories -1)) + 3.0;
	} else {
		numberOfParameters = 11.0 * (this->numberOfRateCategories);
	}
	double penaltyPerParameter = 1.0;
	double optimal_BIC = ( log(sequenceLength) * penaltyPerParameter * numberOfParameters ) - ( 2.0 * this->logLikelihood );	
	return optimal_BIC;
}

void rootedPhylogeny_tree::SetParametersForRateMatrixForNelderMead(double x[], int rateCat){
	for (int par = 0; par < 11; par ++) {
		(*this->parametersPerRateCategory)[rateCat](par,0) = x[par];
	}
	Matrix4f Q;
	// a
	Q(0,1) = x[0];
	// b
	Q(0,2) = x[1];
	// c
	Q(0,3) = x[2];
	// d
	Q(1,0) = x[3];
	// e
	Q(1,2) = x[4];
	// f
	Q(1,3) = x[5];
	// g
	Q(2,0) = x[6];
	// h
	Q(2,1) = x[7];
	// i
	Q(2,3) = x[8];
	// j
	Q(3,0) = x[9];
	// k
	Q(3,1) = x[10];
	// set l to 1
	Q(3,2) = 1.0;
	// D1 
	Q(0,0) = - (Q(0,1) + Q(0,2) + Q(0,3));
	// D2 
	Q(1,1) = - (Q(1,0) + Q(1,2) + Q(1,3));
	// D3 
	Q(2,2) = - (Q(2,0) + Q(2,1) + Q(2,3));
	// D4 
	Q(3,3) = - (Q(3,0) + Q(3,1) + Q(3,2));
	
	(*this->rateMatrixPerRateCategory)[rateCat] = Q;
	float scalingFactor = this->ComputeScalingFactor(Q);
	(*this->scalingFactorForRateCategory)[rateCat] = scalingFactor;
	if (this->root->rateCategory == rateCat){
		MatrixXf stationaryDistribution = this->ComputeStationaryDistribution(Q);
		for (int i = 0; i < 4; i++){
			if (stationaryDistribution(i,0) < 0){
				cout << "stationary distribution has negative entry" << endl;
				cout << "Rate matrix is " << endl << Q << endl;
			}
			this->rootProbability[i] = stationaryDistribution(i,0);
		}
	}
	return;
}

double rootedPhylogeny_tree::GetNegLogLikelihoodForRateParametersForRateCat(double x[], int rateCat){	
	// Set parameters for rate cat
	// Check if constraints are being satisfied
	bool setParameters = 1;
	double valueToReturn;
	for (int par = 0; par < 11; par++){
		if (x[par] < pow(10,-4) or x[par] > pow(10,4)){
			setParameters = 0;
			valueToReturn = pow(10,10);
		}
	}
	if (setParameters){
		this->SetParametersForRateMatrixForNelderMead(x, rateCat);	
		this->ComputeLogLikelihood();	
		valueToReturn = -1 * this->logLikelihood;
	}
	return (valueToReturn);
}


double rootedPhylogeny_tree::GetNegLogLikelihoodForParametersOfRootProb(double x[]){	
	// Set parameters for root prob
	// Check if constraints are being satisfied
	for (int par = 0; par < 3; par++){
		if (x[par] < 0.0 or x[par] > 1.0){
			return (pow(10,10));
		}
	}
	if ((x[0] + x[1] + x[3]) > 1.0){
		return (pow(10,10));
	}
	this->rootProbability[0] = x[0];
	this->rootProbability[1] = x[1];
	this->rootProbability[2] = x[2];
	this->rootProbability[3] = (1 - this->rootProbability[0] - this->rootProbability[1] - this->rootProbability[2]);
	this->ComputeLogLikelihood();	
	// return neg log likelihood
	return (-1 * this->logLikelihood);
}

void rootedPhylogeny_tree::OptimizeModelParametersUsingBFGS(){
	for (int rateCat = 0; rateCat < this->numberOfRateCategories; rateCat++){
		if (this->root->rateCategory == rateCat){
			rootedPhylogeny_vertex * c_l;
			rootedPhylogeny_vertex * c_r;
			c_l = (*this->vertexMap)[this->root->children_id[0]];
			c_r = (*this->vertexMap)[this->root->children_id[1]];
			if (c_l->rateCategory != rateCat and c_r->rateCategory != rateCat){
//				cout << "Case 1" << endl;
				// Construct a modified version of Nelder Mead for the following case
				// no. of vertices that are in the same rate cat as the root equals one		
				this->OptimizeModelParametersForRootProbUsingNelderMead();
			} else {
//				cout << "Case 2" << endl;
				this->ComputeLogLikelihood();
//				cout << endl;
//				cout << "edge length 1 " << endl;				
				this->OptimizeParametersForRateMatrixUsingBFGS(rateCat);
//				cout << "edge length 2 " << endl;
//				cout << this->GetEdgeLength(this->GetVertexId("h_root"),this->GetVertexId("h_1")) << endl;
				this->ComputeLogLikelihood();
//				cout << "updated loglikelihood 1 is " << this->logLikelihood << endl;
				this->OptimizeEdgeLengthsForRateCat(rateCat);				
				this->ComputeLogLikelihood();
//				cout << "updated loglikelihood 2 is " << this->logLikelihood << endl;
			}
		} else {
//			cout << "Case 3" << endl;
			this->ComputeLogLikelihood();
//			cout << "current loglikelihood is " << this->logLikelihood << endl;
			this->OptimizeParametersForRateMatrixUsingBFGS(rateCat);	
			this->ComputeLogLikelihood();
//			cout << "updated loglikelihood 1 is " << this->logLikelihood << endl;
			this->OptimizeEdgeLengthsForRateCat(rateCat);
			this->ComputeLogLikelihood();
//			cout << "updated loglikelihood 2 is " << this->logLikelihood << endl;
		}		
	}
	this->ComputeLogLikelihood();
}


void rootedPhylogeny_tree::OptimizeModelParametersUsingNelderMead(){
	double currentLogLikelihood = 0;
	double updatedLogLikelihood = 0;
	bool convergenceNotReached = 1;
	for (int rateCat = 0; rateCat < this->numberOfRateCategories; rateCat++){
		if (this->root->rateCategory == rateCat){
			rootedPhylogeny_vertex * c_l;
			rootedPhylogeny_vertex * c_r;
			c_l = (*this->vertexMap)[this->root->children_id[0]];
			c_r = (*this->vertexMap)[this->root->children_id[1]];
			if (c_l->rateCategory != rateCat and c_r->rateCategory != rateCat){
//				cout << "Optimizing root prob" << endl;
				// Construct a modified version of Nelder Mead for the following case
				// no. of vertices that are in the same rate cat as the root equals one		
				this->OptimizeModelParametersForRootProbUsingNelderMead();
			} else {				
				
//				cout << endl;
//				cout << "edge length 1 " << endl;				
//				cout << "Starting rate matrix optimization" << endl;				
//				cout << "Rate matrix optimization complete" << endl;
//				cout << "edge length 2 " << endl;
//				cout << this->GetEdgeLength(this->GetVertexId("h_root"),this->GetVertexId("h_1")) << endl;				
//				cout << "updated loglikelihood 1 is " << this->logLikelihood << endl;
//				cout << "Starting edge length optimization" << endl;
				convergenceNotReached = 1;
				this->ComputeLogLikelihood();				
				currentLogLikelihood = this->logLikelihood;
				while (convergenceNotReached){
					this->OptimizeParametersForRateMatrixUsingNelderMead(rateCat);
					this->OptimizeEdgeLengthsForRateCat(rateCat);								
					this->ComputeLogLikelihood();
					updatedLogLikelihood = this->logLikelihood;
					if (updatedLogLikelihood < currentLogLikelihood or abs(updatedLogLikelihood-currentLogLikelihood) < pow(10,-1)){
						convergenceNotReached = 0;
					}
					currentLogLikelihood = updatedLogLikelihood;
//					cout << "updated loglikelihood is " << this->logLikelihood << endl;
				}
				
//				cout << "Edge length optimization complete" << endl;
				
//				cout << "updated loglikelihood 2 is " << this->logLikelihood << endl;
			}
		} else {			
			convergenceNotReached = 1;
			this->ComputeLogLikelihood();				
			currentLogLikelihood = this->logLikelihood;
			while (convergenceNotReached){
				this->OptimizeParametersForRateMatrixUsingNelderMead(rateCat);
				this->OptimizeEdgeLengthsForRateCat(rateCat);								
				this->ComputeLogLikelihood();
				updatedLogLikelihood = this->logLikelihood;
				if (updatedLogLikelihood < currentLogLikelihood or abs(updatedLogLikelihood-currentLogLikelihood) < pow(10,-1)){
					convergenceNotReached = 0;
				}
				currentLogLikelihood = updatedLogLikelihood;
//				cout << "updated loglikelihood is " << this->logLikelihood << endl;
			}
		}		
	}
	this->ComputeLogLikelihood();
}

void rootedPhylogeny_tree::OptimizeModelParametersForRootProbUsingNewtonRaphson(){
	// Compute Jacobian for root prob parameters
	// Compute Hessian for root prob parameters
	// Compute search direction
	// Select step size using Brent's line search // impose bounds using penalty
	// Iterate till convergence
}

void rootedPhylogeny_tree::OptimizeModelParametersForRootProbUsingNelderMead(){
	int i;
	int icount;
	int ifault;
	int kcount;
	int konvge;
	int n;
	int numres;
	double reqmin;
	double *start;
	double *step;
	double *xmin;
	double ynewlo;
	
	n = 3;

	start = new double[n];
	step = new double[n];
	xmin = new double[n];
	
	// Initial estimate of parameters
	for (i = 0; i < n; i ++){
		start[i] = this->rootProbability[i];
	}
		
	reqmin = 1.0E-08;
	double stepSize = pow(10,-2);
    for (i = 0; i < n; i++){
		step[i] = stepSize;
	}
	
	konvge = 10;
	kcount = 500;
	
	ynewlo = this->GetNegLogLikelihoodForParametersOfRootProb( start );
		
	this->NelderMeadForOptimizingParametersOfRootProb( n, start, xmin, &ynewlo, reqmin, step,
	konvge, kcount, &icount, &numres, &ifault );
}

void rootedPhylogeny_tree::OptimizeParametersForRateMatrixUsingBFGS(int rateCat){	
	// Compute gradient
	// Initialize Hessian inverse	
	// Compute search direction vector
	// While norm of search direction vector is greater than pow(10,-4);
		// Select step size // impose bounds using a penalty
		// Update parameters
		// Compute gradient
		// Update Hessian inverse
		// Compute search direction vector
	
//	double convergenceTolerance = pow(10,-4);
//		MatrixXd InverseHessian_current = ArrayXXd::Zero(11,11);
//		MatrixXd InverseHessian_updated;
//		MatrixXd B_current; MatrixXd B_updated;
//		MatrixXd Jacobian_current; MatrixXd Jacobian_updated;
//		MatrixXd s;
//		MatrixXd y;
//		MatrixXd rho = ArrayXXd::Zero(11,11);;
//		MatrixXd IdentityMatrix = ArrayXXd::Zero(11,11);
//		MatrixXd TempMatrix_l = ArrayXXd::Zero(11,11);
//		MatrixXd TempMatrix_r = ArrayXXd::Zero(11,11);
//		int iter;
//		double normOfJacobian;
//		double stepSize;
//		this->rateCategoryForOptimization = rateCat;
//		
//		for (int par_1 = 0; par_1 < 11; par_1 ++){		
//			InverseHessian_current(par_1, par_1) = 1.0;
//			IdentityMatrix(par_1, par_1) = 1.0;
//		}
//		
//		B_current = (*this->parametersPerRateCategory)[rateCat];
//		Jacobian_current = this->GetJacobianForRateCategory(rateCat);
//		normOfJacobian = this->ComputeNormOfJacobian(Jacobian_current);
//		iter = 0;
////		cout << "norm of Jacobian for iteration " << iter << " is " << normOfJacobian << endl;
//		while (normOfJacobian > convergenceTolerance and iter < 100){
//			cout << "B_current is \n" << B_current << endl;
////			cout << "norm of Jacobian for iteration " << iter << " is " << normOfJacobian << endl;
//			iter += 1;
////			cout << "Jacobian current is " << endl;
////			cout << Jacobian_current << endl;
//			this->searchDirection = -1 * (InverseHessian_current * Jacobian_current);
////			cout << "search direction" << endl;
////			cout << this->searchDirection << endl;
////			cout << "B_current is \n" << B_current << endl;
//			stepSize = this->GetOptimalStepSizeUsingLineSearch(B_current);
//			cout << "step size is " << stepSize << endl;			
//			B_updated = B_current + stepSize * this->searchDirection;					
//			this->UpdateParametersForRateCategory(rateCat, B_updated);
//			Jacobian_updated = this->GetJacobianForRateCategory(rateCat);
//			y = Jacobian_updated - Jacobian_current;
//			s = B_updated - B_current;
////			cout << "y" << endl;
////			cout << y << endl;
////			cout << "s" << endl;
////			cout << s << endl;
//			rho = y.transpose()*s;
//			rho = rho.inverse();
////			cout << "rho is " << endl;
////			cout << rho << endl;
////			cout << s * y.transpose() << endl;
//			TempMatrix_l = IdentityMatrix - rho(0,0) * ( s * y.transpose() );
//			TempMatrix_r = IdentityMatrix - rho(0,0) * ( y * s.transpose() );
//			InverseHessian_updated = TempMatrix_l * InverseHessian_current;
//			InverseHessian_updated = InverseHessian_updated * TempMatrix_r;
//			InverseHessian_updated += rho(0,0) * ( s * s.transpose() );
//			B_current = B_updated;
//			InverseHessian_current = InverseHessian_updated;
//			normOfJacobian = this->ComputeNormOfJacobian(Jacobian_updated);
}

void rootedPhylogeny_tree::OptimizeParametersForRateMatrixUsingNelderMead(int rateCat){	
	int icount;
	int ifault;
	int kcount;
	int konvge;
	int n;
	int numres;
	double reqmin;
	double *start;
	double *step;
	double *xmin;
	double ynewlo;
	
	n = 11;

	start = new double[n];
	step = new double[n];
	xmin = new double[n];
	
	// Initial estimate of parameters
	for (int par = 0; par < 11; par ++){
		start[par] = (*this->parametersPerRateCategory)[rateCat](par,0);
	}
		
	reqmin = 1.0E-08;
	double stepSize = pow(10,-2);
    for (int par = 0; par < 11; par++){
		step[par] = stepSize;
	}
	
	konvge = 10;
	kcount = 500;
	
	ynewlo = this->GetNegLogLikelihoodForRateParametersForRateCat( start , rateCat);
	
	this->NelderMeadForOptimizingRateParametersForRateCat( rateCat, n, start, xmin, &ynewlo, reqmin, step,
	konvge, kcount, &icount, &numres, &ifault );
}

void rootedPhylogeny_tree::NelderMeadForOptimizingParametersOfRootProb( int n, double start[], double xmin[], 
  double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
  int *icount, int *numres, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    NELMIN minimizes a function using the Nelder-Mead algorithm.
//
//  Discussion:
//
//    This routine seeks the minimum value of a user-specified function.
//
//    Simplex function minimisation procedure due to Nelder+Mead(1965),
//    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
//    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
//    25, 97) and Hill(1978, 27, 380-2)
//
//    The function to be minimized must be defined by a function of
//    the form
//
//      function fn ( x, f )
//      double fn
//      double x(*)
//
//    and the name of this subroutine must be declared EXTERNAL in the
//    calling routine and passed as the argument FN.
//
//    This routine does not include a termination test using the
//    fitting of a quadratic surface.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by R ONeill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Nelder, Roger Mead,
//    A simplex method for function minimization,
//    Computer Journal,
//    Volume 7, 1965, pages 308-313.
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double FN ( double x[] ), the name of the routine which evaluates
//    the function to be minimized.
//
//    Input, int N, the number of variables.
//
//    Input/output, double START[N].  On input, a starting point
//    for the iteration.  On output, this data may have been overwritten.
//
//    Output, double XMIN[N], the coordinates of the point which
//    is estimated to minimize the function.
//
//    Output, double YNEWLO, the minimum value of the function.
//
//    Input, double REQMIN, the terminating limit for the variance
//    of function values.
//
//    Input, double STEP[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int KONVGE, the convergence check is carried out 
//    every KONVGE iterations.
//
//    Input, int KCOUNT, the maximum number of function 
//    evaluations.
//
//    Output, int *ICOUNT, the number of function evaluations 
//    used.
//
//    Output, int *NUMRES, the number of restarts.
//
//    Output, int *IFAULT, error indicator.
//    0, no errors detected.
//    1, REQMIN, N, or KONVGE has an illegal value.
//    2, iteration terminated because KCOUNT was exceeded without convergence.
//
{
  double ccoeff = 0.5;
  double del;
  double dn;
  double dnn;
  double ecoeff = 2.0;
  double eps = 0.001;
  int i;
  int ihi;
  int ilo;
  int j;
  int jcount;
  int l;
  int nn;
  double *p;
  double *p2star;
  double *pbar;
  double *pstar;
  double rcoeff = 1.0;
  double rq;
  double x;
  double *y;
  double y2star;
  double ylo;
  double ystar;
  double z;
//
//  Check the input parameters.
//
  if ( reqmin <= 0.0 )
  {
    *ifault = 1;
    return;
  }

  if ( n < 1 )
  {
    *ifault = 1;
    return;
  }

  if ( konvge < 1 )
  {
    *ifault = 1;
    return;
  }

  p = new double[n*(n+1)];
  pstar = new double[n];
  p2star = new double[n];
  pbar = new double[n];
  y = new double[n+1];

  *icount = 0;
  *numres = 0;

  jcount = konvge; 
  dn = ( double ) ( n );
  nn = n + 1;
  dnn = ( double ) ( nn );
  del = 1.0;
  rq = reqmin * dn;
//
//  Initial or restarted loop.
//
  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    { 
      p[i+n*n] = start[i];
    }    
	y[n] = this->GetNegLogLikelihoodForParametersOfRootProb( start );
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = start[j];
      start[j] = start[j] + step[j] * del;
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = start[i];
      }
      y[j] = this->GetNegLogLikelihoodForParametersOfRootProb( start );
      *icount = *icount + 1;
      start[j] = x;
    }
//                    
//  The simplex construction is complete.
//                    
//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
//  the vertex of the simplex to be replaced.
//                
    ylo = y[0];
    ilo = 0;

    for ( i = 1; i < nn; i++ )
    {
      if ( y[i] < ylo )
      {
        ylo = y[i];
        ilo = i;
      }
    }
//
//  Inner loop.
//
    for ( ; ; )
    {
      if ( kcount <= *icount )
      {
        break;
      }
      *ynewlo = y[0];
      ihi = 0;

      for ( i = 1; i < nn; i++ )
      {
        if ( *ynewlo < y[i] )
        {
          *ynewlo = y[i];
          ihi = i;
        }
      }
//
//  Calculate PBAR, the centroid of the simplex vertices
//  excepting the vertex with Y value YNEWLO.
//
      for ( i = 0; i < n; i++ )
      {
        z = 0.0;
        for ( j = 0; j < nn; j++ )
        { 
          z = z + p[i+j*n];
        }
        z = z - p[i+ihi*n];  
        pbar[i] = z / dn;
      }
//
//  Reflection through the centroid.
//
      for ( i = 0; i < n; i++ )
      {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
      }
      ystar = this->GetNegLogLikelihoodForParametersOfRootProb( pstar );
      *icount = *icount + 1;
//
//  Successful reflection, so extension.
//
      if ( ystar < ylo )
      {
        for ( i = 0; i < n; i++ )
        {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
        }
		y2star = this->GetNegLogLikelihoodForParametersOfRootProb( p2star );
//        y2star = this->powell ( p2star );
        *icount = *icount + 1;
//
//  Check extension.
//
        if ( ystar < y2star )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Retain extension or contraction.
//
        else
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = p2star[i];
          }
          y[ihi] = y2star;
        }
      }
//
//  No extension.
//
      else
      {
        l = 0;
        for ( i = 0; i < nn; i++ )
        {
          if ( ystar < y[i] )
          {
            l = l + 1;
          }
        }

        if ( 1 < l )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Contraction on the Y(IHI) side of the centroid.
//
        else if ( l == 0 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
          }
		  y2star = this->GetNegLogLikelihoodForParametersOfRootProb( p2star );
//          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Contract the whole simplex.
//
          if ( y[ihi] < y2star )
          {
            for ( j = 0; j < nn; j++ )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                xmin[i] = p[i+j*n];
              }
			  y[j] = this->GetNegLogLikelihoodForParametersOfRootProb( xmin );
//              y[j] = this->powell ( xmin );
              *icount = *icount + 1;
            }
            ylo = y[0];
            ilo = 0;

            for ( i = 1; i < nn; i++ )
            {
              if ( y[i] < ylo )
              {
                ylo = y[i];
                ilo = i;
              }
            }
            continue;
          }
//
//  Retain contraction.
//
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
//
//  Contraction on the reflection side of the centroid.
//
        else if ( l == 1 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
          }
		  y2star = this->GetNegLogLikelihoodForParametersOfRootProb( p2star );
//          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Retain reflection?
//
          if ( y2star <= ystar )
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = pstar[i];
            }
            y[ihi] = ystar;
          }
        }
      }
//
//  Check if YLO improved.
//
      if ( y[ihi] < ylo )
      {
        ylo = y[ihi];
        ilo = ihi;
      }
      jcount = jcount - 1;

      if ( 0 < jcount )
      {
        continue;
      }
//
//  Check to see if minimum reached.
//
      if ( *icount <= kcount )
      {
        jcount = konvge;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + y[i];
        }
        x = z / dnn;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + pow ( y[i] - x, 2 );
        }

        if ( z <= rq )
        {
          break;
        }
      }
    }
//
//  Factorial tests to check that YNEWLO is a local minimum.
//
    for ( i = 0; i < n; i++ )
    {
      xmin[i] = p[i+ilo*n];
    }
    *ynewlo = y[ilo];

    if ( kcount < *icount )
    {
      *ifault = 2;
      break;
    }

    *ifault = 0;

    for ( i = 0; i < n; i++ )
    {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
	  z = this->GetNegLogLikelihoodForParametersOfRootProb( xmin );
//      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
	  z = this->GetNegLogLikelihoodForParametersOfRootProb( xmin );
//      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if ( *ifault == 0 )
    {
      break;
    }
//
//  Restart the procedure.
//
    for ( i = 0; i < n; i++ )
    {
      start[i] = xmin[i];
    }
    del = eps;
    *numres = *numres + 1;
  }
  delete [] p;
  delete [] pstar;
  delete [] p2star;
  delete [] pbar;
  delete [] y;

  return;
}

void rootedPhylogeny_tree::NelderMeadForOptimizingRateParametersForRateCat( int rateCat, int n, double start[], double xmin[], 
  double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
  int *icount, int *numres, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    NELMIN minimizes a function using the Nelder-Mead algorithm.
//
//  Discussion:
//
//    This routine seeks the minimum value of a user-specified function.
//
//    Simplex function minimisation procedure due to Nelder+Mead(1965),
//    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
//    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
//    25, 97) and Hill(1978, 27, 380-2)
//
//    The function to be minimized must be defined by a function of
//    the form
//
//      function fn ( x, f )
//      double fn
//      double x(*)
//
//    and the name of this subroutine must be declared EXTERNAL in the
//    calling routine and passed as the argument FN.
//
//    This routine does not include a termination test using the
//    fitting of a quadratic surface.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by R ONeill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Nelder, Roger Mead,
//    A simplex method for function minimization,
//    Computer Journal,
//    Volume 7, 1965, pages 308-313.
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double FN ( double x[] ), the name of the routine which evaluates
//    the function to be minimized.
//
//    Input, int N, the number of variables.
//
//    Input/output, double START[N].  On input, a starting point
//    for the iteration.  On output, this data may have been overwritten.
//
//    Output, double XMIN[N], the coordinates of the point which
//    is estimated to minimize the function.
//
//    Output, double YNEWLO, the minimum value of the function.
//
//    Input, double REQMIN, the terminating limit for the variance
//    of function values.
//
//    Input, double STEP[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int KONVGE, the convergence check is carried out 
//    every KONVGE iterations.
//
//    Input, int KCOUNT, the maximum number of function 
//    evaluations.
//
//    Output, int *ICOUNT, the number of function evaluations 
//    used.
//
//    Output, int *NUMRES, the number of restarts.
//
//    Output, int *IFAULT, error indicator.
//    0, no errors detected.
//    1, REQMIN, N, or KONVGE has an illegal value.
//    2, iteration terminated because KCOUNT was exceeded without convergence.
//
{
  double ccoeff = 0.5;
  double del;
  double dn;
  double dnn;
  double ecoeff = 2.0;
  double eps = 0.001;
  int i;
  int ihi;
  int ilo;
  int j;
  int jcount;
  int l;
  int nn;
  double *p;
  double *p2star;
  double *pbar;
  double *pstar;
  double rcoeff = 1.0;
  double rq;
  double x;
  double *y;
  double y2star;
  double ylo;
  double ystar;
  double z;
//
//  Check the input parameters.
//
  if ( reqmin <= 0.0 )
  {
    *ifault = 1;
    return;
  }

  if ( n < 1 )
  {
    *ifault = 1;
    return;
  }

  if ( konvge < 1 )
  {
    *ifault = 1;
    return;
  }

  p = new double[n*(n+1)];
  pstar = new double[n];
  p2star = new double[n];
  pbar = new double[n];
  y = new double[n+1];

  *icount = 0;
  *numres = 0;

  jcount = konvge; 
  dn = ( double ) ( n );
  nn = n + 1;
  dnn = ( double ) ( nn );
  del = 1.0;
  rq = reqmin * dn;
//
//  Initial or restarted loop.
//
  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    { 
      p[i+n*n] = start[i];
    }
    y[n] = this->GetNegLogLikelihoodForRateParametersForRateCat( start , rateCat);
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = start[j];
      start[j] = start[j] + step[j] * del;
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = start[i];
      }
      y[j] = this->GetNegLogLikelihoodForRateParametersForRateCat( start , rateCat);
      *icount = *icount + 1;
      start[j] = x;
    }
//                    
//  The simplex construction is complete.
//                    
//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
//  the vertex of the simplex to be replaced.
//                
    ylo = y[0];
    ilo = 0;

    for ( i = 1; i < nn; i++ )
    {
      if ( y[i] < ylo )
      {
        ylo = y[i];
        ilo = i;
      }
    }
//
//  Inner loop.
//
    for ( ; ; )
    {
      if ( kcount <= *icount )
      {
        break;
      }
      *ynewlo = y[0];
      ihi = 0;

      for ( i = 1; i < nn; i++ )
      {
        if ( *ynewlo < y[i] )
        {
          *ynewlo = y[i];
          ihi = i;
        }
      }
//
//  Calculate PBAR, the centroid of the simplex vertices
//  excepting the vertex with Y value YNEWLO.
//
      for ( i = 0; i < n; i++ )
      {
        z = 0.0;
        for ( j = 0; j < nn; j++ )
        { 
          z = z + p[i+j*n];
        }
        z = z - p[i+ihi*n];  
        pbar[i] = z / dn;
      }
//
//  Reflection through the centroid.
//
      for ( i = 0; i < n; i++ )
      {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
      }
      ystar = this->GetNegLogLikelihoodForRateParametersForRateCat( pstar , rateCat);;
      *icount = *icount + 1;
//
//  Successful reflection, so extension.
//
      if ( ystar < ylo )
      {
        for ( i = 0; i < n; i++ )
        {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
        }
		y2star = this->GetNegLogLikelihoodForRateParametersForRateCat( p2star , rateCat);
//        y2star = this->powell ( p2star );
        *icount = *icount + 1;
//
//  Check extension.
//
        if ( ystar < y2star )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Retain extension or contraction.
//
        else
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = p2star[i];
          }
          y[ihi] = y2star;
        }
      }
//
//  No extension.
//
      else
      {
        l = 0;
        for ( i = 0; i < nn; i++ )
        {
          if ( ystar < y[i] )
          {
            l = l + 1;
          }
        }

        if ( 1 < l )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Contraction on the Y(IHI) side of the centroid.
//
        else if ( l == 0 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
          }
		  y2star = this->GetNegLogLikelihoodForRateParametersForRateCat( p2star , rateCat);
//          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Contract the whole simplex.
//
          if ( y[ihi] < y2star )
          {
            for ( j = 0; j < nn; j++ )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                xmin[i] = p[i+j*n];
              }
			  y[j] = this->GetNegLogLikelihoodForRateParametersForRateCat( xmin , rateCat);
//              y[j] = this->powell ( xmin );
              *icount = *icount + 1;
            }
            ylo = y[0];
            ilo = 0;

            for ( i = 1; i < nn; i++ )
            {
              if ( y[i] < ylo )
              {
                ylo = y[i];
                ilo = i;
              }
            }
            continue;
          }
//
//  Retain contraction.
//
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
//
//  Contraction on the reflection side of the centroid.
//
        else if ( l == 1 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
          }
		  y2star = this->GetNegLogLikelihoodForRateParametersForRateCat( p2star , rateCat);
//          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Retain reflection?
//
          if ( y2star <= ystar )
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = pstar[i];
            }
            y[ihi] = ystar;
          }
        }
      }
//
//  Check if YLO improved.
//
      if ( y[ihi] < ylo )
      {
        ylo = y[ihi];
        ilo = ihi;
      }
      jcount = jcount - 1;

      if ( 0 < jcount )
      {
        continue;
      }
//
//  Check to see if minimum reached.
//
      if ( *icount <= kcount )
      {
        jcount = konvge;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + y[i];
        }
        x = z / dnn;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + pow ( y[i] - x, 2 );
        }

        if ( z <= rq )
        {
          break;
        }
      }
    }
//
//  Factorial tests to check that YNEWLO is a local minimum.
//
    for ( i = 0; i < n; i++ )
    {
      xmin[i] = p[i+ilo*n];
    }
    *ynewlo = y[ilo];

    if ( kcount < *icount )
    {
      *ifault = 2;
      break;
    }

    *ifault = 0;

    for ( i = 0; i < n; i++ )
    {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
	  z = this->GetNegLogLikelihoodForRateParametersForRateCat( xmin , rateCat);
//      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
	  z = this->GetNegLogLikelihoodForRateParametersForRateCat( xmin , rateCat);
//      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if ( *ifault == 0 )
    {
      break;
    }
//
//  Restart the procedure.
//
    for ( i = 0; i < n; i++ )
    {
      start[i] = xmin[i];
    }
    del = eps;
    *numres = *numres + 1;
  }
  delete [] p;
  delete [] pstar;
  delete [] p2star;
  delete [] pbar;
  delete [] y;

  return;
}

void rootedPhylogeny_tree::NelderMeadForPowell( int n, double start[], double xmin[], 
  double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
  int *icount, int *numres, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    NELMIN minimizes a function using the Nelder-Mead algorithm.
//
//  Discussion:
//
//    This routine seeks the minimum value of a user-specified function.
//
//    Simplex function minimisation procedure due to Nelder+Mead(1965),
//    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
//    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
//    25, 97) and Hill(1978, 27, 380-2)
//
//    The function to be minimized must be defined by a function of
//    the form
//
//      function fn ( x, f )
//      double fn
//      double x(*)
//
//    and the name of this subroutine must be declared EXTERNAL in the
//    calling routine and passed as the argument FN.
//
//    This routine does not include a termination test using the
//    fitting of a quadratic surface.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by R ONeill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Nelder, Roger Mead,
//    A simplex method for function minimization,
//    Computer Journal,
//    Volume 7, 1965, pages 308-313.
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double FN ( double x[] ), the name of the routine which evaluates
//    the function to be minimized.
//
//    Input, int N, the number of variables.
//
//    Input/output, double START[N].  On input, a starting point
//    for the iteration.  On output, this data may have been overwritten.
//
//    Output, double XMIN[N], the coordinates of the point which
//    is estimated to minimize the function.
//
//    Output, double YNEWLO, the minimum value of the function.
//
//    Input, double REQMIN, the terminating limit for the variance
//    of function values.
//
//    Input, double STEP[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int KONVGE, the convergence check is carried out 
//    every KONVGE iterations.
//
//    Input, int KCOUNT, the maximum number of function 
//    evaluations.
//
//    Output, int *ICOUNT, the number of function evaluations 
//    used.
//
//    Output, int *NUMRES, the number of restarts.
//
//    Output, int *IFAULT, error indicator.
//    0, no errors detected.
//    1, REQMIN, N, or KONVGE has an illegal value.
//    2, iteration terminated because KCOUNT was exceeded without convergence.
//
{
  double ccoeff = 0.5;
  double del;
  double dn;
  double dnn;
  double ecoeff = 2.0;
  double eps = 0.001;
  int i;
  int ihi;
  int ilo;
  int j;
  int jcount;
  int l;
  int nn;
  double *p;
  double *p2star;
  double *pbar;
  double *pstar;
  double rcoeff = 1.0;
  double rq;
  double x;
  double *y;
  double y2star;
  double ylo;
  double ystar;
  double z;
//
//  Check the input parameters.
//
  if ( reqmin <= 0.0 )
  {
    *ifault = 1;
    return;
  }

  if ( n < 1 )
  {
    *ifault = 1;
    return;
  }

  if ( konvge < 1 )
  {
    *ifault = 1;
    return;
  }

  p = new double[n*(n+1)];
  pstar = new double[n];
  p2star = new double[n];
  pbar = new double[n];
  y = new double[n+1];

  *icount = 0;
  *numres = 0;

  jcount = konvge; 
  dn = ( double ) ( n );
  nn = n + 1;
  dnn = ( double ) ( nn );
  del = 1.0;
  rq = reqmin * dn;
//
//  Initial or restarted loop.
//
  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    { 
      p[i+n*n] = start[i];
    }
    y[n] = this->powell( start );
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = start[j];
      start[j] = start[j] + step[j] * del;
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = start[i];
      }
      y[j] = this->powell ( start );
      *icount = *icount + 1;
      start[j] = x;
    }
//                    
//  The simplex construction is complete.
//                    
//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
//  the vertex of the simplex to be replaced.
//                
    ylo = y[0];
    ilo = 0;

    for ( i = 1; i < nn; i++ )
    {
      if ( y[i] < ylo )
      {
        ylo = y[i];
        ilo = i;
      }
    }
//
//  Inner loop.
//
    for ( ; ; )
    {
      if ( kcount <= *icount )
      {
        break;
      }
      *ynewlo = y[0];
      ihi = 0;

      for ( i = 1; i < nn; i++ )
      {
        if ( *ynewlo < y[i] )
        {
          *ynewlo = y[i];
          ihi = i;
        }
      }
//
//  Calculate PBAR, the centroid of the simplex vertices
//  excepting the vertex with Y value YNEWLO.
//
      for ( i = 0; i < n; i++ )
      {
        z = 0.0;
        for ( j = 0; j < nn; j++ )
        { 
          z = z + p[i+j*n];
        }
        z = z - p[i+ihi*n];  
        pbar[i] = z / dn;
      }
//
//  Reflection through the centroid.
//
      for ( i = 0; i < n; i++ )
      {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
      }
      ystar = this->powell ( pstar );
      *icount = *icount + 1;
//
//  Successful reflection, so extension.
//
      if ( ystar < ylo )
      {
        for ( i = 0; i < n; i++ )
        {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
        }
        y2star = this->powell ( p2star );
        *icount = *icount + 1;
//
//  Check extension.
//
        if ( ystar < y2star )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Retain extension or contraction.
//
        else
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = p2star[i];
          }
          y[ihi] = y2star;
        }
      }
//
//  No extension.
//
      else
      {
        l = 0;
        for ( i = 0; i < nn; i++ )
        {
          if ( ystar < y[i] )
          {
            l = l + 1;
          }
        }

        if ( 1 < l )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Contraction on the Y(IHI) side of the centroid.
//
        else if ( l == 0 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
          }
          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Contract the whole simplex.
//
          if ( y[ihi] < y2star )
          {
            for ( j = 0; j < nn; j++ )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                xmin[i] = p[i+j*n];
              }
              y[j] = this->powell ( xmin );
              *icount = *icount + 1;
            }
            ylo = y[0];
            ilo = 0;

            for ( i = 1; i < nn; i++ )
            {
              if ( y[i] < ylo )
              {
                ylo = y[i];
                ilo = i;
              }
            }
            continue;
          }
//
//  Retain contraction.
//
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
//
//  Contraction on the reflection side of the centroid.
//
        else if ( l == 1 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
          }
          y2star = this->powell ( p2star );
          *icount = *icount + 1;
//
//  Retain reflection?
//
          if ( y2star <= ystar )
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = pstar[i];
            }
            y[ihi] = ystar;
          }
        }
      }
//
//  Check if YLO improved.
//
      if ( y[ihi] < ylo )
      {
        ylo = y[ihi];
        ilo = ihi;
      }
      jcount = jcount - 1;

      if ( 0 < jcount )
      {
        continue;
      }
//
//  Check to see if minimum reached.
//
      if ( *icount <= kcount )
      {
        jcount = konvge;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + y[i];
        }
        x = z / dnn;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + pow ( y[i] - x, 2 );
        }

        if ( z <= rq )
        {
          break;
        }
      }
    }
//
//  Factorial tests to check that YNEWLO is a local minimum.
//
    for ( i = 0; i < n; i++ )
    {
      xmin[i] = p[i+ilo*n];
    }
    *ynewlo = y[ilo];

    if ( kcount < *icount )
    {
      *ifault = 2;
      break;
    }

    *ifault = 0;

    for ( i = 0; i < n; i++ )
    {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
      z = this->powell ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if ( *ifault == 0 )
    {
      break;
    }
//
//  Restart the procedure.
//
    for ( i = 0; i < n; i++ )
    {
      start[i] = xmin[i];
    }
    del = eps;
    *numres = *numres + 1;
  }
  delete [] p;
  delete [] pstar;
  delete [] p2star;
  delete [] pbar;
  delete [] y;

  return;
}

void rootedPhylogeny_tree::testNelderMead_1(){
	int i;
	int icount;
	int ifault;
	int kcount;
	int konvge;
	int n;
	int numres;
	double reqmin;
	double *start;
	double *step;
	double *xmin;
	double ynewlo;
	
	n = 4;

	start = new double[n];
	step = new double[n];
	xmin = new double[n];

	cout << "\n";
	cout << "TEST02\n";
	cout << "  Apply NELMIN to POWELL quartic function.\n";

	start[0] =   3.0;
	start[1] = - 1.0;
	start[2] =   0.0;
	start[3] =   1.0;

	reqmin = 1.0E-08;

	step[0] = 1.0;
	step[1] = 1.0;
	step[2] = 1.0;
	step[3] = 1.0;

	konvge = 10;
	kcount = 500;

	cout << "\n";
	cout << "  Starting point X:\n";
	cout << "\n";
	
	for ( i = 0; i < n; i++ )
	{
	cout << "  " << setw(14) << start[i] << "\n";
	}
	
	ynewlo = powell ( start );

	cout << "\n";
	cout << "  F(X) = " << ynewlo << "\n";

	nelmin ( powell, n, start, xmin, &ynewlo, reqmin, step,
	konvge, kcount, &icount, &numres, &ifault );

	cout << "\n";
	cout << "  Return code IFAULT = " << ifault << "\n";
	cout << "\n";
	cout << "  Estimate of minimizing value X*:\n";
	cout << "\n";
	for ( i = 0; i < n; i++ )
	{
	cout << "  " << setw(14) << xmin[i] << "\n";
	}

	cout << "\n";
	cout << "  F(X*) = " << ynewlo << "\n";

	cout << "\n";
	cout << "  Number of iterations = " << icount << "\n";
	cout << "  Number of restarts =   " << numres << "\n";

	delete [] start;
	delete [] step;
	delete [] xmin;
	
	return ;
}

void rootedPhylogeny_tree::testNelderMead_2(){
	int i;
	int icount;
	int ifault;
	int kcount;
	int konvge;
	int n;
	int numres;
	double reqmin;
	double *start;
	double *step;
	double *xmin;
	double ynewlo;
	
	n = 4;

	start = new double[n];
	step = new double[n];
	xmin = new double[n];

	cout << "\n";
	cout << "TEST02\n";
	cout << "  Apply NELMIN to POWELL quartic function.\n";

	start[0] =   3.0;
	start[1] = - 1.0;
	start[2] =   0.0;
	start[3] =   1.0;

	reqmin = 1.0E-08;

	step[0] = 1.0;
	step[1] = 1.0;
	step[2] = 1.0;
	step[3] = 1.0;

	konvge = 10;
	kcount = 500;

	cout << "\n";
	cout << "  Starting point X:\n";
	cout << "\n";
	
	for ( i = 0; i < n; i++ )
	{
	cout << "  " << setw(14) << start[i] << "\n";
	}
	
	ynewlo = this->powell ( start );

	cout << "\n";
	cout << "  F(X) = " << ynewlo << "\n";
	
	this->NelderMeadForPowell( n, start, xmin, &ynewlo, reqmin, step,
	konvge, kcount, &icount, &numres, &ifault );

	cout << "\n";
	cout << "  Return code IFAULT = " << ifault << "\n";
	cout << "\n";
	cout << "  Estimate of minimizing value X*:\n";
	cout << "\n";
	
	for ( i = 0; i < n; i++ )
	{
	cout << "  " << setw(14) << xmin[i] << "\n";
	}

	cout << "\n";
	cout << "  F(X*) = " << ynewlo << "\n";

	cout << "\n";
	cout << "  Number of iterations = " << icount << "\n";
	cout << "  Number of restarts =   " << numres << "\n";

	delete [] start;
	delete [] step;
	delete [] xmin;
	
	return ;
}


double rootedPhylogeny_tree::powell(double x[4]){
	double fx;
	double fx1;
	double fx2;
	double fx3;
	double fx4;

	fx1 = x[0] + 10.0 * x[1];
	fx2 = x[2] - x[3];
	fx3 = x[1] - 2.0 * x[2];
	fx4 = x[0] - x[3];

	fx =            fx1 * fx1
	 +  5.0 * fx2 * fx2
	 +            fx3 * fx3 * fx3 * fx3
	 + 10.0 * fx4 * fx4 * fx4 * fx4;;

	return fx;
}

void rootedPhylogeny_tree::OptimizeEdgeLengthsForRateCat(int rateCat) {
	map <int, array <double, 4>> conditionalLikelihoodMap;
	array <double, 4> conditionalLikelihood;
	map <int, array <double, 4>> firstDerivativeOfConditionalLikelihoodMap;
	map <int, array <double, 4>> secondDerivativeOfConditionalLikelihoodMap;
	array <double, 4> firstDerivativeOfConditionalLikelihood;
	array <double, 4> secondDerivativeOfConditionalLikelihood;
	double firstDerivativeOfpartialLikelihood;
	double secondDerivativeOfpartialLikelihood;
	double partialLikelihood;
	double siteLikelihood;
	double firstDerivativeOfSiteLikelihood;
	double secondDerivativeOfSiteLikelihood;
	double firstDerivativeOfLogLikelihood;
	double secondDerivativeOfLogLikelihood;
	this->logLikelihood = 0;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c = (*this->vertexMap)[this->vertexMap->size()-1];
	rootedPhylogeny_vertex * a;
	rootedPhylogeny_vertex * v;
	rootedPhylogeny_vertex * w;
	float stepSize;
	Matrix4f Q; Matrix4f Q_norm; Matrix4f Q_scaled; float t; Matrix4f P;
	Matrix4f Q_ab; Matrix4f Q_ab_scaled; Matrix4f Q_ab_norm; float t_ab;
	Matrix4f Q_ac; Matrix4f Q_ac_scaled; Matrix4f Q_ac_norm; float t_ac;
	Matrix4f Q_uv; Matrix4f Q_uv_scaled; Matrix4f Q_uv_norm; float t_uv;
	Matrix4f Q_uw; Matrix4f Q_uw_scaled; Matrix4f Q_uw_norm; float t_uw;
	Matrix4f P_uv; Matrix4f P_uw;
	Matrix4f P_ab; Matrix4f P_ac; Matrix4f P_ab_der_1; Matrix4f P_ab_der_2;
	vector <rootedPhylogeny_vertex*> verticesForRateCat;
	for (rootedPhylogeny_vertex* b: *this->verticesForPreOrderTreeTraversal) {
		if (b->id != -1 and b->rateCategory == rateCat) {
			verticesForRateCat.push_back(b);
		}
	}
	bool child_in_path;
	vector <rootedPhylogeny_vertex *> verticesInPathFromRootToANotIncA;
	// Iterate over vertices
	bool convergenceNotReached = 1;
	for (rootedPhylogeny_vertex * b : verticesForRateCat) {
		convergenceNotReached = 1;
		a = (*this->vertexMap)[b->parent_id];
//		cout << "optimizing length for edge " << a->name << " - " << b->name << endl;
//		cout << "current edge length is " << this->GetEdgeLength(a->id,b->id) << endl;
		// Length of edge (a, b) will be optimized
		verticesInPathFromRootToANotIncA.clear();
		while (a->id != -1 and a->parent_id != -1) {
			a = (*this->vertexMap)[a->parent_id];
			verticesInPathFromRootToANotIncA.push_back(a);
		}
		a = (*this->vertexMap)[b->parent_id];
		if (a->id != -1) {
			verticesInPathFromRootToANotIncA.push_back(this->root);
		}
		// Iterate till convergence
		while (convergenceNotReached) {
			// Iterate over sites
			firstDerivativeOfLogLikelihood = 0;
			secondDerivativeOfLogLikelihood = 0;
			for (unsigned int site = 0; site < this->siteWeights.size(); site++){									
				// Compute conditional likelihoods
				conditionalLikelihoodMap.clear();			
				for (pair <int,int> edgeIds : *this->edgesForPostOrderTreeTraversal){					
					p = (*this->vertexMap)[edgeIds.first];
					c = (*this->vertexMap)[edgeIds.second];
					Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];
					t = this->GetEdgeLength(p->id,c->id);
					Q_norm = Q/(*this->scalingFactorForRateCategory)[c->rateCategory];
					Q_scaled = Q_norm * t;
					P = Q_scaled.exp();
					// Initialize conditional likelihood for leaves
					if (c->children_id.size()==0) {
						for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {
							conditionalLikelihood[dna_c] = 0.0;
						}
						conditionalLikelihood[c->compressedSequence[site]] = 1.0;
						conditionalLikelihoodMap.insert(pair<int,array<double,4>>(c->id,conditionalLikelihood));
					}
					// Initialize conditional likelihood for ancestors
					if (conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()) {
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++) {
						conditionalLikelihood[dna_c] = 1.0;
						}				
						conditionalLikelihoodMap.insert(pair<int,array<double,4>>(p->id,conditionalLikelihood));					
					}		
					for (unsigned char dna_p = 0; dna_p < 4; dna_p++) {
						partialLikelihood = 0.0;
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++) {
							partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
						}
						conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;
					}
				}		
				// Compute first and second derivative of conditional likelihood
				firstDerivativeOfConditionalLikelihoodMap.clear();
				secondDerivativeOfConditionalLikelihoodMap.clear();
				// Reset a and c
				a = (*this->vertexMap)[b->parent_id];
				for (int child_id : a->children_id) {
					if (child_id != b->id) {
						c = (*this->vertexMap)[child_id];
					}
				}
				// Initialize first and second derivative
				for (int dna = 0; dna < 4; dna ++) {
					firstDerivativeOfConditionalLikelihood[dna] = 1.0;
					secondDerivativeOfConditionalLikelihood[dna]  = 1.0;
				}
				
				for (int dna_p = 0; dna_p < 4; dna_p ++) {				
					// child b. length of edge (a, b) is being optimized.
					Q_ab = (*this->rateMatrixPerRateCategory)[b->rateCategory];
					Q_ab_norm = Q_ab/(*this->scalingFactorForRateCategory)[b->rateCategory];
					t_ab = this->GetEdgeLength(a->id, b->id);
					Q_ab_scaled = Q_ab_norm * t_ab;
					P_ab = Q_ab_scaled.exp();
					P_ab_der_1 = Q_ab_norm * P_ab;
					P_ab_der_2 = Q_ab_norm * P_ab_der_1;
//					cout << P_ab_der_2 << endl;
					firstDerivativeOfpartialLikelihood = 0;
					secondDerivativeOfpartialLikelihood = 0;
					for (int dna_c = 0; dna_c < 4; dna_c ++) {
						firstDerivativeOfpartialLikelihood += P_ab_der_1(dna_p,dna_c)*conditionalLikelihoodMap[b->id][dna_c];
						secondDerivativeOfpartialLikelihood += P_ab_der_2(dna_p,dna_c)*conditionalLikelihoodMap[b->id][dna_c];
					}
					firstDerivativeOfConditionalLikelihood[dna_p] *= firstDerivativeOfpartialLikelihood;				
					secondDerivativeOfConditionalLikelihood[dna_p] *= secondDerivativeOfpartialLikelihood;				
//					cout << "second der of partial likelihood is " << secondDerivativeOfpartialLikelihood << endl;
					// child c
					Q_ac = (*this->rateMatrixPerRateCategory)[c->rateCategory];
					Q_ac_norm = Q_ac/(*this->scalingFactorForRateCategory)[c->rateCategory];
					t_ac = this->GetEdgeLength(a->id, c->id);
					Q_ac_scaled = Q_ac_norm * t_ac;
					P_ac = Q_ac_scaled.exp();
					partialLikelihood = 0;
					for (int dna_c = 0; dna_c < 4; dna_c ++) {
						partialLikelihood += P_ac(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
					}
					firstDerivativeOfConditionalLikelihood[dna_p] *= partialLikelihood;
					secondDerivativeOfConditionalLikelihood[dna_p] *= partialLikelihood;			
				}
				firstDerivativeOfConditionalLikelihoodMap.insert(pair<int,array<double,4>>(a->id,firstDerivativeOfConditionalLikelihood));
				secondDerivativeOfConditionalLikelihoodMap.insert(pair<int,array<double,4>>(a->id,secondDerivativeOfConditionalLikelihood));
				// Compute the first derivative for ancestors of a
				for (rootedPhylogeny_vertex * u : verticesInPathFromRootToANotIncA) {
					for (int dna = 0; dna < 4; dna ++) {
						firstDerivativeOfConditionalLikelihood[dna] = 1.0;
						secondDerivativeOfConditionalLikelihood[dna] = 1.0;
					}
					for (int dna_p = 0; dna_p < 4; dna_p++) {
						for (int child_id : u->children_id) {
							child_in_path = find(verticesInPathFromRootToANotIncA.begin(),verticesInPathFromRootToANotIncA.end(),(*this->vertexMap)[child_id]) != verticesInPathFromRootToANotIncA.end();
							if (child_id == a->id or child_in_path) {
								v = (*this->vertexMap)[child_id];
								Q_uv = (*this->rateMatrixPerRateCategory)[v->rateCategory];
								Q_uv_norm = Q_uv/(*this->scalingFactorForRateCategory)[v->rateCategory];
								t_uv = this->GetEdgeLength(u->id, v->id);
								Q_uv_scaled = Q_uv_norm * t_uv;
								P_uv = Q_uv_scaled.exp();
								firstDerivativeOfpartialLikelihood = 0;
								secondDerivativeOfpartialLikelihood = 0;
								for (int dna_c = 0; dna_c < 4; dna_c ++) {
									firstDerivativeOfpartialLikelihood += P_uv(dna_p,dna_c)*firstDerivativeOfConditionalLikelihoodMap[v->id][dna_c];
									secondDerivativeOfpartialLikelihood += P_uv(dna_p,dna_c)*secondDerivativeOfConditionalLikelihoodMap[v->id][dna_c];
								}
							} else {
								w = (*this->vertexMap)[child_id];
								Q_uw = (*this->rateMatrixPerRateCategory)[w->rateCategory];
								Q_uw_norm = Q_uw/(*this->scalingFactorForRateCategory)[w->rateCategory];
								t_uw = this->GetEdgeLength(u->id, w->id);
								Q_uw_scaled = Q_uw_norm * t_uw;
								P_uw = Q_uw_scaled.exp();
								partialLikelihood = 0;
								for (int dna_c = 0; dna_c < 4; dna_c ++) {
									partialLikelihood += P_uw(dna_p,dna_c)*conditionalLikelihoodMap[w->id][dna_c];
								}
							}
						}
						firstDerivativeOfConditionalLikelihood[dna_p] *= firstDerivativeOfpartialLikelihood;
						firstDerivativeOfConditionalLikelihood[dna_p] *= partialLikelihood;
						secondDerivativeOfConditionalLikelihood[dna_p] *= secondDerivativeOfpartialLikelihood;
						secondDerivativeOfConditionalLikelihood[dna_p] *= partialLikelihood;
					}
					firstDerivativeOfConditionalLikelihoodMap.insert(pair<int,array<double,4>>(u->id,firstDerivativeOfConditionalLikelihood));
					secondDerivativeOfConditionalLikelihoodMap.insert(pair<int,array<double,4>>(u->id,secondDerivativeOfConditionalLikelihood));
				}
				// Compute site Likelihood
				siteLikelihood = 0;
				firstDerivativeOfSiteLikelihood = 0;
				secondDerivativeOfSiteLikelihood = 0;
				for (int dna = 0; dna < 4; dna ++) {
					siteLikelihood += this->rootProbability[dna] * conditionalLikelihoodMap[-1][dna];
					firstDerivativeOfSiteLikelihood += this->rootProbability[dna] * firstDerivativeOfConditionalLikelihoodMap[-1][dna];
					secondDerivativeOfSiteLikelihood += this->rootProbability[dna] * secondDerivativeOfConditionalLikelihoodMap[-1][dna];
				}
//				cout << "Second der of site likelihood is " << secondDerivativeOfSiteLikelihood << endl;
				// compute first derivative of log likelihood
				firstDerivativeOfLogLikelihood += (firstDerivativeOfSiteLikelihood/siteLikelihood) * this->siteWeights[site];
				// compute second derivative of log likelihood
				secondDerivativeOfLogLikelihood += (secondDerivativeOfSiteLikelihood/siteLikelihood)  * this->siteWeights[site];
				secondDerivativeOfLogLikelihood -= (pow((firstDerivativeOfSiteLikelihood/siteLikelihood),2.0))  * this->siteWeights[site];			
			}
			t_ab = this->GetEdgeLength(a->id, b->id);
			stepSize = 1.0;
//			cout << "Edge length before updating is " << t_ab << "\n";
//			cout << "first derivative is " << firstDerivativeOfLogLikelihood << endl;
//			cout << "second derivative is " << secondDerivativeOfLogLikelihood << endl;
			if (abs(firstDerivativeOfLogLikelihood/secondDerivativeOfLogLikelihood) < pow(10,-4) or t_ab < pow(10,-5)){
				convergenceNotReached = 0;
//				cout << "convergence reached" << endl;
			} else {
				if (t_ab - firstDerivativeOfLogLikelihood/secondDerivativeOfLogLikelihood < 0){
//					cout << "Case 1" << endl;
					while (t_ab - stepSize * (firstDerivativeOfLogLikelihood/secondDerivativeOfLogLikelihood) < 0){
						stepSize /= 2.0;
					}					
					t_ab -= stepSize * (firstDerivativeOfLogLikelihood/secondDerivativeOfLogLikelihood);
				} else {
//					cout << "Case 2" << endl;
					t_ab -= firstDerivativeOfLogLikelihood/secondDerivativeOfLogLikelihood;
				}
//				cout << "Edge length after updating is " << t_ab << "\n";
				this->SetEdgeLength(a->id,b->id,t_ab);
			}					
		}				
	}		
}

double rootedPhylogeny_tree::ComputeNormOfJacobian(MatrixXd Jacobian) {
	double norm = 0;
	for (int par = 0; par < 11; par++) {
		norm += pow(Jacobian(par,0),2);
	}
	return (norm);
}

void rootedPhylogeny_tree::SetParametersForRateMatrixForBFGS(MatrixXd parameters, int rateCat) {
	(*this->parametersPerRateCategory)[rateCat] = parameters;
	Matrix4f Q = this->GetRateMatrixForFreeParameters(parameters);
	(*this->rateMatrixPerRateCategory)[rateCat] = Q;
	(*this->scalingFactorForRateCategory)[rateCat] = 1.0;
}

double rootedPhylogeny_tree::GetOptimalStepSizeUsingLineSearch(MatrixXd B_current) {
	double t = pow(10,-4);
	double lower_limit = 0.0;
	double upper_limit = 100.0;
	// par 0 is pi_1, par 1 is pi_2, par 2 is pi_3
	for (int par = 0; par < 3; par ++) {
		if (this->searchDirection(par,0) < 0) {
			upper_limit = min(upper_limit,(-1.0 * B_current(par,0))/(this->searchDirection(par,0)));
		} else {
			upper_limit = min(upper_limit,(1.0 - B_current(par,0))/(this->searchDirection(par,0)));
		}
	}
	double directionForPi_4 = this->searchDirection(0,0) + this->searchDirection(1,0) + this->searchDirection(2,0);
	double sum_pi_1_pi_2_pi_3 = B_current(0,0) + B_current(1,0) + B_current(2,0);	
	// check for pi_4
	if (directionForPi_4 < 0){
		upper_limit = min(upper_limit, (-1.0 * sum_pi_1_pi_2_pi_3)/directionForPi_4);
	} else {
		upper_limit = min(upper_limit, (1.0 - sum_pi_1_pi_2_pi_3)/directionForPi_4);
	}
	for (int par = 3; par < 10; par++){
		if (this->searchDirection(par,0) < 0){
			upper_limit = min(upper_limit,(-1 * B_current(par,0))/this->searchDirection(par,0));
		}
	}
	cout << "upper limit is " << upper_limit << endl;
	double stepSize;
	double result;
	result = this->BrentLineSearchForSelectingStepSize(lower_limit, upper_limit, t, stepSize);
	bool verbose = 0;
	if (verbose) {
		cout << "Result is " << result << endl;
	}
	return (stepSize);	
	
}

double rootedPhylogeny_tree::BrentLineSearchForSelectingStepSize(double a, double b, double t, double &x) {
	double c;
	double d;
	double e;
	double eps;
	double fu;
	double fv;
	double fw;
	double fx;
	double m;
	double p;
	double q;
	double r;
	double sa;
	double sb;
	double t2;
	double tol;
	double u;
	double v;
	double w;
	
//	const double value = 2.220446049250313E-016;
	//
	//  C is the square of the inverse of the golden ratio.
	//
	c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );

	eps = sqrt ( 2.220446049250313E-016 );

	sa = a;
	sb = b;
	x = sa + c * ( b - a );
	w = x;
	v = w;
	e = 0.0;
	fx = this->GetNegLogLikelihoodForStepSize( x );
//	fx = this->SampleFunctionForMinimization( x );
	fw = fx;
	fv = fw;
	
	for ( ; ; )
  {
    m = 0.5 * ( sa + sb ) ;
    tol = eps * fabs ( x ) + t;
    t2 = 2.0 * tol;
//
//  Check the stopping criterion.
//
    if ( fabs ( x - m ) <= t2 - 0.5 * ( sb - sa ) )
    {
      break;
    }
//
//  Fit a parabola.
//
    r = 0.0;
    q = r;
    p = q;

    if ( tol < fabs ( e ) )
    {
      r = ( x - w ) * ( fx - fv );
      q = ( x - v ) * ( fx - fw );
      p = ( x - v ) * q - ( x - w ) * r;
      q = 2.0 * ( q - r );
      if ( 0.0 < q )
      {
        p = - p;
      }
      q = fabs ( q );
      r = e;
      e = d;
    }

    if ( fabs ( p ) < fabs ( 0.5 * q * r ) && q * ( sa - x ) < p &&  p < q * ( sb - x ) )
    {
//
//  Take the parabolic interpolation step.
//
      d = p / q;
      u = x + d;
//
//  F must not be evaluated too close to A or B.
//
      if ( ( u - sa ) < t2 || ( sb - u ) < t2 )
      {
        if ( x < m )
        {
          d = tol;
        }
        else
        {
          d = - tol;
        }
      }
    }
//
//  A golden-section step.
//
    else
    {
      if ( x < m )
      {
        e = sb - x;
      }
      else
      {
        e = sa - x;
      }
      d = c * e;
    }
//
//  F must not be evaluated too close to X.
//
    if ( tol <= fabs ( d ) )
    {
      u = x + d;
    }
    else if ( 0.0 < d )
    {
      u = x + tol;
    }
    else
    {
      u = x - tol;
    }
	fu = this->GetNegLogLikelihoodForStepSize( u );
//    fu = this->SampleFunctionForMinimization( u );
//
//  Update A, B, V, W, and X.
//
    if ( fu <= fx )
    {
      if ( u < x )
      {
        sb = x;
      }
      else
      {
        sa = x;
      }
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    }
    else
    {
      if ( u < x )
      {
        sa = u;
      }
      else
      {
        sb = u;
      }

      if ( fu <= fw || w == x )
      {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      }
      else if ( fu <= fv || v == x || v == w )
      {
        v = u;
        fv = fu;
      }
    }
  }
  return fx;
}

double rootedPhylogeny_tree::GetNegLogLikelihoodForStepSize(double stepSize){
	int rateCat = this->rateCategoryForOptimization;
	MatrixXd B_current = (*this->parametersPerRateCategory)[rateCat];
	MatrixXd B_updated = B_current + stepSize*this->searchDirection;
	(*this->parametersPerRateCategory)[rateCat] = B_updated;
	(*this->rateMatrixPerRateCategory)[rateCat] = this->GetRateMatrixForFreeParameters(B_updated);
	this->ComputeLogLikelihood();
	return (-1 * this->logLikelihood);
}


Matrix4f rootedPhylogeny_tree::GetRateMatrixForFreeParameters(MatrixXd B){
	Matrix4f Q;
	float pi_1; float pi_2; float pi_3; float pi_4;
	float a; float b; float c; float d; float e; float f;
	float g; float h;float i; float j; float k; float l;
	pi_1 = B(0,0); pi_2 = B(1,0); pi_3 = B(2,0);
	pi_4 = (1 - pi_1 - pi_2 - pi_3);
	// Construct Q
	a = B(3,0); b = B(4,0); c = B(5,0); d = B(6,0); e = B(7,0);
	f = B(8,0); g = B(9,0); h = B(10,0);

	i = (1-(pi_1*(a+b+2*c)+pi_2*(d+e+2*f)+pi_3*(g+h)))/(2*pi_3);
	j = (pi_1*(a+b+c)-pi_2*d-pi_3*g)/pi_4;
	k = (pi_2*(d+e+f)-pi_1*a-pi_3*h)/pi_4;
	l = (1+pi_3*(g+h)-pi_1*(a+2*c+3*b)-pi_2*(d+2*f+3*e))/(2*pi_4);
	
	Q << -(a+b+c), a, b, c,
		  d, -(d+e+f), e, f,
		  g, h, -(g+h+i), i,
		  j, k, l, -(j+k+l);
	
	return (Q);
}

void rootedPhylogeny_tree::SetVerticesForPreOrderTreeTraversal(){
	this->verticesForPreOrderTreeTraversal->clear();
	this->verticesForPreOrderTreeTraversal->push_back(this->root);
	vector <rootedPhylogeny_vertex*> verticesToVisit;
	verticesToVisit.push_back(this->root);
	int numberOfVerticesToVisit = verticesToVisit.size();
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	while (numberOfVerticesToVisit > 0){
		p = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		for (int c_id : p->children_id){
			c = (*this->vertexMap)[c_id];
			this->verticesForPreOrderTreeTraversal->push_back(c);
			if (c->children_id.size() > 0){
				verticesToVisit.push_back(c);				
				numberOfVerticesToVisit += 1;
			}
		}
	}
}

void rootedPhylogeny_tree::SetRateCategories(float baseFreqThreshold){
	this->SetVerticesForPreOrderTreeTraversal();
	// Select change points;
	this->changePoints->clear();
	this->changePoints->push_back(this->root);
	for (pair<rootedPhylogeny_vertex*,float> descendantAndBaseFreqChange : *this->baseFreqChangeForEachDescendant){		
		if (descendantAndBaseFreqChange.second > baseFreqThreshold){
			this->changePoints->push_back(descendantAndBaseFreqChange.first);
		}
	}
	this->numberOfRateCategories = this->changePoints->size();
//	cout << "number of rate categories is " << this->numberOfRateCategories << endl;
	int rateCategory = -1;
	
	for (rootedPhylogeny_vertex * v : (*this->verticesForPreOrderTreeTraversal)){
		vector <rootedPhylogeny_vertex*>::iterator it = find(this->changePoints->begin(),this->changePoints->end(),v);
		if (it == this->changePoints->end()){
			rateCategory = (*this->vertexMap)[v->parent_id]->rateCategory;
			v->rateCategory = rateCategory;			
			
		} else {			
			rateCategory = distance(this->changePoints->begin(),it);
			v->rateCategory = rateCategory;
		}
		if (rateCategory == -1) {
			cout << "check assignment of rate categories" << endl;
		}
	}
}

void rootedPhylogeny_tree::OptimizeModelParametersDep(){	
	// Compute initial estimate of model parameters
	this->InitializeModelParametersForNelderMead();
	MatrixXd B;
	int rateCat;
	Matrix4f Q;
	for (pair<int,MatrixXd> rateCatParamPair : *this->parametersPerRateCategory){
		rateCat = rateCatParamPair.first;
		B = rateCatParamPair.second;		
		Q = this->GetRateMatrixForFreeParameters(B);
		this->rateMatrixPerRateCategory->insert(pair<int,Matrix4f>(rateCat,Q));
	}	
	this->SetMinLengthOfEdges(); // 10^-7;
	//	cout << "here 1" << endl;
	// Iterate over rate category
	for (int rateCat = 0; rateCat < this->numberOfRateCategories; rateCat++){
		this->OptimizeModelParametersForRateCatDep(rateCat);
	}
}

void rootedPhylogeny_tree::OptimizeModelParametersForRateCatDep(int rateCat){
	//	using namespace brent;
	if (1){
		double convergenceTolerance = pow(10,-4);
		MatrixXd InverseHessian_current = ArrayXXd::Zero(11,11);
		MatrixXd InverseHessian_updated;
		MatrixXd B_current; MatrixXd B_updated;
		MatrixXd Jacobian_current; MatrixXd Jacobian_updated;
		MatrixXd s;
		MatrixXd y;
		MatrixXd rho = ArrayXXd::Zero(11,11);;
		MatrixXd IdentityMatrix = ArrayXXd::Zero(11,11);
		MatrixXd TempMatrix_l = ArrayXXd::Zero(11,11);
		MatrixXd TempMatrix_r = ArrayXXd::Zero(11,11);
		int iter;
		double normOfJacobian;
		double stepSize;
		this->rateCategoryForOptimization = rateCat;
		
		for (int par_1 = 0; par_1 < 11; par_1 ++){		
			InverseHessian_current(par_1, par_1) = 1.0;
			IdentityMatrix(par_1, par_1) = 1.0;
		}
		
		B_current = (*this->parametersPerRateCategory)[rateCat];
		Jacobian_current = this->GetJacobianForRateCategory(rateCat);
		normOfJacobian = this->ComputeNormOfJacobian(Jacobian_current);
		iter = 0;
//		cout << "norm of Jacobian for iteration " << iter << " is " << normOfJacobian << endl;
		while (normOfJacobian > convergenceTolerance and iter < 100){
			cout << "B_current is \n" << B_current << endl;
//			cout << "norm of Jacobian for iteration " << iter << " is " << normOfJacobian << endl;
			iter += 1;
//			cout << "Jacobian current is " << endl;
//			cout << Jacobian_current << endl;
			this->searchDirection = -1 * (InverseHessian_current * Jacobian_current);
//			cout << "search direction" << endl;
//			cout << this->searchDirection << endl;
//			cout << "B_current is \n" << B_current << endl;
			stepSize = this->GetOptimalStepSizeUsingLineSearch(B_current);
			cout << "step size is " << stepSize << endl;			
			B_updated = B_current + stepSize * this->searchDirection;					
			this->SetParametersForRateMatrixForBFGS(B_updated,rateCat);
			Jacobian_updated = this->GetJacobianForRateCategory(rateCat);
			y = Jacobian_updated - Jacobian_current;
			s = B_updated - B_current;
//			cout << "y" << endl;
//			cout << y << endl;
//			cout << "s" << endl;
//			cout << s << endl;
			rho = y.transpose()*s;
			rho = rho.inverse();
//			cout << "rho is " << endl;
//			cout << rho << endl;
//			cout << s * y.transpose() << endl;
			TempMatrix_l = IdentityMatrix - rho(0,0) * ( s * y.transpose() );
			TempMatrix_r = IdentityMatrix - rho(0,0) * ( y * s.transpose() );
			InverseHessian_updated = TempMatrix_l * InverseHessian_current;
			InverseHessian_updated = InverseHessian_updated * TempMatrix_r;
			InverseHessian_updated += rho(0,0) * ( s * s.transpose() );
			B_current = B_updated;
			InverseHessian_current = InverseHessian_updated;
			normOfJacobian = this->ComputeNormOfJacobian(Jacobian_updated);
		}
		
	} else {
		double x = 0.1;
		double result;
//		result = this->BrentLineSearchForSelectingStepSize(lower_limit, upper_limit, t, x);
//		result = local_min(lower_limit, upper_limit, t, this->SampleFunctionForMinimization(x), x);
//		result = local_min (lower_limit, upper_limit, t, this->SampleFunctionForMinimization, x);
//		result = local_min (lower_limit, upper_limit, t, this->SampleFunctionForMinimization, x);
		cout << "optimal x is " << x << endl;
		cout << "value of function at minima is " << result << endl;
	}
//	double stepSize;
//	double negLogLikelihood;
//	negLogLikelihood = local_min (lower_limit, upper_limit, t, this->SampleFunctionForMinimization, stepSize);	
}

void rootedPhylogeny_tree::InitializeModelParametersForBFGS(){
	this->rateMatrixPerRateCategory->clear();
	this->scalingFactorForRateCategory->clear();
	this->parametersPerRateCategory->clear();
	// Initialize substitution matrices for each rate category
	map <int, Matrix4f> substitutionMatrixPerRateCategory;
	MatrixXf stationaryDistribution;  
	Matrix4f P;
	// Populate substitution matrices for each rate category
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;
	int rateCat; Matrix4f Q;
	unsigned char dna_p; unsigned char dna_c;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->id != -1){
			if (substitutionMatrixPerRateCategory.find(c->rateCategory) == substitutionMatrixPerRateCategory.end()){
				P = ArrayXXf::Zero(4,4);
				substitutionMatrixPerRateCategory.insert(pair<int,Matrix4f>(c->rateCategory,P));
			}
			p = (*this->vertexMap)[c->parent_id];
			for (unsigned int site = 0; site < this->siteWeights.size(); site++){
				dna_p = p->compressedSequence[site];
				dna_c = c->compressedSequence[site];
				substitutionMatrixPerRateCategory[c->rateCategory](dna_p,dna_c) += this->siteWeights[site];
			}
		}
	}	
	float rowSum; float scalingFactor; MatrixXd parameters;
	for (pair<int,Matrix4f> rateCatSubsMatrixPair : substitutionMatrixPerRateCategory){
		rateCat = rateCatSubsMatrixPair.first;
		P = rateCatSubsMatrixPair.second;
		for (int row = 0; row < 4; row++){
			rowSum = 0;
			for (int col = 0; col < 4; col++){
				rowSum += P(row,col);
			}
			for (int col = 0; col < 4; col++){
				P(row,col) /= rowSum;
			}
		}
				
		Q = P.log();
		
		for (int row = 0; row < 4; row ++){
			for (int col = 0; col < 4; col ++){
				if (row != col){
					if (Q(row,col) < pow(10,-4)){
						Q(row,col) = pow(10,-4);
					}
				}
			}
		}
		
		Q(0,0) = - (Q(0,1) + Q(0,2) + Q(0,3));
		Q(1,1) = - (Q(1,0) + Q(1,2) + Q(1,3));
		Q(2,2) = - (Q(2,0) + Q(2,1) + Q(2,3));
		Q(3,3) = - (Q(3,0) + Q(3,1) + Q(3,2));
		scalingFactor = this->ComputeScalingFactor(Q);
		Q/=scalingFactor;
		this->rateMatrixPerRateCategory->insert(pair<int,Matrix4f>(rateCat,Q));
		scalingFactor = this->ComputeScalingFactor(Q);
		this->scalingFactorForRateCategory->insert(pair<int,float>(rateCat,scalingFactor));
		stationaryDistribution = this->ComputeStationaryDistribution(Q);
		parameters = ArrayXXd::Zero(11,1);
		// pi_1
		parameters(0,0) = stationaryDistribution(0,1); 
		// pi_2
		parameters(1,0) = stationaryDistribution(0,2);
		// pi_3
		parameters(2,0) = stationaryDistribution(0,3);
		// a
		parameters(3,0) = Q(0,1); 
		// b
		parameters(4,0) = Q(0,2);
		// c
		parameters(5,0) = Q(0,3);
		// d
		parameters(6,0) = Q(1,0);
		// e
		parameters(7,0) = Q(1,2);
		// f
		parameters(8,0) = Q(1,3);
		// g
		parameters(9,0) = Q(2,0);
		// h
		parameters(10,0) = Q(2,1);		
		for (int par = 0; par < 11; par++){
			if (parameters(par,0) < 0){
				cout << "Initialized rate matrix has negative entries on non-diagonal elements" << endl;
				cout << Q << endl;
			}
		}
		this->parametersPerRateCategory->insert(pair<int,MatrixXd>(rateCat,parameters));
	}
	
	if (rateMatrixPerRateCategory->find(this->root->rateCategory) == rateMatrixPerRateCategory->end()){
//		cout << "rate cat for root is " << this->root->rateCategory << endl;
//		cout << "rate cat for left child of root is " << (*this->vertexMap)[this->root->children_id[0]]->rateCategory << endl;
//		cout << "rate cat for right child of root is " << (*this->vertexMap)[this->root->children_id[1]]->rateCategory << endl;
		// There is no rate matrix that is associated with the rate category root.cat
//		cout << " There is no rate matrix that is associated with the rate category root.cat" << endl;
		for (int dna = 0; dna < 4; dna ++){
			this->rootProbability[dna] = 0.0;
		}
		unsigned char dna_p;
		float sequenceLength = 0;
		for (unsigned int site = 0; site < this->siteWeights.size(); site++){
			dna_p = this->root->compressedSequence[site];
			this->rootProbability[dna_p] += this->siteWeights[site];
			sequenceLength += this->siteWeights[site];
		}
		for (int dna = 0; dna < 4; dna ++){
			this->rootProbability[dna] /= sequenceLength;
		}
		
	} else {
		// There is one rate matrix that is associated with the rate category root.cat
		Q = (*this->rateMatrixPerRateCategory)[this->root->rateCategory];
		MatrixXf stationaryDistribution = this->ComputeStationaryDistribution(Q);
		this->rootProbability[0] = stationaryDistribution(0,0);
		this->rootProbability[1] = stationaryDistribution(1,0);
		this->rootProbability[2] = stationaryDistribution(2,0);
		this->rootProbability[3] = stationaryDistribution(3,0);
	}
}

void rootedPhylogeny_tree::InitializeModelParametersForNelderMead(){
	this->rateMatrixPerRateCategory->clear();
	this->scalingFactorForRateCategory->clear();
	this->parametersPerRateCategory->clear();
	// Initialize substitution matrices for each rate category
	map <int, Matrix4f> substitutionMatrixPerRateCategory;
	Matrix4f P;
	// Populate substitution matrices for each rate category
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;
	int rateCat; Matrix4f Q;
	unsigned char dna_p; unsigned char dna_c;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->id != -1){
			if (substitutionMatrixPerRateCategory.find(c->rateCategory) == substitutionMatrixPerRateCategory.end()){
				P = ArrayXXf::Zero(4,4);
				substitutionMatrixPerRateCategory.insert(pair<int,Matrix4f>(c->rateCategory,P));
			}
			p = (*this->vertexMap)[c->parent_id];
			for (unsigned int site = 0; site < this->siteWeights.size(); site++){
				dna_p = p->compressedSequence[site];
				dna_c = c->compressedSequence[site];
				substitutionMatrixPerRateCategory[c->rateCategory](dna_p,dna_c) += this->siteWeights[site];
			}
		}
	}	
	float rowSum; float scalingFactor; MatrixXd parameters;
	for (pair<int,Matrix4f> rateCatSubsMatrixPair : substitutionMatrixPerRateCategory){
		rateCat = rateCatSubsMatrixPair.first;
		P = rateCatSubsMatrixPair.second;
		for (int row = 0; row < 4; row++){
			rowSum = 0;
			for (int col = 0; col < 4; col++){
				rowSum += P(row,col);
			}
			for (int col = 0; col < 4; col++){
				P(row,col) /= rowSum;
			}
		}
				
		Q = P.log();
		Q /= Q(3,2);		
		for (int row = 0; row < 4; row ++){
			for (int col = 0; col < 4; col ++){
				if (row != col){
					if (Q(row,col) < 0){
						Q(row,col) *= -1;
					}
				}
			}
		}
		
		for (int row = 0; row < 4; row ++){
			for (int col = 0; col < 4; col ++){
				if (row != col){
					if (Q(row,col) < pow(10,-4)){
						Q(row,col) = pow(10,-4);
					}
				}
			}
		}
		
		Q(0,0) = - (Q(0,1) + Q(0,2) + Q(0,3));
		Q(1,1) = - (Q(1,0) + Q(1,2) + Q(1,3));
		Q(2,2) = - (Q(2,0) + Q(2,1) + Q(2,3));
		Q(3,3) = - (Q(3,0) + Q(3,1) + Q(3,2));
		this->rateMatrixPerRateCategory->insert(pair<int,Matrix4f>(rateCat,Q));
		scalingFactor = this->ComputeScalingFactor(Q);		
		this->scalingFactorForRateCategory->insert(pair<int,float>(rateCat,scalingFactor));
		parameters = ArrayXXd::Zero(11,1);
		// a
		parameters(0,0) = Q(0,1); 
		// b
		parameters(1,0) = Q(0,2);
		// c
		parameters(2,0) = Q(0,3);
		// d
		parameters(3,0) = Q(1,0);
		// e
		parameters(4,0) = Q(1,2);
		// f
		parameters(5,0) = Q(1,3);
		// g
		parameters(6,0) = Q(2,0);
		// h
		parameters(7,0) = Q(2,1);
		// i
		parameters(8,0) = Q(2,3);
		// j
		parameters(9,0) = Q(3,0);	
		// k
		parameters(10,0) = Q(3,1);
		for (int par = 0; par < 11; par++){
			if (parameters(par,0) < 0){
				cout << "Initialized rate matrix has negative entries on non-diagonal elements" << endl;
				cout << Q << endl;
			}
		}
		this->parametersPerRateCategory->insert(pair<int,MatrixXd>(rateCat,parameters));
	}
	
	if (rateMatrixPerRateCategory->find(this->root->rateCategory) == rateMatrixPerRateCategory->end()){
//		cout << "rate cat for root is " << this->root->rateCategory << endl;
//		cout << "rate cat for left child of root is " << (*this->vertexMap)[this->root->children_id[0]]->rateCategory << endl;
//		cout << "rate cat for right child of root is " << (*this->vertexMap)[this->root->children_id[1]]->rateCategory << endl;
		// There is no rate matrix that is associated with the rate category root.cat
//		cout << " There is no rate matrix that is associated with the rate category root.cat" << endl;
		for (int dna = 0; dna < 4; dna ++){
			this->rootProbability[dna] = 0.0;
		}
		unsigned char dna_p;
		float sequenceLength = 0;
		for (unsigned int site = 0; site < this->siteWeights.size(); site++){
			dna_p = this->root->compressedSequence[site];
			this->rootProbability[dna_p] += this->siteWeights[site];
			sequenceLength += this->siteWeights[site];
		}
		for (int dna = 0; dna < 4; dna ++){
			this->rootProbability[dna] /= sequenceLength;
		}
		
	} else {
		// There is one rate matrix that is associated with the rate category root.cat
		Q = (*this->rateMatrixPerRateCategory)[this->root->rateCategory];
		MatrixXf stationaryDistribution = this->ComputeStationaryDistribution(Q);
		this->rootProbability[0] = stationaryDistribution(0,0);
		this->rootProbability[1] = stationaryDistribution(1,0);
		this->rootProbability[2] = stationaryDistribution(2,0);
		this->rootProbability[3] = stationaryDistribution(3,0);
	}
}

float rootedPhylogeny_tree::ComputeAbsDifferenceInBaseFreq(array<float, 4> baseFreq_1, array<float, 4> baseFreq_2){	
	// A is 0, C is 1, G is 2, T is 3
	float abs_nt_diff = 0;
	for (int dna = 0; dna < 4; dna ++){
		abs_nt_diff += abs(baseFreq_1[dna] - baseFreq_2[dna]);
	}	
	return (abs_nt_diff);	
}


void rootedPhylogeny_tree::SetThresholds(){
	this->baseFreqChangeForEachDescendant->clear();	
	this->thresholdList->clear();
	vector <pair<float,rootedPhylogeny_vertex*>> negBaseFreqDiffAndChangePointList;
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;
	float negBaseFreqDiff;
	float baseFreqDiff;
	// goodness-of-clustering vs goodness-of-fit
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
		c = idPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			negBaseFreqDiff = -1 * this->ComputeAbsDifferenceInBaseFreq(p->baseFreq, c->baseFreq);
			this->baseFreqChangeForEachDescendant->insert(pair<rootedPhylogeny_vertex*,float>(c,this->ComputeAbsDifferenceInBaseFreq(p->baseFreq, c->baseFreq)));
			negBaseFreqDiffAndChangePointList.push_back(pair<float,rootedPhylogeny_vertex*>(negBaseFreqDiff,c));
		}		
	}
	sort(negBaseFreqDiffAndChangePointList.begin(), negBaseFreqDiffAndChangePointList.end());
	for (pair<float,rootedPhylogeny_vertex*> negBaseFreqDiffAndChangePoint : negBaseFreqDiffAndChangePointList){
		baseFreqDiff = -1 * negBaseFreqDiffAndChangePoint.first;
		if (find(this->thresholdList->begin(),this->thresholdList->end(),baseFreqDiff) == this->thresholdList->end()){
//			cout << baseFreqDiff << endl;
			this->thresholdList->push_back(baseFreqDiff);
		}
	}
	this->thresholdList->push_back(-1.0);
}

void rootedPhylogeny_tree::ComputeBaseFreq(){
	rootedPhylogeny_vertex * v;
	int dna;
	float sequenceLength = 0;
	for (unsigned int site = 0; site < this->siteWeights.size(); site ++){
		sequenceLength += this->siteWeights[site];
	}	
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
		v = idPtrPair.second;
		for (unsigned int site = 0; site < this->siteWeights.size(); site ++){
			dna = v->compressedSequence[site];
			v->baseFreq[dna] += this->siteWeights[site];
		}
		for (int dna = 0; dna < 4; dna ++){
			v->baseFreq[dna] /= sequenceLength;
		}
	}	
}

double rootedPhylogeny_tree::SampleFunctionForMinimization(double x){	
	return ((x+3)*pow((x-1),2.0));
}

double rootedPhylogeny_tree::ComputeEdgeLogLikelihood(rootedPhylogeny_vertex * p, rootedPhylogeny_vertex * c, Matrix4f P){
	unsigned char dna_p; unsigned char dna_c;
	double edgeLogLik = 0;
	for (unsigned int site = 0; site < this->siteWeights.size(); site++){
		dna_p = p->compressedSequence[site];
		dna_c = c->compressedSequence[site];
		edgeLogLik += log(P(dna_p,dna_c))*this->siteWeights[site];
	}
	return (edgeLogLik);
}

void rootedPhylogeny_tree::AddEdgeLength(int p_id, int c_id, float edgeLength){
	this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(p_id,c_id),edgeLength));
}

void rootedPhylogeny_tree::SetEdgeLength(int p_id, int c_id, float edgeLength){
	if (p_id < c_id){
		(*this->edgeLengthsMap)[pair<int,int>(p_id,c_id)] = edgeLength;
	} else {
		(*this->edgeLengthsMap)[pair<int,int>(c_id,p_id)] = edgeLength;
	}
}

int rootedPhylogeny_tree::GetVertexId(string v_name){
	rootedPhylogeny_vertex * v;
	int idToReturn = -10;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
		v = idPtrPair.second;
		if (v->name == v_name){
			idToReturn = v->id;						
		}
	}
	if (idToReturn == -10){
		cout << "check id to return" << endl;
	}
	return (idToReturn);
}

void rootedPhylogeny_tree::ReadTreeFile(string treeFileName){
	string u_name;
	string v_name;
	int u_id;
	int v_id;
	float edgeLength;
	vector <string> splitLine;
	vector <string> leafNames;
	vector <string> ancestorNames;
	vector <string> nonRootVertexNames;	
	string rootName = "";
	vector <unsigned char> emptySequence;
	v_id = 0;
	ifstream edgeListFile(treeFileName.c_str());
	for (string line; getline(edgeListFile, line);){
		boost::split(splitLine, line, [](char c){return c == '\t';});
		u_name = splitLine[0];
		v_name = splitLine[1];
		if (find(nonRootVertexNames.begin(),nonRootVertexNames.end(),v_name) == nonRootVertexNames.end()){
			nonRootVertexNames.push_back(v_name);
		}		
		if (find(ancestorNames.begin(),ancestorNames.end(),u_name)==ancestorNames.end()){
			ancestorNames.push_back(u_name);
		}
		if (find(leafNames.begin(),leafNames.end(),v_name)==leafNames.end()){
			if(!boost::starts_with(v_name, "h_")){
				leafNames.push_back(v_name);
			}
		}
	}	
	for (string name: leafNames){
		rootedPhylogeny_vertex * v = new rootedPhylogeny_vertex(v_id,name,emptySequence);
		this->vertexMap->insert(pair<int,rootedPhylogeny_vertex*>(v_id,v));
		v_id += 1;
	}
	// Remove root from ancestor names
	for (string name: ancestorNames){
		if (find(nonRootVertexNames.begin(),nonRootVertexNames.end(),name)==nonRootVertexNames.end()){
			rootName = name;
		}		
	}
	// Change root name
	(*this->vertexMap)[-1]->name = rootName;
	ancestorNames.erase(remove(ancestorNames.begin(), ancestorNames.end(), rootName), ancestorNames.end());	
	for (string name: ancestorNames){		
		rootedPhylogeny_vertex * v = new rootedPhylogeny_vertex(v_id,name,emptySequence);
		this->vertexMap->insert(pair<int,rootedPhylogeny_vertex*>(v_id,v));
		v_id += 1;
	}	
	edgeListFile.clear();
	edgeListFile.seekg(0, ios::beg); 	
	for (string line; getline(edgeListFile, line);){
		boost::split(splitLine, line, [](char c){return c == '\t';});
		u_name = splitLine[0];
		v_name = splitLine[1];
		u_id = this->GetVertexId(u_name);
		v_id = this->GetVertexId(v_name);		
		edgeLength = stof(splitLine[2]);
		this->AddEdgeLength(u_id,v_id,edgeLength);
	}
	edgeListFile.close();
}


void rootedPhylogeny_tree::ComputeLogLikelihoodForFullyLabeledTree(){
	this->logLikelihood = 0;
	Matrix4f P; Matrix4f Q; Matrix4f Q_scaled;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	float t;
	this->SetMinLengthOfEdges();
	Q = this->rateMatrix;
	unsigned char dna_p; unsigned char dna_c;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
		c = idPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			t = this->GetEdgeLength(p->id, c->id);
			if (t < pow(10,-7)){
				t = pow(10,-7);
			}
			Q_scaled = Q*t;			
			P = Q_scaled.exp();
			for (unsigned int site = 0; site < this->siteWeights.size(); site++){
				dna_p = p->compressedSequence[site];
				dna_c = c->compressedSequence[site];
				if (t == 0){
					cout << "edgeLength of zero found" << endl;
				}
				if (P(dna_p,dna_c) == 0){
					cout << "edge length " << t << endl;
				}
				this->logLikelihood += log(P(dna_p,dna_c))*this->siteWeights[site];				
			}
		}
	}	
	for (unsigned int site = 0; site < this->siteWeights.size(); site++){
		this->logLikelihood += log(this->rootProbability[(*this->vertexMap)[-1]->compressedSequence[site]])*this->siteWeights[site];
	}	
}

void rootedPhylogeny_tree::SetMinLengthOfEdges(){
	for (pair<pair<int,int>,float> edgeIdsAndLengthsPair: *this->edgeLengthsMap){		
		if (edgeIdsAndLengthsPair.second < pow(10,-7)){
			(*this->edgeLengthsMap)[edgeIdsAndLengthsPair.first]  = pow(10,-7);
		}
	}
}

void rootedPhylogeny_tree::EstimateAncestralSequencesByFittingTheGMM(){
	this->ComputeMPEstimateOfAncestralSequences();	
	map <int,Matrix4f> transitionMatrices;
	map <int,array<double,4>> conditionalLikelihoodMap;
	array <double,4> conditionalLikelihood;
	array <float,4> prior;
	double partialLikelihood;
	double siteLikelihood;
	double currentLogLikelihood = 0;
	double previousLogLikelihood = 0;
	int iter = 0;
	int maxIter = 10;
	pair<int,int> edgeIds;
	unsigned char dna_p; unsigned char dna_c;
	char maxProbState;
	float rowSum;
	float currentProb;
	float maxProb;
	int numberOfVerticesToVisit;
	bool continueEM = 1;
	vector <rootedPhylogeny_vertex*> verticesToVisit;	
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	Matrix4f P;
	this->SetEdgesForPostOrderTreeTraversal();
	// Iterate till convergence of log likelihood
		while (continueEM and iter < maxIter) {
			iter += 1;			
			currentLogLikelihood = 0;
			// Estimate root probablity
			this->ComputeMLEOfRootProbability();		
			// Estimate transition matrices	
//			cout << "here 1" << endl;
			for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
				c = idPtrPair.second;
				if (c->id != -1){
					p = (*this->vertexMap)[c->parent_id];
					P = ArrayXXf::Zero(4,4);
					for (unsigned int site = 0; site < siteWeights.size(); site++){
						dna_p = p->compressedSequence[site];
						dna_c = c->compressedSequence[site];
						P(dna_p,dna_c) += this->siteWeights[site];
					}
					for (unsigned char dna_p = 0; dna_p < 4; dna_p++){
						rowSum = 0;
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
							rowSum += P(dna_p,dna_c);
						}
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
							 P(dna_p,dna_c) /= rowSum;
						}
					}
					transitionMatrices.insert(pair<int,Matrix4f>(c->id,P));
				}
			}
//			cout << "here 2" << endl;
			// Estimate ancestral sequences
			for (pair<int,rootedPhylogeny_vertex*> idPtrPair : (*this->vertexMap)){
				c = idPtrPair.second;
				if (c->children_id.size() > 0){
					c->compressedSequence.clear();
				}
			}		
			// Iterate over sites		
			for (unsigned int site = 0 ; site < this->siteWeights.size(); site++){
				conditionalLikelihoodMap.clear();			
				for (pair<int,int> edgeIds : *this->edgesForPostOrderTreeTraversal){
					p = (*this->vertexMap)[edgeIds.first];
					c = (*this->vertexMap)[edgeIds.second];
					P = transitionMatrices[c->id];					
					// Initialize conditional likelihood for leaves
					if (c->children_id.size()==0){
						for (unsigned char dna_c = 0; dna_c < 4; dna_c ++){
							conditionalLikelihood[dna_c] = 0;
						}
						conditionalLikelihood[c->compressedSequence[site]] = 1;
						conditionalLikelihoodMap.insert(pair<int,array<double,4>>(c->id,conditionalLikelihood));
					}
					// Initialize conditional likelihood for ancestors
					if (conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()){
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
						conditionalLikelihood[dna_c] = 1;
						}				
						conditionalLikelihoodMap.insert(pair<int,array<double,4>>(p->id,conditionalLikelihood));					
					}		
					for (unsigned char dna_p = 0; dna_p < 4; dna_p++){
						partialLikelihood = 0;
						for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
							partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
						}
						conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;
					}
				}
				maxProbState = -1;
				maxProb = 0;
				siteLikelihood = 0; 			
				prior = this->rootProbability;
//				cout << prior << endl;
				for (int dna = 0; dna < 4; dna++){
//					cout << "prior for dna:" << dna << " is " << prior[dna] << endl;
//					cout << "conditional for dna:" << dna << " is " << conditionalLikelihoodMap[-1][dna] << endl;
					currentProb = prior[dna]*conditionalLikelihoodMap[-1][dna];
					siteLikelihood += currentProb;
					if (currentProb > maxProb){
						maxProb = currentProb;
						maxProbState = dna;		
					}
				}				
				currentLogLikelihood += log(siteLikelihood) * this->siteWeights[site];
				if (maxProbState == -1){
					cout << "check state estimation" << endl;
				}			
				(*this->vertexMap)[-1]->compressedSequence.push_back(maxProbState);
				verticesToVisit.clear();			
				for (int c_id: (*this->vertexMap)[-1]->children_id) {
					if ((*this->vertexMap)[c_id]->children_id.size()>0) {
						verticesToVisit.push_back((*this->vertexMap)[c_id]);
					}
				}
				numberOfVerticesToVisit = verticesToVisit.size();
				while (numberOfVerticesToVisit > 0) {
					c = verticesToVisit[numberOfVerticesToVisit-1];				
					verticesToVisit.pop_back();
					numberOfVerticesToVisit -= 1;
					p = (*this->vertexMap)[c->parent_id];
					P = transitionMatrices[c->id];
					dna_p = p->compressedSequence[site];
					maxProbState = -1;
					maxProb = 0;
					for (int dna_c = 0; dna_c < 4; dna_c++){
						currentProb = P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
						if (currentProb > maxProb){
							maxProb = currentProb;
							maxProbState = dna_c;
						}					
					}
					c->compressedSequence.push_back(maxProbState);
					for (int c_id : c->children_id) {
						if ((*this->vertexMap)[c_id]->children_id.size()>0) {
							verticesToVisit.push_back((*this->vertexMap)[c_id]);
							numberOfVerticesToVisit += 1;
						}
					}					
				}
			}
//			cout << "iteration:" << iter << "\t";
//			cout << "previousLogLikelihood:" << previousLogLikelihood << endl;
//			cout << "currentLogLikelihood:" << currentLogLikelihood << endl;			
			if (iter < 2 or currentLogLikelihood > previousLogLikelihood or abs(currentLogLikelihood-previousLogLikelihood) > 0.001){
				continueEM = 1;
				previousLogLikelihood = currentLogLikelihood;
			} else {
				continueEM = 0;
			}
		}
		this->logLikelihood = currentLogLikelihood;
//		cout << "max log likelihood for fitting the GMM is " << currentLogLikelihood << endl;
}

void rootedPhylogeny_tree::ComputeMLEOfEdgeLengths(){
	Matrix4f Q = this->rateMatrix;
	Matrix4f Q_scaled;
	Matrix4f P;
	Matrix4f P_der_1;
	Matrix4f P_der_2;
	float t_new;
	float t_old;
	float proposedChange;
	float s = 1;
	double firstDerivative;
	double secondDerivative;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	double edgeLogLik_current;	
	double edgeLogLik_updated;
	unsigned char dna_p;
	unsigned char dna_c;	
	this->SetMinLengthOfEdges();
	bool continueOptimizingEdgeLength;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			t_old = this->GetEdgeLength(p->id,c->id);
			Q_scaled = Q*t_old;
			P = Q_scaled.exp();		
			continueOptimizingEdgeLength = 1;
			while (continueOptimizingEdgeLength){
				s = 1.0;
				Q_scaled = Q*t_old;
				P = Q_scaled.exp();
				P_der_1 = Q*P;
				P_der_2 = Q*Q;
				P_der_2 = P_der_2*P;
				firstDerivative = 0;
				secondDerivative = 0;
				edgeLogLik_current = this->ComputeEdgeLogLikelihood(p,c,P);
				for (unsigned int site = 0; site < siteWeights.size(); site++){
					dna_p = p->compressedSequence[site];
					dna_c = c->compressedSequence[site];
					firstDerivative += (P_der_1(dna_p,dna_c)/P(dna_p,dna_c))*(this->siteWeights[site]);
					secondDerivative += (P_der_2(dna_p,dna_c)/P(dna_p,dna_c) - pow(P_der_1(dna_p,dna_c)/P(dna_p,dna_c),2.0))*(this->siteWeights[site]);
				}				
				proposedChange = firstDerivative/secondDerivative; 				
				while (t_old - s*proposedChange < 0){
					s/=2.0;					
				}
				t_new = t_old - s*proposedChange;
				Q_scaled = Q*t_new;
				P = Q_scaled.exp();
				edgeLogLik_updated = this->ComputeEdgeLogLikelihood(p,c,P);
				while (edgeLogLik_updated < edgeLogLik_current and abs(edgeLogLik_updated - edgeLogLik_current) > 0.001) {
					s/=2.0;					
					t_new = t_old - s*proposedChange;
					Q_scaled = Q*t_new;
					P = Q_scaled.exp();
					edgeLogLik_updated = this->ComputeEdgeLogLikelihood(p,c,P); 
				}
				if (abs(edgeLogLik_current-edgeLogLik_updated) < 0.001){
					continueOptimizingEdgeLength = 0;
				}
				t_old = t_new;
			}
			if (p->id < c->id) {
				(*this->edgeLengthsMap)[pair<int,int>(p->id,c->id)] = t_old;
			} else {
				(*this->edgeLengthsMap)[pair<int,int>(c->id,p->id)] = t_old;
			}
		}
	}
}

void rootedPhylogeny_tree::ComputeMLEOfRateMatricesForLeafLabeledTrees() {
	
}


void rootedPhylogeny_tree::ScaleEdgeLengths(){ // Modify to include rate categories
	float scalingFactor = this->ComputeScalingFactor(this->rateMatrix);
	for (pair<pair<int,int>,float> edgeIdsAndLengthsPair: *this->edgeLengthsMap){
		(*this->edgeLengthsMap)[edgeIdsAndLengthsPair.first] = edgeIdsAndLengthsPair.second/scalingFactor;
	}
}


void rootedPhylogeny_tree::ComputeMAPEstimateOfAncestralSequences(){
	map <int, array<double,4>> conditionalLikelihoodMap;
	array <double, 4> conditionalLikelihood;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	int maxProbState;
	double maxProb;
	double partialLikelihood;	
	double siteLikelihood;
	double currentProb;
	unsigned char dna_p;
	vector <rootedPhylogeny_vertex *> verticesToVisit;
	// Clear ancestral sequences
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->children_id.size() > 0){
			c->compressedSequence.clear();
		}
	}
	Matrix4f Q = this->rateMatrix;
	Matrix4f Q_scaled;
	Matrix4f P;
	float t;
	int numberOfVerticesToVisit;
	double currentLogLikelihood = 0;
	// Iterate over sites		
	for (unsigned int site = 0 ; site < this->siteWeights.size(); site++){
		conditionalLikelihoodMap.clear();			
		for (pair<int,int> edgeIds : *this->edgesForPostOrderTreeTraversal){
			p = (*this->vertexMap)[edgeIds.first];
			c = (*this->vertexMap)[edgeIds.second];
			t = this->GetEdgeLength(p->id,c->id);
			Q_scaled = Q*t;
			P = Q_scaled.exp();
			// Initialize conditional likelihood for leaves
			if (c->children_id.size()==0){
				for (unsigned char dna_c = 0; dna_c < 4; dna_c ++){
					conditionalLikelihood[dna_c] = 0;
				}
				conditionalLikelihood[c->compressedSequence[site]] = 1;
				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(c->id,conditionalLikelihood));
			}
			// Initialize conditional likelihood for ancestors
			if (conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()){
				for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
				conditionalLikelihood[dna_c] = 1;
				}				
				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(p->id,conditionalLikelihood));					
			}		
			for (unsigned char dna_p = 0; dna_p < 4; dna_p++){
				partialLikelihood = 0;
				for (unsigned char dna_c = 0; dna_c < 4; dna_c++){
					partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
				}
				conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;
			}
		}
		maxProbState = -1;
		maxProb = 0;
		siteLikelihood = 0; 			
		for (int dna = 0; dna < 4; dna++) {
			currentProb = this->rootProbability[dna]*conditionalLikelihoodMap[-1][dna];
			siteLikelihood += currentProb;
			if (currentProb > maxProb){
				maxProb = currentProb;
				maxProbState = dna;		
			}
		}
		currentLogLikelihood += log(siteLikelihood) * this->siteWeights[site];
		if (maxProbState == -1){
			cout << "check state estimation" << endl;
		}			
		(*this->vertexMap)[-1]->compressedSequence.push_back(maxProbState);
		verticesToVisit.clear();			
		for (int c_id: (*this->vertexMap)[-1]->children_id) {
			if((*this->vertexMap)[c_id]->children_id.size()>0){
				verticesToVisit.push_back((*this->vertexMap)[c_id]);
			}			
		}
		numberOfVerticesToVisit = verticesToVisit.size();
		while (numberOfVerticesToVisit > 0) {
			c = verticesToVisit[numberOfVerticesToVisit-1];		
			verticesToVisit.pop_back();
			numberOfVerticesToVisit -= 1;
			p = (*this->vertexMap)[c->parent_id];
			t = this->GetEdgeLength(p->id,c->id);
			Q_scaled = Q*t;
			P = Q_scaled.exp();
			dna_p = p->compressedSequence[site];
			maxProbState = -1;
			maxProb = 0;
			for (int dna_c = 0; dna_c < 4; dna_c++){
				currentProb = P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
				if (currentProb > maxProb){
					maxProb = currentProb;
					maxProbState = dna_c;
				}					
			}
			c->compressedSequence.push_back(maxProbState);
			for (int c_id : c->children_id){
				if ((*this->vertexMap)[c_id]->children_id.size()>0){
					verticesToVisit.push_back((*this->vertexMap)[c_id]);
					numberOfVerticesToVisit += 1;
				}
			}
		}
	}
	this->logLikelihood = currentLogLikelihood;
}


Matrix4f rootedPhylogeny_tree::ComputeFirstDerivativeOfRateMatrix(int rateCat, int par){
	Matrix4f Q_der;	
	MatrixXd B = (*this->parametersPerRateCategory)[rateCat];
	
	float pi_1; float pi_2; float pi_3; float pi_4;
	float a; float b; float c; float d; float e; float f;
	float g; float h;
	
	pi_1 = B(0,0); pi_2 = B(1,0); pi_3 = B(2,0);
	pi_4 = (1 - pi_1 - pi_2 - pi_3);
	
	a = B(3,0); b = B(4,0); c = B(5,0); d = B(6,0); e = B(7,0);
	f = B(8,0); g = B(9,0); h = B(10,0);
		  
	float a_dot = 0; float b_dot = 0; float c_dot = 0; float d_dot = 0; 
	float e_dot = 0; float f_dot = 0; float g_dot = 0; float h_dot = 0; 
	float i_dot = 0; float j_dot = 0; float k_dot = 0; float l_dot = 0; 
	float D1_dot = 0; float D2_dot = 0; float D3_dot = 0; float D4_dot = 0; 	
	
	if (par == 0){
		// par is pi_1
		i_dot = -(a+b+2*c)/(2*pi_3);
		j_dot = (a+b+c)/pi_4 + (pi_1*(a+b+c) - pi_2*d - pi_3*g)/(pow(pi_4,2));
		k_dot = -a/pi_4 + (pi_2*(d+e+f) - pi_1*a - pi_3*h)/pow(pi_4,2);
		l_dot = -(a+2*c+3*b)/(2*pi_4) + (1 + pi_3*(g+h) - pi_1*(a+2*c+3*b) - pi_2*(d+2*f+3*e) )/(2*pow(pi_4,2));
		D1_dot = 0;
		D2_dot = 0;
		D3_dot = (a+b+2*c)/(2*pi_3);
		D4_dot = (a+b)/(2*pi_4) + (pi_1*(a+b) + pi_2*(d+e) + pi_3*(g+h) -1)/(2*pow(pi_4,2));
	} else if (par == 1){
		// par is pi_2
		i_dot = -(d+e+2*f)/(2*pi_3);
		j_dot = -d/pi_4 + (pi_1*(a+b+c) - pi_2*d - pi_3*g)/(pow(pi_4,2));
		k_dot = (d+e+f)/(pi_4) + (pi_2*(d+e+f) - pi_1*a - pi_3*h)/(pow(pi_4,2));
		l_dot = -(d+2*f+3*e)/(2*pi_4) + (1 + pi_3*(g+h) - pi_1*(a+2*c+3*b) - pi_2*(d+2*f+3*e))/(2*pow(pi_4,2));
		D1_dot = 0;
		D2_dot = 0;
		D3_dot = (d+e+2*f)/(2*pi_3);
		D4_dot = (d+e)/(2*pi_4) + (pi_1*(a+b) + pi_2*(d+e) + pi_3*(g+h) -1)/(2*pow(pi_4,2));
	} else if (par == 2){
		// par is pi_3
		i_dot = -(g+h)/(2*pi_3) + (pi_1*(a+b+2*c) + pi_2*(d+e+2*f) + pi_3*(g+h) -1)/(2*pow(pi_3,2));
		j_dot = -g/pi_4 + (pi_1*(a+b+c) - pi_2*d - pi_3*g)/pow(pi_4,2);
		k_dot = -h/pi_4 + (pi_2*(d+e+f) - pi_1*a - pi_3*h)/(pow(pi_4,2));
		l_dot = (g+h)/(2*pi_4) + (1 + pi_3*(g+h) - pi_1*(a+2*c+3*b) - pi_2*(d+2*f+3*e))/(2*pow(pi_4,2));
		D1_dot = 0;
		D2_dot = 0;
		D3_dot = -(g+h)/(2*pi_3) -(pi_1*(a+b+2*c) + pi_2*(d+e+2*f) - pi_3*(g+h) -1)/(2*pow(pi_3,2));
		D4_dot = (g+h)/(2*pi_4) + (pi_1*(a+b) + pi_2*(d+e) + pi_3*(g+h) -1)/(2*pow(pi_4,2));
	} else if (par == 3){
		// par is a
		a_dot = 1.0;
		i_dot = -pi_1/(2*pi_3);
		j_dot = pi_1/pi_4;
		k_dot = -pi_1/pi_4;
		l_dot = -pi_1/(2*pi_4);
		D1_dot = -1.0;
		D2_dot = 0.0;
		D3_dot = pi_1/(2*pi_3);
		D4_dot = pi_1/(2*pi_4);		
	} else if (par == 4){
		// par is b
		b_dot = 1.0;
		i_dot = -pi_1/(2*pi_3);
		j_dot = pi_1/pi_4;
		k_dot = 0;
		l_dot = -(3*pi_1)/(2*pi_4);
		D1_dot = -1.0;
		D2_dot = 0.0;
		D3_dot = pi_1/(2*pi_3);
		D4_dot = pi_1/(2*pi_4);
	} else if (par == 5){
		// par is c
		c_dot = 1.0;
		i_dot = -pi_1/pi_3;
		j_dot = pi_1/pi_4;
		k_dot = 0.0;
		l_dot = -pi_1/pi_4;
		D1_dot = -1.0;
		D2_dot = 0.0;
		D3_dot = pi_1/pi_3;
		D4_dot = 0.0;
	} else if (par == 6){
		// par is d
		d_dot = 1.0;
		i_dot = -pi_2/(2*pi_3);
		j_dot = -pi_2/pi_4;
		k_dot = pi_2/pi_4;
		l_dot = -pi_2/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = -1.0;
		D3_dot = pi_2/(2*pi_3);
		D4_dot = pi_2/(2*pi_4);
	} else if (par == 7){
		// par is e
		e_dot = 1.0;
		i_dot = -pi_2/(2*pi_3);
		j_dot = 0.0;
		k_dot = pi_2/pi_4;
		l_dot = -(3*pi_2)/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = -1.0;
		D3_dot = pi_2/(2*pi_3);
		D4_dot = pi_2/(2*pi_4);	
	} else if (par == 8){
		// par is f
		f_dot = 1.0;
		i_dot = -pi_2/pi_3;
		j_dot = 0;
		k_dot = pi_2/pi_4;
		l_dot = -pi_2/pi_4;
		D1_dot = 0.0;
		D2_dot = -1.0;
		D3_dot = pi_2/pi_3;
		D4_dot = 0.0;
	} else if (par == 9){
		// par is g
		g_dot = 1.0;
		i_dot = -0.5;
		j_dot = -pi_3/pi_4;
		k_dot = 0.0;
		l_dot = pi_3/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = 0.0;
		D3_dot = -0.5;
		D4_dot = pi_3/(2*pi_4);
	} else if (par == 10){
		// par is h
		h_dot = 1.0;
		i_dot = -0.5;
		j_dot = 0.0;
		k_dot = -pi_3/pi_4;
		l_dot = pi_3/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = 0.0;
		D3_dot = -0.5;
		D4_dot = pi_3/(2*pi_4);
	}	
	Q_der << D1_dot, a_dot, b_dot, c_dot,
			 d_dot, D2_dot, e_dot, f_dot,
			 g_dot, h_dot, D3_dot, i_dot,
			 j_dot, k_dot, l_dot, D4_dot;	
	return (Q_der);
}

Matrix4f rootedPhylogeny_tree::ComputeFirstDerivativeOfRateMatrixDep(int par){
	Matrix4f Q_der = ArrayXXf::Zero(4,4);
	if (par == 0){
		// par is a
		Q_der(0,0) = -1;
		Q_der(0,1) = 1;
		
	} else if (par == 1){
		// par is b
		Q_der(0,0) = -1;
		Q_der(0,2) = 1;
		
	} else if (par == 2){
		// par is c
		Q_der(0,0) = -1;
		Q_der(0,3) = 1;
		
	} else if (par == 3){
		// par is d
		Q_der(1,1) = -1;
		Q_der(1,0) = 1;
		
	} else if (par == 4){
		// par is e
		Q_der(1,1) = -1;
		Q_der(1,2) = 1;
		
	} else if (par == 5){
		// par is f
		Q_der(1,1) = -1;
		Q_der(1,3) = 1;
		
	} else if (par == 6){
		// par is g
		Q_der(2,2) = -1;
		Q_der(2,0) = 1;
		
	} else if (par == 7){
		// par is h
		Q_der(2,2) = -1;
		Q_der(2,1) = 1;
		
	} else if (par == 8){
		// par is i
		Q_der(2,2) = -1;
		Q_der(2,3) = 1;
		
	} else if (par == 9){
		// par is j
		Q_der(3,3) = -1;
		Q_der(3,0) = 1;
		
	} else if (par == 10){
		// par is k
		Q_der(3,3) = -1;
		Q_der(3,1) = 1;		
	}
	return Q_der;
}

Matrix4f rootedPhylogeny_tree::ComputeFirstDerivativeOfMatrixExponential(float t, int rateCat, int par){
	Matrix4f Q; Matrix4f Q_der; MatrixXf Q_aug; MatrixXf Q_aug_scaled;
	MatrixXf P_aug; Matrix4f P_der;
	Q = (*this->rateMatrixPerRateCategory)[rateCat];
	Q_der = this->ComputeFirstDerivativeOfRateMatrix(rateCat, par);
	Q_aug = ArrayXXf::Zero(8,8);
	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			Q_aug(row,col) = Q(row,col);
			Q_aug(row+4,col+4) = Q(row,col);
			Q_aug(row,col+4) = Q_der(row,col);
		}
	}
	Q_aug_scaled = Q_aug*t;
	P_aug = Q_aug.exp();
	for (int row = 0; row < 4; row++) {
		for (int col = 0; col < 4; col++) {
			P_der(row,col) = P_aug(row,col+4);
		}
	}
	return P_der;
}

Matrix4f rootedPhylogeny_tree::ComputeSecondDerivativeOfMatrixExponential(float t, int par_1, int par_2){
	Matrix4f Q; Matrix4f Q_der_1; Matrix4f Q_der_2;
	Matrix4f Q_der_1_2; Matrix4f P_der;
	MatrixXf Q_aug; MatrixXf P_aug;
	Q_der_1_2 = ArrayXXf::Zero(4,4);
	Q = this->rateMatrix;
	Q_der_1 = this->ComputeFirstDerivativeOfRateMatrixDep(par_1);
	Q_der_2 = this->ComputeFirstDerivativeOfRateMatrixDep(par_2);
	Q_aug = ArrayXXf::Zero(16,16);
	for (int row = 0; row < 4; row ++){
		for (int col = 0; col < 4; col++){
			Q_aug(row,col) = Q(row,col);
			Q_aug(row+4,col+4) = Q(row,col);
			Q_aug(row+8,col+8) = Q(row,col);
			Q_aug(row+12,col+12) = Q(row,col);
			Q_aug(row,col+4) = Q_der_1(row,col);
			Q_aug(row+8,col+12) = Q_der_1(row,col);
			Q_aug(row,col+8) = Q_der_2(row,col);
			Q_aug(row+4,col+12) = Q_der_2(row,col);
			Q_aug(row,col+12) = Q_der_1_2(row,col);
		}
	}
	Q_aug = Q_aug*t;
	P_aug = Q_aug.exp();
	for (int row = 0; row < 4; row++){
		for (int col = 0; col < 4; col++){
			P_der(row,col) = P_aug(row,col+12);
		}
	}
	return P_der;
}

MatrixXd rootedPhylogeny_tree::GetFreeParametersDep(Matrix4f Q){
	MatrixXd B;
	B = ArrayXXd::Zero(11,1);
	B(0,0) = Q(0,1);
	B(1,0) = Q(0,2);
	B(2,0) = Q(0,3);
	B(3,0) = Q(1,0);
	B(4,0) = Q(1,2);
	B(5,0) = Q(1,3);
	B(6,0) = Q(2,0);
	B(7,0) = Q(2,1);
	B(8,0) = Q(2,3);
	B(9,0) = Q(3,0);
	B(10,0) = Q(3,1);
	return B;
}

void rootedPhylogeny_tree::SetRateMatrixUsingFreeParametersDep(MatrixXd B){
	Matrix4f Q;
	Q(0,1) = B(0,0);
	Q(0,2) = B(1,0);
	Q(0,3) = B(2,0);
	Q(0,0) = -(Q(0,1) + Q(0,2) + Q(0,3));
	Q(1,0) = B(3,0);
	Q(1,2) = B(4,0);
	Q(1,3) = B(5,0);
	Q(1,1) = -(Q(1,0) + Q(1,2) + Q(1,3));
	Q(2,0) = B(6,0);
	Q(2,1) = B(7,0);
	Q(2,3) = B(8,0);
	Q(2,2) = -(Q(2,0) + Q(2,1) + Q(2,3));
	Q(3,0) = B(9,0);
	Q(3,1) = B(10,0);
	Q(3,2) = 1;
	Q(3,3) = -(Q(3,0) + Q(3,1) + 1);
	this->rateMatrix = Q;
}

void rootedPhylogeny_tree::ComputeMLEOfRateMatrices(){
	MatrixXd B_new;
	MatrixXd JacobianDep; MatrixXd HessianDep; MatrixXd proposedChange;
	proposedChange = ArrayXXd::Zero(11,1);
	JacobianDep = ArrayXXd::Zero(11,1);
	HessianDep = ArrayXXd::Zero(11,11);
	B_new = ArrayXXd::Zero(11,1);	
	double logLikelihood_updated;
	double logLikelihood_current;
	bool continueOptimization = 1;
	float s;	
	MatrixXd B_old = this->GetFreeParametersDep(this->rateMatrix);
	while (continueOptimization){
		this->ComputeLogLikelihoodForFullyLabeledTree();
		logLikelihood_current = this->logLikelihood;
		this->ComputeJacobianForFullyLabeledTree();
		this->ComputeHessianForRateMatrixParametersDep();
		proposedChange = this->HessianDep.inverse()*this->JacobianDep;
		// Improve step size selection
		s = 1.0;
		for (int par = 0; par < 11; par++){
			while (B_old(par,0) - proposedChange(par,0)*s < 0){
				s /= 2.0;
			}
		}
		B_new = B_old - proposedChange*s;
		this->SetRateMatrixUsingFreeParametersDep(B_new);
		this->ComputeLogLikelihoodForFullyLabeledTree();
		logLikelihood_updated = this->logLikelihood;
		while (logLikelihood_updated < logLikelihood_current and abs(logLikelihood_updated - logLikelihood_current) > 0.001){
			s /= 2.0;
			B_new = B_old - proposedChange*s;
			this->SetRateMatrixUsingFreeParametersDep(B_new);
			this->ComputeLogLikelihoodForFullyLabeledTree();
			logLikelihood_updated = this->logLikelihood;
		}
		B_new = B_old - proposedChange*s;
		B_old = B_new;
		this->SetRateMatrixUsingFreeParametersDep(B_old);
		if (abs(logLikelihood_updated - logLikelihood_current) < 0.001) {
			continueOptimization = 0;
		}
	}
}

void rootedPhylogeny_tree::ComputeHessianForRateMatrixParametersDep(){
	Matrix4f P; Matrix4f P_der_1; Matrix4f P_der_2; Matrix4f P_der_1_2;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	int rateCat = 0;
	for (int par_1 = 0; par_1 < 11; par_1 ++){
		for (int par_2 = 0; par_2 < 11; par_2++){
			this->HessianDep(par_1,par_2) = 0;
		}
	}
	Matrix4f Q = this->rateMatrix;
	Matrix4f Q_scaled;
	float t;
	unsigned char dna_p;
	unsigned char dna_c;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			t = this->GetEdgeLength(p->id,c->id);
			Q_scaled = Q*t;
			P = Q_scaled.exp();
			for (int par_1 = 0; par_1 < 11; par_1++){
				P_der_1 = this->ComputeFirstDerivativeOfMatrixExponential(t, rateCat, par_1);
				for (int par_2 = 0; par_2 < 11; par_2++){
					P_der_2 = this->ComputeFirstDerivativeOfMatrixExponential(t, rateCat, par_2);
					P_der_1_2 = this->ComputeSecondDerivativeOfMatrixExponential(t, par_1, par_2);
					for (unsigned int site = 0; site < this->siteWeights.size(); site++){
						dna_p = p->compressedSequence[site];
						dna_c = c->compressedSequence[site];
						this->HessianDep(par_1,par_2) += ((P_der_1_2(dna_p,dna_c)/P(dna_p,dna_c)) - (P_der_1(dna_p,dna_c)*P_der_2(dna_p,dna_c))/pow(P(dna_p,dna_c),2))*this->siteWeights[site];
					}
				}
			}
		}
	}	
}

void rootedPhylogeny_tree::ComputeJacobianForFullyLabeledTree(){
	for (int par = 0; par < 11; par ++){
		this->JacobianDep(par,0) = 0;
	}
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;	
	Matrix4f P; Matrix4f Q; Matrix4f Q_scaled; Matrix4f P_der;
	int rateCat = 0;
	Q = this->rateMatrix;	
	float t;
	unsigned char dna_p;
	unsigned char dna_c;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		c = idPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			t = this->GetEdgeLength(p->id,c->id);
			Q_scaled = Q*t;
			P = Q_scaled.exp();
			for (int par = 0; par < 11; par++){
				P_der = this->ComputeFirstDerivativeOfMatrixExponential(t, rateCat, par);					
				for (unsigned int site = 0; site < this->siteWeights.size(); site++){
					dna_p = p->compressedSequence[site];
					dna_c = c->compressedSequence[site];
					this->JacobianDep(par,0) += (this->siteWeights[site]*P_der(dna_p,dna_c))/P(dna_p,dna_c);
				}
			}				
		}
	}
}

void rootedPhylogeny_tree::ComputeMLEOfRootProbability(){
	array <float,4> prob;
	for (int dna = 0; dna < 4; dna++){
		prob[dna] = 0;
	} 
	rootedPhylogeny_vertex * r = (*this->vertexMap)[-1];
	for (unsigned int site = 0; site < this->siteWeights.size(); site++){
		prob[r->compressedSequence[site]] += this->siteWeights[site];		
	}
	float sum = 0;
	for (int dna = 0; dna < 4; dna++){
		sum += prob[dna];
	}
	for (int dna = 0; dna < 4; dna++){
		prob[dna] /= sum;
	}
	this->rootProbability = prob;
}



Matrix4f rootedPhylogeny_tree::ComputeDerivativeOfMatrixExponentialDep(float t, int par){
	VectorXd x(11);
	x = this->freeParametersExcBaseFreq;
	float pi_1; float pi_2; float pi_3; float pi_4;
	float a; float b; float c; float d; float e; float f;
	float g; float h;float i; float j; float k; float l;
	pi_1 = x[0]; pi_2 = x[1]; pi_3 = x[2];
	pi_4 = (1 - pi_1 - pi_2 - pi_3);
	// Construct Q
	a = x[3]; b = x[4]; c = x[5]; d = x[6]; e = x[7];
	f = x[8]; g = x[9]; h = x[10];

	i = (1-(pi_1*(a+b+2*c)+pi_2*(d+e+2*f)+pi_3*(g+h)))/(2*pi_3);
	j = (pi_1*(a+b+c)-pi_2*d-pi_3*g)/pi_4;
	k = (pi_2*(d+e+f)-pi_1*a-pi_3*h)/pi_4;
	l = (1+pi_3*(g+h)-pi_1*(a+2*c+3*b)-pi_2*(d+2*f+3*e))/(2*pi_4);
	Matrix4f Q;
	Q << -(a+b+c), a, b, c,
		  d, -(d+e+f), e, f,
		  g, h, -(g+h+i), i,
		  j, k, l, -(j+k+l);
		  
	float a_dot = 0; float b_dot = 0; float c_dot = 0; float d_dot = 0; 
	float e_dot = 0; float f_dot = 0; float g_dot = 0; float h_dot = 0; 
	float i_dot = 0; float j_dot = 0; float k_dot = 0; float l_dot = 0; 
	float D1_dot = 0; float D2_dot = 0; float D3_dot = 0; float D4_dot = 0; 	
	if (par == 0){
		// par is pi_1
		i_dot = -(a+b+2*c)/(2*pi_3);
		j_dot = (a+b+c)/pi_4 + (pi_1*(a+b+c) - pi_2*d - pi_3*g)/(pow(pi_4,2));
		k_dot = -a/pi_4 + (pi_2*(d+e+f) - pi_1*a - pi_3*h)/pow(pi_4,2);
		l_dot = -(a+2*c+3*b)/(2*pi_4) + (1 + pi_3*(g+h) - pi_1*(a+2*c+3*b) - pi_2*(d+2*f+3*e) )/(2*pow(pi_4,2));
		D1_dot = 0;
		D2_dot = 0;
		D3_dot = (a+b+2*c)/(2*pi_3);
		D4_dot = (a+b)/(2*pi_4) + (pi_1*(a+b) + pi_2*(d+e) + pi_3*(g+h) -1)/(2*pow(pi_4,2));
	} else if (par == 1){
		// par is pi_2
		i_dot = -(d+e+2*f)/(2*pi_3);
		j_dot = -d/pi_4 + (pi_1*(a+b+c) - pi_2*d - pi_3*g)/(pow(pi_4,2));
		k_dot = (d+e+f)/(pi_4) + (pi_2*(d+e+f) - pi_1*a - pi_3*h)/(pow(pi_4,2));
		l_dot = -(d+2*f+3*e)/(2*pi_4) + (1 + pi_3*(g+h) - pi_1*(a+2*c+3*b) - pi_2*(d+2*f+3*e))/(2*pow(pi_4,2));
		D1_dot = 0;
		D2_dot = 0;
		D3_dot = (d+e+2*f)/(2*pi_3);
		D4_dot = (d+e)/(2*pi_4) + (pi_1*(a+b) + pi_2*(d+e) + pi_3*(g+h) -1)/(2*pow(pi_4,2));
	} else if (par == 2){
		// par is pi_3
		i_dot = -(g+h)/(2*pi_3) + (pi_1*(a+b+2*c) + pi_2*(d+e+2*f) + pi_3*(g+h) -1)/(2*pow(pi_3,2));
		j_dot = -g/pi_4 + (pi_1*(a+b+c) - pi_2*d - pi_3*g)/pow(pi_4,2);
		k_dot = -h/pi_4 + (pi_2*(d+e+f) - pi_1*a - pi_3*h)/(pow(pi_4,2));
		l_dot = (g+h)/(2*pi_4) + (1 + pi_3*(g+h) - pi_1*(a+2*c+3*b) - pi_2*(d+2*f+3*e))/(2*pow(pi_4,2));
		D1_dot = 0;
		D2_dot = 0;
		D3_dot = -(g+h)/(2*pi_3) -(pi_1*(a+b+2*c) + pi_2*(d+e+2*f) - pi_3*(g+h) -1)/(2*pow(pi_3,2));
		D4_dot = (g+h)/(2*pi_4) + (pi_1*(a+b) + pi_2*(d+e) + pi_3*(g+h) -1)/(2*pow(pi_4,2));
	} else if (par == 3){
		// par is a
		a_dot = 1.0;
		i_dot = -pi_1/(2*pi_3);
		j_dot = pi_1/pi_4;
		k_dot = -pi_1/pi_4;
		l_dot = -pi_1/(2*pi_4);
		D1_dot = -1.0;
		D2_dot = 0.0;
		D3_dot = pi_1/(2*pi_3);
		D4_dot = pi_1/(2*pi_4);		
	} else if (par == 4){
		// par is b
		b_dot = 1.0;
		i_dot = -pi_1/(2*pi_3);
		j_dot = pi_1/pi_4;
		k_dot = 0;
		l_dot = -(3*pi_1)/(2*pi_4);
		D1_dot = -1.0;
		D2_dot = 0.0;
		D3_dot = pi_1/(2*pi_3);
		D4_dot = pi_1/(2*pi_4);
	} else if (par == 5){
		// par is c
		c_dot = 1.0;
		i_dot = -pi_1/pi_3;
		j_dot = pi_1/pi_4;
		k_dot = 0.0;
		l_dot = -pi_1/pi_4;
		D1_dot = -1.0;
		D2_dot = 0.0;
		D3_dot = pi_1/pi_3;
		D4_dot = 0.0;
	} else if (par == 6){
		// par is d
		d_dot = 1.0;
		i_dot = -pi_2/(2*pi_3);
		j_dot = -pi_2/pi_4;
		k_dot = pi_2/pi_4;
		l_dot = -pi_2/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = -1.0;
		D3_dot = pi_2/(2*pi_3);
		D4_dot = pi_2/(2*pi_4);
	} else if (par == 7){
		// par is e
		e_dot = 1.0;
		i_dot = -pi_2/(2*pi_3);
		j_dot = 0.0;
		k_dot = pi_2/pi_4;
		l_dot = -(3*pi_2)/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = -1.0;
		D3_dot = pi_2/(2*pi_3);
		D4_dot = pi_2/(2*pi_4);	
	} else if (par == 8){
		// par is f
		f_dot = 1.0;
		i_dot = -pi_2/pi_3;
		j_dot = 0;
		k_dot = pi_2/pi_4;
		l_dot = -pi_2/pi_4;
		D1_dot = 0.0;
		D2_dot = -1.0;
		D3_dot = pi_2/pi_3;
		D4_dot = 0.0;
	} else if (par == 9){
		// par is g
		g_dot = 1.0;
		i_dot = -0.5;
		j_dot = -pi_3/pi_4;
		k_dot = 0.0;
		l_dot = pi_3/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = 0.0;
		D3_dot = -0.5;
		D4_dot = pi_3/(2*pi_4);
	} else if (par == 10){
		// par is h
		h_dot = 1.0;
		i_dot = -0.5;
		j_dot = 0.0;
		k_dot = -pi_3/pi_4;
		l_dot = pi_3/(2*pi_4);
		D1_dot = 0.0;
		D2_dot = 0.0;
		D3_dot = -0.5;
		D4_dot = pi_3/(2*pi_4);
	}
	// Construct Q_der
	Matrix4f Q_der;
	Q_der << D1_dot, a_dot, b_dot, c_dot,
			 d_dot, D2_dot, e_dot, f_dot,
			 g_dot, h_dot, D3_dot, i_dot,
			 j_dot, k_dot, l_dot, D4_dot;
	// Construct Q_aug		 
	MatrixXf Q_aug;		
	Q_aug = ArrayXXf::Zero(8,8);			 
	for (int row = 0; row < 4; row++){
		for (int col = 0; col < 4; col++){
			Q_aug(row,col) = Q(row,col);
			Q_aug(row+4,col+4) = Q(row,col);
			Q_aug(row,col+4) = Q_der(row,col);
		}
	}	
	MatrixXf Q_aug_scaled = Q_aug*t;
	MatrixXf P_aug = Q_aug_scaled.exp();
	Matrix4f P_der;
	for (int row = 0; row < 4; row++){
		for (int col = 0; col < 4; col++){
			P_der(row,col) = P_aug(row,col+4);
		}
	}
	return P_der;
}

void rootedPhylogeny_tree::SetParameters(VectorXd x){
	for (int i = 0; i < 11; i++){
		this->freeParametersExcBaseFreq[i] = x[i];
	}	
	float pi_1; float pi_2; float pi_3; float pi_4;
	pi_1 = x[0]; pi_2 = x[1]; pi_3 = x[2]; pi_4 = 1-(pi_1+pi_2+pi_3);
	this->rootProbability[0] = pi_1;
	this->rootProbability[1] = pi_2;
	this->rootProbability[2] = pi_3;
	// check if x[0] + x[1] + x[2] < 1
	if (pi_4 < 0){
		cout << "check free parameters of stationary distribution" << endl;
	}
	this->rootProbability[3] = pi_4;
	float a; float b; float c; float d; float e; float f;
	float g; float h;float i; float j; float k; float l;
	a = x[3]; b = x[4]; c = x[5]; d = x[6]; e = x[7];
	f = x[8]; g = x[9]; h = x[10];
	i = (1-(pi_1*(a+b+2*c)+pi_2*(d+e+2*f)+pi_3*(g+h)))/(2*pi_3);
	j = (pi_1*(a+b+c)-pi_2*d-pi_3*g)/pi_4;
	k = (pi_2*(d+e+f)-pi_1*a-pi_3*h)/pi_4;
	l = (1+pi_3*(g+h)-pi_1*(a+2*c+3*b)-pi_2*(d+2*f+3*e))/(2*pi_4);
	this->rateMatrix << -(a+b+c), a, b, c,
		  d, -(d+e+f), e, f,
		  g, h, -(g+h+i), i,
		  j, k, l, -(j+k+l);
}

float rootedPhylogeny_tree::ComputeScalingFactor(Matrix4f Q){
	MatrixXf Q_aug = ArrayXXf::Zero(4,5);
	for (int row = 0; row < 4; row++){
		for (int col = 0; col < 4; col++){
			Q_aug(row, col) = Q(row, col);
		}
	}
	for (int row = 0; row < 4; row++){
		Q_aug(row, 4) = 1;
	}	
	MatrixXf b = ArrayXXf::Zero(5,1);
	for (int row = 0; row < 4; row++){
		b(row,0) = 0;
	}
	b(4,0) = 1;	
	MatrixXf pi = ArrayXXf::Zero(1,4);
	pi = Q_aug.transpose().colPivHouseholderQr().solve(b).transpose();
	float scalingFactor = 0;
	for (int i = 0; i < 4; i++){
		scalingFactor -= pi(0,i) * Q(i,i);
	}
	return scalingFactor;
}

MatrixXf rootedPhylogeny_tree::ComputeStationaryDistribution(Matrix4f Q){
	MatrixXf Q_aug = ArrayXXf::Zero(4,5);
	for (int row = 0; row < 4; row++){
		for (int col = 0; col < 4; col++){
			Q_aug(row, col) = Q(row, col);
		}
	}
	for (int row = 0; row < 4; row++){
		Q_aug(row, 4) = 1;
	}	
	MatrixXf b = ArrayXXf::Zero(5,1);
	for (int row = 0; row < 4; row++){
		b(row,0) = 0;
	}
	b(4,0) = 1;	
	MatrixXf pi = ArrayXXf::Zero(4,1);
	pi = Q_aug.transpose().colPivHouseholderQr().solve(b);
	return pi;	
}

void rootedPhylogeny_tree::ComputeMPEstimateOfAncestralSequences(){
	vector <int> preOrderVerticesWithoutLeaves_ids;
	vector <int> leafIds;
	preOrderVerticesWithoutLeaves_ids.push_back(-1);
	vector <rootedPhylogeny_vertex*> verticesToVisit;
	rootedPhylogeny_vertex * v;
	verticesToVisit.push_back(this->root);
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0){
		v = verticesToVisit[numberOfVerticesToVisit-1];
		v->compressedSequence.clear();
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		for (int c_id : v->children_id){
			if ((*this->vertexMap)[c_id]->children_id.size() > 0){
				preOrderVerticesWithoutLeaves_ids.push_back(c_id);
				verticesToVisit.push_back((*this->vertexMap)[c_id]);
				numberOfVerticesToVisit += 1;
			} else {
				leafIds.push_back(c_id);
			}
		}
	}
	int p_id; int numberOfPossibleStates; int pos;
	map <int,vector<unsigned char>> VU;
	map <int,unsigned char> V;
	for (unsigned int site = 0; site < siteWeights.size(); site++){
		VU.clear(); V.clear();
		// Set VU and V for leaves;
		for (int v_id : leafIds){
			V[v_id] = (*this->vertexMap)[v_id]->compressedSequence[site];
			VU[v_id].push_back((*this->vertexMap)[v_id]->compressedSequence[site]);
		}
		// Set VU for ancestors
		for (int v_ind = preOrderVerticesWithoutLeaves_ids.size()-1; v_ind > -1; v_ind--){
			p_id = preOrderVerticesWithoutLeaves_ids[v_ind];
			map <unsigned char, int> dnaCount;
			for (unsigned char dna = 0; dna < 4; dna++){
				dnaCount[dna] = 0;
			}
			for (int c_id : (*this->vertexMap)[p_id]->children_id){
				for (unsigned char dna: VU[c_id]){
					dnaCount[dna] += 1;
				}
			}
			int maxCount = 0;
			for (pair<unsigned char, int> dnaCountPair: dnaCount){
				if (dnaCountPair.second > maxCount){
					maxCount = dnaCountPair.second;
				}
			}			
			for (pair<unsigned char, int> dnaCountPair: dnaCount){
				if (dnaCountPair.second == maxCount){
					VU[p_id].push_back(dnaCountPair.first);					
				}
			}			
		}
		// Set V for ancestors
		for (int v_id : preOrderVerticesWithoutLeaves_ids){
			if (v_id == -1){
			// Set V for root
				if (VU[-1].size()==1){
					V[-1] = VU[-1][0];
				} else {
					numberOfPossibleStates = VU[-1].size();
					uniform_int_distribution<int> distribution(0,numberOfPossibleStates-1);
					pos = distribution(generator);
					V[-1] = VU[-1][pos];
				}				
			} else {
				p_id = (*this->vertexMap)[v_id]->parent_id;
				if (find(VU[v_id].begin(),VU[v_id].end(),V[p_id])==VU[v_id].end()){
					numberOfPossibleStates = VU[v_id].size();
					uniform_int_distribution<int> distribution(0,numberOfPossibleStates-1);
					pos = distribution(generator);
					V[v_id] = VU[v_id][pos];					
				} else {
					V[v_id] = V[p_id];
				}				
			}
			// push states to compressedSequence
			(*this->vertexMap)[v_id]->compressedSequence.push_back(V[v_id]);
		}		
	}	
}

array <double, 4> rootedPhylogeny_tree::GetLikelihoodArray(double elem_0, double elem_1, double elem_2, double elem_3){
	array <double, 4> likelihoodArray;		
	likelihoodArray[0] = elem_0;
	likelihoodArray[1] = elem_1;
	likelihoodArray[2] = elem_2;
	likelihoodArray[3] = elem_3;
	return likelihoodArray;
}

void rootedPhylogeny_tree::ComputeInitialEstimateForRateMatrix(){
	Matrix4f P = ArrayXXf::Zero(4, 4);	
	Matrix4f Q;	
	int dna_p; int dna_c;	
	int rateCategory = 1;	
	rootedPhylogeny_vertex * p = (*this->vertexMap)[-1];	    
	rootedPhylogeny_vertex * c;	
	// Select change-point for rateCategory
	vector <rootedPhylogeny_vertex*> verticesToVisit;
	verticesToVisit.push_back(p);
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0){
		p = verticesToVisit[numberOfVerticesToVisit-1];		
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		for (int c_id: p->children_id){
			c = (*this->vertexMap)[c_id];
			if (c->rateCategory == rateCategory or 1){
				verticesToVisit.push_back(c);
				numberOfVerticesToVisit += 1;
				for (unsigned int site = 0; site < siteWeights.size(); site++){
					dna_p = p->compressedSequence[site];
					dna_c = c->compressedSequence[site];
					P(dna_p,dna_c) += this->siteWeights[site];
				}
			} else {
					cout << "rate category mismatch" << endl;
			}
		}
	}
	float rowSum;
	for (int row = 0; row < 4; row++){
		rowSum = 0;
		for (int col = 0; col < 4; col++){
			rowSum += P(row,col);
		}
		for (int col = 0; col < 4; col++){
			P(row,col)/=rowSum;
		}
	}	
	Q = P.log();	
	float l = Q(3,2);
	Q = Q/l;	
	this->rateMatrix = Q;
	
}

void rootedPhylogeny_tree::ComputeInitialEstimateForFreeParametersIncBaseFreq(){
	Matrix4f P = ArrayXXf::Zero(4,4);
	Matrix4f Q;
	this->ComputeMPEstimateOfAncestralSequences();
	float rowSum;
	for (int row = 0; row < 4; row++){
		rowSum = 0;
		for (int col = 0; col < 4; col++){
			rowSum += P(row,col);
		}
		for (int col = 0; col < 4; col++){
			P(row,col)/=rowSum;
		}
	}
	Q = P.log();
	MatrixXf pi;
	pi = this->ComputeStationaryDistribution(Q);
	float mu = 0;
	for (int i=0; i<4; i++){
		mu -= pi(i,0)*Q(i,i);
	}	
	Q = Q/mu;
	this->initialEstimateForFreeParameters << pi(0,0), pi(1,0), pi(2,0),
											   Q(0,1), Q(0,2), Q(0,3),
											   Q(1,0), Q(1,2), Q(1,3),
											   Q(2,0), Q(2,1);
}


void rootedPhylogeny_tree::InitializeConditionalLikelihoods(){
	unsigned char dna;
	rootedPhylogeny_vertex * v;
	for (unsigned int site = 0; site < this->siteWeights.size(); site++){
		this->root->conditionalLikelihood.push_back(this->GetLikelihoodArray(1,1,1,1));		
		for (int v_id = 0; v_id < this->numberOfObservedSequences; v_id++){
			v = (*this->vertexMap)[v_id];
			dna = v->compressedSequence[site];
			v->conditionalLikelihood.push_back(this->GetLikelihoodArray(0,0,0,0));			
			v->conditionalLikelihood[site][dna] = 1;
		}
		for (unsigned int v_id = this->numberOfObservedSequences; v_id < this->vertexMap->size()-1; v_id++){
			v = (*this->vertexMap)[v_id];
			v->conditionalLikelihood.push_back(this->GetLikelihoodArray(1,1,1,1));
		}
	}
}

void rootedPhylogeny_tree::ResetConditionalLikelihoodsForAncestors(){
	rootedPhylogeny_vertex * v;	
	for (unsigned int site = 0; site < this->siteWeights.size(); site++){
		this->root->conditionalLikelihood[site] = this->GetLikelihoodArray(1,1,1,1);		
		for (unsigned int v_id = this->numberOfObservedSequences; v_id < this->vertexMap->size()-1; v_id++){
			v = (*this->vertexMap)[v_id];
			v->conditionalLikelihood[site] = this->GetLikelihoodArray(1,1,1,1);
		}
	}
}

void rootedPhylogeny_tree::SetSiteWeights(vector <int> siteWeightsToSet){
	this->siteWeights = siteWeightsToSet;
}

void rootedPhylogeny_tree::SetEdgesForPostOrderTreeTraversal(){
	this->edgesForPostOrderTreeTraversal->clear();
	vector < rootedPhylogeny_vertex *> verticesToVisit;
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair: *this->vertexMap){
		idPtrPair.second->timesVisited = 0;		
		if (idPtrPair.second->children_id.size()==0){
			verticesToVisit.push_back(idPtrPair.second);
		}
	}	
	int numberOfVerticesToVisit = verticesToVisit.size();
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;
	while (numberOfVerticesToVisit>0){
		c = verticesToVisit[numberOfVerticesToVisit-1];
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];			
			this->edgesForPostOrderTreeTraversal->push_back(pair<int,int>(p->id,c->id));
			p->timesVisited += 1;
			if (p->timesVisited == (int) p->children_id.size()){
				verticesToVisit.push_back(p);
				numberOfVerticesToVisit += 1;
			}
		}
	}	
}

void rootedPhylogeny_tree::ComputeLogLikelihoodUsingStoredConditionals(){
	array <float,4> pi = this->rootProbability;
	this->logLikelihood = 0;
	double siteLikelihood;
	for (unsigned int site=0; site < this->siteWeights.size(); site++){
		siteLikelihood = 0;
		for (unsigned char dna = 0; dna < 4; dna++){
			siteLikelihood += pi[dna]*this->root->conditionalLikelihood[site][dna];
		}
		this->logLikelihood += log(siteLikelihood) * this->siteWeights[site];
	}
}

void rootedPhylogeny_tree::ComputeLogLikelihood(){
	Matrix4f Q;	
	float scalingFactor;
	map <int, array<double,4>> conditionalLikelihoodMap;
	array <double,4> conditionalLikelihood;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	Matrix4f Q_scaled; Matrix4f P;
	float t;
//	this->SetMinLengthOfEdges();
	double partialLikelihood;	
	for (int i = 0; i < 4; i++){
		if (this->rootProbability[i] < 0){
			cout << "root prob is negative" << endl;
			cout << "Rate matrix is " << endl;
			cout << (*this->rateMatrixPerRateCategory)[this->root->rateCategory] << endl;
			break;
		}
//		cout << "root prob for nt " << i << " is " << this->rootProbability[i] << endl;
	}
	this->logLikelihood = 0;
	double siteLikelihood;
	for (unsigned int site = 0; site < this->siteWeights.size(); site++) {
		for (pair<int,int> edge_ids : *this->edgesForPostOrderTreeTraversal) {
			p = (*this->vertexMap)[edge_ids.first];
			c = (*this->vertexMap)[edge_ids.second];			
			// Initialize conditional likelihood for child if child is a leaf
			if (c->id < this->numberOfObservedSequences) {
				for (int dna = 0; dna < 4; dna++) {
					conditionalLikelihood[dna] = 0;
				}				
				conditionalLikelihood[c->compressedSequence[site]] = 1;
				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(c->id,conditionalLikelihood));
			}
			// Initialize conditional likelihood for parent
			if (conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()) {				
				for (int dna = 0; dna < 4; dna++) {
					conditionalLikelihood[dna] = 1;
				}
				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(p->id,conditionalLikelihood));
			}			
			partialLikelihood = 0;
			t = this->GetEdgeLength(p->id, c->id);
			Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];
			scalingFactor = (*this->scalingFactorForRateCategory)[c->rateCategory];
			Q_scaled = Q*(t/scalingFactor);
			P = Q_scaled.exp();			
			for (unsigned char dna_p = 0; dna_p < 4; dna_p ++) {
				partialLikelihood = 0;
				for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {								
					partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
				}
				conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;
			}
			// Erase conditional likelihood for child
			if(conditionalLikelihoodMap.find(c->id) != conditionalLikelihoodMap.end()) {								
				conditionalLikelihoodMap.erase(c->id);				
			} else {
				cout << "Child not present in conditional likelihood map" << endl;
			}			
		}
		if (conditionalLikelihoodMap.find(-1) == conditionalLikelihoodMap.end()) {
			cout << "Conditional likelihood for root is not computed" << endl;
		}		
		siteLikelihood = 0;
		for (unsigned char dna = 0; dna < 4; dna++) {			
			siteLikelihood += this->rootProbability[dna]*conditionalLikelihoodMap[-1][dna];
		}
		if (site == 0 or site == 10 or site == 20 or site == 30 or site == 40 or site == 50) {
//			cout << "Likelihood: sitelikelihood for site " << site << " is " << siteLikelihood << endl;
		}
		if (siteLikelihood == 0) {		
			cout << "Problem with computing siteLikelihood for site " << site << endl;			
			for (int i = 0; i < 4; i++){
				cout << "Root probability for " << i << " is ";
				cout << this->rootProbability[i] << endl;
				cout << "Conditional likelihood for root for " << i << " is ";
				cout << conditionalLikelihoodMap[-1][i] << endl;
				cout << "Conditional likelihood for vertex 17 for " << i << " is ";
				cout << conditionalLikelihoodMap[17][i] << endl;
				cout << "Rate matrix is" << endl;
				cout << (*this->rateMatrixPerRateCategory)[this->root->rateCategory] << endl;
			}
			
//			cout << "root probability for 0 is " << this->rootProbability[0] << endl;
			break;
		}
		this->logLikelihood += log(siteLikelihood) * this->siteWeights[site];
		// Erase conditional likelihood for root
		conditionalLikelihoodMap.erase(-1);
	}
}

MatrixXd rootedPhylogeny_tree::GetJacobianForRateCategory(int rateCategoryForOptimization) {		
	map <int, array <double,4> > conditionalLikelihoodMap;	
	map <int, array <array <double,4>,11> > derivativeOfConditionalLikelihoodMap;
	array <double,4> conditionalLikelihood;
	array <float,4> pi_root;
	array <array <double,4>,11> derivativeOfConditionalLikelihood;
	double partialLikelihood;
	double partialLikelihood_l;
	double partialLikelihood_r;
	double derivativeOfPartialLikelihood_l;
	double derivativeOfPartialLikelihood_r;	
	double derivativeOfSiteLikelihood;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * v;
	Matrix4f Q;	
	Matrix4f Q_scaled; Matrix4f P;
	Matrix4f Q_l ; Matrix4f Q_r ; 
	Matrix4f Q_l_scaled; Matrix4f Q_r_scaled;
	Matrix4f P_l; Matrix4f P_r;
	Matrix4f P_l_der; Matrix4f P_r_der;
	float t; float t_l; float t_r;		
	double siteLikelihood;
	int numberOfVerticesToVisit;
	rootedPhylogeny_vertex * c_l; rootedPhylogeny_vertex * c_r;
	vector <rootedPhylogeny_vertex*> verticesToVisit;
	// Initialize Jacobian matrix
	MatrixXd Jacobian = MatrixXd::Zero(11,1);
	MatrixXd B_rate_cat_root;
	B_rate_cat_root = (*this->parametersPerRateCategory)[this->root->rateCategory];
	pi_root[0] = B_rate_cat_root(0,0);
	pi_root[1] = B_rate_cat_root(1,0);
	pi_root[2] = B_rate_cat_root(2,0);
	pi_root[3] = (1 - pi_root[0] - pi_root[1] - pi_root[2]);
	for (unsigned int site = 0; site < this->siteWeights.size(); site++) {
		// Compute conditional likelihoods
		for (pair<int,int> edge_ids : *this->edgesForPostOrderTreeTraversal) {
			p = (*this->vertexMap)[edge_ids.first];
			c = (*this->vertexMap)[edge_ids.second];			
			// Set conditional likelihood for child if child is a leaf
			if (c->id < this->numberOfObservedSequences) {
				for (int dna = 0; dna < 4; dna++) {
					conditionalLikelihood[dna] = 0;
				}				
				conditionalLikelihood[c->compressedSequence[site]] = 1;
				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(c->id,conditionalLikelihood));
			}
			// Initialize conditional likelihood for parent
			if(conditionalLikelihoodMap.find(p->id) == conditionalLikelihoodMap.end()){				
				for (int dna = 0; dna < 4; dna++){
					conditionalLikelihood[dna] = 1;
				}
				conditionalLikelihoodMap.insert(pair<int,array<double,4>>(p->id,conditionalLikelihood));
			}
			partialLikelihood = 0;
			t = this->GetEdgeLength(p->id, c->id);
			Q = (*this->rateMatrixPerRateCategory)[c->rateCategory];
			Q_scaled = Q*t;
			P = Q_scaled.exp();
			for (unsigned char dna_p = 0; dna_p < 4; dna_p ++) {
				partialLikelihood = 0;
				for (unsigned char dna_c = 0; dna_c < 4; dna_c ++) {								
					partialLikelihood += P(dna_p,dna_c)*conditionalLikelihoodMap[c->id][dna_c];
				}
				conditionalLikelihoodMap[p->id][dna_p] *= partialLikelihood;
			}
		}
		// Set derivative of conditional likelihood for leaves		
		for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap) {
			v = idPtrPair.second;
			if (v->children_id.size()==0) {
				for (int par = 0; par < 11; par ++) {
					for (int dna = 0; dna < 4; dna++) {		
						derivativeOfConditionalLikelihood[par][dna] = 0;
					}
				}
				derivativeOfConditionalLikelihoodMap.insert(pair<int,array<array<double,4>,11>>(v->id,derivativeOfConditionalLikelihood));
			}			
		}
		verticesToVisit.clear();
		for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap) {
			v = idPtrPair.second;
			v->timesVisited = 0;
			if (v->children_id.size()==0) {
				verticesToVisit.push_back(v);
			}
		}
		numberOfVerticesToVisit = verticesToVisit.size();
		while (numberOfVerticesToVisit > 0) {
			c = verticesToVisit[numberOfVerticesToVisit-1];
			numberOfVerticesToVisit -= 1;
			verticesToVisit.pop_back();
			if (c->id != -1) {
				p = (*this->vertexMap)[c->parent_id];
				p->timesVisited += 1;
				if (p->timesVisited == 2) {
					c_l = (*this->vertexMap)[p->children_id[0]];
					c_r = (*this->vertexMap)[p->children_id[1]];
					t_l = this->GetEdgeLength(p->id,c_l->id);
					t_r = this->GetEdgeLength(p->id,c_r->id);
					Q_l = (*this->rateMatrixPerRateCategory)[c_l->rateCategory];
					Q_r = (*this->rateMatrixPerRateCategory)[c_r->rateCategory];
					Q_l_scaled = Q_l * t_l;
					Q_r_scaled = Q_r * t_r;
					P_l = Q_l_scaled.exp();
					P_r = Q_r_scaled.exp();					
					for (int par = 0; par < 11; par++) {
						for (int dna_p = 0; dna_p < 4; dna_p++) {
							derivativeOfConditionalLikelihood[par][dna_p] = 0;
						}
					}
					for (int dna_p = 0; dna_p < 4; dna_p ++) {
						partialLikelihood_l = 0;
						partialLikelihood_r = 0;
						for (int dna_c = 0; dna_c < 4; dna_c ++) {
							partialLikelihood_l += P_l(dna_p,dna_c)*conditionalLikelihoodMap[c_l->id][dna_c];
							partialLikelihood_r += P_r(dna_p,dna_c)*conditionalLikelihoodMap[c_r->id][dna_c];
						}
						for (int par = 0; par < 11; par++) {
							if (c_l->rateCategory == rateCategoryForOptimization) {
								P_l_der = this->ComputeFirstDerivativeOfMatrixExponential(t_l, c_l->rateCategory, par);
							} else {
								P_l_der = MatrixXf::Zero(4,4);
							}
							if (c_r->rateCategory == rateCategoryForOptimization) {
								P_r_der = this->ComputeFirstDerivativeOfMatrixExponential(t_r, c_r->rateCategory, par);
							} else {
								P_r_der = MatrixXf::Zero(4,4);
							}						
							derivativeOfPartialLikelihood_l = 0;
							derivativeOfPartialLikelihood_r = 0;
							for (int dna_c = 0; dna_c < 4; dna_c ++) {
								derivativeOfPartialLikelihood_l += P_l_der(dna_p,dna_c)*conditionalLikelihoodMap[c_l->id][dna_c];
								derivativeOfPartialLikelihood_l += P_l(dna_p,dna_c)*derivativeOfConditionalLikelihoodMap[c_l->id][par][dna_c];
								derivativeOfPartialLikelihood_r += P_r_der(dna_p,dna_c)*conditionalLikelihoodMap[c_r->id][dna_c];
								derivativeOfPartialLikelihood_r += P_r(dna_p,dna_c)*derivativeOfConditionalLikelihoodMap[c_r->id][par][dna_c];
							}
							derivativeOfConditionalLikelihood[par][dna_p] = derivativeOfPartialLikelihood_l*partialLikelihood_r;
							derivativeOfConditionalLikelihood[par][dna_p] += derivativeOfPartialLikelihood_r*partialLikelihood_l;
						}
					}
					derivativeOfConditionalLikelihoodMap.insert(pair<int,array<array<double,4>,11>>(p->id,derivativeOfConditionalLikelihood));					
					verticesToVisit.push_back(p);
					numberOfVerticesToVisit += 1;
				}
			}
		}		
		if (conditionalLikelihoodMap.find(-1) == conditionalLikelihoodMap.end()) {
				cout << "conditional likelihood for root is not computed" << endl;
		}		
		siteLikelihood = 0;
		for (unsigned char dna = 0; dna < 4; dna++) {
			siteLikelihood += pi_root[dna]*conditionalLikelihoodMap[-1][dna];
		}		
		if (derivativeOfConditionalLikelihoodMap.find(-1) == derivativeOfConditionalLikelihoodMap.end()) {
				cout << "derivative of conditional likelihood for root is not computed" << endl;
		}		
		for	(int par = 0; par < 11; par++) {
			derivativeOfSiteLikelihood = 0;			
			if (this->root->rateCategory == rateCategoryForOptimization){
				if (par == 0){
					// par is pi_1
					derivativeOfSiteLikelihood += conditionalLikelihoodMap[-1][0];
					derivativeOfSiteLikelihood -= conditionalLikelihoodMap[-1][3];
				} else if (par == 1){
					// par is pi_2
					derivativeOfSiteLikelihood += conditionalLikelihoodMap[-1][1];
					derivativeOfSiteLikelihood -= conditionalLikelihoodMap[-1][3];
				} else if (par == 2){
					// par is pi_3
					derivativeOfSiteLikelihood += conditionalLikelihoodMap[-1][2];
					derivativeOfSiteLikelihood -= conditionalLikelihoodMap[-1][3];
				}			
			}			
			for (int dna = 0; dna < 4; dna++) {
				derivativeOfSiteLikelihood += this->rootProbability[dna]*derivativeOfConditionalLikelihoodMap[-1][par][dna];
			}						
			Jacobian(par,0) += (derivativeOfSiteLikelihood/siteLikelihood)*this->siteWeights[site];
		}
		conditionalLikelihoodMap.clear();
		derivativeOfConditionalLikelihoodMap.clear();
	}
	return (Jacobian);
}

void rootedPhylogeny_tree::ComputeJacobianOfLogLikelihoodDep() {
	Matrix4f Q;
	Q = this->rateMatrix;
	for (int par = 0; par < 11; par++){
		this->jacobian_dep[par] = 0;
	}	
	map <int, array <long double,4>> conditionalLikelihoodMap;
	map <int, array <array <long double,4>,11>> derivativeOfConditionalLikelihoodMap;
	array <long double,4> conditionalLikelihood;
	array <array <long double,4>,11> derivativeOfConditionalLikelihood;
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	
	
	vector <rootedPhylogeny_vertex *> verticesToVisit;
	for (pair <int,rootedPhylogeny_vertex*> idPtrPair: (*this->vertexMap)){
		idPtrPair.second->timesVisited = 0;
		if (idPtrPair.second->children_id.size() == 0){
			verticesToVisit.push_back(idPtrPair.second);
		}
	}
	vector <rootedPhylogeny_vertex *> postOrderVerticesToVisit;	
	int numberOfVerticesToVisit = verticesToVisit.size();	
	while (numberOfVerticesToVisit > 0){
		c = verticesToVisit[numberOfVerticesToVisit-1];	
		numberOfVerticesToVisit -= 1;
		verticesToVisit.pop_back();
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			p->timesVisited += 1;			
			if (p->timesVisited == 2){
				postOrderVerticesToVisit.push_back(p);
				verticesToVisit.push_back(p);
				numberOfVerticesToVisit += 1;
			}
		}
	}	


	float t_l; float t_r;
	Matrix4f Q_l_scaled; Matrix4f Q_r_scaled;
	Matrix4f P_l; Matrix4f P_r;
	Matrix4f P_l_der; Matrix4f P_r_der;
	P_l_der = ArrayXXf::Zero(4,4);
	P_r_der = ArrayXXf::Zero(4,4);
	long double partialLikelihood_l; long double partialLikelihood_r;
	long double partialDerivativeOfLikelihood_l; long double partialDerivativeOfLikelihood_r; 
	long double siteLikelihood; long double derivativeOfSiteLikelihood;
	bool debug = 0;
	int c_l_id; int c_r_id;
	if (debug) {
		for (int site = 0; site < this->numberOfObservedSequences; site++){
			conditionalLikelihoodMap.clear();
			derivativeOfConditionalLikelihoodMap.clear();
			for (rootedPhylogeny_vertex * p: postOrderVerticesToVisit){
				c_l_id = p->children_id[0];
				c_r_id = p->children_id[1];
				for (int c_id : p->children_id){
					if (c_id < this->numberOfObservedSequences){
						for (int dna = 0; dna < 4; dna++){
							conditionalLikelihood[dna] = 0;
						}
						conditionalLikelihood[(*this->vertexMap)[c_id]->compressedSequence[site]] = 1;
						conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(c_id,conditionalLikelihood));
						for (int par = 0; par < 11; par ++){
							for (int dna = 0; dna < 4; dna++){					
								derivativeOfConditionalLikelihood[par][dna] = 0;
							}
						}
						derivativeOfConditionalLikelihoodMap.insert(pair<int,array<array<long double,4>,11>>(c_id,derivativeOfConditionalLikelihood));
					} else{
						// check if child is present in conditional and derivative of conditional maps
						
					}
				}			
				// Compute conditional likelihood for parent
				t_l = this->GetEdgeLength(p->id,c_l_id);
				t_r = this->GetEdgeLength(p->id,c_r_id);
				Q_l_scaled = Q*t_l;
				Q_r_scaled = Q*t_r;
				P_l = Q_l_scaled.exp();
				P_r = Q_r_scaled.exp();
				partialLikelihood_l = 0; partialLikelihood_r = 0;
				// Initialize conditional likelihood for parent
				for (int dna_p = 0; dna_p < 4; dna_p++){
					conditionalLikelihood[dna_p] = 0;
				}
				for (int dna_p = 0; dna_p < 4; dna_p++){
					for (int dna_c = 0; dna_c < 4; dna_c++){
						partialLikelihood_l += P_l(dna_p,dna_c)*conditionalLikelihoodMap[c_l_id][dna_c];
					}
					for (int dna_c = 0; dna_c < 4; dna_c++){
						partialLikelihood_r += P_r(dna_p,dna_c)*conditionalLikelihoodMap[c_r_id][dna_c];
					}
					conditionalLikelihood[dna_p] = partialLikelihood_l*partialLikelihood_r;
				}			
				conditionalLikelihoodMap.insert(pair<int,array<long double,4>>(p->id,conditionalLikelihood));			
				// Initialize derivative of conditional likelihood for parent			
				for (int par = 0; par < 11; par ++){
					for (int dna = 0; dna < 4; dna++){		
						derivativeOfConditionalLikelihood[par][dna] = 0;
					}
				}
				derivativeOfConditionalLikelihoodMap.insert(pair<int,array<array<long double,4>,11>>(p->id,derivativeOfConditionalLikelihood));
				// Iterate over parameters
				for (int par = 0; par < 11; par++){
					P_l_der = this->ComputeDerivativeOfMatrixExponentialDep(t_l, par);
					P_r_der = this->ComputeDerivativeOfMatrixExponentialDep(t_r, par);
					for (int dna_p = 0; dna_p < 4; dna_p++){
						partialLikelihood_l = 0;
						for (int dna_c = 0; dna_c < 4; dna_c++){
							partialLikelihood_l += P_l(dna_p,dna_c)*conditionalLikelihoodMap[c_l_id][dna_c];
						}
						partialLikelihood_r = 0;
						for (int dna_c = 0; dna_c < 4; dna_c++){
							partialLikelihood_r += P_r(dna_p,dna_c)*conditionalLikelihoodMap[c_r_id][dna_c];
						}
						partialDerivativeOfLikelihood_l = 0; 
						for (int dna_c = 0; dna_c < 4; dna_c++){
							partialDerivativeOfLikelihood_l += P_l_der(dna_p,dna_c)*conditionalLikelihoodMap[c_l_id][dna_c];
							partialDerivativeOfLikelihood_l += P_l(dna_p,dna_c)*derivativeOfConditionalLikelihoodMap[c_l_id][par][dna_c];
						}
						partialDerivativeOfLikelihood_r = 0;
						for (int dna_c = 0; dna_c < 4; dna_c++){
							partialDerivativeOfLikelihood_r += P_r_der(dna_p,dna_c)*conditionalLikelihoodMap[c_r_id][dna_c];
							partialDerivativeOfLikelihood_r += P_r(dna_p,dna_c)*derivativeOfConditionalLikelihoodMap[c_r_id][par][dna_c];
						}
						derivativeOfConditionalLikelihoodMap[p->id][par][dna_p] = partialDerivativeOfLikelihood_l*partialLikelihood_r;
						derivativeOfConditionalLikelihoodMap[p->id][par][dna_p] += partialLikelihood_l*partialDerivativeOfLikelihood_r;					
					}					
				}
				// Remove conditional likelihood for children
				conditionalLikelihoodMap.erase(c_l_id);
				conditionalLikelihoodMap.erase(c_r_id);				
				derivativeOfConditionalLikelihoodMap.erase(c_l_id);
				derivativeOfConditionalLikelihoodMap.erase(c_r_id);
			}
			siteLikelihood = 0;
			for (int dna = 0; dna < 4; dna++){
				siteLikelihood += this->rootProbability[dna]*conditionalLikelihoodMap[-1][dna];
			}
			if (site == 0 or site == 10 or site == 20 or site == 30 or site == 40 or site == 50){
				cout << "Jacobian: sitelikelihood for site " << site << " is " << siteLikelihood << endl;
			}
//			cout << "site likelihood for site " << site << "is " << siteLikelihood << endl;
			for	(int par = 0; par < 11; par++){
				derivativeOfSiteLikelihood = 0;			
				if (par == 0){
					// par is pi_1
					derivativeOfSiteLikelihood += conditionalLikelihoodMap[-1][0];
					derivativeOfSiteLikelihood -= conditionalLikelihoodMap[-1][3];
				} else if (par == 1){
					// par is pi_2
					derivativeOfSiteLikelihood += conditionalLikelihoodMap[-1][1];
					derivativeOfSiteLikelihood -= conditionalLikelihoodMap[-1][3];
				} else if (par == 2){
					// par is pi_3
					derivativeOfSiteLikelihood += conditionalLikelihoodMap[-1][2];
					derivativeOfSiteLikelihood -= conditionalLikelihoodMap[-1][3];
				}			
				for (int dna = 0; dna < 4; dna++){
					derivativeOfSiteLikelihood += this->rootProbability[dna]*derivativeOfConditionalLikelihoodMap[-1][par][dna];
				}
				this->jacobian_dep[par] += (derivativeOfSiteLikelihood/siteLikelihood)*this->siteWeights[site];
			}
			
		}
	}	
}
	
//	double logLikelihood_current = 0;
//	double logLikelihood_updated;
//	int iter = 0;
//	int maxNumberOfIterations = 20;
//	bool continueOptimization = 1;	
//	this->ComputeMPEstimateOfAncestralSequences();
//	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : (*this->vertexMap)){
//		if (idPtrPair.second->compressedSequence.size() != this->siteWeights.size()){
//			cout << "sequence length mismatch" << endl;
//		}
//	}
//	this->ComputeInitialEstimateForRateMatrix();	
//	this->ScaleEdgeLengths();
//	this->SetMinLengthOfEdges();
//	// Initial log likelihood
//	this->ComputeLogLikelihood();
//	//  cout << "initial loglikelihood is " << this->logLikelihood << endl;
//	while (continueOptimization){
//		iter += 1;
//		logLikelihood_current = this->logLikelihood;
//		this->ComputeMLEOfRootProbability();		
//		for (int i = 0; i < 5; i++){
//			// Optimize rate matrices
//			this->ComputeMLEOfRateMatrices();
//			// Optimize edge lengths
//			this->ComputeMLEOfEdgeLengths();
//			// Compute expected states	
//		}		
//		this->ComputeMAPEstimateOfAncestralSequences();
//		logLikelihood_updated = this->logLikelihood;
//		cout << "updated log likelihood is " << logLikelihood_updated << endl;
//	//	cout << "updated loglikelihood is " << this->logLikelihood << endl;
//		for (pair<int,rootedPhylogeny_vertex*> idPtrPair : (*this->vertexMap)){
//			if (idPtrPair.second->compressedSequence.size() != this->siteWeights.size()){
//				cout << "sequence length mismatch" << endl;
//			}
//		}
//		if (iter < 5 or abs(logLikelihood_updated-logLikelihood_current)>0.001){
//			continueOptimization = 1;
//		} else {
//			continueOptimization = 0;
//		}		
//	}	
	
//	this->ComputeLogLikelihoodForFullyLabeledTree();
//	cout << "current logLikelihood is " << endl;
//	cout << this->logLikelihood << endl;
	// ML estimate of rate matrices
//	this->ComputeMLEOfRateMatrices();
//	this->ComputeLogLikelihoodForFullyLabeledTree();
//	cout << "logLikelihood after optimizing rate matrices is " << endl;
//	cout << this->logLikelihood << endl;
	// ML estimate of edge lengths
//	this->ComputeMLEOfEdgeLengths();
//	this->ComputeLogLikelihoodForFullyLabeledTree();
//	cout << "logLikelihood after optimizing edge lengths is " << endl;
	// check convergence of log-likelihood
//	cout << this->logLikelihood << endl;
//	if (currentLogLikelihood < this->logLikelihood){
//		currentLogLikelihood = this->logLikelihood;
//	}

void rootedPhylogeny_tree::AddRateCategoryForVertex(int v_id, int rateCategory){
	(*this->vertexMap)[v_id]->rateCategory = rateCategory;	
}

void rootedPhylogeny_tree::SetChangePointForRateCategory(int v_id, int rateCategory){
	(*this->changePointForRateCategory)[rateCategory] = v_id;
}

void rootedPhylogeny_tree::SetRateMatrixForRateCategory(Matrix4f Q ,int rateCategory){
//	this->parametersPerRateCategory->insert(pair<int,Matrix4f>(rateCategory,Q));	
//	float scalingFactor = this->ComputeScalingFactor(Q);
//	this->scalingFactorForRateCategory->insert(pair<int,float>(rateCategory,scalingFactor));
}

void rootedPhylogeny_tree::AddStationaryDistributionForCategory(MatrixXf  stationaryDistribution ,int rateCategory){
	this->stationaryDistributionForCategory->insert(pair<int,MatrixXf >(rateCategory,stationaryDistribution));
}

void rootedPhylogeny_tree::WriteAncestralSequences(string sequenceFileName){
	ofstream ancestralSequencesFile;
	ancestralSequencesFile.open(sequenceFileName+".ancestralSequences");
	for (pair<int,rootedPhylogeny_vertex*> vIdAndPtr:(*this->vertexMap)){
		if (boost::algorithm::starts_with(vIdAndPtr.second->name , "h_")){			
			ancestralSequencesFile << ">" << vIdAndPtr.second->name << endl;
			ancestralSequencesFile << EncodeAsDNA(vIdAndPtr.second->sequence) << endl;
		} else {			
		}	
	}	
	ancestralSequencesFile.close();
}

bool rootedPhylogeny_tree::IsEdgeContractionFeasbile(int vertex_id_1, int vertex_id_2){	
	bool edgeContractionIsFeasible;
	if (vertex_id_1 > this->numberOfObservedSequences or vertex_id_2 > this->numberOfObservedSequences){
		edgeContractionIsFeasible = 1;
	} else {
		edgeContractionIsFeasible = 0;
	}
	return (edgeContractionIsFeasible);
}

void rootedPhylogeny_tree::ContractEdge(int idOfVertexToKeep, int idOfVertexToRemove){
	rootedPhylogeny_vertex * k = (*this->vertexMap)[idOfVertexToKeep];
	rootedPhylogeny_vertex * r = (*this->vertexMap)[idOfVertexToRemove];
	float edgeLengthToAdd;
	int ind;	
	if (r->parent_id == k->id){
		// Case 1: k is parent of r
		// Remove r from list of children of k
		ind = find(k->children_id.begin(),k->children_id.end(),r->id) - k->children_id.begin();
		k->children_id.erase(k->children_id.begin()+ind);
		for (int childOfR_id: r->children_id){
			// Set parent of childOfR to k
			(*this->vertexMap)[childOfR_id]->parent_id = k->id;			
			// Add edge from k to childOfR
			k->children_id.push_back(childOfR_id);
			edgeLengthToAdd = (*this->edgeLengthsMap)[pair<int,int>(r->id,childOfR_id)];
			this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(k->id,childOfR_id),edgeLengthToAdd));
			// Remove edge from r to child of r
			this->edgeLengthsMap->erase(pair<int,int>(r->id,childOfR_id));
		}		
	} else {
		// Case 2: k is child of r		
		for (int childOfR_id: r->children_id){
			if (childOfR_id != k->id){
				// Set parent of childOfR to k
				(*this->vertexMap)[childOfR_id]->parent_id = k->id;
				// Add edge from k to childOfR
				k->children_id.push_back(childOfR_id);
				edgeLengthToAdd = (*this->edgeLengthsMap)[pair<int,int>(r->id,childOfR_id)];
				this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(k->id,childOfR_id),edgeLengthToAdd));
				// Remove edge from r to child of r
				this->edgeLengthsMap->erase(pair<int,int>(r->id,childOfR_id));
			}
		}
		// Set parent of k to parent of r
		k->parent_id = r->parent_id;
		// If k is not the root then:
		if (k->parent_id != -1){
			// Add edge from parent of k to k
			edgeLengthToAdd = (*this->edgeLengthsMap)[pair<int,int>(r->parent_id,r->id)];
			this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(k->parent_id,k->id),edgeLengthToAdd));
			// Remove edge from parent of r to r
			this->edgeLengthsMap->erase(pair<int,int>(r->parent_id,r->id));	
		}		
	}
	// Remove all children of r
	r->children_id.clear();
	// Set parent of r to r
	r->parent_id = r->id;
}

void rootedPhylogeny_tree::AddNumberOfObservedSequences(int numberOfObservedSequencesToAdd){
	this->numberOfObservedSequences = numberOfObservedSequencesToAdd;
}

void rootedPhylogeny_tree::AddVertex(int idToAdd, string nameToAdd, vector<unsigned char> sequenceToAdd){
	rootedPhylogeny_vertex * v = new rootedPhylogeny_vertex(idToAdd, nameToAdd, sequenceToAdd);
	(*this->vertexMap)[idToAdd] = v;
	
}

void rootedPhylogeny_tree::AddVertex(int idToAdd, string nameToAdd, vector<unsigned char> sequenceToAdd, vector<unsigned char> compressedSequenceToAdd){
	rootedPhylogeny_vertex * v = new rootedPhylogeny_vertex(idToAdd, nameToAdd, sequenceToAdd);
	v->compressedSequence = compressedSequenceToAdd;
	(*this->vertexMap)[idToAdd] = v;
}


void rootedPhylogeny_tree::RemoveEdges(){	
	for (pair<int,rootedPhylogeny_vertex*> idPtrPair : *this->vertexMap){
		idPtrPair.second->parent_id = -1;
		idPtrPair.second->children_id.clear();
	}
}

void rootedPhylogeny_tree::AddDirectedEdges(vector <pair<int,int>> * directedEdgeList_ptr){
	// Remove edge lengths for root
	for (int child_id : this->root->children_id){
		if(this->edgeLengthsMap->find(pair<int,int>(this->root->id,child_id)) != this->edgeLengthsMap->end()){
			this->edgeLengthsMap->erase(pair<int,int>(this->root->id,child_id));
		}
	}
	this->RemoveEdges();	
	for (pair<int,int> pIdAndcIdPair : *directedEdgeList_ptr){		
		if (this->vertexMap->find(pIdAndcIdPair.second) == this->vertexMap->end()){
			cout << "vertex " << pIdAndcIdPair.second << " not in vertex map" << endl;
		}
		this->AddEdge(pIdAndcIdPair.first,pIdAndcIdPair.second);
		if (pIdAndcIdPair.first == -1){			
		}
	}
	int u_id; int v_id;
	// Set edge lengths for root by midpoint rooting	
	u_id = this->root->children_id[0];	
	v_id = this->root->children_id[1];	
	float edgeLength = this->GetEdgeLength(u_id,v_id);	
	this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(this->root->id,u_id),edgeLength/float(2)));
	this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(this->root->id,v_id),edgeLength/float(2)));
	this->SetEdgesForPostOrderTreeTraversal();
	
}

void rootedPhylogeny_tree::AddEdge(int p_id, int c_id){		
	(*this->vertexMap)[p_id]->AddChild(c_id);
	(*this->vertexMap)[c_id]->AddParent(p_id);	
}

void rootedPhylogeny_tree::ComputeAndSetEdgeLength(int u_id, int v_id){
	float edgeLength;	
	rootedPhylogeny_vertex * u = (*this->vertexMap)[u_id];
	rootedPhylogeny_vertex * v = (*this->vertexMap)[v_id];		
	edgeLength = ComputeNormalizedHammingDistance(&u->sequence,&v->sequence);	
	if (u_id < v_id){
		this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(u_id,v_id),edgeLength));
	} else {
		this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(v_id,u_id),edgeLength));
	}
}

float rootedPhylogeny_tree::GetEdgeLength(int u_id, int v_id){
	float edgeLength;
	if (u_id < v_id){
		edgeLength = (*this->edgeLengthsMap)[pair<int,int>(u_id,v_id)];
	} else {
		edgeLength = (*this->edgeLengthsMap)[pair<int,int>(v_id,u_id)];
	}
	return edgeLength;
}
void rootedPhylogeny_tree::ComputeEdgeLengths(){
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	float edgeLength;
	for (pair<int,rootedPhylogeny_vertex *> idPtrPair: *this->vertexMap){
		c = (*this->vertexMap)[idPtrPair.first];
		if (c->parent_id != -1){
			p = (*this->vertexMap)[c->parent_id];
			edgeLength = ComputeNormalizedHammingDistance(&p->sequence,&c->sequence);
			this->edgeLengthsMap->insert(pair<pair<int,int>,float>(pair<int,int>(p->id,c->id),edgeLength));
		}		
	}	
}

void rootedPhylogeny_tree::WriteEdgeList(string sequenceFileName){
	ofstream edgeListFile;
	edgeListFile.open(sequenceFileName + ".edgeList");
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;
	for (pair<int,rootedPhylogeny_vertex*> vIdVPtrPair : (*this->vertexMap)){
		c = vIdVPtrPair.second;
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];
			edgeListFile << p->name;
			edgeListFile << "\t" << c->name;
			edgeListFile << "\t" << this->GetEdgeLength(p->id,c->id) << endl;
			
		}
	}
	edgeListFile.close();
}

void rootedPhylogeny_tree::ComputeNumberOfDescendants(){
	vector <rootedPhylogeny_vertex *> verticesToVisit;
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;
	for (pair<int,rootedPhylogeny_vertex *> vIdVPtrPair : (*this->vertexMap)){
		c = vIdVPtrPair.second;
		if (c->children_id.size() == 0){
			verticesToVisit.push_back(c);
			c->numberOfDescendants = 1;
		}
	}	
	int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0){
		c = verticesToVisit[numberOfVerticesToVisit-1];		
		p = (*this->vertexMap)[c->parent_id];
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;		
		p->timesVisited += 1;
		p->numberOfDescendants += c->numberOfDescendants;
		if (p->timesVisited == int(p->children_id.size()) and p->parent_id != -1){
			verticesToVisit.push_back(p);
			numberOfVerticesToVisit += 1;
		}
	}
}

void rootedPhylogeny_tree::WriteNewickFile(string sequenceFileName){
	vector <rootedPhylogeny_vertex*> verticesToVisit;
	rootedPhylogeny_vertex * c;
	rootedPhylogeny_vertex * p;	
	float edgeLength;
	for (pair<int,rootedPhylogeny_vertex*> idAndVertex: *this->vertexMap){
		idAndVertex.second->timesVisited = 0;
		if (idAndVertex.second->children_id.size() == 0){
			idAndVertex.second->newickLabel = idAndVertex.second->name;
			verticesToVisit.push_back(idAndVertex.second);
		} else {
			idAndVertex.second->newickLabel = "";
		}
	}
	unsigned int numberOfVerticesToVisit = verticesToVisit.size();
	while (numberOfVerticesToVisit > 0){
		c = verticesToVisit[numberOfVerticesToVisit-1];	
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;		
		if (c->id != -1){
			p = (*this->vertexMap)[c->parent_id];			
			p->timesVisited += 1;
			edgeLength = GetEdgeLength(p->id, c->id);
			if (p->timesVisited == int(p->children_id.size())){
				p->newickLabel += "," + c->newickLabel + ":" + to_string(edgeLength) + ")";
				verticesToVisit.push_back(p);
				numberOfVerticesToVisit += 1;
			} else {
				p->newickLabel += "(" + c->newickLabel + ":" + to_string(edgeLength);
			}			
		}		
	}	
	ofstream newickFile;
	newickFile.open(sequenceFileName+".newick");
	newickFile << this->root->newickLabel << ";" << endl;
	newickFile.close();	
}

//void rootedPhylogeny_tree::WriteNewickFile(string sequenceFileName){
//	vector <rootedPhylogeny_vertex*> verticesToVisit;
//	rootedPhylogeny_vertex * c;
//	rootedPhylogeny_vertex * p;	
//	string newickLabelForTree;
//	float edgeLength;
//	for (pair<int,rootedPhylogeny_vertex*> idAndVertex: *this->vertexMap){
//		idAndVertex.second->timesVisited = 0;
//		if (idAndVertex.second->children_id.size() == 0){
//			idAndVertex.second->newickLabel = idAndVertex.second->name;
//			verticesToVisit.push_back(idAndVertex.second);
//		}
//	}
//	unsigned int numberOfVerticesToVisit = verticesToVisit.size();
//	while (numberOfVerticesToVisit > 0){
//		c = verticesToVisit[numberOfVerticesToVisit-1];	
//		p = (*this->vertexMap)[c->parent_id];
//		verticesToVisit.pop_back();
//		numberOfVerticesToVisit -= 1;
//		if (p->parent_id != -1){
//			p->timesVisited += 1;
////			edgeLength = ComputeNormalizedHammingsDistance(&p->sequence, &c->sequence);
//			edgeLength = GetEdgeLength(p->id, c->id);
//			if (p->timesVisited == int(p->children_id.size())){
//				p->newickLabel += "," + c->newickLabel + ":" + to_string(edgeLength) + ")";
//				verticesToVisit.push_back(p);
//				numberOfVerticesToVisit += 1;
//			} else {
//				p->newickLabel += "(" + c->newickLabel + ":" + to_string(edgeLength);
//			}
//		}
//	}	
//	if (p->children_id.size()==3){
//		p->newickLabel += "(";		
//		for (int i=0; i<2; i++){
//			c = (*this->vertexMap)[p->children_id[i]];
//			edgeLength = GetEdgeLength(p->id, c->id);
////			edgeLength = ComputeNormalizedHammingDistance(&p->sequence, &c->sequence);
//			p->newickLabel += c->newickLabel + ":" + to_string(edgeLength) + ",";			
//		}
//		c = (*this->vertexMap)[p->children_id[2]];
//		edgeLength = GetEdgeLength(p->id, c->id);
////		edgeLength = ComputeNormalizedHammingDistance(&p->sequence, &c->sequence);
//		p->newickLabel += c->newickLabel + ":" + to_string(edgeLength) + ");";		
//	} else {
//		c = (*this->vertexMap)[p->children_id[0]];
//		edgeLength = GetEdgeLength(p->id, c->id);
////		edgeLength = ComputeNormalizedHammingDistance(&p->sequence, &c->sequence);
//		p->newickLabel += "(" + c->newickLabel + ":" + to_string(edgeLength);
//		p->newickLabel += "," + p->name + ":0.00);";		
//	}
//	newickLabelForTree = p->newickLabel;
//	ofstream newickFile;
//	newickFile.open(sequenceFileName+".newick");
//	newickFile << newickLabelForTree << endl;
//	newickFile.close();
//}
void rootedPhylogeny_tree::ContractZeroLengthEdges(){
	// Traverse tree from root to leaves and for each visited vertex p.	
	// Check if edge should be contracted, i.e., at least one of p or c is not an observed vertex.
	// Contract edge (p,v) if length of edge for (p,v) is zero.
	rootedPhylogeny_vertex * p;
	rootedPhylogeny_vertex * c;
	for (pair<int,rootedPhylogeny_vertex*> idVPtrPair : (*this->vertexMap)){
		if (idVPtrPair.second->parent_id == -1) {
			p = idVPtrPair.second;
		}
	}
	vector <rootedPhylogeny_vertex*> verticesToVisit;
	verticesToVisit.push_back(p);
	int numberOfVerticesToVisit = verticesToVisit.size();	
	bool edgeContractedInThisRound = 0;
	float edgeLength;
	while (numberOfVerticesToVisit > 0){
		edgeContractedInThisRound = 0;
		p = verticesToVisit[numberOfVerticesToVisit-1];		
		verticesToVisit.pop_back();
		numberOfVerticesToVisit -= 1;		
		vector <int> idsOfChildrenToVisit;		
		idsOfChildrenToVisit = p->children_id;
		for (int child_id: idsOfChildrenToVisit){
			c = (*this->vertexMap)[child_id];
			edgeLength = (*this->edgeLengthsMap)[pair<int,int>(p->id,c->id)];
			if (edgeLength == 0){
				if (IsEdgeContractionFeasbile(p->id,c->id)){
					edgeContractedInThisRound = 1;
					if (p->id < this->numberOfObservedSequences){
						this->ContractEdge(p->id,c->id);
						verticesToVisit.push_back(p);
						numberOfVerticesToVisit += 1;
					} else {
						this->ContractEdge(c->id,p->id);
						verticesToVisit.push_back(c);
						numberOfVerticesToVisit += 1;
					}
				}
			}
			if (edgeContractedInThisRound){
				break;
			}
		}
		if (!edgeContractedInThisRound){
			for (int child_id: idsOfChildrenToVisit){
				c = (*this->vertexMap)[child_id];
				verticesToVisit.push_back(c);
				numberOfVerticesToVisit += 1;				
			}
		}
	}
}

#endif