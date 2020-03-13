#ifndef mstBackbone_H
#define mstBackbone_H

#include <string>
#include "MST.h"
#include "EM.h"
#include "phylogeny_old.h"
#include "rootedPhylogeny_old.h"
#include "optimizer.h"
#include <tuple>
#include <iostream>
#include <stdio.h>
#include <pthread.h>
#include "SEM.h"
using namespace Eigen;
using namespace std;

class MSTBackbone
{
private:
	default_random_engine generator;
	vector <string> sequenceNames;
	map <string,unsigned char> mapDNAtoInteger;	
	int numberOfLargeEdges;
	int numberOfHiddenVertices = 0;
	int edgeWeightThreshold;	
	string MSTFileName;
	// Remove this
	int ComputeHammingDistance(string seq1, string seq2);
//	float ComputeVertexOrderPerturbedDistance(vector<unsigned char> recodedSeq1, vector<unsigned char> recodedSeq2);	
	int ComputeHammingDistance(vector<unsigned char> recodedSeq1, vector<unsigned char> recodedSeq2);
	int GetEdgeIndex (int vertexIndex1, int vertexIndex2, int numberOfVertices);
	MST_tree * MST_ptr;
	phylogeny_tree * P_ptr;
	rootedPhylogeny_tree * RT_ptr;
	SEM * SEM_mlocalPhylogenyger;
	void ComputeMST(string sequenceFileName);	
	void ComputeVMST(string sequenceFileName);
	void WriteOutputFiles();	
	bool debug;
	int numberOfVerticesInSubtree;
	
public:
	void SetThresholds();
	MSTBackbone(string sequenceFileName, int set_numberOfLargeEdges, int set_edgeWeightThreshold){
		auto start_time = chrono::high_resolution_clock::now();
		debug = 0;		
		bool globalSEM = 1;		
		bool localSEM = 0;		
		mapDNAtoInteger["A"] = 0;
		mapDNAtoInteger["C"] = 1;
		mapDNAtoInteger["G"] = 2;
		mapDNAtoInteger["T"] = 3;
		this->MST_ptr = new MST_tree();
		this->RT_ptr = new rootedPhylogeny_tree();
		this->P_ptr = new phylogeny_tree();
		this->SEM_manager = new localPhylogeny();
		this->SEM_manager->sequenceFileName = sequenceFileName;
		
		vector <int> idsOfVerticesToRemove;
		vector <int> idsOfVerticesToKeep;
		vector <int> idsOfExternalVertices;
		vector <int> idsOfAllVerticesForEM;
		vector <tuple <int, string, vector <unsigned char>>> idAndNameAndSeqTupleForVerticesToAdd;
		
		this->edgeWeightThreshold = set_edgeWeightThreshold;
		this->ComputeMST(sequenceFileName);
		int numberOfInputSequences;
		numberOfInputSequences = MST_ptr->GetNumberOfVertices();
		
		if (globalSEM) {
			this->numberOfLargeEdges = numberOfInputSequences;
		} else if (localSEM) {
			this->numberOfLargeEdges = set_numberOfLargeEdges;
		}
		
		MST_ptr->SetEdgeWeightThreshold(this->edgeWeightThreshold);
		MST_ptr->SetNumberOfLargeEdgesThreshold(this->numberOfLargeEdges);		
				
		ofstream mstBackboneLogFile;
		mstBackboneLogFile.open(sequenceFileName + ".mstbackbone_log");
		mstBackboneLogFile << "Tree depth threshold is set at:\t" << this->numberOfLargeEdges << endl;		
		MSTFileName = sequenceFileName + ".mstb_mst";
		
		MST_ptr->SetCompressedSequencesAndSiteWeightsForInputSequences();
		RT_ptr->siteWeights = MST_ptr->siteWeights;
//		auto time_to_computeMST = chrono::high_resolution_clock::now();			
//		mstBackboneLogFile << "MST computed in " << chrono::duration_cast<chrono::milliseconds>(time_to_computeMST-start_time).count() << " milliseconds\n";
		bool subtreeExtractionPossible = 1;		
		vector <string> names;
		vector <vector <unsigned char> > sequences;
		vector <int> sitePatternWeights;
		vector <vector <int> > sitePatternRepetitions;
		
		cout << "Number of input sequences is " << numberOfInputSequences << endl;
		
		RT_ptr->AddNumberOfObservedSequences(numberOfInputSequences);
		int vertex_id = numberOfInputSequences;
		MST_vertex * v_mst;
		vector <unsigned char> seq_u; vector <unsigned char> seq_v;
		vector <unsigned char> compressed_seq_u; vector <unsigned char> compressed_seq_v;
		int u_ind; int v_ind;
		int u_id; int v_id;
		string u_name; string v_name;
		map <int, int> EMVertexIndToPhyloVertexIdMap;
		phylogeny_vertex * v_phylo;		
		bool resetSubtreeSizeThreshold = 1;
		int subtreeSizeThreshold = MST_ptr->numberOfLargeEdgesThreshold;
		while (subtreeExtractionPossible) {
			if (resetSubtreeSizeThreshold) {
				MST_ptr->numberOfLargeEdgesThreshold = subtreeSizeThreshold;
			}
			MST_ptr->ResetVertexAttributesForSubtreeSearch();
		
 			tie (subtreeExtractionPossible, v_mst) = MST_ptr->GetPtrToVertexSubtendingSubtree();
		//	Select sequence set L = V U K where K is the set of vertices in M that is closest to V.						
			if (subtreeExtractionPossible) {
				EMVertexIndToPhyloVertexIdMap.clear();
				idsOfAllVerticesForEM.clear();
				idsOfExternalVertices.clear();
				idsOfVerticesToRemove.clear();
				idsOfVerticesToKeep.clear();
				for (int v_id: v_mst->idsOfVerticesInSubtree) {	
					idsOfAllVerticesForEM.push_back(v_id);
				}
				for (int v_id: MST_ptr->GetIdsOfClosestUnvisitedVertices(v_mst)){
					idsOfExternalVertices.push_back(v_id);
					idsOfAllVerticesForEM.push_back(v_id);
				}
				tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = MST_ptr->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfAllVerticesForEM);														
				numberOfVerticesInSubtree = int(v_mst->idsOfVerticesInSubtree.size());							
				EM_tree * EM_manager = new EM_tree(int(MST_ptr->vertexMap->size()));	
				EM_manager->AddSequences(&sequences);							
				EM_manager->SetNumberOfVerticesInSubtree(numberOfVerticesInSubtree);
				EM_manager->AddSitePatternWeights(&sitePatternWeights);
				EM_manager->ComputeNJTree();
				EM_manager->SelectIndsOfVerticesOfInterestAndEdgesOfInterest();				
				if (EM_manager->indsOfVerticesOfInterest.size() == 0) {
		//	cout << "increasing subtree size threshold" << endl;
					mstBackboneLogFile << "Increasing subtree size threshold" << endl;
					MST_ptr->numberOfLargeEdgesThreshold *= 2;				
					resetSubtreeSizeThreshold = 0;
				} else {
					resetSubtreeSizeThreshold = 1;
					idAndNameAndSeqTupleForVerticesToAdd.clear();
					for (int v_ind = 0; v_ind < numberOfVerticesInSubtree; v_ind ++) {
						if (find(EM_manager->indsOfVerticesToKeepInMST.begin(),EM_manager->indsOfVerticesToKeepInMST.end(),v_ind) == EM_manager->indsOfVerticesToKeepInMST.end()){
							idsOfVerticesToRemove.push_back(idsOfAllVerticesForEM[v_ind]);
						} else {
							idsOfVerticesToKeep.push_back(idsOfAllVerticesForEM[v_ind]);
						}
					}
					int numberOfRemainingVerticesInMST = MST_ptr->vertexMap->size() - idsOfVerticesToRemove.size() + EM_manager->indsOfVerticesOfInterest.size();
					if (numberOfRemainingVerticesInMST < 4) {
						subtreeExtractionPossible = 0;
					} else {
						EM_manager->PerformEM();
						for (pair <int, int> edge: EM_manager->edgesOfInterest) {
							tie (u_ind, v_ind) = edge;
							if (u_ind < EM_manager->numberOfLeaves) {
								compressed_seq_u = sequences[u_ind];				
								u_id = idsOfAllVerticesForEM[u_ind];
								u_name = (*MST_ptr->vertexMap)[u_id]->name;
							} else {				
								compressed_seq_u = (*EM_manager->maximumLikelihoodAncestralSequences)[u_ind - EM_manager->numberOfLeaves];
								if (EMVertexIndToPhyloVertexIdMap.find(u_ind) == EMVertexIndToPhyloVertexIdMap.end()){
									u_id = vertex_id;
									vertex_id += 1;
									EMVertexIndToPhyloVertexIdMap[u_ind] = u_id;
								} else {
									u_id = EMVertexIndToPhyloVertexIdMap[u_ind];
								}
								u_name = "h_" + to_string(u_id - numberOfInputSequences +1);			
							}
							if (v_ind < EM_manager->numberOfLeaves) {
								compressed_seq_v = sequences[v_ind];	
								v_id = idsOfAllVerticesForEM[v_ind];
								v_name = (* MST_ptr->vertexMap)[v_id]->name;
							} else {
								compressed_seq_v = (*EM_manager->maximumLikelihoodAncestralSequences)[v_ind - EM_manager->numberOfLeaves];
								if (EMVertexIndToPhyloVertexIdMap.find(v_ind) == EMVertexIndToPhyloVertexIdMap.end()) {
									v_id = vertex_id;
									vertex_id += 1;
									EMVertexIndToPhyloVertexIdMap[v_ind] = v_id;
								} else {
									v_id = EMVertexIndToPhyloVertexIdMap[v_ind];
								}
								v_name = "h_" + to_string(v_id - numberOfInputSequences +1);
							}
							seq_u = DecompressSequence(&compressed_seq_u,&sitePatternRepetitions);
							seq_v = DecompressSequence(&compressed_seq_v,&sitePatternRepetitions);
							if (!P_ptr->ContainsVertex(u_id)) { 
								P_ptr->AddVertex(u_id, u_name, seq_u);
								if (u_id < numberOfInputSequences) {
									RT_ptr->AddVertex(u_id, u_name, seq_u, (*MST_ptr->vertexMap)[u_id]->globallyCompressedSequence);	
								} else {
									RT_ptr->AddVertex(u_id, u_name, seq_u);
								}
							}
							if (!P_ptr->ContainsVertex(v_id)) {
								P_ptr->AddVertex(v_id, v_name, seq_v);
								if (v_id < numberOfInputSequences) {									
									RT_ptr->AddVertex(v_id, v_name, seq_v, (*MST_ptr->vertexMap)[v_id]->globallyCompressedSequence);										
								} else {
									RT_ptr->AddVertex(v_id, v_name, seq_v);
								}
							}
							P_ptr->AddEdge(u_id, seq_u, v_id, seq_v);
							RT_ptr->ComputeAndSetEdgeLength(u_id,v_id);
						}
						for (int v_ind: EM_manager->indsOfVerticesOfInterest) {				
							compressed_seq_v = (*EM_manager->maximumLikelihoodAncestralSequences)[v_ind - EM_manager->numberOfLeaves];
							seq_v = DecompressSequence(&compressed_seq_v,&sitePatternRepetitions);			
							v_id = EMVertexIndToPhyloVertexIdMap[v_ind];
							v_phylo = (*(P_ptr->vertexMap))[v_id];					
							idAndNameAndSeqTupleForVerticesToAdd.push_back(tuple<int,string,vector<unsigned char>>(v_id,v_phylo->name,seq_v));									
						}
						cout << "Number of vertices in MST is " << MST_ptr->vertexMap->size() << endl;
						mstBackboneLogFile << "Number of vertices in MST is " << MST_ptr->vertexMap->size() << endl;
						MST_ptr->UpdateMSTWithMultipleExternalVertices(idsOfVerticesToKeep, idsOfVerticesToRemove,idAndNameAndSeqTupleForVerticesToAdd,idsOfExternalVertices);									
					}
				}
					delete EM_manager;
			}
		}
		cout << "Number of vertices in MST is " << MST_ptr->vertexMap->size() << endl;
		mstBackboneLogFile << "Number of vertices in MST is " << MST_ptr->vertexMap->size() << endl;		
		cout << "Performing EM for remaining vertices " << endl;
		mstBackboneLogFile << "Performing EM for remaining vertices " << endl;
		
		// Compute NJ tree for remaining sequences
		EMVertexIndToPhyloVertexIdMap.clear();
		idsOfAllVerticesForEM.clear();
		idsOfExternalVertices.clear();
		idsOfVerticesToRemove.clear();
		idsOfVerticesToKeep.clear();
		
		for (pair <int, MST_vertex *> vIdAndPtr : * MST_ptr->vertexMap) {
			idsOfAllVerticesForEM.push_back(vIdAndPtr.first);
		}
		
		tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = MST_ptr->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfAllVerticesForEM);
		cout << "Num of sequences is " << sequences.size() << endl;
		
//		if (1) {
//			//	cout << "Num of sequences is " << sequences.size() << endl;
//			this->debug = 0;
//			if (this->debug) {
//				SEM_manager->TestSEM();		
//			} else {
//				cout << "Adding sequences" << endl;
//				SEM_manager->AddSequences(sequences);
//				cout << "Adding names" << endl;
//				SEM_manager->AddNames(names);
//				SEM_manager->SetNumberOfVerticesInSubtree(numberOfVerticesInSubtree);
//				SEM_manager->AddSitePatternWeights(sitePatternWeights);
//				SEM_manager->OptimizeTopologyAndParameters();
//			}
//			// Write rooted tree to file
//		}
		
		EM_tree * EM_manager = new EM_tree(int(MST_ptr->vertexMap->size()));
		EM_manager->AddSequences(&sequences);
		EM_manager->SetNumberOfVerticesInSubtree(int(idsOfAllVerticesForEM.size()));
		EM_manager->AddSitePatternWeights(&sitePatternWeights);
		EM_manager->ComputeNJTree();
		EM_manager->SelectIndsOfVerticesOfInterestAndEdgesOfInterest();
		EM_manager->PerformEM();
		
		for (pair<int,int> edge: EM_manager->edgesOfInterest) {
			tie (u_ind, v_ind) = edge;
			if (u_ind < EM_manager->numberOfLeaves) {
				compressed_seq_u = sequences[u_ind];				
				u_id = idsOfAllVerticesForEM[u_ind];
				u_name = (*MST_ptr->vertexMap)[u_id]->name;
			} else {				
				compressed_seq_u = (*EM_manager->maximumLikelihoodAncestralSequences)[u_ind - EM_manager->numberOfLeaves];
				if (EMVertexIndToPhyloVertexIdMap.find(u_ind) == EMVertexIndToPhyloVertexIdMap.end()){
					u_id = vertex_id;
					vertex_id +=1;
					EMVertexIndToPhyloVertexIdMap[u_ind] = u_id;
				} else {
					u_id = EMVertexIndToPhyloVertexIdMap[u_ind];
				}
				u_name = "h_" + to_string(u_id - numberOfInputSequences +1);			
			}
			if (v_ind < EM_manager->numberOfLeaves) {
				compressed_seq_v = sequences[v_ind];	
				v_id = idsOfAllVerticesForEM[v_ind];
				v_name = (*MST_ptr->vertexMap)[v_id]->name;
			} else {
				compressed_seq_v = (*EM_manager->maximumLikelihoodAncestralSequences)[v_ind - EM_manager->numberOfLeaves];
				if (EMVertexIndToPhyloVertexIdMap.find(v_ind) == EMVertexIndToPhyloVertexIdMap.end()){
					v_id = vertex_id;
					vertex_id +=1;
					EMVertexIndToPhyloVertexIdMap[v_ind] = v_id;
				} else {
					v_id = EMVertexIndToPhyloVertexIdMap[v_ind];
				}
				v_name = "h_" + to_string(v_id - numberOfInputSequences +1);
			}			
			seq_u = DecompressSequence(&compressed_seq_u,&sitePatternRepetitions);
			seq_v = DecompressSequence(&compressed_seq_v,&sitePatternRepetitions);
			if (!P_ptr->ContainsVertex(u_id)) {
				P_ptr->AddVertex(u_id, u_name, seq_u);
				if (u_id < numberOfInputSequences) {
					RT_ptr->AddVertex(u_id, u_name, seq_u, (*MST_ptr->vertexMap)[u_id]->globallyCompressedSequence);	
				} else {
					RT_ptr->AddVertex(u_id, u_name, seq_u);
				}			
			}
			if (!P_ptr->ContainsVertex(v_id)) {
				P_ptr->AddVertex(v_id, v_name, seq_v);
				if (v_id < numberOfInputSequences) {
					RT_ptr->AddVertex(v_id, v_name, seq_v, (*MST_ptr->vertexMap)[v_id]->globallyCompressedSequence);	
				} else {
					RT_ptr->AddVertex(v_id, v_name, seq_v);
				}
			}
			P_ptr->AddEdge(u_id, seq_u, v_id, seq_v);
			RT_ptr->ComputeAndSetEdgeLength(u_id,v_id);
		}
		P_ptr->sequenceLength = int(seq_u.size());
		EMVertexIndToPhyloVertexIdMap.clear();
		delete EM_manager;
		auto current_time = std::chrono::high_resolution_clock::now();
		cout << "CPU time used for computing unrooted phylogeny is " << chrono::duration_cast<chrono::seconds>(current_time-start_time).count() << " second(s)\n";
		mstBackboneLogFile << "CPU time used for computing unrooted phylogeny is " << chrono::duration_cast<chrono::seconds>(current_time-start_time).count() << " second(s)\n";
//      Root tree by fitting a general Markov model
//		cout << "Selecting edge for rooting\n";
//		mstBackboneLogFile << "Selecting edge for rooting\n";
//		int rootId;
//		bool atLeastOneEdgeChecked = 0;
		pair <int, int> edgeForRooting;
//		double optimalBIC = pow(10,10);
//		double optimalLogLikelihood = pow(10,10);
//		float optimalThreshold = 0;
		string edgeListFileName;
		string scriptFileName;
		string pathForModelSelection = "/TL/euresist_phylodynamics/work/Projects/MSTBasedForests/scripts/modelSelection";
		scriptFileName = sequenceFileName + "_modelSelection.sh";
		ofstream scriptFile;
		scriptFile.open(scriptFileName);
		int edgeInd = 0;
		for (pair <int,int> edge: *P_ptr->edgeList) {
			edgeInd += 1; 									
			P_ptr->RootTreeAtEdge(edge.first, edge.second);			
			RT_ptr->AddDirectedEdges(P_ptr->directedEdgeList);
			edgeListFileName = sequenceFileName + "_rootedAtEdge_" + to_string(edgeInd);
			RT_ptr->WriteEdgeList(edgeListFileName);
			scriptFile << pathForModelSelection << " " << sequenceFileName << " " << edgeListFileName << ".edgeList" << endl;
//			RT_ptr->WriteNewickFile(edgeListFileName);
		}
		scriptFile.close();
		// write script file for model selection
//		cout << "Rooting tree under the general Markov model" << endl;		
		// sort edges using ML score?
		// store edge for rooting, and optimal threshold 
		// pass edge for rooting and receive BIC, optimal threshold
//		cout << "Performing model selection by fitting " << endl;
//		int numberOfEdgesTried = 0;		
//		for (pair<int,int> edge: *P_ptr->edgeList){			
//			cout << "BIC for rooting tree at edge " << edge.first << "\t" << edge.second << " is ";
//			numberOfEdgesTried += 1;			
//			P_ptr->RootTreeAtEdge(edge.first, edge.second);
//			// Add directed edges to rooted tree
//			RT_ptr->AddDirectedEdges(P_ptr->directedEdgeList);	
//			// Perform Model Selection
//			RT_ptr->PerformModelSelectionUsingNelderMead();					
//			if (optimalBIC > RT_ptr->optimal_BIC or numberOfEdgesTried ==1){
//				optimalBIC = RT_ptr->optimal_BIC;				
//				optimalThreshold = RT_ptr->optimalThreshold;
//				edgeForRooting = edge;
//			}			
//			cout << RT_ptr->optimal_BIC << endl;
//		}
//		
//		P_ptr->RootTreeAtEdge(edgeForRooting.first, edgeForRooting.second);
//		RT_ptr->AddDirectedEdges(P_ptr->directedEdgeList);
//		RT_ptr->OptimizeModelParametersForAGivenThresholdUsingNelderMead(optimalThreshold);
//		RT_ptr->ComputeLogLikelihood();
//		RT_ptr->WriteRateCategoryPerVertexAndModelParameters(sequenceFileName);
//		cout << "Log-likelihood of optimal rooted tree is:\t" << setprecision(10) << RT_ptr->logLikelihood << endl;
//		mstBackboneLogFile << "Log-likelihood of rooted tree is:\t" << setprecision(10) << RT_ptr->logLikelihood << endl;
//		cout << "BIC of optimal rooted tree is:\t" << setprecision(10) << optimalBIC << endl;
//		mstBackboneLogFile << "BIC of rooted tree is:\t" << setprecision(10) << optimalBIC << endl;
//		RT_ptr->WriteEdgeList(sequenceFileName);
////		cout << "Log-likelihood of rooted tree is:\t" << setprecision(10) << P_ptr->maxLogLikelihood << endl;
////		mstBackboneLogFile << "Log-likelihood of rooted tree is:\t" << setprecision(10) << P_ptr->maxLogLikelihood << endl;
//		cout << "Writing rooted tree in newick format\n";
//		mstBackboneLogFile << "Writing rooted tree in newick format\n";	
//		RT_ptr->WriteNewickFile(sequenceFileName);
//		cout << "Contracting zero-length edges\n";
//		mstBackboneLogFile << "Contracting zero-length edges\n";
//		RT_ptr->ContractZeroLengthEdges();
//		cout << "Contracting short edges such that BIC is minimized\n";
//		mstBackboneLogFile << "Contracting short edges such that BIC is minimized\n";
//		cout << "Writing generally labeled rooted tree in edge list format\n";
//		mstBackboneLogFile << "Writing generally labeled rooted tree in edge list format\n";
//		RT_ptr->WriteEdgeList(sequenceFileName);
//		cout << "Writing estimated ancestral sequences\n";
//		mstBackboneLogFile << "Writing estimated ancestral sequences\n";
//		RT_ptr->WriteAncestralSequences(sequenceFileName);	
		auto end_time = std::chrono::high_resolution_clock::now();
		cout << "Total CPU time used is " << chrono::duration_cast<chrono::seconds>(end_time-start_time).count() << " second(s)\n";
		mstBackboneLogFile << "Total CPU time used is " << chrono::duration_cast<chrono::seconds>(end_time-start_time).count() << " second(s)\n";
		mstBackboneLogFile.close();
	}
	~MSTBackbone(){
		delete this->MST_ptr;
		delete this->P_ptr;
		delete this->RT_ptr;
		delete this->SEM_manager;
//		delete this->baseFreqForInputSequences;
	}
};

void MSTBackbone::SetThresholds() {
	
}


int MSTBackbone::GetEdgeIndex (int vertexIndex1, int vertexIndex2, int numberOfVertices){
	int edgeIndex;
	edgeIndex = numberOfVertices*(numberOfVertices-1)/2;
	edgeIndex -= (numberOfVertices-vertexIndex1)*(numberOfVertices-vertexIndex1-1)/2;
	edgeIndex += vertexIndex2 - vertexIndex1 - 1;
	return edgeIndex;
}

int MSTBackbone::ComputeHammingDistance(string seq1, string seq2) {
	int hammingDistance = 0;
	for (unsigned int i=0;i<seq1.length();i++){
		if (seq1[i] != seq2[i]){
			hammingDistance+=1;
		}		
	}
	return (hammingDistance);
};

int MSTBackbone::ComputeHammingDistance(vector<unsigned char> recodedSeq1, vector<unsigned char> recodedSeq2) {
	int hammingDistance = 0;
	for (unsigned int i=0;i<recodedSeq1.size();i++) {
		if (recodedSeq1[i] != recodedSeq2[i]) {
			hammingDistance+=1;
		}		
	}
	return (hammingDistance);
};

//float MSTBackbone::ComputeVertexOrderPerturbedDistance(vector<unsigned char> recodedSeq1, vector<unsigned char> recodedSeq2) {
//	float hammingDistance = 0;
//	for (unsigned int i=0;i<recodedSeq1.size();i++){
//		if (recodedSeq1[i] != recodedSeq2[i]){
//			hammingDistance+=1;
//		}		
//	}
//	return (hammingDistance);
//};

void MSTBackbone::ComputeMST(string sequenceFileName) {	
//	cout << *this->baseFreqForInputSequences << endl;
//	for (int dna = 0; dna < 4; dna++){
//		(*this->baseFreqForInputSequences)[dna] = 0;
//	}
	vector <unsigned char> recodedSequence;
	
	unsigned int site = 0;
	ifstream inputFile(sequenceFileName.c_str());
	string seqName;
	string seq = "";	
	int vertex_ind = 0;
	for(string line; getline(inputFile, line );) {
		if (line[0]=='>') {
			if (seq != "") {
				sequenceNames.push_back(seqName);
				for (char const dna: seq) {
					recodedSequence.push_back(mapDNAtoInteger[string(1,toupper(dna))]);					
					site += 1;
					}
				MST_ptr->AddVertex(vertex_ind,seqName,recodedSequence);
				vertex_ind += 1;
				recodedSequence.clear();
			} 
			seqName = line.substr(1,line.length());
			seq = "";
			site = 0;			
		}
		else {
			seq += line ;
		}		
	}		
	for (char const dna: seq) {
		recodedSequence.push_back(mapDNAtoInteger[string(1,toupper(dna))]);		
		site += 1;
	}	
	MST_ptr->AddVertex(vertex_ind,seqName,recodedSequence);
	recodedSequence.clear();
	sequenceNames.push_back(seqName);
	inputFile.close();
	
	int numberOfVertices = vertex_ind+1;		
	const int numberOfEdges = numberOfVertices*(numberOfVertices-1)/2;	
	
//	for (unsigned char dna = 0; dna < 4; dna++){
//		(*this->baseFreqForInputSequences)[dna] /= numberOfVertices*sequenceLength;			
//	}
//	for (pair<string,unsigned char> p : mapDNAtoInteger){
//		cout << "pi(" << p.first << ") = " << (*this->baseFreqForInputSequences)[p.second] << endl;
//	}
	
	int * weights;
	weights = new int [numberOfEdges];
		
	int edgeIndex = 0;
	for (int i=0; i<numberOfVertices; i++) {
		for (int j=i+1; j<numberOfVertices; j++) {			
			weights[edgeIndex] = ComputeHammingDistance((*MST_ptr->vertexMap)[i]->sequence,(*MST_ptr->vertexMap)[j]->sequence);
			edgeIndex += 1;
		}
	}
	typedef pair <int,int > E;

	E * edges;
	edges = new E [numberOfEdges];
	edgeIndex = 0;
	for (int i=0; i<numberOfVertices; i++) {
		for (int j=i+1; j<numberOfVertices; j++) {
			edges[edgeIndex] = E(i,j);
			edgeIndex += 1;
		}
	}
	typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS, boost::property<boost::vertex_distance_t, int>, boost::property < boost::edge_weight_t, int> > Graph;
	Graph g(edges, edges + numberOfEdges, weights, numberOfVertices);
	cout << "Number of vertices in MST are " << num_vertices(g) << endl;
	vector < boost::graph_traits < Graph >::vertex_descriptor >  p(num_vertices(g));
	prim_minimum_spanning_tree(g, &p[0]);
	delete[] edges;		
	int edgeCount = 0;
	ofstream MSTFile;
	MSTFile.open(sequenceFileName+".mst");
	for (size_t u = 0; u != p.size(); u++) {
		if (p[u] != u) {
			edgeCount += 1;
			if (u < p[u]) {
				edgeIndex = GetEdgeIndex(u,p[u],numberOfVertices);
			} else {
				edgeIndex = GetEdgeIndex(p[u],u,numberOfVertices);
			}
			MST_ptr->AddEdge(u, p[u], weights[edgeIndex]);
			MSTFile << (*MST_ptr->vertexMap)[u]->name << "\t" << (*MST_ptr->vertexMap)[p[u]]->name << "\t" << weights[edgeIndex] << endl;
		}
	}
	MSTFile.close();
	delete[] weights;
};

void MSTBackbone::ComputeVMST(string sequenceFileName) {
	vector <unsigned char> recodedSequence;
	ifstream inputFile(sequenceFileName.c_str());
	string seqName;
	string seq = "";	
	string seq2;
	int vertex_ind = 0;
	for(string line; getline(inputFile, line );) {
		if (line[0]=='>') {
			if (seq != "") {
				sequenceNames.push_back(seqName);
				for (char const dna: seq) {
					recodedSequence.push_back(mapDNAtoInteger[string(1,toupper(dna))]);
					}
				this->MST_ptr->AddVertex(vertex_ind,seqName,recodedSequence);
				vertex_ind += 1;
				recodedSequence.clear();
			} 
			seqName = line.substr(1,line.length());
			seq = "";			
		}
		else {
			seq += line ;
		}		
	}		
	for (char const dna: seq) {				
		recodedSequence.push_back(mapDNAtoInteger[string(1,toupper(dna))]);	
		}
	MST_ptr->AddVertex(vertex_ind,seqName,recodedSequence);
	recodedSequence.clear();
	sequenceNames.push_back(seqName);
	inputFile.close();
	
	int numberOfVertices = vertex_ind+1;		
	const int numberOfEdges = numberOfVertices*(numberOfVertices-1)/2;
	
	int * weights;
	weights = new int [numberOfEdges];
		
	int edgeIndex = 0;
	for (int i=0; i<numberOfVertices; i++) {
		for (int j=i+1; j<numberOfVertices; j++) {
			weights[edgeIndex] = ComputeHammingDistance((*MST_ptr->vertexMap)[i]->sequence,(*MST_ptr->vertexMap)[j]->sequence);
			edgeIndex += 1;
		}
	}
	typedef pair <int,int > E;

	E * edges;
	edges = new E [numberOfEdges];
	edgeIndex = 0;
	for (int i=0; i<numberOfVertices; i++) {
		for (int j=i+1; j<numberOfVertices; j++) {
			edges[edgeIndex] = E(i,j);
			edgeIndex += 1;
		}
	}
	typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS, boost::property<boost::vertex_distance_t, int>, boost::property < boost::edge_weight_t, int> > Graph;
	Graph g(edges, edges + numberOfEdges, weights, numberOfVertices);
	vector < boost::graph_traits < Graph >::vertex_descriptor >  p(num_vertices(g));
	prim_minimum_spanning_tree(g, &p[0]);
	delete[] edges;		
	int edgeCount = 0;
	ofstream MSTFile;
	MSTFile.open(sequenceFileName+".mst");
	for (size_t u = 0; u != p.size(); ++u) {
		if (p[u] != u){
			edgeCount += 1;
			if (u < p[u]){
				edgeIndex = GetEdgeIndex(u,p[u],numberOfVertices);
			} else {
				edgeIndex = GetEdgeIndex(p[u],u,numberOfVertices);
			}
			MST_ptr->AddEdge(u, p[u], weights[edgeIndex]);
			MSTFile << u << "\t" << p[u] << "\t" << weights[edgeIndex] << endl;
		}
	}
	MSTFile.close();
	delete[] weights;
};
#endif
