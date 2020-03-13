#ifndef optimizer_H
#define optimizer_H
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include "cppoptlib/meta.h"
#include "cppoptlib/problem.h"
#include "cppoptlib/solver/bfgssolver.h"
#include "rootedPhylogeny_old.h"

using namespace std;
using namespace Eigen;
using namespace cppoptlib;
using Eigen::VectorXd;

class Optimizer : public Problem<double> {
  public:
    using typename cppoptlib::Problem<double>::Scalar;
    using typename cppoptlib::Problem<double>::TVector;
	rootedPhylogeny_tree * RT_ptr;	
	void SetPtrToRootedPhylogeny(rootedPhylogeny_tree* ptrToSet) {
		RT_ptr = ptrToSet;
	}
    double value(const TVector &x) {	
//		RT_ptr->SetRateMatrixForRateCategory(Q, RT_ptr->rateCategoryForOptimization);
		RT_ptr->SetParameters(x);
		RT_ptr->ComputeLogLikelihood();
		cout << x << endl;
		cout << RT_ptr->logLikelihood << endl;
        return (-1 * RT_ptr->logLikelihood);
    }
    void gradient(const TVector &x, TVector &grad) {
//		RT_ptr->SetRateMatrixForRateCategory(Q, RT_ptr->rateCategoryForOptimization);
		RT_ptr->SetParameters(x);
//		RT_ptr->GetJacobianForRateCategory(rateCat);
		// Set gradient of conditional likelihoods
		for (int par = 0; par < 11; par ++){
			grad[par]  = -1 * RT_ptr->JacobianDep(par,0);
		}
    }
};


#endif