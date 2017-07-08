
#include "ceres/ceres.h"
#include "glog/logging.h"
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

// Consider the problem of finding the minimum of the function (1/2)(10-x)^2
// The first step is to write a functor
struct CostFunctor {
   template <typename T>
   bool operator()(const T* const x, T* residual) const {
     residual[0] = T(10.0) - x[0]; // f(x) = 10 - x;
     return true;
   }
};
// The important thing to note here is that operator() is a templated method, which assumes that all its inputs and outputs are of some type T. The use of templating here allows Ceres to call CostFunctor::operator<T>(), with T=double when just the value of the residual is needed, and with a special type T=Jet when the Jacobians are needed.

int main(int argc, char** argv) {
	//google::InitGoogleLogging(argv[0]);

	// The variable to solve for with its initial value.
  	double initial_x = 5.0;
  	double x = initial_x;

  	// Build the problem.
  	Problem problem;

  	// Set up the only cost function (also known as residual). This uses
  	// auto-differentiation to obtain the derivative (jacobian).
  	CostFunction* cost_function =
      	new AutoDiffCostFunction<CostFunctor, 1, 1>(new CostFunctor);
  	problem.AddResidualBlock(cost_function, NULL, &x);

  	// Run the solver!
  	Solver::Options options;
  	options.linear_solver_type = ceres::DENSE_QR;
  	options.minimizer_progress_to_stdout = true;
  	Solver::Summary summary;
  	Solve(options, &problem, &summary);

  	std::cout << summary.BriefReport() << "\n";
  	std::cout << "x : " << initial_x
            	<< " -> " << x << "\n";
  	return 0;
}
