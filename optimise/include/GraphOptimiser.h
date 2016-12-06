/**
 * @file GraphOptimiser.h
 * @brief Contains functions for multiple views pairwise geometry estimation.
 * @detail Implements functions needed for robust nonlinear least-squares batch SLAM.
 * @copyright Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
 * @authors Tariq Abuhashim, Nicolo' Genesio
 * @email t.abuhashim@gmail.com, nicogene@hotmail.it
 * @date Nov 2016
 * @acknowledgement This research has received funding from the European Union’s 
 * Seventh Framework Programme for research, technological development and demonstration 
 * under grant agreement No. 611909(KoroiBot).
 * @license Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 */

#ifndef GRAPHOPTIMISER_H
#define GRAPHOPTIMISER_H

#include "RecoverMoments.h"
#include <unistd.h> // sleep command
#include "eigen3/Eigen/Sparse"
#include "eigen3/Eigen/Dense"
typedef Eigen::SparseMatrix<double> SpMat ;
typedef Eigen::Triplet<double> T ;

class GraphOptimiser {

public:

	// constructor
	GraphOptimiser();

	// destructor
	~GraphOptimiser();

	void initialise_a_constraint(
			const std::vector<double>& edge,
			const std::vector<double>& z,
			const std::vector<double>& R);

	void compute_gate(
			double *gate,
			double *x,
			double *s,
			double *xs,
			long unsigned int ncols);

private:

	struct constraint {
		const std::vector<double> edge; /* edge by connected nodes */
		const std::vector<double> z; /* measurements */
		const std::vector<double> R; /* measurements noise covariance */
		constraint( ); // Do not use this default constructor (structure will be empty)
		constraint(
				std::vector<double> a,
				std::vector<double> b,
				std::vector<double> c):
					edge(a), z(b), R(c){} // Member initialization in a constructor
	};

	std::vector<constraint> constraints;

	void constraint_model(
			std::vector<double>& zs,
			double *xs,
			int i,
			int j);

	void constraint_jacobian_nodepair_model(
			Eigen::SparseMatrix<double> &H,
			double *xs,
			int i,
			int j);

	void compute_rotation_angles(
			std::vector<double> &R,
			std::vector<double> &x);

	void compute_rotation_matrix(
			std::vector<double> &R,
			std::vector<double> &x);

	void compute_rotation_matrix_derivatives(
			std::vector<double> &R,
			std::vector<double> &dRa,
			std::vector<double> &dRb,
			std::vector<double> &dRc,
			std::vector<double> &x);

	void derivative_R1tR2(
			std::vector<double> &dRRdaij,
			std::vector<double> &x1,
			std::vector<double> &x2);

	void derivative_R2w(
			std::vector<double> &dadRR,
			std::vector<double> &R);

	void compute_constraint_model_derivatives(
			std::vector<double> &h,
			std::vector<double> &x1,
			std::vector<double> &x2);

	void multiply(
			std::vector<double> &C,
			std::vector<double> &A,
			std::vector<double> &B,
			int nrows,
			int ncomm,
			int ncols);

	void transpose(
			std::vector<double> &At,
			std::vector<double> &A,
			int nrows, int ncols);

	void display_matrix(
			std::vector<double> &At,
			int nrows, int ncols);

	void pi_to_pi(
			double &x);

};

#endif /* GRAPHOPTIMISER_H */

/*
 * Best wishes;
 * Tariq Abuhashim, for iCub Facility
 * September, 2016
 * Thanks to Tim Bailey and Lorenzo Natale
 */
