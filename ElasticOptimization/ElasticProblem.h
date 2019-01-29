/*!
 * 	\file ElasticProblem.h
 * 	\brief A problem structure for the optimization of a thin shell's elastic energy
 *
 * 	A problem structure for the optimization of a thin shell's elastic energy using
 * 	The problem structure is equally applicable to problems in the classical nonlinear
 * 	elasticity of thin shells and to the incompatible elasticity of non-Euclidean shells.
 *
 * 	\author Dillon Cislo
 * 	\date 01/11/2018
 *
 */

#ifndef _ELASTIC_PROBLEM_H_
#define _ELASTIC_PROBLEM_H_

#include <iostream>
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "StretchOperator.h"
#include "HingeOperator.h"
#include "BendOperator.h"
#include "FixedPointOperator.h"
#include "../ElasticMesh/ElasticUpdater.h"

///
/// A problem structure for the optimization of a thin shell's elastic energy using
/// gradient based methods ( FIRE, L-BFGS, ... )
///
class ElasticProblem {

	private:

		typedef typename CGAL::Simple_cartesian<double> 		Kernel;
		typedef typename CGAL::Polyhedron_3<Kernel, ElasticItems> 	Polyhedron;

		typedef typename Polyhedron::Vertex_handle 	Vertex_handle;
		typedef typename Polyhedron::Face_handle 	Face_handle;
		typedef typename Polyhedron::Halfedge_handle 	Halfedge_handle;
		typedef typename Polyhedron::Vertex_iterator	Vertex_iterator;

		typedef typename Eigen::VectorXi 		VectorXi;
		typedef typename Eigen::VectorXd 		VectorXd;
		typedef typename Eigen::MatrixXd 		MatrixXd;

		typedef typename Eigen::SparseMatrix<double> 	SparseMatrix;

	protected:

		///
		/// The polyhedral mesh representing the elastic shell
		///
		Polyhedron m_P;

		///
		/// A boolean determining whether any vertices are given target locations
		///
		bool anyFixed;

		///
		/// The stretch operator used to calculate the stretching energy
		///
		StretchOperator m_SO;

		///
		/// The hinge operator used to construct hinge-based derivative objects
		///
		HingeOperator m_HO;

		///
		/// The bend operator used to calculate the bending energy
		///
		BendOperator m_BO;

		///
		/// The fixed point operator used to calculate the target vertex
		/// correspondence energy
		///
		FixedPointOperator m_FPO;

		///
		/// The elastic updater used to update the polyhedron's physical geometry
		/// between optimization iterations
		///
		ElasticUpdater m_EU;


	public:

		///
		/// Null constructor
		///
		ElasticProblem() {};

		///
		/// Constructor in the case that no vertices have target locations
		/// Note that we expect the polyhedron object to have a fully updated
		/// target geometry prior to construction of the problem structure
		///
		ElasticProblem( Polyhedron &P, double h, double nu ) : m_P( P ) {

				this->anyFixed = false;

				this->m_SO = StretchOperator( nu );
				this->m_BO = BendOperator( h, nu );

		};

		///
		/// Constructor in the case that there are user supplied target locations
		/// for a subset of the mesh vertices.  Note that we expect the polyhedron
		/// object to have a fully updated target geometry prior to the construction
		/// of the problem structure
		ElasticProblem( Polyhedron &P, double h, double nu, double alpha,
				const VectorXi &target_ID );

		///
		/// Evaluate the energy
		///
		double operator()( const VectorXd &x );
		

		///
		/// Evaulate the energy and the energy gradient
		///
		double operator()( const VectorXd &x, VectorXd &grad );

		///
		/// Evaulate the energy, the energy gradient, and the energy Hessian
		///
		double operator()( const VectorXd &x, VectorXd &grad, SparseMatrix &hess );


};

///
/// Constructor in the case that there are user supplied target locations for a subset of
/// the mesh vertices.  Note that we expect the polyhedron object to have a fully updated
/// target geometry prior to the construction of the problem structure
///
ElasticProblem::ElasticProblem( Polyhedron &P, double h, double nu, double alpha,
	       const VectorXi &target_ID ) : m_P( P ) {

	this->anyFixed = true;

	this->m_SO = StretchOperator( nu );
	this->m_BO = BendOperator( h, nu );

	int Nv = this->m_P.size_of_vertices();
	this->m_FPO = FixedPointOperator( Nv, alpha, target_ID );

};

///
/// Evaluate the energy
///
double ElasticProblem::operator()( const VectorXd &x ) {

	// Update the current mesh geometry
	this->m_EU.updateCurrentGeometry( this->m_P, x );

	// Calculate the stretching energy
	double Estretch = this->m_SO( this->m_P );

	// Calculate the bending energy
	double Ebend = this->m_BO( this->m_P );

	double Etotal = Estretch + Ebend;

	// OPTIONAL: Calculate the target vertex correspondence energy
	if ( this->anyFixed ) {

		Etotal += this->m_FPO( this->m_P );

	}

	return Etotal;

};

///
/// Evaluate the energy and the energy gradient
///
double ElasticProblem::operator()( const VectorXd &x, VectorXd &grad ) {

	// Update the current mesh geometry
	this->m_EU.updateCurrentGeometry( this->m_P, x );

	// Update hinge stencils for bending energy calculation
	this->m_HO.constructGradient( this->m_P );

	// Reset the global gradient vector
	grad = VectorXd::Zero( grad.size() );

	// Calculate stretching energy
	double Estretch = this->m_SO( this->m_P, grad );

	// std::cout << "Sgrad = " << grad.norm();

	// Calculate bending energy
	double Ebend = this->m_BO( this->m_P, grad );

	// std::cout << " + Bgrad = " << grad.norm() << std::endl;

	double Etotal = Estretch + Ebend;

	// OPTIONAL: Calculate the target vertex correspondence energy
	if ( this->anyFixed ) {

		Etotal += this->m_FPO( this->m_P, grad );

	}

	

	return Etotal;

};

///
/// Evaluate the energy, the energy gradient, and the energy Hessian.
///
double ElasticProblem::operator()( const VectorXd &x, VectorXd &grad, SparseMatrix &hess ) {

	// Update the current mesh geometry
	this->m_EU.updateCurrentGeometry( this->m_P, x );

	// Update hinge stencils for bending energy calculation
	this->m_HO.constructGradientAndHessian( this->m_P );

	// Reset the global gradient vector
	grad = VectorXd::Zero( grad.size() );
	
	// Create stand-in Hessian matrix
	SparseMatrix spHess( grad.size(), grad.size() );

	// Calculate stretching energy
	double Estretch = this->m_SO( this->m_P, grad, spHess );

	// Calculate bending energy
	double Ebend = this->m_BO( this->m_P, grad, spHess );

	double Etotal = Estretch + Ebend;

	// OPTIONAL: Calculate the target vertex correspondence energy
	if ( this->anyFixed ) {

		Etotal += this->m_FPO( this->m_P, grad, spHess );

	}

	// Set the global Hessian matrix from stand-in
	hess = spHess;

	return Etotal;

};

#endif
