/*!
 * 	\file FixedPointOperator.h
 * 	\brief A class to calculate the target vertex correspondence energy and its derivatives
 *
 * 	\author Dillon Cislo
 * 	\date 01/11/2019
 *
 */

#ifndef _FIXED_POINT_OPERATOR_H_
#define _FIXED_POINT_OPERATOR_H_

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include "../ElasticMesh/ElasticItems.h"

///
/// A class to calculate the target vertex correspondence energy and its derivatives
///
class FixedPointOperator {

	private:

		typedef typename CGAL::Simple_cartesian<double> 		Kernel;
		typedef typename CGAL::Polyhedron_3<Kernel, ElasticItems> 	Polyhedron;

		typedef typename Polyhedron::Vertex_iterator 		Vertex_iterator;
		typedef typename Polyhedron::Vertex_handle 		Vertex_handle;

		typedef typename Eigen::RowVector3d 			RowVector3d;
		typedef typename Eigen::Vector3d 			Vector3d;
		typedef typename Eigen::VectorXd 			VectorXd;
		typedef typename Eigen::MatrixXd 			MatrixXd;

		typedef typename Eigen::SparseMatrix<double> 		SparseMatrix;
		typedef typename Eigen::Triplet<double> 		Triplet;

	protected:

		///
		/// The sparse constand energy Hessian matrix
		///
		SparseMatrix m_fixedHessian;

		///
		/// A vector of pointers to the fixed vertices
		///
		std::vector<Vertex_handle> m_fixedVertices;

		///
		/// The target locations of the fixed vertices
		///
		MatrixXd m_targetPositions;

		///
		/// Fixed point constraint tolerance coefficient
		/// 
		double m_alpha;

	public:
		///
		/// Null constructor
		///
		FixedPointOperator() {};

		///
		/// Default constructor
		///
		FixedPointOperator( int Nv, double alpha,
				std::vector<Vertex_handle> fV, MatrixXd tP );

		///
		/// An overloaded function to calculate the target vertex correspondence
		/// energy.  For use when it it only necessary to calculate the energy
		/// and not any of its derivatives
		///
		double operator()( Polyhedron &P );

		///
		/// An overloaded function to evaluate the target vertex correspondence
		/// energy and its gradient. For use with gradient-based methods ( FIRE, 
		/// L-BFGS, ... )
		///
		double operator()( Polyhedron &P, VectorXd &grad );

		///
		/// An overloaded function to evaluate the target vertex correspondence
		/// energy, its gradient, and its Hessian matrix.
		/// For use with a fully implemented Newton method.
		///
		double operator()( Polyhedron &P, VectorXd &grad, SparseMatrix &hess );

};

///
/// Default constructor
///
FixedPointOperator::FixedPointOperator( int Nv, double alpha,
		std::vector<Vertex_handle> fV, MatrixXd tP ) :
	m_alpha( alpha ), m_fixedVertices( fV ), m_targetPositions( tP ) {
	
	// Populate the sparse Hessian matrix
	// Note: this matrix will remain constant between optimization iterations
	SparseMatrix fH( 3*Nv, 3*Nv );

	std::vector<Triplet> tripletList;
	tripletList.reserve( 3 * fV.size() );
	for( int k = 0; k < fV.size(); k++ ) {

		int vID = fV[k]->id();

		tripletList.push_back( Triplet( vID, vID, 2.0 * alpha ) );
		tripletList.push_back( Triplet( vID+Nv, vID+Nv, 2.0 * alpha ) );
		tripletList.push_back( Triplet( vID+(2*Nv), vID+(2*Nv), 2.0 * alpha ) );

	}

	fH.setFromTriplets( tripletList.begin(), tripletList.end() );

	this->m_fixedHessian = fH;

};

///
/// An overloaded function to calculate the target vertex correspondence energy.  For use when
/// it is only necessary to calculate the energy and not any of its derivatives
///
double FixedPointOperator::operator()( Polyhedron &P ) {

	// NOTE: CURRENT GEOMETRY OF THE POLYHEDRON SHOULD BE UP TO DATE
	double EFP = 0.0;

	for( int k = 0; k < this->m_fixedVertices.size(); k++ ) {

		RowVector3d x0 = this->m_targetPositions.row( k );
		RowVector3d x = this->m_fixedVertices[k]->v().transpose();
		RowVector3d dx = x - x0;

		// The single vertex contribution to the energy
		EFP += this->m_alpha * dx.squaredNorm();

	}

	return EFP;

};

///
/// An overloaded function to evaluate the target vertex correspondence energy and its gradient
/// For use with gradient-based methods ( FIRE, L-BFGS, ... )
///
double FixedPointOperator::operator()( Polyhedron &P, VectorXd &grad ) {

	// NOTE: CURRENT GEOMETRY OF THE POLYHEDRON SHOULD BE UP TO DATE
	double EFP = 0.0;

	int Nv = P.size_of_vertices();
	for( int k = 0; k < this->m_fixedVertices.size(); k++ ) {

		int vID = this->m_fixedVertices[k]->id();

		RowVector3d x0 = this->m_targetPositions.row( k );
		RowVector3d x = this->m_fixedVertices[k]->v().transpose();
		RowVector3d dx = x - x0;

		// The single vertex contribution to the energy
		EFP += this->m_alpha * dx.squaredNorm();

		// The single vertex contribution to the energy gradient
		RowVector3d gradEFP = 2.0 * m_alpha * dx;

		grad( vID ) += gradEFP(0);
		grad( vID+Nv ) += gradEFP(1);
		grad( vID+(2*Nv) ) += gradEFP(2);

	}

	return EFP;

};

///
/// An overloaded function to evaulate the target vertex correspondence energy, its gradient
/// and its Hessian matrix.  For use with a fully implemented Newton method
///
double FixedPointOperator::operator()( Polyhedron &P, VectorXd &grad, SparseMatrix &hess ) {

	// NOTE: CURRENT GEOMETRY OF THE POLYHEDRON SHOULD BE UP TO DATE
	hess += this->m_fixedHessian;

	double EFP = 0.0;

	int Nv = P.size_of_vertices();
	for( int k = 0; k < this->m_fixedVertices.size(); k++ ) {

		int vID = this->m_fixedVertices[k]->id();

		RowVector3d x0 = this->m_targetPositions.row( k );
		RowVector3d x = this->m_fixedVertices[k]->v().transpose();
		RowVector3d dx = x - x0;

		// The single vertex contribution to the energy
		EFP += this->m_alpha * dx.squaredNorm();

		// The single vertex contribution to the energy gradient
		RowVector3d gradEFP = 2.0 * m_alpha * dx;

		grad( vID ) += gradEFP(0);
		grad( vID+Nv ) += gradEFP(1);
		grad( vID+(2*Nv) ) += gradEFP(2);

	}

	return EFP;

};

#endif
