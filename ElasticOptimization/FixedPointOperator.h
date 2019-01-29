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
		typedef typename Eigen::VectorXi 			VectorXi;
		typedef typename Eigen::VectorXd 			VectorXd;
		typedef typename Eigen::MatrixXd 			MatrixXd;

		typedef typename Eigen::SparseMatrix<double> 		SparseMatrix;
		typedef typename Eigen::Triplet<double> 		Triplet;

	protected:

		///
		/// The number of vertices in the mesh
		///
		int m_Nv;

		///
		/// The sparse constand energy Hessian matrix
		///
		SparseMatrix m_fixedHessian;

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
		FixedPointOperator( int Nv, double alpha, const VectorXi &target_ID );

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
FixedPointOperator::FixedPointOperator( int Nv, double alpha, const VectorXi &target_ID ) :
	m_Nv( Nv ), m_alpha( alpha ) {
	
	// Populate the sparse Hessian matrix
	// Note: this matrix will remain constant between optimization iterations
	SparseMatrix fH( 3*Nv, 3*Nv );

	std::vector<Triplet> tripletList;
	tripletList.reserve( 3 * target_ID.size() );
	for( int k = 0; k < target_ID.size(); k++ ) {

		int vID = target_ID(k);

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

	Vertex_iterator v;
	for( v = P.vertices_begin(); v != P.vertices_end(); v++ ) {

		if ( v->isTarget() ) {

			Vector3d dx = v->v() - v->tarV();

			// The single vertex contribution to the energy
			EFP += dx.squaredNorm();

		}
	}

	EFP *= this->m_alpha;

	return EFP;

};

///
/// An overloaded function to evaluate the target vertex correspondence energy and its gradient
/// For use with gradient-based methods ( FIRE, L-BFGS, ... )
///
double FixedPointOperator::operator()( Polyhedron &P, VectorXd &grad ) {

	// NOTE: CURRENT GEOMETRY OF THE POLYHEDRON SHOULD BE UP TO DATE 
	double EFP = 0.0;

	Vertex_iterator v;
	for( v = P.vertices_begin(); v != P.vertices_end(); v++ ) {

		if ( v->isTarget() ) {

			int vID = v->id();
			Vector3d dx = v->v() - v->tarV();

			// The single vertex contribution to the energy
			EFP += dx.squaredNorm();

			// The single vertex contribution to the energy gradient
			Vector3d gradEFP = 2.0 * this->m_alpha * dx;

			grad( vID ) += gradEFP(0);
			grad( vID+m_Nv ) += gradEFP(1);
			grad( vID+(2*m_Nv) ) += gradEFP(2);

		}
	}

	EFP *= this->m_alpha;

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

	Vertex_iterator v;
	for( v = P.vertices_begin(); v != P.vertices_end(); v++ ) {

		if ( v->isTarget() ) {

			int vID = v->id();
			Vector3d dx = v->v() - v->tarV();

			// The single vertex contribution to the energy
			EFP += dx.squaredNorm();

			// The single vertex contribution to the energy gradient
			Vector3d gradEFP = 2.0 * this->m_alpha * dx;

			grad( vID ) += gradEFP(0);
			grad( vID+m_Nv ) += gradEFP(1);
			grad( vID+(2*m_Nv) ) += gradEFP(2);

		}
	}

	EFP *= this->m_alpha;

	return EFP;

};

#endif
