/*!
 * 	\file BendOperator.h
 * 	\brief A function class to calculate the bending energy and its derivatives
 *
 * 	\author Dillon Cislo
 * 	\date 01/10/2019
 *
 */

#ifndef _BEND_OPERATOR_H_
#define _BEND_OPERATOR_H_

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include "../ElasticMesh/ElasticItems.h"
#include "FlapMap.h"

///
/// A function class to calculate the bending energy and its derivatives
///
class BendOperator {

	private:

		typedef typename CGAL::Simple_cartesian<double> 		Kernel;
		typedef typename CGAL::Polyhedron_3<Kernel, ElasticItems> 	Polyhedron;

		typedef typename Polyhedron::Facet_iterator 			Facet_iterator;

		typedef typename Polyhedron::Halfedge_handle 			Halfedge_handle;
		typedef typename Polyhedron::Vertex_handle 			Vertex_handle;
		typedef typename Polyhedron::Face_handle 			Face_handle;

		typedef typename Eigen::Vector3d 		Vector3d;
		typedef typename Eigen::VectorXd 		VectorXd;
		typedef typename Eigen::Matrix3d 		Matrix3d;
		typedef typename Eigen::MatrixXd 		MatrixXd;
		typedef typename Eigen::Matrix<double,1,12> 	HingeGrad;
		typedef typename Eigen::Matrix<double,12,12> 	HingeHess;
		typedef typename Eigen::Matrix<double,1,18> 	FlapGrad;
		typedef typename Eigen::Matrix<double,18,18> 	FlapHess;

		typedef typename Eigen::SparseMatrix<double> 	SparseMatrix;
		typedef typename Eigen::Triplet<double> 	Triplet;

	protected:

		///
		/// The thickness of the elastic shell
		///
		double m_h;

		///
		/// Poisson's ratio
		///
		double m_nu;

	public:

		///
		/// Null constructor
		///
		BendOperator() {};

		///
		/// Default constructor
		///
		BendOperator( double h, double nu ) : m_h( h ), m_nu( nu ) {};

		///
		/// An overloaded function to evaluate the bending energy.
		/// For use when it is only necessary to calculate the energy and
		/// not any of its derivatives
		///
		double operator()( Polyhedron &P );

		///
		/// An overloaded function to evaluate the bending energy and its gradient
		/// For use with gradient based-methods ( FIRE, L-BFGS, ... )
		///
		double operator()( Polyhedron &P, VectorXd &grad );

		///
		/// An overloaded function to evaluate the bending energy, its gradient,
		/// and its Hessian matrix.  For use with a fully implemented Newton method.
		///
		double operator()( Polyhedron &P, VectorXd &grad, SparseMatrix &hess );

	protected:

		///
		/// An overloaded function that maps local gradient quantities to the
		/// global gradient vector
		///
		void mapLocalToGlobal( Face_handle f,
				const FlapGrad &LGrad, VectorXd &GGrad );

		///
		/// An overloaded function that maps local Hessian quantities to the vector
		/// of Eigen-style which will be used to construct the global Hessian
		/// matrix of the bending energy
		///
		void mapLocalToGlobal( Face_handle f, int Nv,
				const FlapHess &LHess, std::vector<Triplet> &GTrip );

};

///
/// An overlaoded function that maps local gradient quantities to the global gradient vector
///
void BendOperator::mapLocalToGlobal( Face_handle f,
		const FlapGrad &LGrad, VectorXd &GGrad ) {

	// NOTE: VERTEX IDS SHOULD ALREADY BE UP TO DATE!
	int vi, vj, vk;
	int va, vb, vc;

	vi = f->halfedge()->next()->vertex()->id();
	vj = f->halfedge()->prev()->vertex()->id();
	vk = f->halfedge()->vertex()->id();

	va = f->halfedge()->opposite()->next()->vertex()->id();
	vb = f->halfedge()->next()->opposite()->next()->vertex()->id();
	vc = f->halfedge()->prev()->opposite()->next()->vertex()->id();

	if ( ( GGrad.size() % 3 ) != 0 ) {
		std::cerr << "Gradient vector is improperly sized!" << std::endl;
		return;
	}

	int Nv = GGrad.size() / 3;

	// Vi ---------------------------------------------------
	GGrad( vi ) += LGrad(0);
	GGrad( vi + Nv ) += LGrad(1);
	GGrad( vi + (2*Nv) ) += LGrad(2);

	// Vj ---------------------------------------------------
	GGrad( vj ) += LGrad(3);
	GGrad( vj + Nv ) += LGrad(4);
	GGrad( vj + (2*Nv) ) += LGrad(5);

	// Vk ---------------------------------------------------
	GGrad( vk ) += LGrad(6);
	GGrad( vk + Nv ) += LGrad(7);
	GGrad( vk + (2*Nv) ) += LGrad(8);

	// Va ----------------------------------------------------
	GGrad( va ) += LGrad(9);
	GGrad( va + Nv ) += LGrad(10);
	GGrad( va + (2*Nv) ) += LGrad(11);

	// Vb ----------------------------------------------------
	GGrad( vb ) += LGrad(12);
	GGrad( vb + Nv ) += LGrad(13);
	GGrad( vb + (2*Nv) ) += LGrad(14);

	// Vc -----------------------------------------------------
	GGrad( vc ) += LGrad(15);
	GGrad( vc + Nv ) += LGrad(16);
	GGrad( vc + (2*Nv) ) += LGrad(17);

};

///
/// An overloaded function that maps local Hessian quantities to the vector of
/// Eigen-style triplets which will be used to construct the global Hessian
/// matrix of the bending energy
///
void BendOperator::mapLocalToGlobal( Face_handle f, int Nv,
		const FlapHess &LHess, std::vector<Triplet> &GTrip ) {

	// NOTE: VERTEX IDS SHOULD ALREADY BE UP TO DATE!
	int vi, vj, vk;
	int va, vb, vc;

	vi = f->halfedge()->next()->vertex()->id();
	vj = f->halfedge()->prev()->vertex()->id();
	vk = f->halfedge()->vertex()->id();

	va = f->halfedge()->opposite()->next()->vertex()->id();
	vb = f->halfedge()->next()->opposite()->next()->vertex()->id();
	vc = f->halfedge()->prev()->opposite()->next()->vertex()->id();

	std::vector<int> vID = { vi, vj, vk, va, vb, vc };
	for( int m = 0; m < 6; m++ ) {

		int mID = vID[m];
		int M = 3 * m;

		for ( int n = 0; n < 6; n++ ) {

			int nID = vID[n];
			int N = 3 * n;

			// X-coordinate --------------------------------------------
			GTrip.push_back( Triplet( mID, nID, LHess(M,N) ) );
			GTrip.push_back( Triplet( mID, nID+Nv, LHess(M,N+1) ) );
			GTrip.push_back( Triplet( mID, nID+(2*Nv), LHess(M,N+2) ) );

			// Y-coordinate --------------------------------------------
			GTrip.push_back( Triplet( mID+Nv, nID, LHess(M+1,N) ) );
			GTrip.push_back( Triplet( mID+Nv, nID+Nv, LHess(M+1,N+1) ) );
			GTrip.push_back( Triplet( mID+Nv, nID+(2*Nv), LHess(M+1,N+2) ) );
			
			// Z-coordinate --------------------------------------------
			GTrip.push_back( Triplet( mID+(2*Nv), nID, LHess(M+2,N) ) );
			GTrip.push_back( Triplet( mID+(2*Nv), nID+Nv, LHess(M+2,N+1) ) );
			GTrip.push_back( Triplet( mID+(2*Nv), nID+(2*Nv), LHess(M+2,N+2) ) );

		}

	}

};

///
/// An overloaded function to calculate the bending energy
/// For use when it is only necessary to calculate the energy and
/// not any of its derivatives
///
double BendOperator::operator()( Polyhedron &P ) {

	// NOTE: CURRENT GEOMETRY OF POLYHEDRON AND HINGE DERIVATIVES SHOULD BE UP TO DATE
	double Ebend = 0.0; // The global bending energy

	Facet_iterator f;
	for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

		double tarA = f->tarFaceArea();
		Vector3d hBar = f->hBar();
		Matrix3d xi = f->xi();

		// Re-format face quantities for energy calculation
		Vector3d Phi;
		Halfedge_handle he = f->halfedge();
		for( int k = 0; k < 3; k++ ) {

			Phi(k) = ( he->Phi() ) - ( he->tarPhi() );
			he = he->next();

		}

		// Construct the traces and trace gradients
		double trB = 0.0;
		double trB2 = 0.0;
		for( int i = 0; i < 3; i++ ) {

			trB += Phi(i) / hBar(i);
			trB2 += xi(i,i) * Phi(i) * Phi(i);

			for( int j = (i+1); j < 3; j++ ) {

				trB2 += xi(i,j) * Phi(i) * Phi(j);

			}

		}

		// Single stencil contribution to the energy
		Ebend += m_h * m_h * tarA * ( m_nu * trB * trB + (1-m_nu) * trB2 ) / 24.0;

	}

	return Ebend;

};



///
/// An overloaded function to evaluate the bending energy and its gradient
/// For use with gradient-based methods ( FIRE, L-BFGS, ... )
///
double BendOperator::operator()( Polyhedron &P, VectorXd &grad ) {

	// NOTE: CURRENT GEOMETRY OF POLYHEDRON AND HINGE DERIVATIVES SHOULD BE UP TO DATE
	double Ebend = 0.0; // The global bending energy

	Facet_iterator f;
	for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

		double tarA = f->tarFaceArea();
		Vector3d hBar = f->hBar();
		Matrix3d xi = f->xi();

		// Re-format face quantities for energy calculation
		Vector3d Phi;
		Halfedge_handle he = f->halfedge();
		for( int k = 0; k < 3; k++ ) {

			Phi(k) = ( he->Phi() ) - ( he->tarPhi() );
			he = he->next();

		}

		std::vector<FlapGrad> gradPhi = FlapMap::hingeGradients( f );

		// Construct the traces and trace gradients
		double trB = 0.0;
		double trB2 = 0.0;
		FlapGrad gradTrB = FlapGrad::Zero();
		FlapGrad gradTrB2 = FlapGrad::Zero();

		for( int i = 0; i < 3; i++ ) {

			trB += Phi(i) / hBar(i);
			trB2 += xi(i,i) * Phi(i) * Phi(i);

			gradTrB += gradPhi[i] / hBar(i);
			gradTrB2 += 2.0 * xi(i,i) * Phi(i) * gradPhi[i];

			for( int j = (i+1); j < 3; j++ ) {

				trB2 += xi(i,j) * Phi(i) * Phi(j);

				gradTrB2 += xi(i,j)*( Phi(i)*gradPhi[j] + Phi(j)*gradPhi[i] );

			}

		}

		// Single stencil contribution to the energy
		Ebend += m_h * m_h * tarA * ( m_nu * trB * trB + (1-m_nu) * trB2 ) / 24.0;

		// Singe stencil contribution to the energy gradient
		FlapGrad gradEB = m_h * m_h * tarA * ( 2.0 * m_nu * trB * gradTrB
				+ (1-m_nu) * gradTrB2 ) / 24.0;
		this->mapLocalToGlobal( f, gradEB, grad );

	}

	return Ebend;

};

///
/// An overloaded function to evaluate the bending energy, its gradients, and
/// its Hessian matrix.  For use with a fully implemented Newton method.
///
double BendOperator::operator()( Polyhedron &P, VectorXd &grad, SparseMatrix &hess ) {

	// NOTE: CURRENT GEOMETRY OF POLYHEDRON AND HINGE DERIVATIVES SHOULD BE UP TO DATE
	int Nv = P.size_of_vertices();
	if( hess.size() != ( 9 * Nv * Nv ) ) {
		std::runtime_error( "Hessian matrix is improperly sized!" );
	}

	double Ebend = 0.0; // The global bending energy

	// The list of triplets that will be used to construct the
	// contributions to the global Hessian matrix of the bending energy
	std::vector<Triplet> tripletList;
	// ///////////////////////////////////////////////////////////////////////////////////
	// TO DO: FIGURE OUT HOW TO RESERVE THE PROPER AMOUNT OF MEMORY FOR THIS VECTOR!!!!!!!
	// ///////////////////////////////////////////////////////////////////////////////////

	Facet_iterator f;
	for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

		double tarA = f->tarFaceArea();
		Vector3d hBar = f->hBar();
		Matrix3d xi = f->xi();

		// Re-format face quantities for energy calculation
		Vector3d Phi;
		Halfedge_handle he = f->halfedge();
		for( int k = 0; k < 3; k++ ) {

			Phi(k) = ( he->Phi() ) - ( he->tarPhi() );
			he = he->next();

		}

		std::vector<FlapGrad> gradPhi = FlapMap::hingeGradients( f );
		std::vector<FlapHess> hessPhi = FlapMap::hingeHessians( f );

		// Construct the traces, trace gradients, and trace Hessians
		double trB = 0.0;
		double trB2 = 0.0;

		FlapGrad gradTrB = FlapGrad::Zero();
		FlapGrad gradTrB2 = FlapGrad::Zero();

		FlapHess hessTrB = FlapHess::Zero();
		FlapHess hessTrB2 = FlapHess::Zero();

		for( int i = 0; i < 3; i++ ) {

			trB += Phi(i) / hBar(i);
			trB2 += xi(i,i) * Phi(i) * Phi(i);

			gradTrB += gradPhi[i] / hBar(i);
			gradTrB2 += 2.0 * xi(i,i) * Phi(i) * gradPhi[i];

			hessTrB += hessPhi[i] / hBar(i);
			hessTrB2 += 2.0 * xi(i,i) * (
					gradPhi[i].transpose() * gradPhi[i] +
					Phi(i) * hessPhi[i] );

			for( int j = (i+1); j < 3; j++ ) {

				trB2 += xi(i,j) * Phi(i) * Phi(j);

				gradTrB2 += xi(i,j)*( Phi(i)*gradPhi[j] + Phi(j)*gradPhi[i] );

				hessTrB2 += xi(i,j) * (
						Phi(i) * hessPhi[i] +
						Phi(j) * hessPhi[j] +
						gradPhi[i].transpose() * gradPhi[j] +
						gradPhi[j].transpose() * gradPhi[i] );
			}

		}

		// Single stencil contribution to the energy
		Ebend += m_h * m_h * tarA * ( m_nu * trB * trB + (1-m_nu) * trB2 ) / 24.0;

		// Single stencil contribution to the energy gradient
		FlapGrad gradEB = m_h * m_h * tarA * ( 2.0 * m_nu * trB * gradTrB
				+ (1-m_nu) * gradTrB2 ) / 24.0;
		this->mapLocalToGlobal( f, gradEB, grad );

		// Single stencil contribution to the energy Hessian
		FlapHess hessEB = m_h * m_h * tarA * (
				2.0 * m_nu * ( gradTrB.transpose() * gradTrB + trB * hessTrB ) +
				(1.0 - m_nu) * hessTrB2 ) / 24.0;
		this->mapLocalToGlobal( f, Nv, hessEB, tripletList );

	}

	// Compile the triplet list into the global Hessian matrix
	SparseMatrix GHess( 3*Nv, 3*Nv );
	GHess.setFromTriplets( tripletList.begin(), tripletList.end() );
	hess += GHess;

	return Ebend;

};

#endif
