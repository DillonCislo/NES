/*!
 * 	\file StretchOperator.h
 * 	\brief A class to calculate the stretching energy and its derivatives
 *
 * 	\author Dillon Cislo
 * 	\date 01/08/2019
 *
 */

#ifndef _STRETCH_OPERATOR_H_
#define _STRETCH_OPERATOR_H_

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include "../ElasticMesh/ElasticItems.h"

///
/// A class to calculate the stretching energy and its derivatives
///
class StretchOperator {

	private:

		typedef typename CGAL::Simple_cartesian<double> 		Kernel;
		typedef typename CGAL::Polyhedron_3<Kernel, ElasticItems>	Polyhedron;

		typedef typename Polyhedron::Vertex_iterator 		Vertex_iterator;
		typedef typename Polyhedron::Facet_iterator 		Facet_iterator;
		typedef typename Polyhedron::Edge_iterator 		Edge_iterator;
		typedef typename Polyhedron::Halfedge_iterator 		Halfedge_iterator;

		typedef typename Polyhedron::Halfedge_handle 		Halfedge_handle;
		typedef typename Polyhedron::Face_handle 		Face_handle;
		typedef typename Polyhedron::Vertex_handle 		Vertex_handle;

		typedef typename Polyhedron::Halfedge_around_vertex_circulator 	HV_circulator;
		typedef typename Polyhedron::Halfedge_around_facet_circulator 	HF_circulator;

		typedef typename Eigen::RowVector3d 			RowVector3d;
		typedef typename Eigen::Vector3d 			Vector3d;
		typedef typename Eigen::VectorXd 			VectorXd;
		typedef typename Eigen::VectorXi 			VectorXi;
		typedef typename Eigen::Matrix3d 			Matrix3d;
		typedef typename Eigen::MatrixXd 			MatrixXd;
		typedef typename Eigen::Matrix<double, 1, 9> 		RowGradS;
		typedef typename Eigen::Matrix<double, 9, 9> 		Matrix9d;

		typedef typename Eigen::SparseMatrix<double> 		SparseMatrix;
		typedef typename Eigen::Triplet<double> 		Triplet;

	protected:

		///
		/// A vector containing the constant edge vector gradient matrices
		///
		std::vector< Eigen::Matrix<double, 3, 9> > m_gradE;

		///
		/// A vector containing the constant Hessian matrices of the squared
		/// edge lengths
		///
		std::vector< Matrix9d > m_hessL2;

		///
		/// Poisson's ratio
		///
		double m_nu;

	public:
		///
		/// Null constructor
		///
		StretchOperator() {};

		///
		/// Default constructor
		///
		StretchOperator( double nu );

		///
		/// An overloaded function to calculate the stretching energy.
		/// For use when it is only necessary to calculate the energy and not
		/// any of its derivatives
		///
		double operator()( Polyhedron &P );

		///
		/// An overloaded function to evaluate the stretching energy and its gradient
		/// For use with gradient-based methods ( FIRE, L-BFGS, ... )
		///
		double operator()( Polyhedron &P, VectorXd &grad );

		///
		/// An overloaded function to evaluate the stretching energy, its gradient,
		/// and its Hessian matrix.  For use with a fully implemented Newton method.
		///
		double operator()( Polyhedron &P, VectorXd &grad, SparseMatrix &hess );

	protected:

		///
		/// An overloaded function that maps local gradient quantities to the global
		/// gradient vector
		///
		void mapLocalToGlobal( Face_handle f, 
				const RowGradS &LGrad, VectorXd &GGrad );

		///
		/// An overloaded function that maps local Hessian quantities to the vector
		/// of Eigen-style triplets which will be used to construct the global Hessian
		/// matrix of the stretching energy
		///
		void mapLocalToGlobal( Face_handle f, int Nv,
				const Matrix9d &LHess, std::vector<Triplet> &GTrip );

};

///
/// Default constructor
///
StretchOperator::StretchOperator( double nu ) : m_nu( nu ) {


	Matrix3d Z3 = Matrix3d::Zero();
	Matrix3d I3 = Matrix3d::Identity();

	Eigen::Matrix<double, 3, 9> gradEI;
	gradEI << Z3, -I3, I3;

	Eigen::Matrix<double, 3, 9> gradEJ;
	gradEJ << I3, Z3, -I3;

	Eigen::Matrix<double, 3, 9> gradEK;
	gradEK << -I3, I3, Z3;

	std::vector< Eigen::Matrix<double, 3, 9> > gradE{ gradEI, gradEJ, gradEK };
	this->m_gradE = gradE;

	MatrixXd Z9 = Matrix9d::Zero();

	Matrix9d hessL2I;
	hessL2I << Z9, -gradEI, gradEI;

	Matrix9d hessL2J;
	hessL2J << gradEJ, Z9, -gradEJ;

	Matrix9d hessL2K;
	hessL2K << -gradEK, gradEK, Z9;

	std::vector< Matrix9d > hessL2{ hessL2I, hessL2J, hessL2K };
	this->m_hessL2 = hessL2;

};


///
/// An overloaded function that maps local gradient quantities to the global gradient vector
///
void StretchOperator::mapLocalToGlobal( Face_handle f,
	       	const RowGradS &LGrad, VectorXd &GGrad ) {

	// NOTE: VERTEX IDs SHOULD ALREADY BE UP TO DATE
	int vi, vj, vk;

	vi = f->halfedge()->next()->vertex()->id();
	vj = f->halfedge()->prev()->vertex()->id();
	vk = f->halfedge()->vertex()->id();

	if ( (GGrad.size() % 3) != 0 ) {
		std::runtime_error( "Gradient vector is improperly sized!" );
	}

	int Nv = GGrad.size() / 3;

	// Vi --------------------------------------------------------------
	GGrad( vi ) += LGrad(0);
	GGrad( vi + Nv ) += LGrad(1);
	GGrad( vi + (2*Nv) ) += LGrad(2);

	// Vj --------------------------------------------------------------
	GGrad( vj ) += LGrad(3);
	GGrad( vj + Nv ) += LGrad(4);
	GGrad( vj + (2*Nv) ) += LGrad(5);

	// Vk --------------------------------------------------------------
	GGrad( vk ) += LGrad(6);
	GGrad( vk + Nv ) += LGrad(7);
	GGrad( vk + (2*Nv) ) += LGrad(8);

};
	
///
/// An overloaded function that maps local Hessian quantities to the vector of
/// Eigen-style triplets, which will be used to construct the global Hessian 
/// matrix of the stretching energy
///
void StretchOperator::mapLocalToGlobal( Face_handle f, int Nv,
		const Matrix9d &LHess, std::vector<Triplet> &GTrip ) {

	// NOTE: VERTEX IDs SHOULD ALREADY BE UP TO DATE
	int vi, vj, vk;

	vi = f->halfedge()->next()->vertex()->id();
	vj = f->halfedge()->prev()->vertex()->id();
	vk = f->halfedge()->vertex()->id();
	
	std::vector<int> vID = { vi, vj, vk };
	for( int m = 0; m < 3; m++ ) {

		int mID = vID[m];
		int M = 3 * m;

		for( int n = 0; n < 3; n++ ) {

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
/// An overloaded function to calculate the stretching energy.
/// For use when it is only necessary to calculate the energy and not any of its derivatives
///
double StretchOperator::operator()( Polyhedron &P ) {

	// NOTE: CURRENT GEOMETRY OF POLYHEDRON SHOULD BE UP TO DATE
	double Estretch = 0.0; // The global stretching energy

	Facet_iterator f;
	for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

		double tarA = f->tarFaceArea();
		Vector3d C = f->C();
		Matrix3d zeta = f->zeta();

		// Re-format face quantities for energy calculation
		Vector3d str;
		Halfedge_handle he = f->halfedge();
		for( int k = 0; k < 3; k++ ) {

			str(k) = he->edgeStrain();
			he = he->next();

		}

		// Construct the traces
		double trE = 0.0;
		double trE2 = 0.0;
		for( int i = 0; i < 3; i++ ) {

			trE += C(i) * str(i);
			trE2 += zeta(i,i) * str(i) * str(i);

			for( int j = (i+1); j < 3; j++ ) {

				trE2 += zeta(i,j) * str(i) * str(j);

			}

		}

		// Single stencil contribution to the energy
		Estretch += tarA * ( m_nu * trE * trE + (1.0-m_nu) * trE2 ) / 8.0;
		
	}

	return Estretch;

};


///
/// An overloaded function to evaluate the stretching energy and its gradient
/// For use with gradient-based methods ( FIRE, L-BFGS, ... )
///
double StretchOperator::operator()( Polyhedron &P, VectorXd &grad ) {

	// NOTE: CURRENT GEOMETRY OF POLYHEDRON SHOULD BE UP TO DATE
	double Estretch = 0.0; // The global stretching energy

	Facet_iterator f;
	for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

		double tarA = f->tarFaceArea();
		Vector3d C = f->C();
		Matrix3d zeta = f->zeta();

		// Re-format face quantities for energy calculation
		Vector3d str;
		Halfedge_handle he = f->halfedge();
		for( int k = 0; k < 3; k++ ) {

			str(k) = he->edgeStrain();
			he = he->next();

		}

		// Assemble edge-length gradients
		RowVector3d ei = 2.0 * f->halfedge()->edgeVector().transpose();
		RowVector3d ej = 2.0 * f->halfedge()->next()->edgeVector().transpose();
		RowVector3d ek = 2.0 * f->halfedge()->prev()->edgeVector().transpose();

		RowGradS gradLI2, gradLJ2, gradLK2;

		gradLI2 << RowVector3d::Zero(), -ei, ei;
		gradLJ2 << ej, RowVector3d::Zero(), -ej;
		gradLK2 << -ek, ek, RowVector3d::Zero();

		std::vector<RowGradS> gradL2{ gradLI2, gradLJ2, gradLK2 };

		// Construct the traces and trace gradients
		double trE = 0.0;
		double trE2 = 0.0;
		RowGradS gradTrE = RowGradS::Zero();
		RowGradS gradTrE2 = RowGradS::Zero();
		for( int i = 0; i < 3; i++ ) {

			trE += C(i) * str(i);
			trE2 += zeta(i,i) * str(i) * str(i);

			gradTrE += C(i) * gradL2[i];
			gradTrE2 += 2.0 * zeta(i,i) * str(i) * gradL2[i];

			for( int j = (i+1); j < 3; j++ ) {

				trE2 += zeta(i,j) * str(i) * str(j);
				
				gradTrE2 += zeta(i,j) * ( str(i)*gradL2[j] + str(j)*gradL2[i] );

			}

		}

		// Single stencil contribution to the energy
		Estretch += tarA * ( m_nu * trE * trE + (1.0-m_nu) * trE2 ) / 8.0;

		// Single stencil constribution to the energy gradient
		RowGradS gradES = tarA * ( 2.0 * m_nu * trE * gradTrE +
				(1.0 - m_nu) * gradTrE2 ) / 8.0;
		this->mapLocalToGlobal( f, gradES, grad );
		
	}

	return Estretch;

};

///
/// An overloaded function to evaulate the stretching energy, its gradient, and
/// its Hessian matrix.  For use with a fully implemented Newton method.
///
double StretchOperator::operator()( Polyhedron &P, VectorXd &grad, SparseMatrix &hess ) {

	// NOTE: CURRENT GEOMETRY OF POLYHEDRON SHOULD BE UP TO DATE
	int Nv = P.size_of_vertices();
	if( hess.size() != ( 9 * Nv * Nv ) ) {
		std::runtime_error( "Hessian matrix is improperly sized!" );
	}

	double Estretch = 0.0; // The global stretching energy

	// The list of triplets that will be used to construct the
	// contributions to the global Hessian matrix of the stretching energy
	std::vector<Triplet> tripletList;
	// /////////////////////////////////////////////////////////////////////////////////
	// TO DO: FIGURE OUT HOW TO RESERVE THE PROPER AMOUNT OF MEMORY FOR THIS VECTOR!!!!!
	// /////////////////////////////////////////////////////////////////////////////////

	Facet_iterator f;
	for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

		double tarA = f->tarFaceArea();
		Vector3d C = f->C();
		Matrix3d zeta = f->zeta();
		Matrix9d hessTrE = f->hessTrE();

		// Re-format face quantities for energy calculation
		Vector3d str;
		Halfedge_handle he = f->halfedge();
		for( int k = 0; k < 3; k++ ) {

			str(k) = he->edgeStrain();
			he = he->next();

		}

		// Assemble edge-length gradients
		RowVector3d ei = 2.0 * f->halfedge()->edgeVector().transpose();
		RowVector3d ej = 2.0 * f->halfedge()->next()->edgeVector().transpose();
		RowVector3d ek = 2.0 * f->halfedge()->prev()->edgeVector().transpose();

		RowGradS gradLI2, gradLJ2, gradLK2;

		gradLI2 << RowVector3d::Zero(), -ei, ei;
		gradLJ2 << ej, RowVector3d::Zero(), -ej;
		gradLK2 << -ek, ek, RowVector3d::Zero();

		std::vector<RowGradS> gradL2{ gradLI2, gradLJ2, gradLK2 };

		// Construct the traces, trace gradients, and trace Hessians
		double trE = 0.0;
		double trE2 = 0.0;
		RowGradS gradTrE = RowGradS::Zero();
		RowGradS gradTrE2 = RowGradS::Zero();
		Matrix9d hessTrE2 = Matrix9d::Zero();
		for( int i = 0; i < 3; i++ ) {

			trE += C(i) * str(i);
			trE2 += zeta(i,i) * str(i) * str(i);

			gradTrE += C(i) * gradL2[i];
			gradTrE2 += 2.0 * zeta(i,i) * str(i) * gradL2[i];

			hessTrE2 += 2.0 * zeta(i,i) * (
				gradL2[i].transpose() * gradL2[i] +
				str(i) * m_hessL2[i] );

			for( int j = (i+1); j < 3; j++ ) {

				trE2 += zeta(i,j) * str(i) * str(j);
				
				gradTrE2 += zeta(i,j) * ( str(i)*gradL2[j] + str(j)*gradL2[i] );

				hessTrE2 += zeta(i,j) * (
					str(i) * m_hessL2[j] +
					str(j) * m_hessL2[i] +
					gradL2[i].transpose() * gradL2[j] +
					gradL2[j].transpose() * gradL2[i] );


			}

		}

		// Single stencil contribution to the energy
		Estretch += tarA * ( m_nu * trE * trE + (1.0-m_nu) * trE2 ) / 8.0;

		// Single stencil constribution to the energy gradient
		RowGradS gradES = tarA * ( 2.0 * m_nu * trE * gradTrE +
				(1.0 - m_nu) * gradTrE2 ) / 8.0;
		this->mapLocalToGlobal( f, gradES, grad );

		// Single stencil contribution to the energy Hessian
		Matrix9d hessES = tarA * (
				2.0 * m_nu * ( gradTrE.transpose() * gradTrE + trE * hessTrE ) +
				(1.0-m_nu) * hessTrE2 ) / 8.0;
		this->mapLocalToGlobal( f, Nv, hessES, tripletList );

		
	}

	// Compile the triplet list into the global Hessian matrix
	SparseMatrix GHess( 3*Nv, 3*Nv );
	GHess.setFromTriplets( tripletList.begin(), tripletList.end() );
	hess += GHess;

	return Estretch;

};


#endif
