/*!
 * 	\file HingeOperator.h
 * 	\brief A function class to calculate the hinge function and its derivatives
 *
 * 	\author Dillon Cislo
 * 	\date 01/10/2019
 *
 */

#ifndef _HINGE_OPERATOR_H_
#define _HINGE_OPERATOR_H_

#include <Eigen/Core>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include "../ElasticMesh/ElasticItems.h"

///
/// A function class to calculate the hinge function and its derivatives
///
class HingeOperator {

	private:

		typedef typename CGAL::Simple_cartesian<double> 		Kernel;
		typedef typename CGAL::Polyhedron_3<Kernel, ElasticItems> 	Polyhedron;

		typedef typename Polyhedron::Edge_iterator 	Edge_iterator;

		typedef typename Polyhedron::Halfedge_handle 	Halfedge_handle;
		typedef typename Polyhedron::Face_handle 	Face_handle;
		typedef typename Polyhedron::Vertex_handle 	Vertex_handle;

		typedef typename Eigen::RowVector3d 		RowVector3d;
		typedef typename Eigen::Vector3d 		Vector3d;
		typedef typename Eigen::VectorXd 		VectorXd;
		typedef typename Eigen::Matrix3d 		Matrix3d;
		typedef typename Eigen::MatrixXd 		MatrixXd;
		typedef typename Eigen::Matrix<double,1,12> 	HingeGrad;
		typedef typename Eigen::Matrix<double,12,12> 	HingeHess;

	protected:

		///
		/// Returns the symmetrized form of an Eigen-style
		/// (3x3)-matrix.  Mainly for internal use
		///
		inline Matrix3d sym3( const Matrix3d &A ) {
			return ( A + A.transpose() );
		};

		///
		/// Convenience function for the M-Blocks of the Hessian matrix construction
		/// Mainly for internal use
		///
		inline Matrix3d Mijk( double cosAlphaI, double hi, double hj,
				const Vector3d &mj, const Vector3d &nk ) {

			return ( cosAlphaI * mj * nk.transpose() / ( hi * hj ) );

		};

		///
		/// Convenience function for the N-Blocks of the Hessian matrix construction
		/// Mainly for internal use
		///
		inline Matrix3d Nij( double h0i, double hj,
				const Vector3d &ni, const Vector3d &mj ) {

			return ( ni * mj.transpose() / ( h0i * hj ) );

		};

		///
		/// Convenience function for the B-Blocks of the Hessian matrix construction
		/// Mainly for internal use
		///
		inline Matrix3d Bi( double l0, const Vector3d &ni, const Vector3d &m0i ) {

			return ( ni * m0i.transpose() / ( l0 * l0 ) );

		};

	public:

		///
		/// Evaluate the hinge function and its gradient
		/// For use with gradient based methods ( FIRE, L-BFGS, ... )
		///
		void constructGradient( Polyhedron &P );

		///
		/// Evaluate the hinge function, its gradient, and its Hessian
		/// matrix.  For use with a fully implemented Newton method
		///
		void constructGradientAndHessian( Polyhedron &P );


};

/* BASIC METHOD
///
/// Evaluate the hinge function and its gradient
/// For use with gradient based methods ( FIRE, L-BFGS, ... )
///
void HingeOperator::constructGradient( Polyhedron &P ) {

	Edge_iterator e;
	for( e = P.edges_begin(); e != P.edges_end(); e++ ) {

		// Skip incomplete boundary hinges
		if ( e->is_border() || e->opposite()->is_border() ) continue;

		// /////////////////////////////////////////////////////////////////////////////
		// Construct geometric hinge quantities
		// /////////////////////////////////////////////////////////////////////////////
		
		Vector3d n1 = e->face()->faceNormal();
		Vector3d n2 = e->opposite()->face()->faceNormal();

		double A1 = e->face()->faceArea();
		double A2 = e->opposite()->face()->faceArea();

		double l0 = e->length();

		Vector3d e0Hat = e->edgeUnitVector();
		Vector3d e1 = e->prev()->opposite()->edgeVector();
		Vector3d e2 = e->opposite()->next()->edgeVector();
		Vector3d e3 = e->next()->edgeVector();
		Vector3d e4 = e->opposite()->prev()->opposite()->edgeVector();

		// /////////////////////////////////////////////////////////////////////////////
		// Gradient construction
		// /////////////////////////////////////////////////////////////////////////////

		// Construct bend angle gradient components
		RowVector3d gradTheta0;
		gradTheta0 = -( e0Hat.dot(e3) * n1.transpose() ) / ( 2.0 * A1 );
		gradTheta0 -= ( e0Hat.dot(e4) * n2.transpose() ) / ( 2.0 * A2 );

		RowVector3d gradTheta1;
		gradTheta1 =  ( e0Hat.dot(e1) * n1.transpose() ) / ( 2.0 * A1 );
		gradTheta1 += ( e0Hat.dot(e2) * n2.transpose() ) / ( 2.0 * A2 );

		RowVector3d gradTheta2 = -( l0 * n1.transpose() ) / ( 2.0 * A1 );

		RowVector3d gradTheta3 = -( l0 * n2.transpose() ) / ( 2.0 * A2 );

		// Construct hinge function gradient
		HingeGrad gradTheta;
		gradTheta << gradTheta0, gradTheta1, gradTheta2, gradTheta3;

		HingeGrad gradPhi = 2.0 * gradTheta / ( 1.0 + n1.transpose() * n2 );

		e->setGradPhi( gradPhi );
		e->opposite()->setGradPhi( gradPhi );

	}

};
*/

///
/// Evaluate the hinge function and its gradient
/// For use with gradient based methods ( FIRE, L-BFGS, ... )
///
void HingeOperator::constructGradient( Polyhedron &P ) {

	Edge_iterator e0;
	for( e0 = P.edges_begin(); e0 != P.edges_end(); e0++ ) {

		Halfedge_handle e00 = e0->opposite();

		// Skip incomplete boundary hinges
		if ( e0->is_border() || e00->is_border() ) continue;

		// /////////////////////////////////////////////////////////////////////////////
		// Construct geometric hinge quantities
		// /////////////////////////////////////////////////////////////////////////////
		
		Vector3d n1 = e0->face()->faceNormal();
		Vector3d n2 = e00->face()->faceNormal();

		double A1 = e0->face()->faceArea();
		double A2 = e00->face()->faceArea();

		Halfedge_handle e1 = e0->prev()->opposite();
		Halfedge_handle e2 = e00->next();
		Halfedge_handle e3 = e0->next();
		Halfedge_handle e4 = e00->prev()->opposite();

		double l0 = e0->length();
		double l1 = e1->length();
		double l2 = e2->length();
		double l3 = e3->length();
		double l4 = e4->length();

		double h01 = 2.0 * A1 / l0;
		double h02 = 2.0 * A2 / l0;
		double h1 = 2.0 * A1 / l1;
		double h2 = 2.0 * A2 / l2;
		double h3 = 2.0 * A1 / l3;
		double h4 = 2.0 * A2 / l4;

		double cosAlpha1 = e0->edgeUnitVector().dot( e1->edgeUnitVector() );
		double cosAlpha2 = e0->edgeUnitVector().dot( e2->edgeUnitVector() );
		double cosAlpha3 = e00->edgeUnitVector().dot( e3->edgeUnitVector() );
		double cosAlpha4 = e00->edgeUnitVector().dot( e4->edgeUnitVector() );

		// /////////////////////////////////////////////////////////////////////////////
		// Gradient construction
		// /////////////////////////////////////////////////////////////////////////////
		
		// Construct bend angle gradient components
		RowVector3d gradTheta0 = cosAlpha3 * n1.transpose() / h3;
		gradTheta0 += cosAlpha4 * n2.transpose() / h4;

		RowVector3d gradTheta1 = cosAlpha1 * n1.transpose() / h1;
		gradTheta1 += cosAlpha2 * n2.transpose() / h2;

		RowVector3d gradTheta2 = -n1.transpose() / h01;

		RowVector3d gradTheta3 = -n2.transpose() / h02;

		// Assemble bend angle gradient
		HingeGrad gradTheta;
		gradTheta << gradTheta0, gradTheta1, gradTheta2, gradTheta3;

		// Assemble hinge function gradient
		double cosTheta = n1.transpose() * n2;
		HingeGrad gradPhi = 2.0 * gradTheta / ( 1.0 + cosTheta );

		e0->setGradPhi( gradPhi );
		e00->setGradPhi( gradPhi );

	}

};

/* MY ORIGINAL METHOD
///
/// Evaluate the hinge function and its gradient
/// For use with gradient based methods ( FIRE, L-BFGS, ... )
///
void HingeOperator::constructGradient( Polyhedron &P ) {

	Edge_iterator e0;
	for( e0 = P.edges_begin(); e0 != P.edges_end(); e0++ ) {

		Halfedge_handle e00 = e0->opposite();

		// Skip incomplete boundary hinges
		if ( e0->is_border() || e00->is_border() ) continue;

		// /////////////////////////////////////////////////////////////////////////////
		// Construct geometric hinge quantities
		// /////////////////////////////////////////////////////////////////////////////
		
		// Face normals
		Vector3d n1 = e0->face()->faceNormal();
		Vector3d n2 = e00->face()->faceNormal();

		// Face areas
		double A1 = e0->face()->faceArea();
		double A2 = e00->face()->faceArea();

		// Hinge edge handles
		Halfedge_handle e1 = e0->prev()->opposite();
		Halfedge_handle e2 = e00->next();
		Halfedge_handle e3 = e0->next();
		Halfedge_handle e4 = e00->prev()->opposite();

		// Hinge edge lengths
		double l0 = e0->length();
		double l1 = e1->length();
		double l2 = e2->length();
		double l3 = e3->length();
		double l4 = e4->length();

		// Hinge edge altitudes
		double h01 = 2.0 * A1 / l0;
		double h02 = 2.0 * A2 / l0;
		double h1 = 2.0 * A1 / l1;
		double h2 = 2.0 * A2 / l2;
		double h3 = 2.0 * A1 / l3;
		double h4 = 2.0 * A2 / l4;

		// Re-scaled in-plane edge unit normals
		Vector3d mh01 = e0->edgeNormal() / h01;
		Vector3d mh02 = e00->edgeNormal() / h02;
		Vector3d mh1 = e1->opposite()->edgeNormal() / h1;
		Vector3d mh2 = e2->edgeNormal() / h2;
		Vector3d mh3 = e3->edgeNormal() / h3;
		Vector3d mh4 = e4->opposite()->edgeNormal() / h4;
		
		// Hinge angle constants
		double cCosTheta = 1.0 + n1.dot(n2); // ( 1 + cosTheta )
		double sinTheta = n1.cross(n2).dot(e0->edgeUnitVector());

		// /////////////////////////////////////////////////////////////////////////////
		// Gradient construction
		// /////////////////////////////////////////////////////////////////////////////
		
		// Form inner product matrices
		Matrix3d n1n1 = n1 * n1.transpose();
		Matrix3d n1n2 = n1 * n2.transpose();
		Matrix3d n2n1 = n2 * n1.transpose();
		Matrix3d n2n2 = n2 * n2.transpose();

		// Form the gradient of cosTheta
		HingeGrad cosGrad;
		RowVector3d cosGrad0, cosGrad1, cosGrad2, cosGrad3;

		cosGrad0 << mh3.transpose() * n2n1 + mh4.transpose() * n1n2;
		cosGrad1 << mh1.transpose() * n2n1 + mh2.transpose() * n1n2;
		cosGrad2 << mh01.transpose() * n2n1;
		cosGrad3 << mh02.transpose() * n1n2;

		cosGrad << cosGrad0, cosGrad1, cosGrad2, cosGrad3;

		// Form the gradient of sinTheta
		HingeGrad sinGrad;
		RowVector3d sinGrad0, sinGrad1, sinGrad2, sinGrad3;

		sinGrad0 << mh3.transpose() * n1n1 + mh4.transpose() * n2n2;
		sinGrad1 << mh1.transpose() * n1n1 + mh2.transpose() * n2n2;
		sinGrad2 << mh01.transpose() * n1n1;
		sinGrad3 << mh02.transpose() * n2n2;

		sinGrad << sinGrad0, sinGrad1, sinGrad2, sinGrad3;

		// Form the gradient of the hinge function
		HingeGrad gradPhi;

		gradPhi << 2.0 * ( cCosTheta * sinGrad - sinTheta * cosGrad );
		gradPhi = gradPhi / ( cCosTheta * cCosTheta );

		e0->setGradPhi( gradPhi );
		e00->opposite()->setGradPhi( gradPhi );

	}

};
*/

///
/// Evaluate the hinge function, its gradient, and its Hessian matrix.
/// For use with a fully implemented Newton method
///
void HingeOperator::constructGradientAndHessian( Polyhedron &P ) {

	Edge_iterator e0;
	for( e0 = P.edges_begin(); e0 != P.edges_end(); e0++ ) {

		Halfedge_handle e00 = e0->opposite();

		// Skip incomplete boundary hinges
		if ( e0->is_border() || e00->is_border() ) continue;

		// /////////////////////////////////////////////////////////////////////////////
		// Construct geometric hinge quantities
		// /////////////////////////////////////////////////////////////////////////////
		
		Vector3d n1 = e0->face()->faceNormal();
		Vector3d n2 = e00->face()->faceNormal();

		double A1 = e0->face()->faceArea();
		double A2 = e00->face()->faceArea();

		Halfedge_handle e1 = e0->prev()->opposite();
		Halfedge_handle e2 = e00->next();
		Halfedge_handle e3 = e0->next();
		Halfedge_handle e4 = e00->prev()->opposite();

		double l0 = e0->length();
		double l1 = e1->length();
		double l2 = e2->length();
		double l3 = e3->length();
		double l4 = e4->length();

		double h01 = 2.0 * A1 / l0;
		double h02 = 2.0 * A2 / l0;
		double h1 = 2.0 * A1 / l1;
		double h2 = 2.0 * A2 / l2;
		double h3 = 2.0 * A1 / l3;
		double h4 = 2.0 * A2 / l4;

		Vector3d m01 = e0->edgeNormal();
		Vector3d m02 = e00->edgeNormal();
		Vector3d m1 = e1->opposite()->edgeNormal();
		Vector3d m2 = e2->edgeNormal();
		Vector3d m3 = e3->edgeNormal();
		Vector3d m4 = e4->opposite()->edgeNormal();

		double cosAlpha1 = e0->edgeUnitVector().dot( e1->edgeUnitVector() );
		double cosAlpha2 = e0->edgeUnitVector().dot( e2->edgeUnitVector() );
		double cosAlpha3 = e00->edgeUnitVector().dot( e3->edgeUnitVector() );
		double cosAlpha4 = e00->edgeUnitVector().dot( e4->edgeUnitVector() );

		// /////////////////////////////////////////////////////////////////////////////
		// Gradient construction
		// /////////////////////////////////////////////////////////////////////////////
		
		// Construct bend angle gradient components
		RowVector3d gradTheta0 = cosAlpha3 * n1.transpose() / h3;
		gradTheta0 += cosAlpha4 * n2.transpose() / h4;

		RowVector3d gradTheta1 = cosAlpha1 * n1.transpose() / h1;
		gradTheta1 += cosAlpha2 * n2.transpose() / h2;

		RowVector3d gradTheta2 = -n1.transpose() / h01;

		RowVector3d gradTheta3 = -n2.transpose() / h02;

		// Assemble bend angle gradient
		HingeGrad gradTheta;
		gradTheta << gradTheta0, gradTheta1, gradTheta2, gradTheta3;

		// Assemble hinge function gradient
		double cosTheta = n1.transpose() * n2;
		HingeGrad gradPhi = 2.0 * gradTheta / ( 1.0 + cosTheta );

		e0->setGradPhi( gradPhi );
		e00->setGradPhi( gradPhi );

		// /////////////////////////////////////////////////////////////////////////////
		// Hessian construction
		// /////////////////////////////////////////////////////////////////////////////
		
		// Construct basic bend angle Hessian sub-block constituents
		Matrix3d M331 = Mijk( cosAlpha3, h3, h3, m3, n1 );
		Matrix3d M442 = Mijk( cosAlpha4, h4, h4, m4, n2 );
		Matrix3d M311 = Mijk( cosAlpha3, h3, h1, m1, n1 );
		Matrix3d M131 = Mijk( cosAlpha1, h1, h3, m3, n1 );
		Matrix3d M422 = Mijk( cosAlpha4, h4, h2, m2, n2 );
		Matrix3d M242 = Mijk( cosAlpha2, h2, h4, m4, n2 );
		Matrix3d M111 = Mijk( cosAlpha1, h1, h1, m1, n1 );
		Matrix3d M222 = Mijk( cosAlpha2, h2, h2, m2, n2 );
		
		Matrix3d M3011 = Mijk( cosAlpha3, h3, h01, m01, n1 );
		Matrix3d M4022 = Mijk( cosAlpha4, h4, h02, m02, n2 );
		Matrix3d M1011 = Mijk( cosAlpha1, h1, h01, m01, n1 );
		Matrix3d M2022 = Mijk( cosAlpha2, h2, h02, m02, n2 );

		Matrix3d N13 = Nij( h01, h3, n1, m3 );
		Matrix3d N24 = Nij( h02, h4, n2, m4 );
		Matrix3d N11 = Nij( h01, h1, n1, m1 );
		Matrix3d N22 = Nij( h02, h2, n2, m2 );

		Matrix3d N101 = Nij( h01, h01, n1, m01 );
		Matrix3d N202 = Nij( h02, h02, n2, m02 );

		Matrix3d B1 = Bi( l0, n1, m01 );
		Matrix3d B2 = Bi( l0, n2, m02 );

		// Construct bend angle Hessian sub-blocks
		Matrix3d H00 = sym3(M331) - B1 + sym3(M442) - B2;

		Matrix3d H01 = M311 + M131.transpose() + B1 + M422 + M242.transpose() + B2;
		Matrix3d H10 = H01.transpose();

		Matrix3d H02 = M3011 - N13;
		Matrix3d H20 = H02.transpose();

		Matrix3d H03 = M4022 - N24;
		Matrix3d H30 = H03.transpose();

		Matrix3d H11 = sym3(M111) - B1 + sym3(M222) - B2;

		Matrix3d H12 = M1011 - N11;
		Matrix3d H21 = H12.transpose();

		Matrix3d H13 = M2022 - N22;
		Matrix3d H31 = H13.transpose();

		Matrix3d H22 = -sym3(N101);
		
		Matrix3d H23 = Matrix3d::Zero();
		Matrix3d H32 = Matrix3d::Zero();

		Matrix3d H33 = -sym3(N202);

		// Assemble bend angle Hessian matrix
		HingeHess hessTheta;
		hessTheta << H00, H01, H02, H03,
			     H10, H11, H12, H13,
			     H20, H21, H22, H23,
			     H30, H31, H32, H33;

		// Assemble hinge function Hessian matrix
		double sinTheta = n1.cross(n2).dot( e0->edgeUnitVector() );
		
		HingeHess hessPhi;
		hessPhi = 2.0 * sinTheta * gradTheta.transpose() * gradTheta;
		hessPhi = hessPhi / ( ( 1.0 + cosTheta ) * ( 1.0 + cosTheta ) );
		hessPhi += 2.0 * hessTheta / ( 1.0 + cosTheta );

		e0->setHessianPhi( hessPhi );
		e00->setHessianPhi( hessPhi );

	}

};

#endif
