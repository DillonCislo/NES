/*!
 *  \file RestrictGrowthOperator.h
 *  \brief A class to calculate the growth restriction energy and its derivatives
 *
 *  \author Dillon Cislo
 *  \date 12/25/2020
 *
 */

#ifndef _RESTRICT_GROWTH_OPERATOR_H_
#define _RESTRICT_GROWTH_OPERATOR_H_


#include <vector>
#include <cassert>
#include <cmath>
#include <limits>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include "../ElasticMesh/ElasticItems.h"

// Branchless, type-safe sign function (used to calculate hinge angle)
template <class T>
inline int sgn( T val ) {
	return ( T(0) < val ) - ( val < T(0) );
};

///
/// A class to calculate the growth restriction energy and its derivatives
///
class RestrictGrowthOperator {

  private:

    typedef typename CGAL::Simple_cartesian<double>             Kernel;
    typedef typename CGAL::Polyhedron_3<Kernel, ElasticItems>   Polyhedron;

    typedef typename Polyhedron::Facet_iterator                 Facet_iterator;
    typedef typename Polyhedron::Halfedge_handle                Halfedge_handle;
    typedef typename Polyhedron::Face_handle                    Face_handle;
    typedef typename Polyhedron::Vertex_handle                  Vertex_handle;

    typedef typename Eigen::RowVector3d                         RowVector3d;
    typedef typename Eigen::Vector3d                            Vector3d;
    typedef typename Eigen::VectorXd                            VectorXd;
    typedef typename Eigen::VectorXi                            VectorXi;
    typedef typename Eigen::Matrix3d                            Matrix3d;
    typedef typename Eigen::MatrixXd                            MatrixXd;
    typedef typename Eigen::Matrix<double, 1, 9>                RowGradS;
    typedef typename Eigen::Matrix<double, 3, 9>                RowGradVec;
    typedef typename Eigen::Matrix<double, 9, 9>                Matrix9d;

    typedef typename Eigen::SparseMatrix<double>                SparseMatrix;
    typedef typename Eigen::Triplet<double>                     Triplet;

  protected:

		///
		/// A vector containing the constant edge vector gradient matrices
		///
		std::vector<RowGradVec, Eigen::aligned_allocator<RowGradVec> > m_gradE;

    ///
    /// The growth restriction coefficient
    ///
    double m_mu;

  public:

    ///
    /// Null constructor
    ///
    RestrictGrowthOperator() {};

    ///
    /// Default constructor
    ///
     RestrictGrowthOperator( double mu );

    ///
    /// An overloaded function calculate the growth restriction energy. For use when it
    /// is only necessary to calculate the energy and not any of its derivatives.
    ///
    double operator()( Polyhedron &P );

    ///
    /// An overloaded function to evaluate the growth restriction energy and its gradient.
    /// For use with gradient-based optimization methods ( FIRE, L-BFGS, ... )
    ///
    double operator()( Polyhedron &P, VectorXd &grad );

    ///
    /// An overloaded function to evaluate the growth restriction energy, its gradient
    /// and its Hessian matrix.  For use with a fully implemented Newton method.
    ///
    double operator()( Polyhedron &P, VectorXd &grad, SparseMatrix &hess );

  protected:

    ///
    /// An overloaded function that maps local gradient quantities to the
    /// global gradient vector
    ///
    void mapLocalToGlobal( Face_handle f, const RowGradS &LGrad, VectorXd &GGrad );

    ///
    /// An overloaded function that maps local Hessian quantities to the
    /// vector of Eigen-style triplets that will be used to construct the global
    /// Hessian matrix of the fixed volume energy
    ///
    void mapLocalToGlobal( Face_handle f, int Nv,
        const Matrix9d &LHess, std::vector<Triplet> &GTrip );

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

};

///
/// Default constructor
///
RestrictGrowthOperator::RestrictGrowthOperator( double mu ) : m_mu( mu ) {

	Matrix3d Z3 = Matrix3d::Zero();
	Matrix3d I3 = Matrix3d::Identity();

	RowGradVec gradEI;
	gradEI << Z3, -I3, I3;

	RowGradVec gradEJ;
	gradEJ << I3, Z3, -I3;

	RowGradVec gradEK;
	gradEK << -I3, I3, Z3;

	std::vector<RowGradVec,
    Eigen::aligned_allocator<RowGradVec> >
      gradE{ gradEI, gradEJ, gradEK };

	this->m_gradE = gradE;

};

///
/// An overloaded function that maps local gradient quantities to the
/// global gradient vector
///
void RestrictGrowthOperator::mapLocalToGlobal( Face_handle f,
    const RowGradS &LGrad, VectorXd &GGrad ) {

	// NOTE: VERTEX IDs SHOULD ALREADY BE UP TO DATE
	int vi, vj, vk;

	vi = f->halfedge()->next()->vertex()->id();
	vj = f->halfedge()->prev()->vertex()->id();
	vk = f->halfedge()->vertex()->id();

  /*
	if ( (GGrad.size() % 3) != 0 ) {
		std::runtime_error( "Gradient vector is improperly sized!" );
	}
  */

  assert( ((GGrad.size() % 3) == 0) &&
      "Gradient vector is improperly sized!" );

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
/// matrix of the fixed volume energy
///
void RestrictGrowthOperator::mapLocalToGlobal( Face_handle f, int Nv,
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
/// An overloaded function calculate the growth restriction energy. For use when it
/// is only necessary to calculate the energy and not any of its derivatives.
///
double RestrictGrowthOperator::operator()( Polyhedron &P ) {

  // NOTE: CURRENT GEOMETRY OF POLYHEDRON SHOULDBE UP TO DATE
  double Erg = 0.0; // The global growth restriction energy

  bool constraintViolation = false;

  Facet_iterator f;
  for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

    double gi = f->halfedge()->initEdgeGrowth();
    double gj = f->halfedge()->next()->initEdgeGrowth();
    double gk = f->halfedge()->prev()->initEdgeGrowth();

    Vector3d ei = f->halfedge()->edgeVector();
    Vector3d ej = f->halfedge()->next()->edgeVector();
    Vector3d ek = f->halfedge()->prev()->edgeVector();

    double Li = f->halfedge()->restrictedLength();
    double Lj = f->halfedge()->next()->restrictedLength();
    double Lk = f->halfedge()->prev()->restrictedLength();

    // Map v0 into the deformed face
    Vector3d vHat = gj * ek - gk * ej;
    vHat = vHat / vHat.norm();

    double ci = Li - std::abs(ei.dot(vHat));
    double cj = Lj - std::abs(ej.dot(vHat));
    double ck = Lk - std::abs(ek.dot(vHat));

    if ( ci > 0.0 ) {
      Erg += std::log( ci );
    } else {
      constraintViolation = true;
      break;
    }

    if ( cj > 0.0 ) {
      Erg += std::log( cj );
    } else {
      constraintViolation = true;
      break;
    }

    if ( ck > 0.0 ) {
      Erg += std::log( ck );
    } else {
      constraintViolation = true;
      break;
    }

  }

  if ( constraintViolation ) {
    Erg = std::numeric_limits<double>::infinity();
  } else {
    Erg = -Erg / this->m_mu;
  }

  return Erg;

};

///
/// An overloaded function to evaluate the growth restriction energy and its gradient.
/// For use with gradient-based optimization methods ( FIRE, L-BFGS, ... )
///
double RestrictGrowthOperator::operator()( Polyhedron &P, VectorXd &grad ) {

  // NOTE: CURRENT GEOMETRY OF POLYHEDRON SHOULDBE UP TO DATE
  double Erg = 0.0; // The global growth restriction energy

  bool constraintViolation = false;

  Facet_iterator f;
  for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

    double gi = f->halfedge()->initEdgeGrowth();
    double gj = f->halfedge()->next()->initEdgeGrowth();
    double gk = f->halfedge()->prev()->initEdgeGrowth();

    Vector3d ei = f->halfedge()->edgeVector();
    Vector3d ej = f->halfedge()->next()->edgeVector();
    Vector3d ek = f->halfedge()->prev()->edgeVector();

    double Li = f->halfedge()->restrictedLength();
    double Lj = f->halfedge()->next()->restrictedLength();
    double Lk = f->halfedge()->prev()->restrictedLength();

    RowGradVec gradEI = this->m_gradE[0];
    RowGradVec gradEJ = this->m_gradE[1];
    RowGradVec gradEK = this->m_gradE[2];

    // Map v0 into the deformed face
    Vector3d v = gj * ek - gk * ej;
    double normv = v.norm();
    Vector3d vHat = v / normv;

    double dotiv = ei.dot(vHat);
    double dotjv = ej.dot(vHat);
    double dotkv = ek.dot(vHat);

    double ci = Li - std::abs(dotiv);
    double cj = Lj - std::abs(dotjv);
    double ck = Lk - std::abs(dotkv);

    // Single stencil contribution to the energy
    if ( ci > 0.0 ) {
      Erg += std::log( ci );
    } else {
      constraintViolation = true;
      break;
    }

    if ( cj > 0.0 ) {
      Erg += std::log( cj );
    } else {
      constraintViolation = true;
      break;
    }

    if ( ck > 0.0 ) {
      Erg += std::log( ck );
    } else {
      constraintViolation = true;
      break;
    }

    Matrix3d I3 = Matrix3d::Identity();
    RowGradVec gradVHat = gj * gradEK - gk * gradEJ;
    gradVHat = ((I3 - vHat * vHat.transpose()) * gradVHat) / normv;

    RowGradS gradEGI, gradEGJ, gradEGK, gradE;

    gradEGI << v.transpose() * gradEI + ei.transpose() * gradVHat;
    gradEGI = sgn( dotiv ) * gradEGI / ( this->m_mu * std::max(ci, 1.0e-10) );

    gradEGJ << v.transpose() * gradEJ + ej.transpose() * gradVHat;
    gradEGJ = sgn( dotjv ) * gradEGJ / ( this->m_mu * std::max(cj, 1.0e-10) );

    gradEGK << v.transpose() * gradEK + ek.transpose() * gradVHat;
    gradEGK = sgn( dotkv ) * gradEGK / ( this->m_mu * std::max(ck, 1.0e-10) );

    // Single stencil contribution to the energy gradient
    gradE = gradEGI + gradEGJ + gradEGK;
    this->mapLocalToGlobal( f, gradE, grad );

  }

  if ( constraintViolation ) {
    Erg = std::numeric_limits<double>::infinity();
  } else {
    Erg = -Erg / this->m_mu;
  }

  return Erg;

};

///
/// An overloaded function to evaluate the growth restriction energy, its gradient
/// and its Hessian matrix.  For use with a fully implemented Newton method.
///
double RestrictGrowthOperator::operator()( Polyhedron &P,
    VectorXd &grad, SparseMatrix &hess ) {

  assert( false && "This functionality has not yet been implemented!" );

  return -1.0;

};

#endif
