/*!
 *  \file FixedVolumeOperator.h
 *  \brief A class to calculate the fixed volume energy and its derivatives
 *
 *  \author Dillon Cislo
 *  \date 03/05/2020
 *
 */

#ifndef _FIXED_VOLUME_OPERATOR_H_
#define _FIXED_VOLUME_OPERATOR_H_

#include <vector>
#include <cassert>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include "../ElasticMesh/ElasticItems.h"

///
/// A class to calculate the fixed volume energy and its derivatives
///
class FixedVolumeOperator {

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
    typedef typename Eigen::Matrix<double, 9, 9>                Matrix9d;

    typedef typename Eigen::SparseMatrix<double>                SparseMatrix;
    typedef typename Eigen::Triplet<double>                     Triplet;

  protected:

    ///
    /// The target volume
    ///
    double m_tarVol;

    ///
    /// The fixed volume constraint coefficient
    ///
    double m_beta;

  public:

    ///
    /// Null constructor
    ///
    FixedVolumeOperator() {};

    ///
    /// Default constructor
    ///
    FixedVolumeOperator( double tarVol, double beta );

    ///
    /// An overloaded function calculate the fixed volume energy. For use when it
    /// is only necessary to calculate the energy and not any of its derivatives.
    ///
    double operator()( Polyhedron &P );

    ///
    /// An overloaded function to evaluate the fixed volume energy and its gradient.
    /// For use with gradient-based optimization methods ( FIRE, L-BFGS, ... )
    ///
    double operator()( Polyhedron &P, VectorXd &grad );

    ///
    /// An overloaded function to evaluate the fixed volume energy, its gradient
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
FixedVolumeOperator::FixedVolumeOperator( double tarVol, double beta ) :
  m_tarVol( tarVol ), m_beta( beta ) { };

///
/// An overloaded function that maps local gradient quantities to the
/// global gradient vector
///
void FixedVolumeOperator::mapLocalToGlobal( Face_handle f,
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
void FixedVolumeOperator::mapLocalToGlobal( Face_handle f, int Nv,
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
/// An overloaded function to calculate the fixed volume energy.
/// For use when it is only necessary to calculate the energy and not any of its derivatives
///
double FixedVolumeOperator::operator()( Polyhedron &P ) {

  // NOTE: CURRENT GEOMETRY OF THE POLYHEDRON SHOULD BE UP TO DATE

  double V = 0.0; // The total enclosed volume of the surface

  Facet_iterator f;
  for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

    // Extract face quantities
    Vector3d n = f->faceAreaWeightedNormal();
    Vector3d vi = f->halfedge()->next()->vertex()->v();
    Vector3d vj = f->halfedge()->prev()->vertex()->v();
    Vector3d vk = f->halfedge()->vertex()->v();

    // Calculate the facet centroid
    Vector3d R = ( vi + vj + vk ) / 3.0;

    // Calculate contribution to total volume
    V += R.dot(n);

  }

  V = V / 6.0;

  // Calculate the fixed volume energy
  double EFV = this->m_beta * ( V - this->m_tarVol ) * ( V - this->m_tarVol );
  EFV = EFV / ( this->m_tarVol * this->m_tarVol );

  return EFV;

};

///
/// An overloaded function to evaluate the fixed volume energy and
/// its gradient. For use with gradient-based optimization methods ( FIRE, L-BFGS, ... )
///
double FixedVolumeOperator::operator()( Polyhedron &P, VectorXd &grad ) {

  // NOTE: CURRENT GEOMETRY OF THE POLYHEDRON SHOULD BE UP TO DATE

  double V = 0.0; // The total enclosed volume of the surface

  // A temporary vector used to calculate the volume gradient
  VectorXd volGrad = VectorXd::Zero(grad.size());

  Facet_iterator f;
  for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

    // Extract face quantities
    Vector3d n = f->faceAreaWeightedNormal();

    Vector3d vi = f->halfedge()->next()->vertex()->v();
    Vector3d vj = f->halfedge()->prev()->vertex()->v();
    Vector3d vk = f->halfedge()->vertex()->v();

    Vector3d ei = f->halfedge()->edgeVector();
    Vector3d ej = f->halfedge()->next()->edgeVector();
    Vector3d ek = f->halfedge()->prev()->edgeVector();

    // Calculate the facet centroid
    Vector3d R = ( vi + vj + vk ) / 3.0;

    // Calculate the cross products of the position vector with the edge vectors
    RowVector3d Ri = R.cross(ei).transpose();
    RowVector3d Rj = R.cross(ej).transpose();
    RowVector3d Rk = R.cross(ek).transpose();

    // Calculate contribution to total volume
    V += R.dot(n);

    // Construct local contribution to energy gradient
    RowGradS gradEFV1, gradEFV2, gradEFV;

    gradEFV1 << n.transpose(), n.transpose(), n.transpose();
    gradEFV2 << Ri, Rj, Rk;

    gradEFV = (gradEFV1 / 3.0) + gradEFV2;

    // Map local gradient to global gradient vector
    this->mapLocalToGlobal( f, gradEFV, volGrad );

  }

  V = V / 6.0;

  double EFV = this->m_beta * ( V - this->m_tarVol ) * ( V - this->m_tarVol );
  EFV = EFV / ( this->m_tarVol * this->m_tarVol );

  volGrad = ( this->m_beta * ( V - this->m_tarVol ) / 3.0 ) * volGrad;
  volGrad = volGrad / ( this->m_tarVol * this->m_tarVol );
  grad += volGrad;

  return EFV;

};

///
/// An overloaded function to evaluate the fixed volume energy, its gradient,
/// and its Hessian matrix.  For use with a fully implemented Newton method.
///
double FixedVolumeOperator::operator()( Polyhedron &P, VectorXd &grad, SparseMatrix &hess ) {

  assert( false && "This functionality has not yet been implemented!" );

  return -1.0;

};

#endif
