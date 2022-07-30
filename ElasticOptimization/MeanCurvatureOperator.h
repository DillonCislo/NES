/*!
 * \file MeanCurvatureOperator.h
 * \brief A class to calculate an energy that smooths the mean
 * curvature of the embedding and its derivatives
 *
 * \author Dillon Cislo
 * \date 02/04/2022
 *
 */

#ifndef _MEAN_CURVATURE_OPERATOR_H_
#define _MEAN_CURVATURE_OPERATOR_H_

#include <vector>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include <igl/cotmatrix_intrinsic.h>
#include <igl/grad_intrinsic.h>

#include "../ElasticMesh/ElasticItems.h"

///
/// A class to calculate the energy that smooths the mean curvature of
/// the embedding
///
class MeanCurvatureOperator {

	private:

		typedef typename CGAL::Simple_cartesian<double> 		        Kernel;
		typedef typename CGAL::Polyhedron_3<Kernel, ElasticItems> 	Polyhedron;

		typedef typename Polyhedron::Vertex_iterator 		  Vertex_iterator;
		typedef typename Polyhedron::Facet_iterator 		  Facet_iterator;
		typedef typename Polyhedron::Edge_iterator 		    Edge_iterator;
		typedef typename Polyhedron::Halfedge_iterator 		Halfedge_iterator;

		typedef typename Polyhedron::Halfedge_handle 		Halfedge_handle;
		typedef typename Polyhedron::Face_handle 		    Face_handle;
		typedef typename Polyhedron::Vertex_handle 		  Vertex_handle;

		typedef typename Polyhedron::Halfedge_around_vertex_circulator 	HV_circulator;
		typedef typename Polyhedron::Halfedge_around_facet_circulator 	HF_circulator;

		typedef typename Eigen::RowVector3d 	RowVector3d;
		typedef typename Eigen::Vector3d 			Vector3d;
		typedef typename Eigen::VectorXi 			VectorXi;
		typedef typename Eigen::VectorXd 			VectorXd;
		typedef typename Eigen::MatrixXd 			MatrixXd;
    typedef typename Eigen::MatrixXi      MatrixXi;

		typedef typename Eigen::SparseMatrix<double> 		SparseMatrix;
		typedef typename Eigen::Triplet<double> 		    Triplet;

  protected:

    ///
    /// The total area of non-boundary faces in the mesh
    /// defined by the target edge lengths
    ///
    double m_totalAF = 0.0;

    ///
    /// The scalar weight of the mean curvature energy
    ///
    double m_kappa = 0.0;

    /// A length scale for the physical system. Used to normalize the energy
    double m_systemSize = 1.0;

  public:

    ///
    /// Null constructor
    ///
    MeanCurvatureOperator() {};

    ///
    /// Default constructor
    ///
    MeanCurvatureOperator( Polyhedron &P, double kappa );

    ///
    /// An overloaded function to calculate the mean curvature
    /// energy. For use when it is only necessary to calculate the
    /// energy and not any of its derivatives
    ///
    double operator()( Polyhedron &P );

    ///
    /// An overloaded function to evaluate the mean curvature
    /// energy and its gradient. For use with gradient based methods
    /// ( FIRE, L-BFGS, ... )
    ///
    double operator()( Polyhedron &P, VectorXd &grad );

    ///
    /// An overloaded function to evaluate the mean curvature energy,
    /// its gradient, and its Hessian matrix. For use with a
    /// fully implemented Newton method.
    ///
    double operator()( Polyhedron &P, VectorXd &grad, SparseMatrix &hess );

  protected:

    ///
    /// Branchless sign function
    ///
    inline double sgn( double val ) { return ( (0.0 < val) - (val < 0.0) ); };

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

};

///
/// Default constructor
///
MeanCurvatureOperator::MeanCurvatureOperator(
    Polyhedron &P, double kappa ) : m_kappa( kappa ) {

  // Accrue total mesh area excluding boundary faces
  double totalAF = 0.0;
  Facet_iterator f;
  for( f = P.facets_begin(); f != P.facets_end(); f++ )
    if (!f->isBoundary())
      totalAF += f->tarFaceArea();

  this->m_totalAF = totalAF;

};

///
/// An overloaded function to calculate the mean curvature
/// energy. For use when it is only necessary to calculate the
/// energy and not any of its derivatives
///
double MeanCurvatureOperator::operator()( Polyhedron &P ) {

  // NOTE: CURRENT MESH GEOMETRY SHOULD BE UP TO DATE
  // Update current vertex geometry
  Vertex_iterator v;
  for( v = P.vertices_begin(); v != P.vertices_end(); v++ ) {

    Vector3d vn = Vector3d::Zero();
    Vector3d Hn = Vector3d::Zero();

    HV_circulator he = v->vertex_begin();
    do {

      if ( !he->is_border() ) {

        // Accumulate the current face normal weighted by the
        // incident angle
        vn += he->prev()->targetAngle() * he->face()->faceNormal();

      }

      // Accumulate the current edge's contribution to the mean
      // curvature normal
      Hn += he->laplaceBeltramiWeight() *
        (he->opposite()->vertex()->v() - v->v());

      he++;

    } while ( he != v->vertex_begin() );

    // Normalize the vertex normal by the total incident angle sum
    vn = vn / v->targetAngleSum();
    vn = vn / std::sqrt( vn.dot(vn) );

    // Normalize the mean curvature normal by the vertex Voronoi area
    Hn = Hn / (2.0 * v->targetArea());

    double absH = std::sqrt( Hn.dot(Hn) );
    double sgnH = this->sgn( Hn.dot(vn) );
    double H = sgnH * absH;

    Hn = Hn / absH;

    if ( sgnH == 0.0 )
      throw std::runtime_error("Mean curvature normal is perpendicular to vertex normal");

    v->setVertexNormal( vn );
    v->setMeanCurvatureNormal( Hn );
    v->setAbsH( absH );
    v->setSgnH( sgnH );
    v->setH( H );

  }

  // Calculate mean curvature smoothness energy
  double EH = 0.0;
  Facet_iterator f;
  for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

    // Ignore the contruibution of boundary faces
    if ( f->isBoundary() )
      continue;

    double Li2 = f->halfedge()->targetLength(); Li2 *= Li2;
    double Lj2 = f->halfedge()->next()->targetLength(); Lj2 *= Lj2;
    double Lk2 = f->halfedge()->prev()->targetLength(); Lk2 *= Lk2;

    double Hi = f->halfedge()->next()->vertex()->H();
    double Hj = f->halfedge()->prev()->vertex()->H();
    double Hk = f->halfedge()->vertex()->H();

    double curEH = 0.0;
    curEH += (Hk - Hi) * (Hj - Hi) * Li2;
    curEH += (Hi - Hj) * (Hk - Hj) * Lj2;
    curEH += (Hj - Hk) * (Hi - Hk) * Lk2;
    curEH = curEH / f->tarFaceArea();

    EH += curEH;

  }

  EH = this->m_kappa * EH / ( 4.0 * this->m_totalAF );

  return EH;

};

///
/// An overloaded function to evaluate the mean curvature
/// energy and its gradient. For use with gradient based methods
/// ( FIRE, L-BFGS, ... )
///
double MeanCurvatureOperator::operator()( Polyhedron &P, VectorXd &grad ) {

  int Nv = P.size_of_vertices();

  // NOTE: CURRENT MESH GEOMETRY SHOULD BE UP TO DATE
  // Update current vertex geometry
  Vertex_iterator v;
  for( v = P.vertices_begin(); v != P.vertices_end(); v++ ) {

    Vector3d vn = Vector3d::Zero();
    Vector3d Hn = Vector3d::Zero();

    HV_circulator he = v->vertex_begin();
    do {

      if ( !he->is_border() ) {

        // Accumulate the current face normal weighted by the
        // incident angle
        vn += he->prev()->targetAngle() * he->face()->faceNormal();

      }

      // Accumulate the current edge's contribution to the mean
      // curvature normal
      Hn += he->laplaceBeltramiWeight() *
        (he->opposite()->vertex()->v() - v->v());

      he++;

    } while ( he != v->vertex_begin() );

    // Normalize the vertex normal by the total incident angle sum
    vn = vn / v->targetAngleSum();
    vn = vn / std::sqrt( vn.dot(vn) );

    // Normalize the mean curvature normal by the vertex Voronoi area
    Hn = Hn / (2.0 * v->targetArea());

    double absH = std::sqrt( Hn.dot(Hn) );
    double sgnH = this->sgn( Hn.dot(vn) );
    double H = sgnH * absH;

    Hn = Hn / absH;

    if ( sgnH == 0.0 )
      throw std::runtime_error("Mean curvature normal is perpendicular to vertex normal");

    v->setVertexNormal( vn );
    v->setMeanCurvatureNormal( Hn );
    v->setAbsH( absH );
    v->setSgnH( sgnH );
    v->setH( H );

  }

  // Calculate mean curvature smoothness energy and its gradient
  double EH = 0.0;
  Facet_iterator f;
  for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

    // Ignore the contruibution of boundary faces
    if ( f->isBoundary() )
      continue;

    double Li2 = f->halfedge()->targetLength(); Li2 *= Li2;
    double Lj2 = f->halfedge()->next()->targetLength(); Lj2 *= Lj2;
    double Lk2 = f->halfedge()->prev()->targetLength(); Lk2 *= Lk2;

    Vertex_handle vi = f->halfedge()->next()->vertex();
    Vertex_handle vj = f->halfedge()->prev()->vertex();
    Vertex_handle vk = f->halfedge()->vertex();

    double Hi = vi->H();
    double Hj = vj->H();
    double Hk = vk->H();

    double sgnHi = vi->sgnH();
    double sgnHj = vj->sgnH();
    double sgnHk = vk->sgnH();

    Vector3d Hni = vi->Hn();
    Vector3d Hnj = vj->Hn();
    Vector3d Hnk = vk->Hn();

    double Ai = vi->targetArea();
    double Aj = vj->targetArea();
    double Ak = vk->targetArea();

    // Energy calculation ---------------------------------------------------------------

    double curEH = 0.0;
    curEH += (Hk - Hi) * (Hj - Hi) * Li2;
    curEH += (Hi - Hj) * (Hk - Hj) * Lj2;
    curEH += (Hj - Hk) * (Hi - Hk) * Lk2;
    curEH = curEH / f->tarFaceArea();

    EH += curEH;

    // Gradient calculation -------------------------------------------------------------

    // Gradient coefficients
    double gradHICoeff, gradHJCoeff, gradHKCoeff;
    gradHICoeff = (2.0 * Hi - Hj - Hk) * Li2 + (Hj - Hk) * (Lk2 - Lj2);
    gradHJCoeff = (2.0 * Hj - Hk - Hi) * Lj2 + (Hk - Hi) * (Li2 - Lk2);
    gradHKCoeff = (2.0 * Hk - Hi - Hj) * Lk2 + (Hi - Hj) * (Lj2 - Li2);

    gradHICoeff *= this->m_kappa / ( 4.0 * this->m_totalAF * f->tarFaceArea() );
    gradHJCoeff *= this->m_kappa / ( 4.0 * this->m_totalAF * f->tarFaceArea() );
    gradHKCoeff *= this->m_kappa / ( 4.0 * this->m_totalAF * f->tarFaceArea() );

    // Accumulate the contributions from the first vertex
    Vector3d dHIdXI = -sgnHi * vi->laplaceBeltramiSum() * Hni / (2.0 * Ai);
    grad( vi->id() ) += gradHICoeff * dHIdXI(0);
    grad( vi->id() + Nv ) += gradHICoeff * dHIdXI(1);
    grad( vi->id() + 2*Nv ) += gradHICoeff * dHIdXI(2);

    HV_circulator he = vi->vertex_begin();
    do {

      if ( !he->is_border() ) {

        Vertex_handle vm = he->opposite()->vertex();
        Vector3d dHIdXM = sgnHi * he->laplaceBeltramiWeight() * Hni / (2.0 * Ai);

        grad( vm->id() ) += gradHICoeff * dHIdXM(0);
        grad( vm->id() + Nv ) += gradHICoeff * dHIdXM(1);
        grad( vm->id() + 2*Nv ) += gradHICoeff * dHIdXM(2);

      }

      he++;

    } while ( he != vi->vertex_begin() );

    // Accumulate the contributions from the second vertex
    Vector3d dHJdXJ = -sgnHj * vj->laplaceBeltramiSum() * Hnj / (2.0 * Aj);
    grad( vj->id() ) += gradHJCoeff * dHJdXJ(0);
    grad( vj->id() + Nv ) += gradHJCoeff * dHJdXJ(1);
    grad( vj->id() + 2*Nv ) += gradHJCoeff * dHJdXJ(2);

    he = vj->vertex_begin();
    do {

      if ( !he->is_border() ) {

        Vertex_handle vm = he->opposite()->vertex();
        Vector3d dHJdXM = sgnHj * he->laplaceBeltramiWeight() * Hnj / (2.0 * Aj);

        grad( vm->id() ) += gradHJCoeff * dHJdXM(0);
        grad( vm->id() + Nv ) += gradHJCoeff * dHJdXM(1);
        grad( vm->id() + 2*Nv ) += gradHJCoeff * dHJdXM(2);

      }

      he++;

    } while ( he != vj->vertex_begin() );

    // Accumulate the contributions from the third vertex
    Vector3d dHKdXK = -sgnHk * vk->laplaceBeltramiSum() * Hnk / (2.0 * Ak);
    grad( vk->id() ) += gradHKCoeff * dHKdXK(0);
    grad( vk->id() + Nv ) += gradHKCoeff * dHKdXK(1);
    grad( vk->id() + 2*Nv ) += gradHKCoeff * dHKdXK(2);

    he = vk->vertex_begin();
    do {

      if ( !he->is_border() ) {

        Vertex_handle vm = he->opposite()->vertex();
        Vector3d dHKdXM = sgnHk * he->laplaceBeltramiWeight() * Hnk / (2.0 * Ak);

        grad( vm->id() ) += gradHKCoeff * dHKdXM(0);
        grad( vm->id() + Nv ) += gradHKCoeff * dHKdXM(1);
        grad( vm->id() + 2*Nv ) += gradHKCoeff * dHKdXM(2);

      }

      he++;

    } while ( he != vk->vertex_begin() );

  }

  EH = this->m_kappa * EH / ( 4.0 * this->m_totalAF );

  return EH;

};

///
/// An overloaded function to evaluate the mean curvature energy,
/// its gradient, and its Hessian matrix. For use with a
/// fully implemented Newton method.
///
double MeanCurvatureOperator::operator()(
    Polyhedron &P, VectorXd &grad, SparseMatrix &hess ) {

  throw std::runtime_error("This method is not yet supported");

  return EXIT_FAILURE;

};


#endif
