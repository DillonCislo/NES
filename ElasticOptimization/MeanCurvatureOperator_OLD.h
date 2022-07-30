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
    /// The number of vertices in the mesh
    ///
    int m_Nv;

    ///
    /// The sparse mean curvature gradient operator
    ///
    SparseMatrix m_MGO;

    ///
    /// The total area of non-boundary faces in the mesh
    /// defined by the target edge lengths
    ///
    double m_totalAF;

    ///
    /// The scalar weight of the mean curvature energy
    ///
    double m_kappa;

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
    /// Construct the FEM gradient operator
    ///
    inline void build_grad( const MatrixXd V, const MatrixXi &F, SparseMatrix &G );

    ///
    /// An overloaded function to calculate the mean curvature
    /// energy. For use when it is only necessary to calculate the
    /// energy and not any of its derivatives
    ///
    double operator()( const VectorXd &x );

    ///
    /// An overloaded function to evaluate the mean curvature
    /// energy and its gradient. For use with gradient based methods
    /// ( FIRE, L-BFGS, ... )
    ///
    double operator()( const VectorXd &x, VectorXd &grad );

    ///
    /// An overloaded function to evaluate the mean curvature energy,
    /// its gradient, and its Hessian matrix. For use with a
    /// fully implemented Newton method.
    ///
    double operator()( const VectorXd &x, VectorXd &grad, SparseMatrix &hess );

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

};

///
/// Default constructor
///
MeanCurvatureOperator::MeanCurvatureOperator(
    Polyhedron &P, double kappa ) : m_kappa( kappa ) {

  std::cout << "MCO Check 0" << std::endl;

  // The total number of vertices in the mesh
  int Nv = P.size_of_vertices();

  // The total number of faces in the mesh
  int Nf = P.size_of_facets();

  // Accrue target vertex area matrix ---------------------------------------------------

  SparseMatrix invAV(3*Nv, 3*Nv);
  invAV.reserve(3*Nv);

  std::vector<Triplet > IJV_AV;
  IJV_AV.reserve(3*Nv);

  Vertex_iterator v;
  for( v = P.vertices_begin(); v != P.vertices_end(); v++ ) {

    int vID = v->id();
    double AV = 0.0;

    HV_circulator he = v->vertex_begin();
    do {

      if ( !he->is_border() )
        AV += he->face()->tarFaceArea();

      he++;

    } while ( he != v->vertex_begin() );

    if (AV == 0.0)
      throw std::runtime_error("Target vertex area equals zero");

    IJV_AV.push_back( Triplet( vID, vID, 1.0 / (AV / 3.0) ) );
    IJV_AV.push_back( Triplet( vID+Nv, vID+Nv, 1.0 / (AV / 3.0) ) );
    IJV_AV.push_back( Triplet( vID+2*Nv, vID+2*Nv, 1.0 / (AV / 3.0) ) );

  }

  invAV.setFromTriplets( IJV_AV.begin(), IJV_AV.end() );

  std::cout << "MCO Check 1" << std::endl;

  // Accrue target face area weighting matrix -------------------------------------------

  SparseMatrix AF(6*Nf, 6*Nf);
  AF.reserve(6*Nf);

  std::vector<Triplet > IJV_AF;
  IJV_AF.reserve(6*Nf);

  Facet_iterator f;
  int fID = 0;
  double totalAF = 0.0;
  for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

    // Ignore faces on the boundary
    bool is_border_face = f->halfedge()->is_border_edge() ||
      f->halfedge()->next()->is_border_edge() ||
      f->halfedge()->prev()->is_border_edge();


    if ( !is_border_face ) {

      IJV_AF.push_back( Triplet( fID, fID, f->tarFaceArea() ) );
      IJV_AF.push_back( Triplet( fID+Nf, fID+Nf, f->tarFaceArea() ) );
      IJV_AF.push_back( Triplet( fID+2*Nf, fID+2*Nf, f->tarFaceArea() ) );
      IJV_AF.push_back( Triplet( fID+3*Nf, fID+3*Nf, f->tarFaceArea() ) );
      IJV_AF.push_back( Triplet( fID+4*Nf, fID+4*Nf, f->tarFaceArea() ) );
      IJV_AF.push_back( Triplet( fID+5*Nf, fID+5*Nf, f->tarFaceArea() ) );

      totalAF += f->tarFaceArea();

    }

    fID++;

  }

  AF.setFromTriplets( IJV_AF.begin(), IJV_AF.end() );
  AF = kappa * AF / totalAF;
  m_totalAF = totalAF;

  std::cout << "MCO Check 2" << std::endl;

  // Generate cotan Laplacian matrix from target edge lengths ---------------------------

  SparseMatrix L(3*Nv, 3*Nv);
  L.reserve(30*Nv);

  // Convenience variable for sparse matrix construction
  MatrixXi edges(3,2);
  edges <<
    1,2,
    2,0,
    0,1;

  // Gather cotangents, edge lengths, and face connectivity list
  // ( This definitely does not exploit the power of the )
  // ( DCEL, but it only needs to be called once         )
  MatrixXd C = MatrixXd::Zero(Nf, 3);
  MatrixXi F = MatrixXi::Zero(Nf, 3);
  MatrixXd L_F = MatrixXd::Zero(Nf, 3);
  fID = 0;
  for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

    double tarA = f->tarFaceArea();

    double Li = f->halfedge()->targetLength();
    double Lj = f->halfedge()->next()->targetLength();
    double Lk = f->halfedge()->prev()->targetLength();

    L_F(fID, 0) = Li;
    L_F(fID, 1) = Lj;
    L_F(fID, 2) = Lk;

    double Li2 = Li * Li;
    double Lj2 = Lj * Lj;
    double Lk2 = Lk * Lk;

    C(fID, 0) = ( Lj2 + Lk2 - Li2 ) / tarA / 8.0;
    C(fID, 1) = ( Lk2 + Li2 - Lj2 ) / tarA / 8.0;
    C(fID, 2) = ( Li2 + Lj2 - Lk2 ) / tarA / 8.0;

    int vi = f->halfedge()->next()->vertex()->id();
    int vj = f->halfedge()->prev()->vertex()->id();
    int vk = f->halfedge()->vertex()->id();

    F(fID, 0) = vi;
    F(fID, 1) = vj;
    F(fID, 2) = vk;

    fID++;

  }

  // Construct basic Laplacian
  std::vector<Triplet > IJV;
  IJV.reserve( 36 * Nf );
  for( int i = 0; i < Nf; i++ ) {
    for( int e = 0; e < 3; e++ ) {

      int source = F(i, edges(e,0));
      int dest = F(i, edges(e,1));

      IJV.push_back( Triplet( source, dest, C(i,e) ) );
      IJV.push_back( Triplet( dest, source, C(i,e) ) );
      IJV.push_back( Triplet( source, source, -C(i,e) ) );
      IJV.push_back( Triplet( dest, dest, -C(i,e) ) );

      IJV.push_back( Triplet( source+Nv, dest+Nv, C(i,e) ) );
      IJV.push_back( Triplet( dest+Nv, source+Nv, C(i,e) ) );
      IJV.push_back( Triplet( source+Nv, source+Nv, -C(i,e) ) );
      IJV.push_back( Triplet( dest+Nv, dest+Nv, -C(i,e) ) );

      IJV.push_back( Triplet( source+2*Nv, dest+2*Nv, C(i,e) ) );
      IJV.push_back( Triplet( dest+2*Nv, source+2*Nv, C(i,e) ) );
      IJV.push_back( Triplet( source+2*Nv, source+2*Nv, -C(i,e) ) );
      IJV.push_back( Triplet( dest+2*Nv, dest+2*Nv, -C(i,e) ) );

    }
  }

  L.setFromTriplets( IJV.begin(), IJV.end() );

  SparseMatrix LCheck(Nv, Nv);
  igl::cotmatrix_intrinsic(L_F, F, LCheck);
  double Lerr = 0.0;
  for( int i = 0; i < LCheck.rows(); i++ )
    for( int j = 0; j < LCheck.cols(); j++ )
      Lerr = std::max( Lerr, std::abs(LCheck.coeff(i,j)-L.coeff(i,j)) );

  std::cout << "Maximum L Error = " << Lerr << std::endl;

  // Normalize by inverse vertex areas
  L = invAV * L;

  std::cout << "MCO Check 3" << std::endl;

  // Generate intrinsic gradient operator -----------------------------------------------

  // Generate phantom vertex locations
  MatrixXd V2 = MatrixXd::Zero(3*Nf, 2);

  // Place 3rd vertex at [L_F(:,1), 0]
  V2.block(2*Nf, 0, Nf, 1) = L_F.col(0);

  // Place 2nd vertex at [0 0];

  // Place 1st vertex at [x,y]
  V2.block(0, 0, Nf, 1) =
    (L_F.col(1).cwiseAbs2()-L_F.col(0).cwiseAbs2()-L_F.col(2).cwiseAbs2()).array() /
    ( -2.0 * L_F.col(0)).array();
  V2.block(0, 1, Nf, 1) =
    (L_F.col(2).cwiseAbs2() - V2.block(0, 0, Nf, 1).cwiseAbs2()).array().sqrt();

  // Construct phantom face list and correspondence operator
  MatrixXi F2(Nf, 3);
  std::vector<Triplet > PPijv;
  PPijv.reserve(F.size());
  for( int i = 0; i < Nf; i++ ) {
    for( int j = 0; j < 3; j++ ) {

      F2(i,j) = i + j*Nf;
      PPijv.emplace_back( F2(i,j), F(i,j), 1 );

    }
  }

  SparseMatrix PP(3*Nf, Nv);
  PP.setFromTriplets( PPijv.begin(), PPijv.end() );

  std::cout << "MCO Check 4" << std::endl;

  // Construct explicit gradient operator from phantom face/vertex lists
  SparseMatrix G2;
  this->build_grad(V2, F2, G2);

  std::cout << "MCO Check 5" << std::endl;

  // Construct complete #(2*F)x#V intrinsic gradient matrix
  SparseMatrix G1 = G2 * PP;

  // Convert to #(6*F)x#(3*V) block diagonal matrix
  SparseMatrix G(6*Nf, 3*Nv);
  G.reserve(6*Nf*Nv);

  std::vector<Triplet > Gijv;
  Gijv.reserve(6*Nf*Nv);
  for( int i = 0; i < G1.rows(); i++ ) {
    for( int j = 0; j < G1.cols(); j++ ) {

      double Gij = G1.coeff(i,j);

      if (Gij != 0.0) {

        Gijv.push_back( Triplet( i, j, Gij ) );
        Gijv.push_back( Triplet( i+2*Nf, j+Nv, Gij ) );
        Gijv.push_back( Triplet( i+4*Nf, j+2*Nv, Gij ) );

      }

    }
  }

  G.setFromTriplets( Gijv.begin(), Gijv.end() );

  SparseMatrix GCheck(2*Nf, Nv);
  igl::grad_intrinsic(L_F, F, GCheck);
  double Gerr = 0.0;
  for( int i = 0; i < GCheck.rows(); i++ )
    for( int j = 0; j < GCheck.cols(); j++ )
      Gerr = std::max( Gerr, std::abs(GCheck.coeff(i,j)-G.coeff(i,j)) );

  std::cout << "Maximum G error = " << Gerr << std::endl;

  std::cout << "MCO Check 6" << std::endl;

  // Generate mean curvature gradient operator ------------------------------------------

  m_MGO = L.transpose() * G.transpose() * AF * G * L;

  std::cout << "MCO Check 7" << std::endl;

};


///
/// Build FEM gradient operator
///
inline void MeanCurvatureOperator::build_grad(
    const MatrixXd V, const MatrixXi &F, SparseMatrix &G ) {

  const int m = F.rows(); // Number of faces
  const int nv = V.rows(); // Number of vertices
  const int dims = V.cols(); // Number of dimensions

  Eigen::Matrix<double, Eigen::Dynamic, 3>
    eperp21(m,3), eperp13(m,3);

  for( int i = 0; i < m; i++ ) {

    // Renaming indices of vertices of triangles for convenience
    int i1 = F(i,0);
    int i2 = F(i,1);
    int i3 = F(i,2);

    // Triangle edge vectors, named after opposite vertices
    RowVector3d v32 = RowVector3d::Zero(1,3);
    RowVector3d v13 = RowVector3d::Zero(1,3);
    RowVector3d v21 = RowVector3d::Zero(1,3);
    v32.head(dims) = V.row(i3) - V.row(i2);
    v13.head(dims) = V.row(i1) - V.row(i3);
    v21.head(dims) = V.row(i2) - V.row(i1);

    // Face normal vectors
    RowVector3d n = v32.cross(v13);
    double dblA = std::sqrt(n.dot(n));
    Vector3d u = n / dblA;

    // Rotate each vector 90 degrees around normal
    double norm21 = std::sqrt(v21.dot(v21));
    double norm13 = std::sqrt(v13.dot(v13));

    eperp21.row(i) = u.cross(v21);
    eperp21.row(i) = eperp21.row(i) / std::sqrt(eperp21.row(i).dot(eperp21.row(i)));
    eperp21.row(i) *= norm21 / dblA;

    eperp13.row(i) = u.cross(v13);
    eperp13.row(i) = eperp13.row(i) / std::sqrt(eperp13.row(i).dot(eperp13.row(i)));
    eperp13.row(i) *= norm13 / dblA;

  }

  // Create sparse gradient operator matrix
  G.resize(dims*m, nv);
  std::vector<Triplet > Gijv;
  Gijv.reserve(4*dims*m);
  for( int f = 0; f < m; f++ ) {
    for( int d = 0; d < dims; d++ ) {

      Gijv.emplace_back( f+d*m, F(f,1),  eperp13(f,d) );
      Gijv.emplace_back( f+d*m, F(f,0), -eperp13(f,d) );
      Gijv.emplace_back( f+d*m, F(f,2),  eperp21(f,d) );
      Gijv.emplace_back( f+d*m, F(f,0), -eperp21(f,d) );

    }
  }

  G.setFromTriplets( Gijv.begin(), Gijv.end() );

};

///
/// An overloaded function to calculate the mean curvature
/// energy. For use when it is only necessary to calculate the
/// energy and not any of its derivatives
///
double MeanCurvatureOperator::operator()( const VectorXd &x ) {

  double EH = x.transpose() * m_MGO * x;

  return EH;

};

///
/// An overloaded function to evaluate the mean curvature
/// energy and its gradient. For use with gradient based methods
/// ( FIRE, L-BFGS, ... )
///
double MeanCurvatureOperator::operator()( const VectorXd &x, VectorXd &grad ) {

  VectorXd gradEH = m_MGO * x;

  double EH = x.transpose() * gradEH;

  grad = ( grad.array() + (2.0 * gradEH.array()) ).matrix();

  return EH;

};

///
/// An overloaded function to evaluate the mean curvature energy,
/// its gradient, and its Hessian matrix. For use with a
/// fully implemented Newton method.
///
double MeanCurvatureOperator::operator()(
    const VectorXd &x, VectorXd &grad, SparseMatrix &hess ) {

  hess += 2.0 * m_MGO;

  VectorXd gradEH = m_MGO * x;

  double EH = x.transpose() * gradEH;

  grad = ( grad.array() + (2.0 * gradEH.array()) ).matrix();

  return EH;

};


#endif
