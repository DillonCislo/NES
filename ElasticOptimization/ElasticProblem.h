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
#include "FixedVolumeOperator.h"
#include "RestrictGrowthOperator.h"
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
    /// The phantom polyhedral mesh used to calculated enclosed volumes
    /// for surfaces with holes
    ///
    Polyhedron m_PP;

		///
		/// A boolean determining whether any vertices are given target locations
		///
		bool anyFixed;

    ///
    /// A boolean determining whether the enclosed volume of the elastic body
    /// is fixed
    ///
    bool fixedVolume;

    ///
    /// A boolean determining whether the enclosed volume of the elastic body
    /// is calculated using the basic or the phantom triangulations
    ///
    bool usePhantom;

    ///
    /// A boolean determining whether growth should be restricted along a
    /// particular direction
    ///
    bool restrictGrowth;

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
    /// The fixed volume operator used to calculate the fixed volume energy
    ///
    FixedVolumeOperator m_FVO;

    ///
    /// The growth restriction operator used to calculate the growth
    /// restriction energy
    ///
    RestrictGrowthOperator m_RGO;

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
		/// Basic constructor.
		/// Note that we expect the polyhedron object to have a fully updated
		/// target geometry prior to construction of the problem structure
		///
		ElasticProblem( Polyhedron &P, double nu, double mu,
        double beta, double targetVolume, bool usePP, Polyhedron &PP,
        double alpha, const VectorXi &target_ID ) :
      m_P( P ), m_PP( PP ), usePhantom( usePP ) {

      // Create elastic energy operators
			this->m_SO = StretchOperator( nu );
			this->m_BO = BendOperator( nu );

      // Fixed point handling
      if ( alpha > 0.0 ) {

        this->anyFixed = true;
        int Nv = this->m_P.size_of_vertices();
        this->m_FPO = FixedPointOperator( Nv, alpha, target_ID );

      } else {

        this->anyFixed = false;

      }

      // Fixed volume handling
      if ( beta > 0.0 ) {

        this->fixedVolume = true;
        this->m_FVO = FixedVolumeOperator( targetVolume, beta );

      } else {

        this->fixedVolume = false;

      }

      // Growth restriction handling
      if ( mu > 0.0 ) {

        this->restrictGrowth = true;
        this->m_RGO = RestrictGrowthOperator( mu );

      } else {

        this->restrictGrowth = false;

      }


		};

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

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;


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
	if ( this->anyFixed ) { Etotal += this->m_FPO( this->m_P ); }

  // OPTIONAL: Calculate the fixed volume energy
  if ( this->fixedVolume ) {

    if ( usePhantom ) {

      this->m_EU.updateCurrentGeometry( this->m_PP, x );
      Etotal += this->m_FVO( this->m_PP );

    } else {

      Etotal += this->m_FVO( this->m_P );

    }

  }

  // OPTIONAL: Calculate the growth restriction energy
  if ( this->restrictGrowth ) { Etotal += this->m_RGO( this->m_P ); }

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

	// Calculate bending energy
	double Ebend = this->m_BO( this->m_P, grad );

	double Etotal = Estretch + Ebend;

	// OPTIONAL: Calculate the target vertex correspondence energy
	if ( this->anyFixed ) { Etotal += this->m_FPO( this->m_P, grad ); }

  // OPTIONAL: Calculate the fixed volume energy
  if ( this->fixedVolume ) {

    if ( usePhantom ) {

      this->m_EU.updateCurrentGeometry( this->m_PP, x );
      Etotal += this->m_FVO( this->m_PP, grad );

    } else {

      Etotal += this->m_FVO( this->m_P, grad );

    }

  }

  // OPTIONAL: Calculate the growth restriction energy
  if ( this->restrictGrowth ) { Etotal += this->m_RGO( this->m_P, grad ); }

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
	if ( this->anyFixed ) { Etotal += this->m_FPO( this->m_P, grad, spHess ); }

  // OPTIONAL: Calculate the fixed volume energy
  if ( this->fixedVolume ) {

    if ( usePhantom ) {

      this->m_EU.updateCurrentGeometry( this->m_PP, x );
      Etotal += this->m_FVO( this->m_PP, grad, spHess );

    } else {

      Etotal += this->m_FVO( this->m_P, grad, spHess );

    }

  }

  // OPTIONAL: Calculate the growth restriction energy
  if ( this->restrictGrowth ) { Etotal += this->m_RGO( this->m_P, grad, spHess ); }

	// Set the global Hessian matrix from stand-in
	hess = spHess;

	return Etotal;

};

#endif
