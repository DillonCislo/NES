/*!
 * 	\file ElasticHalfEdge.h
 * 	\brief Extension of CGAL base halfedge class
 *
 * 	This subclass extends the CGAL base halfedge class to include functionalities
 * 	needed to calculate the elastic energy/gradients of a Discrete Non-Euclidean
 * 	Koiter Surface
 *
 * 	\author Dillon Cislo
 * 	\date 11/01/2018
 *
 */

#ifndef _ELASTIC_HALFEDGE_H_
#define _ELASTIC_HALFEDGE_H_

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/HalfedgeDS_halfedge_max_base_with_id.h>
#include <Eigen/Core>
#include <cmath>

//! The extension of the base CGAL facet class to the elastic mesh
template <class Refs>
struct ElasticHalfEdge : public CGAL::HalfedgeDS_halfedge_max_base_with_id<Refs, std::size_t> {

	public:

		typedef typename Eigen::Vector3d	 	          Vector3d;
		typedef typename Eigen::Matrix<double,1,12>	  HingeGrad;
		typedef typename Eigen::Matrix<double,12,12> 	HingeHess;

	protected:

		//! The directed edge vector associated with the halfedge
		Vector3d m_edgeVector = Vector3d::Zero();

		//! The directed unit edge vector associated with the halfedge
		Vector3d m_edgeUnitVector = Vector3d::Zero();

    //! The a precomputed quantity used to calculate the restriction energy
    double m_initEdgeGrowth = 0.0;

		//! The in-plane mid-edge unit normal
		Vector3d m_edgeNormal = Vector3d::Zero();

		//! The length of the halfedge
		double m_length = 0.0;

		//! The target length of the halfedge
		double m_targetLength = 0.0;

    /// The maximum length of the projection of the halfedge
    /// along its respective restriction vector
    double m_restrictLength = 0.0;

		//! The strain component of associated with the edge (see documentation)
		double m_edgeStrain = 0.0;

		//! The hinge angle function of the edge (see documentation)
		double m_Phi = 0.0;

		//! The target hinge angle function of the edge (see documentation)
		double m_tarPhi = 0.0;

		//! The gradient of the hinge angle function defined on a hinge stencil
		HingeGrad m_gradPhi = HingeGrad::Zero();

		//! The Hessian of the hinge angle function defined on a hinge stencil
		HingeHess m_hessPhi = HingeHess::Zero();

		//! Determines whether this halfedge is the major edge of a hinge
		bool m_isMajor = false;

    //! The edge weight of the target cotan Laplace-Beltrami operator
    double m_laplaceBeltramiWeight = 0.0;

    //! Target internal angle opposite the halfedge (defined on [0, 2*pi))
    double m_targetAngle = 0.0;

	public:

		/*******************************************************************************
		 * SETTERS
		 ******************************************************************************/

		//! Calculate all quantites associated with directed edge
		void calculateEdgeQuantities();

    //! Set the inital edge growth
    void setInitialEdgeGrowth( double initEdgeGrowth ) {
      this->m_initEdgeGrowth = initEdgeGrowth;
    };

		//! Set the target edge length
		void setTargetEdgeLength( double tarL ) {
			this->m_targetLength = tarL;
		};

    //! Set the restricted edge length
    void setRestrictedLength( double restrictL ) {
      this->m_restrictLength = restrictL;
    };

		//! Set the target hinge function
		void setTargetHingeFunction( double tarPhi ) {
			this->m_tarPhi = tarPhi;
		};

		//! Set the hinge function
		void setPhi( double Phi ) {
			this->m_Phi = Phi;
		};

		//! Set the gradient of the hinge angle function
		void setGradPhi( const HingeGrad &gradPhi ) {
			this->m_gradPhi = gradPhi;
		};

		//! Set the Hessian of the hinge angle function
		void setHessianPhi( const HingeHess &hessPhi ) {
			this->m_hessPhi = hessPhi;
		};

		//! Set the major edge indicator
		void setMajor( bool isMajor ) {
			this->m_isMajor = isMajor;
		};

    //! Set the edge weight of the target cotan Laplace-Beltrami operator
    void setLaplaceBeltramiWeight( double LBW ) {
      this->m_laplaceBeltramiWeight = LBW;
    };

    //! Set the internal angle opposite the halfedge
    void setTargetAngle( double ang ) {
      this->m_targetAngle = ang;
    };

		/*******************************************************************************
		 * GETTERS
		 ******************************************************************************/

		//! Get the directed edge vector
		Vector3d edgeVector() { return this->m_edgeVector; };

		//! Get the directed edge unit vector
		Vector3d edgeUnitVector() { return this->m_edgeUnitVector; };

    //! Get the initial edge growth
    double initEdgeGrowth() { return this->m_initEdgeGrowth; };

		//! Get the in-plane mid-edge unit normal vector
		Vector3d edgeNormal() { return this->m_edgeNormal; };

		//! Get the edge length
		double length() { return this->m_length; };

		//! Get the target edge length
		double targetLength() { return this->m_targetLength; };

    //! Get the restricted length
    double restrictedLength() { return this->m_restrictLength; };

		//! Get the target hinge function
		double tarPhi() { return this->m_tarPhi; };

		//! Get the edge strain
		double edgeStrain() { return this->m_edgeStrain; };

		//! Get the hinge angle function
		double Phi() { return this->m_Phi; };

		//! Get the gradient of the hinge angle function
		HingeGrad gradPhi() { return this->m_gradPhi; };

		//! Get the Hessian of the hinge angle function
		HingeHess hessPhi() { return this->m_hessPhi; };

		//! Get the major edge indicator
		bool isMajor() { return this->m_isMajor; };

    //! Get the edge weight of the target cotan Laplace-Beltrami operator
    double laplaceBeltramiWeight() { return this->m_laplaceBeltramiWeight; };

    //! Get the internal angle opposite the halfedge
    double targetAngle() { return this->m_targetAngle; };

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

};

//! Calculate all quantities associate with a directed edge
template <class Refs>
void ElasticHalfEdge<Refs>::calculateEdgeQuantities() {

	// NOTE: ALL VERTEX POSITIONS SHOULD BE UP TO DATE
	// NOTE: ALL FACE AREAS SHOULD BE UP TO DATE
	this->m_edgeVector = this->vertex()->v() - this->prev()->vertex()->v();
	this->m_length = this->m_edgeVector.norm();

	Vector3d eHat = this->m_edgeVector / this->m_length;
	this->m_edgeUnitVector = eHat;

	this->m_edgeStrain = ( this->m_length ) * ( this->m_length ) -
		( this->m_targetLength ) * ( this->m_targetLength );

	if( this->face() == NULL ) return;
	Vector3d n = this->face()->faceNormal();
	this->m_edgeNormal = eHat.cross( n );

};

#endif
