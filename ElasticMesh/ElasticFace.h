/*!
 * 	\file ElasticFace.h
 * 	\brief Extenstion of CGAL base facet class
 *
 * 	This sub class extends the CGAL base facet class to include functionalities
 * 	needed to calculated the elastic energy/gradients of a Discrete Non-Euclidean
 * 	Koiter Surface
 *
 * 	\author Dillon Cislo
 * 	\date 12/20/2018
 *
 */

#ifndef _ELASTIC_FACE_H_
#define _ELASTIC_FACE_H_

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/HalfedgeDS_face_max_base_with_id.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

typedef CGAL::Simple_cartesian<double>	 		Kernel;
typedef Kernel::Plane_3 		 		Plane;

//! The extension of the base CGAL facet class to the elastic mesh
template <class Refs>
struct ElasticFace : public CGAL::HalfedgeDS_face_max_base_with_id<Refs, Plane, std::size_t> {

	public:

		typedef typename Refs::Halfedge_handle 		    Halfedge_handle;
		typedef typename Refs::Halfedge_const_handle 	Halfedge_const_handle;

		typedef typename Eigen::Vector3d		        Vector3d;
		typedef typename Eigen::Matrix3d		        Matrix3d;
		typedef typename Eigen::Matrix<double,9,9> 	Matrix9d;

	protected:

		//! The area of the face
		double m_faceArea = 0.0;

		//! The target area of the face
		double m_tarFaceArea = 0.0;

    //! The thickness of the face
    double m_thickness = 0.0;

		//! The unit normal vector of the face
		Vector3d m_faceNormal = Vector3d::Zero();

    //! The double area weighted normal vector of the face
    Vector3d m_faceAreaWeightedNormal = Vector3d::Zero();

		//! The stretching trace constants
		Vector3d m_C = Vector3d::Zero();

		//! The target edge heights
		Vector3d m_hBar = Vector3d::Zero();

		//! The stretching energy constant matrix
		Matrix3d m_zeta = Matrix3d::Zero();

		//! The bending energy constant matrix
		Matrix3d m_xi = Matrix3d::Zero();

		//! The Hessian of the trace of the strain tensor
		Matrix9d m_hessTrE = Matrix9d::Zero();

    //! A unit tangent vector in the plane of the face along which to restrict growth
    Vector3d m_restrictVector = Vector3d::Zero();

    //! A boolean specifying whether the face lies on the mesh boundary
    bool m_isBoundary = false;

	public:

		/*******************************************************************************
		 * SETTERS
		 ******************************************************************************/

		//! Calculate face area and unit normal vector
		void calculateFaceAreaAndNormal();

		//! Set target face area
		void setTargetFaceArea( const double &tarFaceArea ) {
			this->m_tarFaceArea = tarFaceArea;
		};

    //! Set the face thickness
    void setThickness( const double &thickness ) {
      this->m_thickness = thickness;
    };

		//! Set the stretching trace constants
		void setC( const Vector3d &C ) {
			this->m_C = C;
		};

		//! Set the target edge heights
		void setHBar( const Vector3d &hBar ) {
			this->m_hBar = hBar;
		};

		//! Set the stretching energy constant matrix
		void setZeta( const Matrix3d &zeta ) {
			this->m_zeta = zeta;
		}

		//! Set the bending energy constant matrix
		void setXi( const Matrix3d &xi ) {
			this->m_xi = xi;
		};

		//! Set the Hessian of the trace of the strain tensor
		void setHessTrE( const Matrix9d &hE ) {
			this->m_hessTrE = hE;
		};

    //! Set the restriction vector
    void setRestrictionVector( const Vector3d &restrictVector ) {
      this->m_restrictVector = restrictVector.normalized();
    };

    //! Set whether the face is on the mesh boundary
    void setBoundaryTag( bool isBoundary ) {
      this->m_isBoundary = isBoundary;
    };

		/*******************************************************************************
		 * GETTERS
		 ******************************************************************************/

		//! Get face area
		double faceArea() { return this->m_faceArea; };

    //! Get face thickness
    double thickness() { return this->m_thickness; };

		//! Get face normal
		Vector3d faceNormal() { return this->m_faceNormal; };

    //! Get double area weighted face normal
    Vector3d faceAreaWeightedNormal() { return this->m_faceAreaWeightedNormal; };

		//! Get target face area
		double tarFaceArea() { return this->m_tarFaceArea; };

		//! Get the stretching trace constants
		Vector3d C() { return this->m_C; };

		//! Get the target edge heights
		Vector3d hBar() { return this->m_hBar; };

		//! Get the stretching energy constant matrix
		Matrix3d zeta() { return this->m_zeta; };

		//! Get the bending energy constant matrix
		Matrix3d xi() { return this->m_xi; };

		//! Get the Hessian of the trace of the strain tensor
		Matrix9d hessTrE() { return this->m_hessTrE; };

    //! Get the restriction vector
    Vector3d restrictionVector() { return this->m_restrictVector; };

    //! Get whether the face is on the boundary
    bool isBoundary() { return this->m_isBoundary; };

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

};

//! Calculate the face area and unit normal vector
template <class Refs>
void ElasticFace<Refs>::calculateFaceAreaAndNormal() {

	// NOTE: THE VERTEX COORDINATES SHOULD ALREADE BE UP TO DATE!

	Halfedge_handle h1 = this->halfedge();

	Vector3d v1 = h1->vertex()->v() - h1->prev()->vertex()->v();
	Vector3d v2 = h1->next()->vertex()->v() - h1->vertex()->v();

	Vector3d fN = v1.cross(v2);

  this->m_faceAreaWeightedNormal = fN;
	this->m_faceArea =  fN.norm() / 2.0;
	this->m_faceNormal = fN / fN.norm();

};


#endif
