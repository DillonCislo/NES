/*!
 * 	\file ElasticVertex.h
 * 	\brief Extension of CGAL base vertex class
 *
 * 	This subclass extends the CGAL base vertex class to include functionalities
 * 	needed to calculate the elastic energy/gradients of a Discrete Non-Euclidean
 * 	Koiter Surface
 *
 * 	\author Dillon Cislo
 * 	\date 12/20/2018
 *
 */

#ifndef _ELASTIC_VERTEX_H_
#define _ELASTIC_VERTEX_H_

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/HalfedgeDS_vertex_max_base_with_id.h>

typedef CGAL::Simple_cartesian<double> 		Kernel;
typedef Kernel::Point_3 			Point;

//! The extension of the base CGAL vertex class to the elastic mesh
template <class Refs>
struct ElasticVertex : public CGAL::HalfedgeDS_vertex_max_base_with_id<Refs, Point, std::size_t> {

	public:

		typedef typename CGAL::Simple_cartesian<double>	Kernel;
		typedef typename Kernel::Point_3 		Point;

		typedef typename Eigen::Vector3d 		Vector3d;

	public:

		//! Default constructor
		ElasticVertex() :
		       	CGAL::HalfedgeDS_vertex_max_base_with_id<Refs, Point, std::size_t>() {};

		//! Initialized with point p
		ElasticVertex( const Point &p ) :
			CGAL::HalfedgeDS_vertex_max_base_with_id<Refs, Point, std::size_t>( p ) {

			Vector3d v(3);
			v << p[0], p[1], p[2];
			this->m_v = v;

		};

		//! Initialized with point p and ID
		ElasticVertex( const Point &p, std::size_t i ) :
			CGAL::HalfedgeDS_vertex_max_base_with_id<Refs, Point, std::size_t>( p, i ) {

			Vector3d v(3);
			v << p[0], p[1], p[2];
			this->m_v = v;

		};

	protected:

		//! The vertex coordinates
		Vector3d m_v = Vector3d::Zero();

		//! The target vertex coordinates
		Vector3d  m_tarV = Vector3d::Zero();

		//! Will be true if this vertex has a user supplied target location
		bool m_targetCheck = false;

    //! The (face-averaged) vertex unit normal
    Vector3d m_vertexNormal = Vector3d::Zero();

    //! The mean curvature unit normal
    Vector3d m_Hn = Vector3d::Zero();

    //! The absolute mean curvature
    double m_absH = 0.0;

    // The sign of the mean curvature
    double m_sgnH = 0.0;

    // The signed mean curvature
    double m_H = 0.0;

    // The target (circumcentric) vertex area
    double m_targetArea = 0.0;

    // The sum of all of the target internal angles incident on the vertex
    double m_targetAngleSum = 0.0;

    // The sum of all the Laplace-Beltrami edge weights around the vertex
    double m_laplaceBeltramiSum = 0.0;

	public:

    /************************************************************************************
     * SETTERS
     ***********************************************************************************/

		//! Set the vertex coordinates
		void setV( const Vector3d &v ) { this->m_v = v; };

		//! Set the target vertex coordinates
		void setTarV( const Vector3d &tarV ) { this->m_tarV = tarV; };

		//! Set the target check
		void setTarget( bool isTarget ) { this->m_targetCheck = isTarget; };

    //! Set the (face-averaged) vertex unit normal
    void setVertexNormal( const Vector3d &vn ) { this->m_vertexNormal = vn; };

    //! Set the mean curvature unit normal
    void setMeanCurvatureNormal( const Vector3d &Hn ) { this->m_Hn = Hn; };

    //! Set the absolute mean curvature
    void setAbsH( double absH ) { this->m_absH = absH; };

    //! Set the sign of the mean curvature
    void setSgnH( double sgnH ) { this->m_sgnH = sgnH; };

    //! Set the signed mean curvature
    void setH( double H ) { this->m_H = H; };

    //! Set the target (circumcentric) vertex area
    void setTargetVertexArea( double vA ) { this->m_targetArea = vA; };

    //! Set the sum of all of the target internal angles incident on the vertex
    void setTargetAngleSum( double angSum ) { this->m_targetAngleSum = angSum; };

    //! Set the sum of the target Laplace-Beltrami edge weights incident on the vertex
    void setLaplaceBeltramiSum( double LBWSum ) { this->m_laplaceBeltramiSum = LBWSum; };

    /************************************************************************************
     * GETTERS
     ***********************************************************************************/

		//! Get the vertex coordinates
		Vector3d v() { return this->m_v; };

		//! Get the target vertex coordinates
		Vector3d tarV() { return this->m_tarV; };

		//! Get the target check
		bool isTarget() { return this->m_targetCheck; };

    //! Get the (face-averaged) vertex unit normal
    Vector3d normal() { return this->m_vertexNormal; };

    //! Get the mean curvature unit normal
    Vector3d Hn() { return this->m_Hn; };

    //! Get the absolute mean curvature
    double absH() { return this->m_absH; };

    //! Get the sign of the mean curvature
    double sgnH() { return this->m_sgnH; };

    //! Get the signed mean curvature
    double H() { return this->m_H; };

    //! Get the target (circumcentric) vertex area
    double targetArea() { return this->m_targetArea; };

    //! Get the sum of all of the target internal angles incident on the vertex
    double targetAngleSum() { return this->m_targetAngleSum; };

    //! Get the sum of the target Laplace-Beltrami edge weighs incident on the vertex
    double laplaceBeltramiSum() { return this->m_laplaceBeltramiSum; };

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

};

#endif
