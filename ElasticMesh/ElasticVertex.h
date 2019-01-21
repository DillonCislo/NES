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

	public:

		//! Set the vertex coordinates
		void setV( const Vector3d &v ) { this->m_v = v; };

		//! Set the target vertex coordinates
		void setTarV( const Vector3d &tarV ) { this->m_tarV = tarV; };

		//! Set the target check
		void setTarget( bool isTarget ) { this->m_targetCheck = isTarget; };

		//! Get the vertex coordinates
		Vector3d v() { return this->m_v; };

		//! Get the target vertex coordinates
		Vector3d tarV() { return this->m_tarV; };

		//! Get the target check
		bool isTarget() { return this->m_targetCheck; };

};

#endif
