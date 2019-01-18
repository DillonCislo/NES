/*!
 * 	\file PolyhedronSoupToPolyhedronMesh.h
 * 	\brief Builds a Polyhedron_3 object from a 'polyhedron soup' structure
 *
 * 	\author Dillon Cislo
 * 	\date 01/11/2018
 *
 */

#ifndef _POLYHEDRON_SOUP_TO_POLYHEDRON_MESH_H_
#define _POLYHEDRON_SOUP_TO_POLYHEDRON_MESH_H_

#include <CGAL/basic.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

#include <cstddef>
#include <iostream>

#include <Eigen/Core>

///
/// The class containing the methods used to scan in the mesh geometry/topology from the user
/// supplied 'polyhedron soup' structure
///
template <class HDS, class Scalar>
class PolyhedronSoupToPolyhedronMesh : public CGAL::Modifier_base<HDS> {

	public:

		typedef HDS 				Halfedge_data_structure;

		typedef typename HDS::Traits 		Traits;
		typedef typename Traits::Point_3 	Point;

		typedef typename Eigen::MatrixXi 					MatrixXi;
		typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> 	MatrixXx;


	protected:
		
		///
		/// The list of 3-D vertex coordinates with size (Nv x 3)
		///
		const MatrixXx &m_points;

		///
		/// The N-vertex face connectivity list with size (Nf x N)
		/// For example a triangulation has size (Nf x 3)
		///
		const MatrixXi &m_faces;

	public:

		///
		/// The constructor for the modifer object
		/// \param points the vertex coordinates of the polygon soup
		/// \param faces the face connectivity list of the polygon soup
		///
		PolyhedronSoupToPolyhedronMesh( const MatrixXx &points, const MatrixXi &faces ) :
       				m_points(points), m_faces(faces) {};

		///
		/// Build the polyhedral mesh
		///
		void operator()( HDS &target );

};

///
/// Build the polyhedral mesh
///
template <class HDS, class Scalar>
void PolyhedronSoupToPolyhedronMesh<HDS, Scalar>::operator()( HDS &target ) {

	// Postcondition: HDS is a valid polyhedral surface
	CGAL::Polyhedron_incremental_builder_3<HDS> B( target, true );
	B.begin_surface( m_points.rows(), m_faces.rows(), 0 );

	// Read in all vertices
	for( std::size_t i = 0; i < m_points.rows(); i++ ) {

		Eigen::RowVector3d vec = m_points.row(i);
		Point newPnt( vec(0), vec(1), vec(2) );
		B.add_vertex( newPnt );
	
	}

	if( B.error() ) {
		B.rollback();
		return;
	}

	// Read in all facets
	for ( std::size_t i = 0; i < m_faces.rows(); i++ ) {

		B.begin_facet();

		std::size_t no = m_faces.cols();
		if( B.error() || no < 3 ) {
			std::cerr << " " << std::endl;
			std::cerr << "PolyhedronBuildFromArray<HDS>::" << std::endl;
			std::cerr << "buildFromEigenArray: input error: facet" << i
				<< " has less than 3 vertices." << std::endl;
	
			B.rollback();
			return;
		}

		for( std::size_t j = 0; j < no; j++ ) {

			B.add_vertex_to_facet( m_faces(i,j) );

		}

		B.end_facet();
	}

	if( B.error() ) {
		B.rollback();
		return;
	}
	
	// Remove unconnected vertices
	if( B.check_unconnected_vertices() ) {
		if( !B.remove_unconnected_vertices() ) {
			std::cerr << " " << std::endl;
			std::cerr << "PolyhedronBuildFromArray<HDS>::" << std::endl;
			std::cerr << "buildFromEigenArray : input error: cannot "
				"successfully remove isolated vertices."
				<< std::endl;

			B.rollback();
			return;

		}
	}

	B.end_surface();

};

#endif
