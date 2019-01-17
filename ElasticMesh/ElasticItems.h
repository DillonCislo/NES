/*!
 * 	\file ElasticItems.h
 * 	\brief Extention of CGAL Polyhedron_3 Items class
 *
 * 	This class extends the CGAL Polyhedron_3 Items class to include the geometric objects
 * 	( ElasticVertex, ElasticFace, ElasticHalfEdge ) needed to calculated the elastic
 * 	energy/gradients of a Discrete Non-Euclidean Koiter Surface
 *
 * 	\author Dillon Cislo
 * 	\date 11/01/2018
 *
 */

#ifndef _ELASTIC_ITEMS_H_
#define _ELASTIC_ITEMS_H_

#include <CGAL/Polyhedron_items_with_id_3.h>
#include "ElasticVertex.h"
#include "ElasticFace.h"
#include "ElasticHalfEdge.h"

//! The extension of the base CGAL items class to the elastic mesh
struct ElasticItems : public CGAL::Polyhedron_items_with_id_3 {

	template <class Refs, class Traits>
	struct Vertex_wrapper { typedef ElasticVertex<Refs> Vertex; };

	template <class Refs, class Traits>
	struct Face_wrapper { typedef ElasticFace<Refs> Face; };

	template <class Refs, class Traits>
	struct Halfedge_wrapper { typedef ElasticHalfEdge<Refs> Halfedge; };

};

#endif
