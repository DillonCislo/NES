/*!
 * 	\file ElasticUpdater.h
 * 	\brief A functor to update the elastic mesh during equilibration
 *
 * 	\author Dillon Cislo
 * 	\date 12/28/2018
 *
 */

#ifndef _ELASTIC_UPDATER_H_
#define _ELASTIC_UPDATER_H_

#include <Eigen/Core>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include "ElasticItems.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <cassert>

/* INCLUDED IN "RestrictGrowthOperator.h"
// Branchless, type-safe sign function (used to calculate hinge angle)
template <class T>
inline int sgn( T val ) {
	return ( T(0) < val ) - ( val < T(0) );
};
*/

///
/// Class to update the elements of the elastic mesh before and during equilibration
///
class ElasticUpdater {

	private:

		typedef typename CGAL::Simple_cartesian<double> 		Kernel;
		typedef typename CGAL::Polyhedron_3<Kernel, ElasticItems> 	Polyhedron;

		typedef typename Polyhedron::Vertex_iterator 		Vertex_iterator;
		typedef typename Polyhedron::Facet_iterator 		Facet_iterator;
		typedef typename Polyhedron::Edge_iterator		Edge_iterator;
		typedef typename Polyhedron::Halfedge_iterator 		Halfedge_iterator;
		typedef typename Polyhedron::Vertex_handle 		Vertex_handle;

		typedef typename Polyhedron::Halfedge_around_vertex_circulator 	HV_circulator;
		typedef typename Polyhedron::Halfedge_around_facet_circulator 	HF_circulator;

		typedef typename Eigen::Vector3d 			Vector3d;
		typedef typename Eigen::VectorXd 			VectorXd;
		typedef typename Eigen::VectorXi			VectorXi;
		typedef typename Eigen::MatrixXd 			MatrixXd;
		typedef typename Eigen::Matrix3d			Matrix3d;
    typedef typename Eigen::MatrixXi      MatrixXi;
    typedef typename Eigen::Matrix<double,3,9>    RowGradVec;
		typedef typename Eigen::Matrix<double,9,9> 		Matrix9d;

	protected:

		///
		/// A vector containing the constant edge vector gradient matrices
		///
		std::vector<RowGradVec, Eigen::aligned_allocator<RowGradVec> > m_gradE;

		///
		/// A vector containing the constant Hessian matrices of the squared
		/// edge lengths
		///
		std::vector<Matrix9d, Eigen::aligned_allocator<Matrix9d> > m_hessL2;

	public:

		///
		/// Default constructor
		///
		ElasticUpdater();

		// /////////////////////////////////////////////////////////////////////////////
		// ASSIGN TARGET VERTEX INFORMATION
		// /////////////////////////////////////////////////////////////////////////////

		//! Assign vertex IDs
		void assignVertexID( Polyhedron &P );

		//! Assign target vertex information
		std::vector<Vertex_handle> assignTargetVertex( Polyhedron &P,
				const VectorXi &tarV_ID, const MatrixXd &tarV );

    // /////////////////////////////////////////////////////////////////////////////
    // SET BOUNDARY TAGS
    // /////////////////////////////////////////////////////////////////////////////

    //! Assign boundary tags to mesh faces
    void assignFacetBoundaryTags( Polyhedron &P );

    // /////////////////////////////////////////////////////////////////////////////
    // ASSIGN MATERIAL PARAMETERS
    // /////////////////////////////////////////////////////////////////////////////

    //! Assign material parameters
    void assignMaterialParameters( Polyhedron &P,
        const VectorXd &thickness );

		// /////////////////////////////////////////////////////////////////////////////
		// ASSIGN MAJOR EDGES
		// /////////////////////////////////////////////////////////////////////////////

		//! Assign major edges
		void assignMajorEdges( Polyhedron &P );

		// /////////////////////////////////////////////////////////////////////////////
		// UPDATE TARGET GEOMETRY
		// /////////////////////////////////////////////////////////////////////////////

		//! Complete assignment of target geometry and constants
		void updateTargetGeometry( Polyhedron &P,
				const VectorXi &v1_ID, const VectorXi &v2_ID,
				const VectorXd &tarL, const VectorXd &tarTheta );

		//! Assign target edge lengths and hinge functions
		void assignTargetLengthsAndAngles( Polyhedron &P,
				const VectorXi &v1_ID, const VectorXi &v2_ID,
				const VectorXd &tarL, const VectorXd &tarTheta );

		// /////////////////////////////////////////////////////////////////////////////
		// UPDATE CURRENT GEOMETRY
		// /////////////////////////////////////////////////////////////////////////////

		//! Complete update of current geometry from a vector of vertex coordinates
		void updateCurrentGeometry( Polyhedron &P, const VectorXd &x );

		//! Update the current vertex positions
		void updateCurrentVertices( Polyhedron &P, const VectorXd &x );

		//! Update the current faces
		void updateCurrentFaces( Polyhedron &P );

		//! Update current edges
		void updateCurrentEdges( Polyhedron &P );

    // /////////////////////////////////////////////////////////////////////////////
    // SET INITIAL GEOMETRY
    // /////////////////////////////////////////////////////////////////////////////

    //! Set the initial geometry
    void setInitialGeometry( Polyhedron &P );

    // /////////////////////////////////////////////////////////////////////////////
    // SET RESTRICTION GEOMETRY
    // /////////////////////////////////////////////////////////////////////////////

    //! Set the restriction geometry
    void setRestrictionGeometry( Polyhedron &P, const MatrixXd &restrictVectors,
        const MatrixXd &restrictedLengths );

  public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

};

///
/// Default constructor
///
ElasticUpdater::ElasticUpdater() {

	Matrix3d Z3 = Matrix3d::Zero();
	Matrix3d I3 = Matrix3d::Identity();

	RowGradVec gradEI;
	gradEI << Z3, -I3, I3;

	RowGradVec gradEJ;
	gradEJ << I3, Z3, -I3;

	RowGradVec gradEK;
	gradEK << -I3, I3, Z3;

	std::vector<RowGradVec,
    Eigen::aligned_allocator<RowGradVec> >
      gradE{ gradEI, gradEJ, gradEK };

	this->m_gradE = gradE;

	RowGradVec Z9 = RowGradVec::Zero();

	Matrix9d hessL2I;
	hessL2I << Z9, -gradEI, gradEI;

	Matrix9d hessL2J;
	hessL2J << gradEJ, Z9, -gradEJ;

	Matrix9d hessL2K;
	hessL2K << -gradEK, gradEK, Z9;

	std::vector<Matrix9d,
    Eigen::aligned_allocator<Matrix9d> >
      hessL2{ hessL2I, hessL2J, hessL2K };

	this->m_hessL2 = hessL2;

};

// /////////////////////////////////////////////////////////////////////////////////////////////
// SET BOUNDARY TAGS
// /////////////////////////////////////////////////////////////////////////////////////////////

//! Assign boundary tags to mesh faces
void ElasticUpdater::assignFacetBoundaryTags( Polyhedron &P ) {

  // A facet is considered a boundary facet if ANY of its vertices
  // lie on the mesh boundary
  Facet_iterator f;
  for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

    bool onBdy = false;

    // Determine if the first vertex is on the boundary
    HV_circulator he = f->halfedge()->next()->vertex()->vertex_begin();
    do {

      onBdy = onBdy || he->is_border_edge();
      he++;

    } while ( he != f->halfedge()->next()->vertex()->vertex_begin() );

    // Determine if the second vertex is on the boundary
    he = f->halfedge()->prev()->vertex()->vertex_begin();
    do {

      onBdy = onBdy || he->is_border_edge();
      he++;

    } while ( he != f->halfedge()->prev()->vertex()->vertex_begin() );

    // Determine if the third vertex is on the boundary
    he = f->halfedge()->vertex()->vertex_begin();
    do {

      onBdy = onBdy || he->is_border_edge();
      he++;

    } while ( he != f->halfedge()->vertex()->vertex_begin() );

    f->setBoundaryTag( onBdy );

  }

};


// /////////////////////////////////////////////////////////////////////////////////////////////
// ASSIGN TARGET VERTEX INFORMATION
// /////////////////////////////////////////////////////////////////////////////////////////////

//! Assign vertex IDs
void ElasticUpdater::assignVertexID( Polyhedron &P ) {

	int index = 0;
	Vertex_iterator v;
	for( v = P.vertices_begin();  v != P.vertices_end(); v++ ) {
		v->id() = index;
		index++;
	}

};


//! Assign target vertex information
std::vector< CGAL::Polyhedron_3<CGAL::Simple_cartesian<double>, ElasticItems>::Vertex_handle >
ElasticUpdater::assignTargetVertex( Polyhedron &P,
		const VectorXi &tarV_ID, const MatrixXd &tarV ) {

	// NOTE: VERTEX IDs SHOULD ALREADY BE ASSIGNED
	// Validate input argmuents
	if( tarV_ID.size() != tarV.rows() ) {
		throw std::invalid_argument("Target vertex list is improperly sized");
	}

	// Assign target vertices and target locations
	std::vector<Vertex_handle> targetVertices;
	targetVertices.reserve( tarV_ID.size() );

	Vertex_iterator v;
	for( int i = 0; i < tarV_ID.size(); i++ ) {

		v = P.vertices_begin();
		bool vtxSet = false;
		while( ( v != P.vertices_end() ) && !vtxSet ) {

			if( v->id() == tarV_ID(i) ) {

				v->setTarV( tarV.row(i).transpose() );
				v->setTarget( true );
				targetVertices.push_back( v );
				vtxSet = true;

			}

			v++;

		}

		if( !vtxSet )
			throw std::runtime_error("Target vertex not found");

	}

	return targetVertices;

};

// /////////////////////////////////////////////////////////////////////////////////////////////
// ASSIGN MAJOR EDGES
// /////////////////////////////////////////////////////////////////////////////////////////////

//! Assign major edges
void ElasticUpdater::assignMajorEdges( Polyhedron &P ) {

	Edge_iterator e;
	for( e = P.edges_begin(); e != P.edges_end(); e++ ) { e->setMajor( true ); }

};

// /////////////////////////////////////////////////////////////////////////////////////////////
// ASSIGN MATERIAL PARAMETERS
// /////////////////////////////////////////////////////////////////////////////////////////////

//! Assign material parameters
void ElasticUpdater::assignMaterialParameters( Polyhedron &P,
    const VectorXd &thickness ) {

  int fID = 0;
  Facet_iterator f;
  for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

    f->setThickness( thickness(fID) );
    fID++;

  }

};

// /////////////////////////////////////////////////////////////////////////////////////////////
// UPDATE TARGET GEOMETRY
// /////////////////////////////////////////////////////////////////////////////////////////////

//! Assign target edge lengths and hinge functions
void ElasticUpdater::assignTargetLengthsAndAngles( Polyhedron &P,
		const VectorXi &v1_ID, const VectorXi &v2_ID,
		const VectorXd &tarL, const VectorXd &tarTheta ) {

	// NOTE: VERTEX IDs SHOULD ALREADY BE ASSIGNED
	Vertex_iterator v;
	HV_circulator he;

	bool edgeSet;
	bool allEdgesSet = true;

	double theta;
	double phiBar;

	for( int i = 0; i < v1_ID.size(); i++ ) {

		v = P.vertices_begin();
		edgeSet = false;
		while( ( v != P.vertices_end() ) && !edgeSet ) {

			if( v->id() == v1_ID(i) ) {

				he = v->vertex_begin();
				do {
					if( he->opposite()->vertex()->id() == v2_ID(i) ) {

						he->setTargetEdgeLength( tarL(i) );
						he->opposite()->setTargetEdgeLength( tarL(i) );

						theta = tarTheta(i);
						phiBar = 2.0 * std::tan( theta / 2.0 );
						he->setTargetHingeFunction( phiBar );
						he->opposite()->setTargetHingeFunction( phiBar );

						edgeSet = true;

					}

					he++;

				} while( ( he != v->vertex_begin() ) && !edgeSet );

			}

			v++;

		}

		if (!edgeSet) allEdgesSet = false;

	}

	if (!allEdgesSet)
    throw std::runtime_error("Not all edges properly set.");

};

void ElasticUpdater::updateTargetGeometry( Polyhedron &P,
		const VectorXi &v1_ID, const VectorXi &v2_ID,
		const VectorXd &tarL, const VectorXd &tarTheta ) {

	// NOTE: VERTEX IDs SHOULD ALREADY BE ASSIGNED
	this->assignTargetLengthsAndAngles( P, v1_ID, v2_ID, tarL, tarTheta );

  // Set face-based target quantities and target internal angles
	Facet_iterator f;
	for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

		double Li = f->halfedge()->targetLength();
		double Lj = f->halfedge()->next()->targetLength();
		double Lk = f->halfedge()->prev()->targetLength();
		Vector3d L;
		L << Li, Lj, Lk;

    // Calculate target face areas
		double s = ( Li + Lj + Lk ) / 2.0;
		double tarArea = std::sqrt( s * (s-Li) * (s-Lj) * (s-Lk) );
		f->setTargetFaceArea( tarArea );

    // Calculate target internal angles
    double thetaI = std::acos( ( Lk*Lk + Lj*Lj - Li*Li ) / ( 2.0 * Lj * Lk ) );
    double thetaJ = std::acos( ( Li*Li + Lk*Lk - Lj*Lj ) / ( 2.0 * Lk * Li ) );
    double thetaK = std::acos( ( Lj*Lj + Li*Li - Lk*Lk ) / ( 2.0 * Li * Lj ) );
    f->halfedge()->setTargetAngle( thetaI );
    f->halfedge()->next()->setTargetAngle( thetaJ );
    f->halfedge()->prev()->setTargetAngle( thetaK );

		Vector3d C = Vector3d::Zero();		    // Stretching constants Ci
		Vector3d hBar = Vector3d::Zero();	    // Target edge heights
		Matrix3d zeta = Matrix3d::Zero();	    // Stretching constant matrix
		Matrix3d xi = Matrix3d::Zero(); 	    // Bending constant matrix
		Matrix9d hessTrE = Matrix9d::Zero(); 	// The Hessian of the trace of the strain
							                            // tensor

		for( int i = 0; i < 3; i++ ) {

			// Stretching constants Ci
			C(i) = -( L[i]*L[i] - L[(i+1)%3]*L[(i+1)%3] - L[(i+2)%3]*L[(i+2)%3] );
			C(i) = C(i) / ( 8.0 * tarArea * tarArea );

			// Target edge heights
			hBar(i) = 2.0 * tarArea / L[i];

			// Diagonal elements of zeta matrix
			zeta(i,i) = C(i)*C(i) + 1 / ( 8.0 * tarArea * tarArea );

			// Diagonal elementa of xi matrix
			xi(i,i) = 1 / ( hBar(i) * hBar(i) );

			// The Hessian of the trace of the strain tensor
			hessTrE += C(i) * this->m_hessL2[i];

		}

		for( int i = 0; i < 3; i++ ) {
			for( int j = (i+1); j < 3; j++ ) {

				int k = ( 2 * ( i + j ) ) % 3;

				// Off-diagonal elements of the zeta matrix
				zeta(i,j) = 1.0 -2.0 * L[k]*L[k] * C(k);
				zeta(i,j) /= 4.0 * tarArea * tarArea;

				// Off-diagonal elements of the xi matrix
				xi(i,j) = L[k]*L[k] - L[i]*L[i] - L[j]*L[j];
				xi(i,j) = xi(i,j) * xi(i,j);
				xi(i,j) = xi(i,j) / ( 8.0 * tarArea * tarArea * L[i] * L[j] );

			}
		}

		f->setC( C );
		f->setHBar( hBar );
		f->setZeta( zeta );
		f->setXi( xi );
		f->setHessTrE( hessTrE );

	}

  // Set Laplace-Beltrami edge weights
  Edge_iterator e;
	for( e = P.edges_begin(); e != P.edges_end(); e++ ) {

    double LBW = 0.0;

    if ( !e->is_border() )
      LBW += std::cos(e->targetAngle()) / std::sin(e->targetAngle());

    if ( !e->opposite()->is_border() ) {

      LBW += std::cos(e->opposite()->targetAngle()) /
        std::sin(e->opposite()->targetAngle());

    }

    LBW = LBW / 2.0;

    e->setLaplaceBeltramiWeight( LBW );
    e->opposite()->setLaplaceBeltramiWeight( LBW );

  }

  // Set the vertex target areas, angle sums, and Laplace-Beltrami edge weight sums
  Vertex_iterator v;
  for( v = P.vertices_begin(); v != P.vertices_end(); v++ ) {

    double angSum = 0.0;
    double AV = 0.0;
    double LBWSum = 0.0;

    HV_circulator he = v->vertex_begin();
    do {

      if ( !he->is_border() ) {

        // Accumulate the incident angle
        angSum += he->prev()->targetAngle();

        // Accumulate the incident Laplace-Beltrami edge weight
        LBWSum += he->laplaceBeltramiWeight();

        double AF = he->face()->tarFaceArea();

        double Li2 = he->prev()->targetLength(); Li2 *= Li2;
        double Lj2 = he->targetLength(); Lj2 *= Lj2;
        double Lk2 = he->next()->targetLength(); Lk2 *= Lk2;

        // double cotI = ( Lj2 + Lk2 - Li2 ) / ( 4.0 * AF );
        double cotJ = ( Lk2 + Li2 - Lj2 ) / ( 4.0 * AF );
        double cotK = ( Li2 + Lj2 - Lk2 ) / ( 4.0 * AF );

        // Accumulate the contribution to the vertex area
        AV += Lj2 * cotJ + Lk2 * cotK;

      }

      he++;

    } while ( he != v->vertex_begin() );

    v->setTargetVertexArea( AV / 8.0 );
    v->setTargetAngleSum( angSum );
    v->setLaplaceBeltramiSum( LBWSum );

  }

};

// /////////////////////////////////////////////////////////////////////////////////////////////
// UPDATE CURRENT GEOMETRY
// /////////////////////////////////////////////////////////////////////////////////////////////


//! Update the current vertex positions
void ElasticUpdater::updateCurrentVertices( Polyhedron &P, const VectorXd &x ) {

	int Nv = P.size_of_vertices();

	if ( x.size() != (3*Nv) )
		throw std::runtime_error( "Vertex vector is improperly sized!" );

	int i = 0;
	Vertex_iterator v;
	for( v = P.vertices_begin(); v != P.vertices_end(); v++ ) {

		Vector3d newPnt;
		newPnt << x(i), x(i+Nv), x(i+(2*Nv));

		v->setV( newPnt );

		i++;

	}

};

//! Update the current faces
void ElasticUpdater::updateCurrentFaces( Polyhedron &P ) {

	// NOTE: CURRENT VERTEX POSITIONS SHOULD BE UP TO DATE
	Facet_iterator f;
	for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

		f->calculateFaceAreaAndNormal();

	}

};


//! Update the current edges
void ElasticUpdater::updateCurrentEdges( Polyhedron &P ) {

	// NOTE:: CURRENT VERTEX POSITIONS AND FACE QUANTITIES SHOULD BE UP TO DATE
	Halfedge_iterator he;
	for( he = P.halfedges_begin(); he != P.halfedges_end(); he++ ) {

		he->calculateEdgeQuantities();

	}

	Edge_iterator e;
	for( e = P.edges_begin(); e != P.edges_end(); e++ ) {

		// Skip boundary edges
		if ( e->is_border() || e->opposite()->is_border() ) continue;

		Vector3d n1 = e->face()->faceNormal();
		Vector3d n2 = e->opposite()->face()->faceNormal();
		Vector3d e0Hat = e->edgeUnitVector();

		Vector3d z = n1.cross(n2);
		double phi = 2.0 * z.dot(e0Hat) / ( 1.0 + n1.dot(n2) );

		e->setPhi( phi );
		e->opposite()->setPhi( phi );

	}

};

//! Complete update of current geometry from a vector of vertex coordinates
void ElasticUpdater::updateCurrentGeometry( Polyhedron &P, const VectorXd &x ) {

	this->updateCurrentVertices( P, x );
	this->updateCurrentFaces( P );
	this->updateCurrentEdges( P );

};

// /////////////////////////////////////////////////////////////////////////////////////////////
// SET INITIAL GEOMETRY
// /////////////////////////////////////////////////////////////////////////////////////////////

//! Set the initial geometry
void ElasticUpdater::setInitialGeometry( Polyhedron &P ) {

  // NOTE: CURRENT AND RESTRICTION GEOMETRIES SHOULD BOTH ALREADY BE UP TO DATE
  // This function just copies the current geometry into
  // the initial geometry parameters

  // Set initial halfedge growth scalars
	Facet_iterator f;
	for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

    double dblA0 = 2.0 * f->faceArea();
    Vector3d n0 = f->faceNormal();
    Vector3d v0 = f->restrictionVector();

    Vector3d ei0 = f->halfedge()->edgeVector();
    Vector3d ej0 = f->halfedge()->next()->edgeVector();
    Vector3d ek0 = f->halfedge()->prev()->edgeVector();

    double gi = (v0.dot(n0.cross(ei0))) / dblA0;
    double gj = (v0.dot(n0.cross(ej0))) / dblA0;
    double gk = (v0.dot(n0.cross(ek0))) / dblA0;

    f->halfedge()->setInitialEdgeGrowth( gi );
    f->halfedge()->next()->setInitialEdgeGrowth( gj );
    f->halfedge()->prev()->setInitialEdgeGrowth( gk );

	}


};


// /////////////////////////////////////////////////////////////////////////////////////////////
// SET RESTRICTION GEOMETRY
// /////////////////////////////////////////////////////////////////////////////////////////////

//! Set the initial restriction geometry
void ElasticUpdater::setRestrictionGeometry( Polyhedron &P,
    const MatrixXd &restrictVectors, const MatrixXd &restrictedLengths ) {

  int i = 0;
  Facet_iterator f;
  for( f = P.facets_begin(); f != P.facets_end(); f++ ) {

    f->setRestrictionVector( restrictVectors.row(i).transpose() );

    // By default, the incident halfedge to the face points at the first
    // vertex in that face. The input restricted edge length list is
    // supplied such that the entry corresponding to a particular vertex
    // sets the length for the edge OPPOSITE that vertex, rather than
    // incident to it. The correspondence is:
    // [ 0, 1, 2 ] --> [ 1, 2, 0 ]

    // Incident to v0 - Opposite to v1
    f->halfedge()->setRestrictedLength( restrictedLengths(i,1) );

    // Incident to v1 - Opposite to v2
    f->halfedge()->next()->setRestrictedLength( restrictedLengths(i,2) );

    // Incident to v2 - Opposite to v0
    f->halfedge()->prev()->setRestrictedLength( restrictedLengths(i,0) );

    i++;

  }

};

#endif
