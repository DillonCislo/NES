/* =============================================================================================
 *
 * 	calculate_edge_angles.cpp
 *
 * 	Calculates the bending angles for each edge of a triangular mesh
 *
 * 	by Dillon Cislo
 * 	05/02/2021
 *
 * 	This is a MEX-file for MATLAB
 *
 * ============================================================================================*/

#include "mex.h" // for MATLAB

#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cassert>

#include <Eigen/Core>

#include "../../ElasticOptimization/ElasticProblem.h"
#include "../../ElasticMesh/PolyhedronSoupToPolyhedronMesh.h"

typedef CGAL::Simple_cartesian<double>	 		      Kernel;
typedef CGAL::Polyhedron_3<Kernel, ElasticItems> 	Polyhedron;

typedef Polyhedron::Vertex_handle 		Vertex_handle;
typedef Polyhedron::Vertex_iterator 	Vertex_iterator;
typedef Polyhedron::Facet_iterator    Facet_iterator;
typedef Polyhedron::Edge_iterator     Edge_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::HalfedgeDS 				HalfedgeDS;

typedef Eigen::Vector3d     Vector3d;
typedef Eigen::VectorXd 		VectorXd;
typedef Eigen::MatrixXd 		MatrixXd;
typedef Eigen::VectorXi 		VectorXi;
typedef Eigen::MatrixXi 		MatrixXi;

// Brief main function to call computational functionalities
void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] ) {

	// -------------------------------------------------------------------------------------
	// INPUT PROCESSING
	// -------------------------------------------------------------------------------------
	
	// VARIABLE DECLARATIONS ---------------------------------------------------------------
	
	MatrixXi faces; 			      // The face connectivity list

	MatrixXd vertex; 		        // The initial vertex positions

	VectorXi v1_ID;	 		        // The start vertex ID of each edge

	VectorXi v2_ID; 			      // The end vertex ID of each edge

	int Nf = 0; 	// Number of faces
	int Nv = 0; 	// Number of vertices
	int Ne = 0;	  // Number of edges

	// EXTRACT INPUT PARAMETERS ------------------------------------------------------------

	// Check for proper number of arguments
	if (nrhs != 4 ) {
		mexErrMsgIdAndTxt("MATLAB:calculate_edge_angles:nargin",
				"CALCULATE_EDGE_ANGLES requires four input arguments.");
	} else if ( nlhs != 1 ) {
		mexErrMsgIdAndTxt("MATLAB:calculate_edge_angles:nargout",
				"MINIMIZE_EDGE_ANGLES requires one output argument.");
	}

	double *m_faces = mxGetPr( prhs[0] );
	Nf = (int) mxGetM( prhs[0] );
  MatrixXd Fd = Eigen::Map<MatrixXd>( m_faces, Nf, 3 );
  faces = Fd.cast <int> ();
  faces = (faces.array() - 1).matrix(); // Subtract 1 to account for MATLAB indexing

	double *m_vertex = mxGetPr( prhs[1] );
	Nv = (int) mxGetM( prhs[1] );
	vertex = Eigen::Map<MatrixXd>( m_vertex, Nv, 3 );

	double *m_v1_ID = mxGetPr( prhs[2] );
	Ne = (int) mxGetM( prhs[2] );
  VectorXd v1_IDd = Eigen::Map<VectorXd>( m_v1_ID, Ne );
	v1_ID = v1_IDd.cast <int> ();
  v1_ID = (v1_ID.array() - 1).matrix(); // Subtract 1 to account for MATLAB indexing

	double *m_v2_ID = mxGetPr( prhs[3] );
	VectorXd v2_IDd = Eigen::Map<VectorXd>( m_v2_ID, Ne );
  v2_ID = v2_IDd.cast <int> ();
  v2_ID = (v2_ID.array() - 1).matrix(); // Subtract 1 to account for MATLAB indexing

  // PREPARE POLYHEDRAL MESH -------------------------------------------------------------
  
  Polyhedron P;
  
  // Read in basic mesh topology and vertex locations
  PolyhedronSoupToPolyhedronMesh<HalfedgeDS, double> inPoly( vertex, faces );
  P.delegate( inPoly );

  assert( P.is_valid() && "Failed to properly construct polyhedral mesh!" );

  // Assign vertex IDs
  ElasticUpdater EU;
  EU.assignVertexID( P );

  // Update current geometry
  MatrixXd vertexCopy = vertex;
  VectorXd x = Eigen::Map<Eigen::VectorXd>( vertexCopy.data(), vertexCopy.size() );
  EU.updateCurrentGeometry( P, x );
  
  // -------------------------------------------------------------------------------------
  // CALCULATE EDGE ANGLES
  // -------------------------------------------------------------------------------------
  
  VectorXd edgeAngles = VectorXd::Zero( Ne );
  bool allEdgesSet = true;

  Edge_iterator e;
  for( e = P.edges_begin(); e != P.edges_end(); e++ ) {

    // Skip boundary edges
    if ( e->is_border() || e->opposite()->is_border() ) continue;

    Vector3d n1 = e->face()->faceNormal();
    Vector3d n2 = e->opposite()->face()->faceNormal();
    Vector3d e0Hat = e->edgeUnitVector();

    Vector3d z = n1.cross(n2);

    double ang = 2.0 * std::atan2( z.dot(e0Hat), (1.0 + n1.dot(n2)) );

    // Assign to the corresponding edge in the output list
    int v1 = e->vertex()->id();
    int v2 = e->opposite()->vertex()->id();

    bool edgeSet = false;
    int edgeCounter = 0;
    do {

      bool setCurEdge = (v1_ID(edgeCounter) == v1) && (v2_ID(edgeCounter) == v2);
      setCurEdge = setCurEdge ||
        ( (v1_ID(edgeCounter) == v2) && (v2_ID(edgeCounter) == v1) );

      if( setCurEdge ) {

        edgeAngles(edgeCounter) = ang;
        edgeSet = true;

      }

      edgeCounter++;

    } while( (edgeCounter < Ne) && !edgeSet );

    if (!edgeSet) allEdgesSet = false;

  }

  assert( allEdgesSet && "Not all edges angles properly calculated!" );

	// -------------------------------------------------------------------------------------
	// OUTPUT PROCESSING
	// -------------------------------------------------------------------------------------
  
  plhs[0] = mxCreateDoubleMatrix( Ne, 1, mxREAL );
  Eigen::Map<MatrixXd>( mxGetPr( plhs[0] ), Ne, 1 ) = edgeAngles;

  return;

};

    


