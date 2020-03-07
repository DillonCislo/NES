/* =============================================================================================
 *
 * 	minimize_elastic_energy.cpp
 *
 * 	Minimizes the elastic energy of a non-Euclidean shell in order to find the
 * 	equilibrium configuration associated with a given intrinsic geometry.
 *
 * 	by Dillon Cislo
 * 	01/16/2019
 *
 * 	This is a MEX-file for MATLAB
 *
 * ============================================================================================*/

#include "mex.h" // for MATLAB

#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>

#include <Eigen/Core>
#include <LBFGSpp/LBFGS.h>

#include "../ElasticOptimization/ElasticProblem.h"
#include "../ElasticMesh/PolyhedronSoupToPolyhedronMesh.h"

typedef CGAL::Simple_cartesian<double>	 		Kernel;
typedef CGAL::Polyhedron_3<Kernel, ElasticItems> 	Polyhedron;

typedef Polyhedron::Vertex_handle 			Vertex_handle;
typedef Polyhedron::Vertex_iterator 			Vertex_iterator;
typedef Polyhedron::HalfedgeDS 				HalfedgeDS;

typedef Eigen::VectorXd 		VectorXd;
typedef Eigen::MatrixXd 		MatrixXd;
typedef Eigen::VectorXi 		VectorXi;
typedef Eigen::MatrixXi 		MatrixXi;


// Prepare a polyhedral mesh object for the elastic minimization procedure
void preparePolyhedron( Polyhedron &P, const MatrixXi &faces, const MatrixXd &vertex,
    const VectorXi &v1_ID, const VectorXi &v2_ID,
		const VectorXd &targetLengths, const VectorXd &targetAngles,
    double alpha, const VectorXi &target_ID, const MatrixXd &targetLocations ) {


	// Read in basic mesh topology and vertex locations
	PolyhedronSoupToPolyhedronMesh<HalfedgeDS, double> inPoly( vertex, faces );
	P.delegate( inPoly );

	if ( !P.is_valid() ) {
		std::runtime_error( "Failed to properly construct polyhedral mesh." );
	}

	// Update target geometry
	ElasticUpdater EU;

	EU.assignVertexID( P );

	if ( alpha > 0.0 ) {
		std::vector<Vertex_handle> tV;
		tV = EU.assignTargetVertex( P, target_ID, targetLocations );
	}

	EU.assignMajorEdges( P );
	EU.updateTargetGeometry( P, v1_ID, v2_ID, targetLengths, targetAngles );

	// Update current geometry
	MatrixXd vertexCopy = vertex;
	VectorXd x = Eigen::Map<Eigen::VectorXd>( vertexCopy.data(), vertexCopy.size() );

	EU.updateCurrentGeometry( P, x );

	return;

};

// Run elastic minimization procedure
VectorXd minimizeElasticEnergy( const MatrixXi &faces, const MatrixXd &vertex,
		const VectorXi &v1_ID, const VectorXi &v2_ID,
		const VectorXd &targetLengths, const VectorXd &targetAngles,
		const LBFGSpp::LBFGSParam<double> param,
		double h, double nu, double alpha,
		const VectorXi &target_ID, const MatrixXd &targetLocations,
    double beta, double targetVolume ) {

	Polyhedron P;
	preparePolyhedron( P, faces, vertex,
			v1_ID, v2_ID, targetLengths, targetAngles,
			alpha, target_ID, targetLocations );

	// Initialize the problem structure
	ElasticProblem f(P, h, nu, beta, targetVolume, alpha, target_ID );

	// Initial guess
	MatrixXd vertexCopy = vertex;
	VectorXd x = Eigen::Map<VectorXd>( vertexCopy.data(), vertexCopy.size() );

	// Initialize solver
	LBFGSpp::LBFGSSolver<double> solver(param);

	// Run solver
	double fx;
	int niter = solver.minimize( f, x, fx );

	return x;

};



// Brief main function to call computational functionalities
void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] ) {

	// -------------------------------------------------------------------------------------
	// INPUT PROCESSING
	// -------------------------------------------------------------------------------------
	
	// VARIABLE DECLARATIONS ---------------------------------------------------------------
	
	int *m_faces; 			        // The face connectivity list
	MatrixXi faces;

	double *m_vertex; 		      // The initial vertex positions
	MatrixXd vertex;

	double *m_targetLengths; 	  // The list of target edge lengths
	VectorXd targetLengths;

	double *m_targetAngles; 	  // The list of target bend angles
	VectorXd targetAngles;
	
	int *m_v1_ID;	 		          // The start vertex ID of each edge
	VectorXi v1_ID;

	int *m_v2_ID; 			        // The end vertex ID of each edge
	VectorXi v2_ID;

	int *m_target_ID; 		      // The vertex IDs of those vertices with a target
	VectorXi target_ID;		      // location

	double *m_targetLocations;	// The target locations of those vertices
	MatrixXd targetLocations;

	double h; 			      // The thickness of the elastic sheet
	double nu; 			      // Poisson's ratio
	double alpha; 	      // The coefficient of target vertex correspondence
  double beta;          // The coefficient of the fixed folume energy
  double targetVolume;  // The target volume for the fixed volume energy

	// L-BFGS parameter objects
	LBFGSpp::LBFGSParam<double> param;

	int Nf = 0; 	// Number of faces
	int Nv = 0; 	// Number of vertices
	int Ne = 0;	  // Number of edges
	int Nt = 0; 	// Number of target vertices

	// EXTRACT INPUT PARAMETERS ------------------------------------------------------------

	// Check for proper number of arguments
	if (nrhs != 14 ) {
		mexErrMsgIdAndTxt("MATLAB:minimize_elastic_energy:nargin",
				"MINIMIZE_ELASTIC_ENERGY requires fourteen input arguments.");
	} else if ( nlhs != 1 ) {
		mexErrMsgIdAndTxt("MATLAB:minimize_elastic_energy:nargout",
				"MINIMIZE_ELASTIC_ENERGY requires one output argument.");
	}

	m_faces = (int*) mxGetData( prhs[0] );
	Nf = (int) mxGetM( prhs[0] );
	faces = Eigen::Map<MatrixXi>( m_faces, Nf, 3 );

	m_vertex = mxGetPr( prhs[1] );
	Nv = (int) mxGetM( prhs[1] );
	vertex = Eigen::Map<MatrixXd>( m_vertex, Nv, 3 );

	m_v1_ID = (int*) mxGetData( prhs[2] );
	Ne = (int) mxGetM( prhs[2] );
	v1_ID = Eigen::Map<VectorXi>( m_v1_ID, Ne );

	m_v2_ID = (int*) mxGetData( prhs[3] );
	v2_ID = Eigen::Map<VectorXi>( m_v2_ID, Ne );

	m_targetLengths = mxGetPr( prhs[4] );
	targetLengths = Eigen::Map<VectorXd>( m_targetLengths, Ne );

	m_targetAngles = mxGetPr( prhs[5] );
	targetAngles = Eigen::Map<VectorXd>( m_targetAngles, Ne );

	// SET MINIMIZATION PARAMETERS ---------------------------------------------------------
	
	int idx;

	// Number of LBFGS history vectors retained
	if ( (idx = mxGetFieldNumber( prhs[6], "m")) == -1 ) {
		mexErrMsgTxt("No LBFGS history field!");
	} else {
		param.m = (int) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
	}

	// Tolerance for convergence test
	if ( (idx = mxGetFieldNumber( prhs[6], "epsilon")) == -1 ) {
		mexErrMsgTxt("No epsilon field!");
	} else {
		param.epsilon = (double) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
	}

	// Distance for delta-based convergence test
	if ( (idx = mxGetFieldNumber( prhs[6], "past")) == -1 ) {
		mexErrMsgTxt("No past field!");
	} else {
		param.past = (int) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
	}

	// Delta for convergence test
	if ( (idx = mxGetFieldNumber( prhs[6], "delta")) == -1 ) {
		mexErrMsgTxt("No delta field!");
	} else {
		param.delta = (double) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
	}

	// Maximum number of iterations
	if ( (idx = mxGetFieldNumber( prhs[6], "max_iterations")) == -1 ) {
		mexErrMsgTxt("No max_iterations field!");
	} else {
		param.max_iterations = (int) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
	}

	// Maximum number of trials for the line search
	if ( (idx = mxGetFieldNumber( prhs[6], "max_linesearch")) == -1 ) {
		mexErrMsgTxt("No max_linesearch field!");
	} else {
		param.max_linesearch = (int) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
	}

	// Minimum step in the line search
	if ( (idx = mxGetFieldNumber( prhs[6], "min_step")) == -1 ) {
		mexErrMsgTxt("No min_step field!");
	} else {
		param.min_step = (double) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
	}

	// Maximum step in the line search
	if ( (idx = mxGetFieldNumber( prhs[6], "max_step")) == -1 ) {
		mexErrMsgTxt("No max step field!");
	} else {
		param.max_step = (double) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
	}

	// Accuracy of the line search routine
	if ( (idx = mxGetFieldNumber( prhs[6], "ftol")) == -1 ) {
		mexErrMsgTxt("No ftol field!");
	} else {
		param.ftol = (double) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
	}

	// A coefficient for the Wolfe condition
	if ( (idx = mxGetFieldNumber( prhs[6], "wolfe")) == -1 ) {
		mexErrMsgTxt("No wolfe field!");
	} else {
		param.wolfe = (double) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
	}

	// Iterative display option
	int tmpLS;
	if ( (idx = mxGetFieldNumber( prhs[6], "iterDisp")) == -1 ) {
		mexErrMsgTxt("No iterDisp field!");
	} else {
		tmpLS = (int) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
	}

	if (tmpLS == 1) {
		param.iterDisp = true;
	} else { param.iterDisp = false; };

	// Linesearch procedure
	if ( (idx = mxGetFieldNumber( prhs[6], "linesearch")) == -1 ) {
		mexErrMsgTxt("No linesearch field!");
	} else {
		tmpLS = (int) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
	}

	if (tmpLS == 1) {
		param.linesearch = LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
	} else if (tmpLS == 2) {
		param.linesearch = LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
	} else if (tmpLS == 3) {
		param.linesearch = LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
	} else {
		mexErrMsgTxt("Linesearch field invalid!");
	};

	// GET MATERIAL PARAMETERS -------------------------------------------------------------
	h = *mxGetPr( prhs[7] );
	nu = *mxGetPr( prhs[8] );

	// GET TARGET VERTEX CORRESPONDENCE PARAMETERS -----------------------------------------
	alpha = *mxGetPr( prhs[9] );

	m_target_ID = (int*) mxGetData( prhs[10] );
	Nt = (int) mxGetM( prhs[10] );
	target_ID = Eigen::Map<VectorXi>( m_target_ID, Nt );

	m_targetLocations = mxGetPr( prhs[11] );
	targetLocations = Eigen::Map<MatrixXd>( m_targetLocations, Nt, 3 );

  // GET FIXED VOLUME PARAMETERS ----------------------------------------------------------
  beta = *mxGetPr( prhs[12] );
  targetVolume = *mxGetPr( prhs[13] );

	// -------------------------------------------------------------------------------------
	// RUN MINIMIZATION
	// -------------------------------------------------------------------------------------
	
	VectorXd x = minimizeElasticEnergy( faces, vertex,
			v1_ID, v2_ID, targetLengths, targetAngles,
			param, h, nu,
			alpha, target_ID, targetLocations,
      beta, targetVolume );

	// -------------------------------------------------------------------------------------
	// OUTPUT PROCESSING
	// -------------------------------------------------------------------------------------
	
	plhs[0] = mxCreateDoubleMatrix( Nv, 3, mxREAL );
	Eigen::Map<MatrixXd>( mxGetPr( plhs[0] ), Nv, 3 ) = x;

	return;

};
