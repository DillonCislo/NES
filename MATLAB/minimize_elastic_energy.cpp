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
#include <cassert>

#include <Eigen/Core>

#include "../external/LBFGSpp/include/LBFGS.h"
#include "../ElasticOptimization/ElasticProblem.h"
#include "../ElasticMesh/PolyhedronSoupToPolyhedronMesh.h"

typedef CGAL::Simple_cartesian<double>	 		      Kernel;
typedef CGAL::Polyhedron_3<Kernel, ElasticItems> 	Polyhedron;

typedef Polyhedron::Vertex_handle 		Vertex_handle;
typedef Polyhedron::Vertex_iterator 	Vertex_iterator;
typedef Polyhedron::Facet_iterator    Facet_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::HalfedgeDS 				HalfedgeDS;

typedef Eigen::VectorXd 		VectorXd;
typedef Eigen::MatrixXd 		MatrixXd;
typedef Eigen::VectorXi 		VectorXi;
typedef Eigen::MatrixXi 		MatrixXi;

// Enumeration of line search methods
enum LINE_SEARCH_METHOD {

  // Line search via backtracking
  LBFGS_LINESEARCH_METHOD_BACKTRACKING = 1,

  // Line search via backtracking + bracketing
  LBFGS_LINESEARCH_METHOD_BRACKETING = 2,
  
  // Line search method by Nocedal and Wright for the
  // strong Wolfe conditions
  LBFGS_LINESEARCH_METHOD_NOCEDAL = 3

};

// Prepare a polyhedral mesh object for the elastic minimization procedure
void preparePolyhedron( Polyhedron &P, const MatrixXi &faces, const MatrixXd &vertex,
    const VectorXi &v1_ID, const VectorXi &v2_ID,
		const VectorXd &targetLengths, const VectorXd &targetAngles,
    const VectorXd &thickness, double alpha,
    const VectorXi &target_ID, const MatrixXd &targetLocations,
    double mu, const MatrixXd &restrictVectors, const MatrixXd &restrictedLengths ) {


	// Read in basic mesh topology and vertex locations
	PolyhedronSoupToPolyhedronMesh<HalfedgeDS, double> inPoly( vertex, faces );
	P.delegate( inPoly );

  assert( P.is_valid() && "Failed to properly construct polyhedral mesh!" );

	// Update target geometry
	ElasticUpdater EU;

	EU.assignVertexID( P );

	if ( alpha > 0.0 ) {
		std::vector<Vertex_handle> tV;
		tV = EU.assignTargetVertex( P, target_ID, targetLocations );
	}

	EU.assignMajorEdges( P );
	EU.updateTargetGeometry( P, v1_ID, v2_ID, targetLengths, targetAngles );

  // Assign material parameters
  EU.assignMaterialParameters( P, thickness );

	// Update current geometry
	MatrixXd vertexCopy = vertex;
	VectorXd x = Eigen::Map<Eigen::VectorXd>( vertexCopy.data(), vertexCopy.size() );

	EU.updateCurrentGeometry( P, x );

  // Set initial/growth restriction geometry
  if ( mu > 0.0 ) {
    EU.setRestrictionGeometry( P, restrictVectors, restrictedLengths );
    EU.setInitialGeometry( P );
  }

	return;

};

// Run elastic minimization procedure
VectorXd minimizeElasticEnergy( const MatrixXi &faces, const MatrixXd &vertex,
		const VectorXi &v1_ID, const VectorXi &v2_ID,
		const VectorXd &targetLengths, const VectorXd &targetAngles,
		const LBFGSpp::LBFGSParam<double> &param, int linesearch_method,
		const VectorXd &thickness, double nu, double alpha,
		const VectorXi &target_ID, const MatrixXd &targetLocations,
    double beta, double targetVolume, bool usePhantom, const MatrixXi &phantomFaces,
    double mu, const MatrixXd &restrictVectors, const MatrixXd &restrictedLengths ) {

	Polyhedron P;
	preparePolyhedron( P, faces, vertex,
			v1_ID, v2_ID, targetLengths, targetAngles,
			thickness, alpha, target_ID, targetLocations,
      mu, restrictVectors, restrictedLengths );

  Polyhedron PP;
  if ( (beta > 0.0) && usePhantom ) {

    // NOTE: The phantom mesh is initiated with zero thickness
    // This is fine since the phantom mesh is only used in the
    // calculation of the target volume energy
    VectorXd zeroThickness = VectorXd::Zero( phantomFaces.rows(), 1 );

    preparePolyhedron( PP, phantomFaces, vertex,
        v1_ID, v2_ID, targetLengths, targetAngles,
        zeroThickness, alpha, target_ID, targetLocations,
        mu, restrictVectors, restrictedLengths );

  }

	// Initialize the problem structure
	ElasticProblem f(P, nu, mu, beta, targetVolume,
      usePhantom, PP, alpha, target_ID, targetLocations );

	// Initial guess
	MatrixXd vertexCopy = vertex;
	VectorXd x = Eigen::Map<VectorXd>( vertexCopy.data(), vertexCopy.size() );

	// Initialize solver - this is hideous there has to be a better way to do this
  LBFGSpp::LBFGSSolver<double, LBFGSpp::LineSearchBacktracking> solverBT(param);
  LBFGSpp::LBFGSSolver<double, LBFGSpp::LineSearchBracketing> solverBR(param);
  LBFGSpp::LBFGSSolver<double, LBFGSpp::LineSearchNocedalWright> solverNW(param);
  
	// Run solver
	double fx;
	int niter;

  try {

    switch (linesearch_method) {

      case LBFGS_LINESEARCH_METHOD_BACKTRACKING :
        niter = solverBT.minimize( f, x, fx );
        break;

      case LBFGS_LINESEARCH_METHOD_BRACKETING :
        niter = solverBR.minimize( f, x, fx );
        break;

      case LBFGS_LINESEARCH_METHOD_NOCEDAL :
        niter = solverNW.minimize( f, x, fx );
        break;

    }

  } catch ( const std::runtime_error &ere ) {

    mexWarnMsgTxt(ere.what());

  }

  /*
  // Test face ordering
  Facet_iterator fID;
  for( fID = P.facets_begin(); fID != P.facets_end(); fID++ ) {

    int vj = fID->halfedge()->next()->vertex()->id();
    int vk = fID->halfedge()->prev()->vertex()->id();
    int vi = fID->halfedge()->vertex()->id();

    vi++;
    vj++;
    vk++;

    std::cout << vi << " " << vj << " " << vk << std::endl;

  }
  */

	return x;

};

// Brief main function to call computational functionalities
void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] ) {

	// -------------------------------------------------------------------------------------
	// INPUT PROCESSING
	// -------------------------------------------------------------------------------------
	
	// VARIABLE DECLARATIONS ---------------------------------------------------------------
	
	int *m_faces; 			          // The face connectivity list
	MatrixXi faces;

	double *m_vertex; 		        // The initial vertex positions
	MatrixXd vertex;

	double *m_targetLengths; 	    // The list of target edge lengths
	VectorXd targetLengths;

	double *m_targetAngles; 	    // The list of target bend angles
	VectorXd targetAngles;
	
	int *m_v1_ID;	 		            // The start vertex ID of each edge
	VectorXi v1_ID;

	int *m_v2_ID; 			          // The end vertex ID of each edge
	VectorXi v2_ID;

	int *m_target_ID; 		        // The vertex IDs of those vertices with a target
	VectorXi target_ID;		        // location

	double *m_targetLocations;	  // The target locations of those vertices
	MatrixXd targetLocations;

  int *m_phantomFaces;          // The phantom face connectivity list used to
  MatrixXi phantomFaces;        // calculate the volume of a surface with holes

  double *m_restrictVectors;    // Tangent vectors in the plane of each face
  MatrixXd restrictVectors;     // along which to restrict growth

  double *m_restrictedLengths;  // The maximum projected length of each edge along
  MatrixXd restrictedLengths;   // the respective restriction vector

  double *m_thickness;          // The thickness of the elastic sheet
  VectorXd thickness;

	double nu; 			      // Poisson's ratio
	double alpha; 	      // The coefficient of target vertex correspondence
  double beta;          // The coefficient of the fixed folume energy
  double targetVolume;  // The target volume for the fixed volume energy
  double mu;            // The coefficient of the growth restriction energy

  bool usePhantom;      // If true, the fixed volume energy will be calculated
                        // using the phantom face connectivity list

  int linesearch_method;  // The line search method to employ

	// L-BFGS parameter objects
	LBFGSpp::LBFGSParam<double> param;

	int Nf = 0; 	// Number of faces
	int Nv = 0; 	// Number of vertices
	int Ne = 0;	  // Number of edges
	int Nt = 0; 	// Number of target vertices
  int Npf = 0;  // Number of phantom faces

	// EXTRACT INPUT PARAMETERS ------------------------------------------------------------

	// Check for proper number of arguments
	if (nrhs != 19 ) {
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

	// Absolute tolerance for convergence test
	if ( (idx = mxGetFieldNumber( prhs[6], "epsilon")) == -1 ) {
		mexErrMsgTxt("No absolute epsilon field!");
	} else {
		param.epsilon = (double) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
	}

  // Relative tolerance for convergence test
  if ( (idx = mxGetFieldNumber( prhs[6], "epsilon_rel")) == -1 ) {
    mexErrMsgTxt("No relative epsilon field!");
  } else {
    param.epsilon_rel = (double) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
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

	// Linesearch termination condition
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

  // Linesearch method
  if ( (idx = mxGetFieldNumber( prhs[6], "linesearch_method")) == -1 ) {
    mexErrMsgTxt("No linesearch method field!");
  } else {
    tmpLS = (int) *mxGetPr(mxGetFieldByNumber( prhs[6], 0, idx ));
  }

  if (tmpLS == 1) {
    linesearch_method = LBFGS_LINESEARCH_METHOD_BACKTRACKING;
  } else if (tmpLS == 2) {
    linesearch_method = LBFGS_LINESEARCH_METHOD_BRACKETING;
  } else if (tmpLS == 3) {
    linesearch_method = LBFGS_LINESEARCH_METHOD_NOCEDAL;
  } else {
    mexErrMsgTxt("Linesearch method field invalid!");
  };

	// GET MATERIAL PARAMETERS -------------------------------------------------------------
	m_thickness = mxGetPr( prhs[7] );
  thickness = Eigen::Map<VectorXd>( m_thickness, Nf, 1 );

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

  usePhantom = *mxGetLogicals( prhs[14] );

  m_phantomFaces = (int*) mxGetData( prhs[15] );
  Npf = (int) mxGetM( prhs[15] );
  phantomFaces = Eigen::Map<MatrixXi>( m_phantomFaces, Npf, 3 );

  // GET GROWTH RESTRICTION PARAMETERS ---------------------------------------------------
  mu = *mxGetPr( prhs[16] );

  m_restrictVectors = mxGetPr( prhs[17] );
  restrictVectors = Eigen::Map<MatrixXd>( m_restrictVectors, Nf, 3 );

  m_restrictedLengths = mxGetPr( prhs[18] );
  restrictedLengths = Eigen::Map<MatrixXd>( m_restrictedLengths, Nf, 3 );

	// -------------------------------------------------------------------------------------
	// RUN MINIMIZATION
	// -------------------------------------------------------------------------------------
  
  VectorXd x( 3*Nv, 1 );
	
  try {

	  x = minimizeElasticEnergy( faces, vertex,
			  v1_ID, v2_ID, targetLengths, targetAngles,
		  	param, linesearch_method, thickness, nu,
			  alpha, target_ID, targetLocations,
        beta, targetVolume, usePhantom, phantomFaces,
        mu, restrictVectors, restrictedLengths );

  } catch ( const std::invalid_argument &eia ) {

    mexErrMsgTxt(eia.what());

  } catch ( const std::logic_error &ele ) {

    mexErrMsgTxt(ele.what());

  } catch ( ... ) {

    mexErrMsgTxt("Unknown exception caught. Optimization terminating abnormally");

  }

	// -------------------------------------------------------------------------------------
	// OUTPUT PROCESSING
	// -------------------------------------------------------------------------------------
	
	plhs[0] = mxCreateDoubleMatrix( 3*Nv, 1, mxREAL );
	Eigen::Map<MatrixXd>( mxGetPr( plhs[0] ), 3*Nv, 1 ) = x;

	return;

};
