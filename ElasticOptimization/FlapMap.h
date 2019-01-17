/*!
 * 	\file FlapMap.h
 * 	\brief A class that maps derivative quantities from hinges to a local flap stencil
 *
 * 	\author Dillon Cislo
 * 	\date 01/11/2019
 *
 */

#ifndef _FLAP_MAP_H_
#define _FLAP_MAP_H_

#include <vector>

#include <Eigen/Core>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include "../ElasticMesh/ElasticItems.h"

///
/// A class that maps gradient and Hessian quantities from hinges to a local flap stencil
///
class FlapMap {

	private:

		typedef typename CGAL::Simple_cartesian<double> 		Kernel;
		typedef typename CGAL::Polyhedron_3<Kernel, ElasticItems>	Polyhedron;

		typedef typename Polyhedron::Face_handle 		Face_handle;
		typedef typename Polyhedron::Halfedge_handle 		Halfedge_handle;

		typedef typename Eigen::RowVector3d 			RowVector3d;
		typedef typename Eigen::RowVectorXd 			RowVectorXd;
		typedef typename Eigen::Vector3d 			Vector3d;
		typedef typename Eigen::VectorXd 			VectorXd;
		typedef typename Eigen::Matrix3d 			Matrix3d;
		typedef typename Eigen::MatrixXd 			MatrixXd;

		typedef typename Eigen::Matrix<double,1,12> 		HingeGrad;
		typedef typename Eigen::Matrix<double,12,12> 		HingeHess;
		typedef typename Eigen::Matrix<double,1,18> 		FlapGrad;
		typedef typename Eigen::Matrix<double,18,18> 		FlapHess;

	public:

		///
		///  Map all hinge gradients of a particular face to the local flap stencil
		///
		static std::vector<FlapGrad> hingeGradients( Face_handle f );

		///
		/// Map the hinge gradient of a particular hinge to the local flap stencil
		///
		static FlapGrad mapHingeGradToFlap( Halfedge_handle he, int i );

		///
		/// Map all hinge Hessians of a particular face to the local flap stencil
		///
		static std::vector<FlapHess> hingeHessians( Face_handle f );

		///
		/// Map the hinge Hessian of a particular hinge to the local flap stencil
		///
		static FlapHess mapHingeHessToFlap( Halfedge_handle he, int i );

};

///
/// Map all hinge gradients of a particular face to the local flap stencil
///
std::vector< Eigen::Matrix<double,1,18> >
FlapMap::hingeGradients( Face_handle f ) {

	std::vector<FlapGrad> fGS;
	fGS.reserve( 3 );

	Halfedge_handle he = f->halfedge();
	for( int i = 0; i < 3; i++ ) {

		if ( he->is_border() || he->opposite()->is_border() ) {

			fGS.push_back( FlapGrad::Zero() );

		} else {

			fGS.push_back( mapHingeGradToFlap( he, i ) );

		}

		he = he->next();

	}

	return fGS;

};

///
/// Map all hinge Hessians of a particular face to the local flap stencil
///
std::vector< Eigen::Matrix<double,18,18> >
FlapMap::hingeHessians( Face_handle f ) {

	std::vector<FlapHess> fHS;
	fHS.reserve( 3 );

	Halfedge_handle he = f->halfedge();
	for( int i = 0; i < 3; i++ ) {

		if ( he->is_border() || he->opposite()->is_border() ) {

			fHS.push_back( FlapHess::Zero() );

		} else {

			fHS.push_back( mapHingeHessToFlap( he, i ) );

		}

		he = he->next();

	}

	return fHS;

};

///
/// Map the hinge gradient of a particular hinge to the local flap stencil
///
Eigen::Matrix<double,1,18>
FlapMap::mapHingeGradToFlap( Halfedge_handle he, int i ) {

	HingeGrad hG = he->gradPhi();

	RowVector3d hG0 = hG.segment<3>(0);
	RowVector3d hG1 = hG.segment<3>(3);
	RowVector3d hG2 = hG.segment<3>(6);
	RowVector3d hG3 = hG.segment<3>(9);
	RowVector3d Z3 = RowVector3d::Zero();

	// Choose appropriate mapping based on relative order of the hinge in the local flap
	FlapGrad fG;
	switch( i ) {

		case 0:

			if ( he->isMajor() ) {
			       	fG << hG2, hG0, hG1, hG3, Z3, Z3;
			} else {
				fG << hG3, hG1, hG0, hG2, Z3, Z3;
			}
			
			break;

		case 1:

			if ( he->isMajor() ) {
				fG << hG1, hG2, hG0, Z3, hG3, Z3;
			} else {
				fG << hG0, hG3, hG1, Z3, hG2, Z3;
			}

			break;

		case 2:
			if ( he->isMajor() ) {
				fG << hG0, hG1, hG2, Z3, Z3, hG3;
			} else {
				fG << hG1, hG0, hG3, Z3, Z3, hG2;
			}

			break;

	}

	return fG;

};

///
/// Map the hinge Hessian of a particular hinge to the local flap stencil
///
Eigen::Matrix<double,18,18>
FlapMap::mapHingeHessToFlap( Halfedge_handle he, int i ) {

	HingeHess hH = he->hessPhi();

	Matrix3d H00 = hH.block<3,3>(0,0);
	Matrix3d H01 = hH.block<3,3>(0,3);
	Matrix3d H02 = hH.block<3,3>(0,6);
	Matrix3d H03 = hH.block<3,3>(0,9);

	Matrix3d H10 = hH.block<3,3>(3,0);
	Matrix3d H11 = hH.block<3,3>(3,3);
	Matrix3d H12 = hH.block<3,3>(3,6);
	Matrix3d H13 = hH.block<3,3>(3,9);

	Matrix3d H20 = hH.block<3,3>(6,0);
	Matrix3d H21 = hH.block<3,3>(6,3);
	Matrix3d H22 = hH.block<3,3>(6,6);
	Matrix3d H23 = hH.block<3,3>(6,9);

	Matrix3d H30 = hH.block<3,3>(9,0);
	Matrix3d H31 = hH.block<3,3>(9,3);
	Matrix3d H32 = hH.block<3,3>(9,6);
	Matrix3d H33 = hH.block<3,3>(9,9);

	Matrix3d Z3 = Matrix3d::Zero();

	// Choose appropriate mapping based on relative order of the hinge in the local flap
	FlapHess fH;
	switch( i ) {

		case 0:

			if ( he->isMajor() ) {

				fH << H22, H20, H21, H23, Z3, Z3,
				      H02, H00, H01, H03, Z3, Z3,
				      H12, H10, H11, H13, Z3, Z3,
				      H32, H30, H31, H33, Z3, Z3,
				      Z3,  Z3,  Z3,  Z3,  Z3, Z3,
				      Z3,  Z3,  Z3,  Z3,  Z3, Z3;
			
			} else {

				fH << H33, H31, H30, H32, Z3, Z3,
				      H13, H11, H10, H12, Z3, Z3,
				      H03, H01, H00, H02, Z3, Z3,
				      H23, H21, H20, H22, Z3, Z3,
				      Z3,  Z3,  Z3,  Z3,  Z3, Z3,
				      Z3,  Z3,  Z3,  Z3,  Z3, Z3;

			}

			break;

		case 1:

			if ( he->isMajor() ) {

				fH << H11, H12, H10, Z3, H13, Z3,
				      H21, H22, H20, Z3, H23, Z3,
				      H01, H02, H00, Z3, H03, Z3,
				      Z3,  Z3,  Z3,  Z3, Z3,  Z3,
				      H31, H32, H30, Z3, H33, Z3,
				      Z3,  Z3,  Z3,  Z3, Z3,  Z3;

			} else {

				fH << H00, H03, H01, Z3, H02, Z3,
				      H30, H33, H31, Z3, H32, Z3,
				      H10, H13, H11, Z3, H12, Z3,
				      Z3,  Z3,  Z3,  Z3, Z3,  Z3,
				      H20, H23, H21, Z3, H22, Z3,
				      Z3,  Z3,  Z3,  Z3, Z3,  Z3;

			}

			break;

		case 2:

			if ( he->isMajor() ) {

				fH << H00, H01, H02, Z3, Z3, H03,
				      H10, H11, H12, Z3, Z3, H13,
				      H20, H21, H22, Z3, Z3, H23,
				      Z3,  Z3,  Z3,  Z3, Z3, Z3,
				      Z3,  Z3,  Z3,  Z3, Z3, Z3,
				      H30, H31, H32, Z3, Z3, H33;

			} else {

				fH << H11, H10, H13, Z3, Z3, H12,
				      H01, H00, H03, Z3, Z3, H02,
				      H31, H30, H33, Z3, Z3, H32,
				      Z3,  Z3,  Z3,  Z3, Z3, Z3,
				      Z3,  Z3,  Z3,  Z3, Z3, Z3,
				      H21, H20, H23, Z3, Z3, H22;

			}

			break;

	}

	return fH;

};

#endif
