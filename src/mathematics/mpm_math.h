/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */


#ifndef MPM_MATH_H_
#define MPM_MATH_H_

#include <singular_value_decompose.h>
#include <inverse.h>
#include <determinant.h>


#include <iostream>

namespace MPM_Math {

  /*
 * deviator of a tensor
 */
static inline Matrix3d Deviator(const Matrix3d M) {
	Matrix3d eye = Matrix3d::identity();
	eye *= M.trace() / 3.0;
	return M - eye;
}

static inline Matrix2d Deviator(const Matrix2d M) {
	Matrix2d eye = Matrix2d::identity();
	eye *= 0.5 * M.trace();
	return M - eye;
}
/*
 * Polar Decomposition M = R * T
 * where R is a rotation and T a pure translation/stretch matrix.
 *
 * The decomposition is achieved using SVD, i.e. M = U S V^T,
 * where U = R V and S is diagonal.
 *
 *
 * For any physically admissible deformation gradient, the determinant of R must equal +1.
 * However, scenerios can arise, where the particles interpenetrate and cause inversion, leading to a determinant of R equal to -1.
 * In this case, the inversion direction is heuristically identified with the eigenvector of the smallest entry of S, which should work for most cases.
 * The sign of this corresponding eigenvalue is flipped, the original matrix M is recomputed using the flipped S, and the rotation and translation matrices are
 * obtained again from an SVD. The rotation should proper now, i.e., det(R) = +1.
 */

static inline bool PolDec(Matrix3d M, Matrix3d &R, Matrix3d &T, bool scaleF) {
	Matrix3d S = M;
	const std::pair<Matrix3d, Matrix3d> &uv = singular_value_decompose(S); // SVD(A) = U S V*
	Matrix3d U = uv.first;
	Matrix3d V = uv.second;
	Matrix3d eye = Matrix3d::identity();

	// now do polar decomposition into M = R * T, where R is rotation
	// and T is translation matrix
	R = U * V.transpose();
	T = V * S * V.transpose();

	if (determinant(R) < 0.0) { // this is an improper rotation
		// identify the smallest entry in S and flip its sign
		int imin = S(0, 0) < S(1, 1) && S(0, 0) < S(2, 2)? 0:
							 S(1, 1) < S(0, 0) && S(1, 1) < S(2, 2)? 1:
																											 2;
		S(imin, imin) *= -1.0;

		R = M * V * inverse(S) * V.transpose(); // recompute R using flipped stretch eigenvalues
	}

	/*
	 * scale S to avoid small principal strains
	 */

	if (scaleF) {
		double min = 0.3; // 0.3^2 = 0.09, should suffice for most problems
		double max = 2.0;
		for (int i = 0; i < 3; i++) {
			if (S(i, i) < min) {
				S(i, i) = min;
			} else if (S(i, i) > max) {
				S(i, i) = max;
			}
		}
		T = V * S * V.transpose();
	}

	if (determinant(R) > 0.0) {
		return true;
	} else {
		return false;
	}
}

static inline bool PolDec(Matrix3d M, Matrix3d &R) {
	Matrix3d S = M;
	const std::pair<Matrix3d, Matrix3d> &uv = singular_value_decompose(S); // SVD(A) = U S V*
	Matrix3d U = uv.first;
	Matrix3d V = uv.second;
	Matrix3d eye = Matrix3d::identity();

  // now do polar decomposition into M = R * T, where R is rotation
  // and T is translation matrix
  R = U * V.transpose();

  if (determinant(R) < 0.0) { // this is an improper rotation
    // identify the smallest entry in S and flip its sign
		int imin = S(0, 0) < S(1, 1) && S(0, 0) < S(2, 2)? 0:
							 S(1, 1) < S(0, 0) && S(1, 1) < S(2, 2)? 1:
																											 2;
		S(imin, imin) *= -1.0;

    R = M * V * inverse(S) *
        V.transpose(); // recompute R using flipped stretch eigenvalues
  }

  if (determinant(R) > 0.0) {
    return true;
  } else {
    return false;
  }
}

}

#endif /* SMD_MATH_H_ */
