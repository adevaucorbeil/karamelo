/* -*- c++ -*- ----------------------------------------------------------*/


#ifndef MPM_MATH_H_
#define MPM_MATH_H_

#include <Eigen/Eigen>
#include <iostream>
using namespace Eigen;
using namespace std;

namespace MPM_Math {

  /*
 * deviator of a tensor
 */
static inline Matrix3d Deviator(const Matrix3d M) {
	Matrix3d eye;
	eye.setIdentity();
	eye *= M.trace() / 3.0;
	return M - eye;
}

static inline Matrix2d Deviator(const Matrix2d M) {
	Matrix2d eye;
	eye.setIdentity();
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

	JacobiSVD<Matrix3d> svd(M, ComputeFullU | ComputeFullV); // SVD(A) = U S V*
	Vector3d S_eigenvalues = svd.singularValues();
	Matrix3d S = svd.singularValues().asDiagonal();
	Matrix3d U = svd.matrixU();
	Matrix3d V = svd.matrixV();
	Matrix3d eye;
	eye.setIdentity();

	// now do polar decomposition into M = R * T, where R is rotation
	// and T is translation matrix
	R = U * V.transpose();
	T = V * S * V.transpose();

	if (R.determinant() < 0.0) { // this is an improper rotation
		// identify the smallest entry in S and flip its sign
		int imin;
		S_eigenvalues.minCoeff(&imin);
		S(imin, imin) *= -1.0;

		R = M * V * S.inverse() * V.transpose(); // recompute R using flipped stretch eigenvalues
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

	if (R.determinant() > 0.0) {
		return true;
	} else {
		return false;
	}
}

/*
 * Pseudo-inverse via SVD
 */

static inline void pseudo_inverse_SVD(Matrix3d &M) {

	//JacobiSVD < Matrix3d > svd(M, ComputeFullU | ComputeFullV);
	JacobiSVD<Matrix3d> svd(M, ComputeFullU); // one Eigevector base is sufficient because matrix is square and symmetric

	Vector3d singularValuesInv;
	Vector3d singularValues = svd.singularValues();

//cout << "Here is the matrix V:" << endl << V * singularValues.asDiagonal() * U << endl;
//cout << "Its singular values are:" << endl << singularValues << endl;

	double pinvtoler = 1.0e-16; // 2d machining example goes unstable if this value is increased (1.0e-16).
	for (int row = 0; row < 3; row++) {
		if (singularValues(row) > pinvtoler) {
			singularValuesInv(row) = 1.0 / singularValues(row);
		} else {
			singularValuesInv(row) = 0.0;
		}
	}

	M = svd.matrixU() * singularValuesInv.asDiagonal() * svd.matrixU().transpose();

	/*
	 * commented version below is for non-symmetric matrices
	 */

//	JacobiSVD < Matrix3d > svd(M, ComputeFullU | ComputeFullV);
//
//	Vector3d singularValuesInv;
//	Vector3d singularValues = svd.singularValues();
//
//	//cout << "Here is the matrix V:" << endl << V * singularValues.asDiagonal() * U << endl;
//	//cout << "Its singular values are:" << endl << singularValues << endl;
//
//	double pinvtoler = 1.0e-16; // 2d machining example goes unstable if this value is increased (1.0e-16).
//	for (int row = 0; row < 3; row++) {
//		if (singularValues(row) > pinvtoler) {
//			singularValuesInv(row) = 1.0 / singularValues(row);
//		} else {
//			singularValuesInv(row) = 0.0;
//		}
//	}
//
//	M = svd.matrixU() * singularValuesInv.asDiagonal() * svd.matrixV().transpose();
}

/*
 * specialized SVD for symmetric input matrices.
 * limits the eigenvalues to the supplied limit value.
 */

static inline void pseudo_inverse_SVD_limit_eigenvalue(Matrix3d &M, const double limit) {

	JacobiSVD<Matrix3d> svd(M, ComputeFullU); // one Eigevector base is sufficient because matrix is square and symmetric

	Vector3d singularValuesInv;
	Vector3d singularValues = svd.singularValues();

	double pinvtoler = 1.0e-16; // 2d machining example goes unstable if this value is increased (1.0e-16).
	for (int row = 0; row < 3; row++) {
		if (singularValues(row) > pinvtoler) {
			singularValuesInv(row) = 1.0 / singularValues(row);
			if (singularValuesInv(row) > limit) {
				printf("eigenvalue is %f\n", singularValuesInv(row));
				singularValuesInv(row) = limit;
			}
		} else {
			singularValuesInv(row) = 0.0;
		}
	}

	M = svd.matrixU() * singularValuesInv.asDiagonal() * svd.matrixU().transpose();
}

/*
 * test if two matrices are equal
 */
static inline double TestMatricesEqual(Matrix3d A, Matrix3d B, double eps) {
	Matrix3d diff;
	diff = A - B;
	double norm = diff.norm();
	if (norm > eps) {
		printf("Matrices A and B are not equal! The L2-norm difference is: %g\n", norm);
		cout << "Here is matrix A:" << endl << A << endl;
		cout << "Here is matrix B:" << endl << B << endl;
	}
	return norm;
}

/* ----------------------------------------------------------------------
 Limit eigenvalues of a matrix to upper and lower bounds.
 ------------------------------------------------------------------------- */

static inline Matrix3d LimitEigenvalues(Matrix3d S, double limitEigenvalue) {

	/*
	 * compute Eigenvalues of matrix S
	 */
	SelfAdjointEigenSolver<Matrix3d> es;
	es.compute(S);

	double max_eigenvalue = es.eigenvalues().maxCoeff();
	double min_eigenvalue = es.eigenvalues().minCoeff();
	double amax_eigenvalue = fabs(max_eigenvalue);
	double amin_eigenvalue = fabs(min_eigenvalue);

	if ((amax_eigenvalue > limitEigenvalue) || (amin_eigenvalue > limitEigenvalue)) {
		if (amax_eigenvalue > amin_eigenvalue) { // need to scale with max_eigenvalue
			double scale = amax_eigenvalue / limitEigenvalue;
			Matrix3d V = es.eigenvectors();
			Matrix3d S_diag = V.inverse() * S * V; // diagonalized input matrix
			S_diag /= scale;
			Matrix3d S_scaled = V * S_diag * V.inverse(); // undiagonalize matrix
			return S_scaled;
		} else { // need to scale using min_eigenvalue
			double scale = amin_eigenvalue / limitEigenvalue;
			Matrix3d V = es.eigenvectors();
			Matrix3d S_diag = V.inverse() * S * V; // diagonalized input matrix
			S_diag /= scale;
			Matrix3d S_scaled = V * S_diag * V.inverse(); // undiagonalize matrix
			return S_scaled;
		}
	} else { // limiting does not apply
		return S;
	}
}

static inline bool LimitMinMaxEigenvalues(Matrix3d &S, double min, double max) {

	/*
	 * compute Eigenvalues of matrix S
	 */
	SelfAdjointEigenSolver<Matrix3d> es;
	es.compute(S);

	if ((es.eigenvalues().maxCoeff() > max) || (es.eigenvalues().minCoeff() < min)) {
		Matrix3d S_diag = es.eigenvalues().asDiagonal();
		Matrix3d V = es.eigenvectors();
		for (int i = 0; i < 3; i++) {
			if (S_diag(i, i) < min) {
				//printf("limiting eigenvalue %f --> %f\n", S_diag(i, i), min);
				//printf("these are the eigenvalues of U: %f %f %f\n", es.eigenvalues()(0), es.eigenvalues()(1), es.eigenvalues()(2));
				S_diag(i, i) = min;
			} else if (S_diag(i, i) > max) {
				//printf("limiting eigenvalue %f --> %f\n", S_diag(i, i), max);
				S_diag(i, i) = max;
			}
		}
		S = V * S_diag * V.inverse(); // undiagonalize matrix
		return true;
	} else {
		return false;
	}
}

static inline void reconstruct_rank_deficient_shape_matrix(Matrix3d &K) {

	JacobiSVD<Matrix3d> svd(K, ComputeFullU | ComputeFullV);
	Vector3d singularValues = svd.singularValues();

	for (int i = 0; i < 3; i++) {
		if (singularValues(i) < 1.0e-8) {
			singularValues(i) = 1.0;
		}
	}

//		int imin;
//		double minev = singularValues.minCoeff(&imin);
//
//		printf("min eigenvalue=%f has index %d\n", minev, imin);
//		Vector3d singularVec = U.col(0).cross(U.col(1));
//		cout << "the eigenvalues are " << endl << singularValues << endl;
//		cout << "the singular vector is " << endl << singularVec << endl;
//
//		// reconstruct original K
//
//		singularValues(2) = 1.0;

	K = svd.matrixU() * singularValues.asDiagonal() * svd.matrixV().transpose();
	//cout << "the reconstructed K is " << endl << K << endl;
	//exit(1);
}

/* ----------------------------------------------------------------------
 helper functions for crack_exclude
 ------------------------------------------------------------------------- */
static inline bool IsOnSegment(double xi, double yi, double xj, double yj, double xk, double yk) {
	return (xi <= xk || xj <= xk) && (xk <= xi || xk <= xj) && (yi <= yk || yj <= yk) && (yk <= yi || yk <= yj);
}

static inline char ComputeDirection(double xi, double yi, double xj, double yj, double xk, double yk) {
	double a = (xk - xi) * (yj - yi);
	double b = (xj - xi) * (yk - yi);
	return a < b ? -1.0 : a > b ? 1.0 : 0;
}

/** Do line segments (x1, y1)--(x2, y2) and (x3, y3)--(x4, y4) intersect? */
static inline bool DoLineSegmentsIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
	char d1 = ComputeDirection(x3, y3, x4, y4, x1, y1);
	char d2 = ComputeDirection(x3, y3, x4, y4, x2, y2);
	char d3 = ComputeDirection(x1, y1, x2, y2, x3, y3);
	char d4 = ComputeDirection(x1, y1, x2, y2, x4, y4);
	return (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0)))
			|| (d1 == 0 && IsOnSegment(x3, y3, x4, y4, x1, y1)) || (d2 == 0 && IsOnSegment(x3, y3, x4, y4, x2, y2))
			|| (d3 == 0 && IsOnSegment(x1, y1, x2, y2, x3, y3)) || (d4 == 0 && IsOnSegment(x1, y1, x2, y2, x4, y4));
}

/*
 * Pseudo-inverse via SVD for a 4d matrix
 */

static inline void pinv4(Matrix4d &M) {

	JacobiSVD < Matrix4d > svd(M, ComputeFullU | ComputeFullV);
	//JacobiSVD<Matrix4d> svd(M, ComputeFullU); // one Eigevector base is sufficient because matrix is square and symmetric

	Vector4d singularValuesInv;
	Vector4d singularValues = svd.singularValues();

//cout << "Here is the matrix V:" << endl << svd.MatrixV * singularValues.asDiagonal() * U << endl;
cout << "Its singular values are:" << endl << singularValues << endl;

	double pinvtoler = 1.0e-3; // 2d machining example goes unstable if this value is increased (1.0e-16).
	for (int row = 0; row < 4; row++) {
		if (singularValues(row) > pinvtoler) {
			singularValuesInv(row) = 1.0 / singularValues(row);
		} else {
			singularValuesInv(row) = 0.0;
		}
	}

	M = svd.matrixV() * singularValuesInv.asDiagonal() * svd.matrixU().transpose();
}

/* circle_circle_intersection() *
 * Determine the points where 2 circles in a common plane intersect.
 *
 * int circle_circle_intersection(
 *                                // center and radius of 1st circle
 *                                double x0, double y0, double r0,
 *                                // center and radius of 2nd circle
 *                                double x1, double y1, double r1,
 *                                // 1st intersection point
 *                                double *xi, double *yi,
 *                                // 2nd intersection point
 *                                double *xi_prime, double *yi_prime)
 *
 * This is a public domain work. 3/26/2005 Tim Voght
 *
 */
//#include <stdio.h>
//#include <math.h>

static int circle_circle_intersection(double x0, double y0, double r0,
                               double x1, double y1, double r1,
                               double &xi, double &yi,
                               double &xi_prime, double &yi_prime)
{
  double a, dx, dy, d, h, rx, ry;
  double x2, y2;

  /* dx and dy are the vertical and horizontal distances between
   * the circle centers.
   */
  dx = x1 - x0;
  dy = y1 - y0;

  /* Determine the straight-line distance between the centers. */
  //d = sqrt((dy*dy) + (dx*dx));
  d = hypot(dx,dy); // Suggested by Keith Briggs

  /* Check for solvability. */
  if (d > (r0 + r1))
  {
    /* no solution. circles do not intersect. */
    return 0;
  }
  if (d < fabs(r0 - r1))
  {
    /* no solution. one circle is contained in the other */
    return 0;
  }

  /* 'point 2' is the point where the line through the circle
   * intersection points crosses the line between the circle
   * centers.
   */

  /* Determine the distance from point 0 to point 2. */
  a = ((r0*r0) - (r1*r1) + (d*d)) / (2.0 * d) ;

  /* Determine the coordinates of point 2. */
  x2 = x0 + (dx * a/d);
  y2 = y0 + (dy * a/d);

  /* Determine the distance from point 2 to either of the
   * intersection points.
   */
  h = sqrt((r0*r0) - (a*a));

  /* Now determine the offsets of the intersection points from
   * point 2.
   */
  rx = -dy * (h/d);
  ry = dx * (h/d);

  /* Determine the absolute intersection points. */
  xi = x2 + rx;
  xi_prime = x2 - rx;
  yi = y2 + ry;
  yi_prime = y2 - ry;

  return 1;
}

static inline double clamp(double x, double lowerlimit, double upperlimit)
{
    if (x < lowerlimit) x = lowerlimit;
    if (x > upperlimit) x = upperlimit;
    return x;
}

static inline double smoothstep(double edge0, double edge1, double x)
{
    // Scale, bias and saturate x to 0..1 range
    x = clamp((x - edge0)/(edge1 - edge0), 0.0, 1.0);
    // Evaluate polynomial
    return x*x*(3 - 2*x);
}

}

#endif /* SMD_MATH_H_ */
