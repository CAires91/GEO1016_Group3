/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"


using namespace easy3d;


bool Calibration::calibration(
    const std::vector<Vector3D> &points_3d, /// input: An array of 3D points.
    const std::vector<Vector2D> &points_2d, /// input: An array of 2D image points.
    double &fx, /// output: focal length (i.e., K[0][0]).
    double &fy, /// output: focal length (i.e., K[1][1]).
    double &cx, /// output: x component of the principal point (i.e., K[0][2]).
    double &cy, /// output: y component of the principal point (i.e., K[1][2]).
    double &s, /// output: skew factor (i.e., K[0][1]), which is s = -alpha * cot(theta).
    Matrix33 &R, /// output: the 3x3 rotation matrix encoding camera rotation.
    Vector3D &t) /// output：a 3D vector encoding camera translation. {
{

    // ----------------------- check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)

    int num_3d_points = points_3d.size();
    int num_2d_points = points_2d.size();

    if (num_3d_points < 6 || num_2d_points < 6 || num_3d_points != num_2d_points) {
        return false;
    }

    // ----------------------- construct the P matrix (so P * m = 0).
    int num_rows = num_3d_points * 2;

    Matrix P(num_rows, 12, 0.0);
    for (int i = 0; i < num_3d_points; i++) {
        Vector3D P1 = points_3d[i];
        Vector2D P2 = points_2d[i];

        P.set_row(2 * i, {
                      P1.x(), P1.y(), P1.z(), 1, 0, 0, 0, 0, -P2.x() * P1.x(), -P2.x() * P1.y(), -P2.x() * P1.z(), -P2.x()
                  });
        P.set_row(2 * i + 1, {
                      0, 0, 0, 0, P1.x(), P1.y(), P1.z(), 1, -P2.y() * P1.x(), -P2.y() * P1.y(), -P2.y() * P1.z(), -P2.y()
                  });;
    }

    // ----------------------- solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.

    const int m = 2 * num_3d_points, n = 12;

    Matrix U(m, m, 0.0); // initialized with 0s
    Matrix S(m, n, 0.0); // initialized with 0s
    Matrix V(n, n, 0.0); // initialized with 0s

    // compute the SVD decomposition of A
    svd_decompose(P, U, S, V);

    /// get the last column of the V matrix

    Vector m_vector = V.get_column(V.cols() - 1);

    /// build small m or cursive M matrix (m_matrix)

    Matrix m_matrix(3, 4, 0.0);

    int idx = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            m_matrix[i][j] = m_vector[idx++];
        }
    }

    // get a3 from m_matrix:

    double a3_1 = m_matrix(2, 0);

    double a3_2 = m_matrix(2, 1);

    double a3_3 = m_matrix(2, 2);

    Vector3D a3 = Vector3D(a3_1, a3_2, a3_3);

    double rho = 1 / norm(a3);

    Matrix M = rho * m_matrix;

    // ---------------------- extract intrinsic parameters from M.

    // calculate a1 using vector3d
    double a1_1 = m_matrix(0, 0);
    double a1_2 = m_matrix(0, 1);
    double a1_3 = m_matrix(0, 2);
    Vector3D a1 = Vector3D(a1_1, a1_2, a1_3);

    // calculate a2 using vector3d
    double a2_1 = m_matrix(1, 0);
    double a2_2 = m_matrix(1, 1);
    double a2_3 = m_matrix(1, 2);
    Vector3D a2 = Vector3D(a2_1, a2_2, a2_3);

    // calculate cx = rho^2 *(a1 * a3)
    cx = (rho * rho) * dot(a1, a3);

    // calculate cy = rho^2 *(a2 * a3)
    cy = (rho * rho) * dot(a2, a3);

    // calculate cos(ø) = ((a1*a3)*(a2*a3))/(||a1*a3||*||a2*a3||)
    Vector3D cross_a1_a3 = cross(a1, a3);
    Vector3D cross_a2_a3 = cross(a2, a3);

    double norm_cross_a1_a3 = norm(cross_a1_a3);
    double norm_cross_a2_a3 = norm(cross_a2_a3);

    double cos_theta = -1 * dot(cross_a1_a3, cross_a2_a3) / (norm_cross_a1_a3 * norm_cross_a2_a3);

    // Clamp cos_theta to the range [-1, 1]
    if (cos_theta > 1.0) {
        cos_theta = 1.0;
    } else if (cos_theta < -1.0) {
        cos_theta = -1.0;
    }

    // calculate alpha= rho^2 * ||a1*a3|| * sin(ø)

    double sin_theta = sqrt(1 - cos_theta * cos_theta);

    double alpha = (rho * rho) * norm_cross_a1_a3 * sin_theta;


    // calculate beta = rho^2 * ||a2*a3|| * sin(ø)

    double beta = (rho * rho) * norm_cross_a2_a3 * sin_theta;

    // calculate fx = alpha
    fx = alpha;

    // calculate fy = beta / sin(ø)

    fy = beta / sin_theta;

    // calculate s = -alpha * cot(ø)
    s = -alpha / tan(acos(cos_theta));


    // ----------------------  extract extrinsic parameters from M.

    // extract b
    double b1 = m_matrix(0, 3);
    double b2 = m_matrix(1, 3);
    double b3 = m_matrix(2, 3);
    Vector3D b = Vector3D(b1, b2, b3);

    Matrix33 K(fx, s, cx,
               0, fy, cy,
               0, 0, 1);

    // calculate t= rho*K-1*b
    t = rho * inverse(K) * b;

    // if t(2) < 0, then rho = -rho
    if (t.z() < 0) {
        rho = -rho;
        t = rho * inverse(K) * b;
    }

    // calculate r1 = (a2*a3)/||a2*a3||
    Vector3D r1 = cross_a2_a3 / norm_cross_a2_a3;

    // calculate r3= rho * a3
    Vector3D r3 = rho * a3;

    // calculate r2 = r3 x r1
    Vector3D r2 = cross(r3, r1);

    // calculate R
    R = Matrix33(r1.x(), r1.y(), r1.z(),
                 r2.x(), r2.y(), r2.z(),
                 r3.x(), r3.y(), r3.z());


    return true;
}
                                                                                                                                               