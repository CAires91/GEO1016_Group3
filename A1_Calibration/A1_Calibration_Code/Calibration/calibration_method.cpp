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


/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters are returned by fx, fy, cx, cy, skew, R, and t).
 */
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

    std::cout << "[Liangliang]:\n"
            "\tIn this assignment, two essential data structures, 'Matrix' and 'Vector', are provided for the\n"
            "\tmanipulation and storage of matrices and vectors. These data structures are defined in:\n"
            "\t    - Calibration/matrix.h: handles matrices of arbitrary dimensions and related functions.\n"
            "\t    - Calibration/vector.h: manages vectors of arbitrary sizes and related functions.\n"
            "\tCamera calibration requires computing the SVD and inverse of matrices. These functions, along\n"
            "\twith several other relevant ones, are provided in:\n"
            "\t    - Calibration/matrix_algo.h: contains functions for determinant, inverse, SVD, linear least-squares...\n"
            "\tIn the 'Calibration::calibration(...)' function, code snippets are provided for your reference.\n"
            "\tFor more details about these data structures and a complete list of related functions, please\n"
            "\trefer to the header files mentioned above.\n\n"
            "\tFor your final submission, adhere to the following guidelines:\n"
            "\t    - submit ONLY the 'Calibration/calibration_method.cpp' file.\n"
            "\t    - remove ALL unrelated test code, debugging code, and comments.\n"
            "\t    - ensure that your code compiles and can reproduce your results WITHOUT ANY modification.\n\n" <<
            std::flush;

    std::cout << "\n[Liangliang]:\n"
            "\tThis function takes two arrays as input parameters:\n"
            "\t\t- points_3d: An array of 3D points representing the scene\n"
            "\t\t- points_2d: An array of 2D image points corresponding to the 3D points\n"
            "\tThe function should return either 'true' upon successful calibration or 'false' otherwise.\n"
            "\tUpon success, the following parameters must be stored in the specified variables:\n"
            "\t\t- fx and fy: focal lengths along the x and y axes, respectively\n"
            "\t\t- cx and cy: coordinates of the principal point\n"
            "\t\t- s: the skew factor, i.e., s = -alpha * cot(theta)\n"
            "\t\t- R: the 3x3 rotation matrix encoding camera orientation\n"
            "\t\t- t: a 3D vector encoding camera location.\n"
            "\tIMPORTANT: don't forget to write your recovered parameters to the above variables." << std::endl;

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
    std::cout << "P: \n" << P << std::endl;

    // ----------------------- solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.

    const int m = 2 * num_3d_points, n = 12;

    Matrix U(m, m, 0.0); // initialized with 0s
    Matrix S(m, n, 0.0); // initialized with 0s
    Matrix V(n, n, 0.0); // initialized with 0s

    // Compute the SVD decomposition of A
    svd_decompose(P, U, S, V);

    // Check 1: U is orthogonal, so U * U^T must be identity
    std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;

    // Check 2: V is orthogonal, so V * V^T must be identity
    std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;

    // Check 3: S must be a diagonal matrix
    std::cout << "S: \n" << S << std::endl;

    // Check 4: according to the definition, P = U * S * V^T
    std::cout << "M - U * S * V^T: \n" << P - U * S * transpose(V) << std::endl;

    /// get the last column of the V matrix

    std::cout << "V: \n" << V << std::endl;

    Vector m_vector = V.get_column(V.cols() - 1);

    /// build m matrix

    Matrix m_matrix(3, 4, 0.0);

    int idx = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            m_matrix[i][j] = m_vector[idx++];
        }
    }

    std::cout << "m_matrix: \n" << m_matrix << std::endl;


    // get a3 from m_matrix:

    double a3_1 = m_matrix(2, 0);

    double a3_2 = m_matrix(2, 1);

    double a3_3 = m_matrix(2, 2);

    Vector3D a3 = Vector3D(a3_1, a3_2, a3_3);

    std::cout << "a3: \n" << a3 << std::endl;

    double rho = 1 / norm(a3);

    Matrix M = rho * m_matrix;

    std::cout << "M: \n" << M << std::endl;

    // ---------------------- Optional: check if your M is correct by applying M on the 3D points. If correct, the projected point
    // ----------------------- should be very close to your input images points.
    std::cout << "Verifying the intermediate results: applying M on the 3D points" << std::endl;

    for (size_t i = 0; i < points_3d.size(); ++i) {
        const Vector3D &point = points_3d[i];
        const double x_input = points_2d[i][0];  // Access corresponding 2D point x
        const double y_input = points_2d[i][1];  // Access corresponding 2D point y

        Matrix point_test_matrix(4, 1, 1.0);

        for (int i = 0; i < 3; i++) {
            point_test_matrix(i, 0) = point[i];
        }

        Matrix s_pi = M * point_test_matrix;

        double su = s_pi(0, 0);
        double sv = s_pi(1, 0);
        double sc = s_pi(2, 0);

        double x = su / sc;
        double y = sv / sc;


        std::cout << "Input 2D Point: (" << x_input << ", " << y_input << ")" << std::endl;
        std::cout << "Computed 2D Point: (" << x << ", " << y << ")\n" << std::endl;

    }

    // ---------------------- extract intrinsic parameters from M.

    // calculate a1 using vector3d
    double a1_1 = m_matrix(0, 0);
    double a1_2 = m_matrix(0, 1);
    double a1_3 = m_matrix(0, 2);
    Vector3D a1 = Vector3D(a1_1, a1_2, a1_3);
    std::cout << "a1: " << a1 << std::endl;

    // calculate a2 using vector3d
    double a2_1 = m_matrix(1, 0);
    double a2_2 = m_matrix(1, 1);
    double a2_3 = m_matrix(1, 2);
    Vector3D a2 = Vector3D(a2_1, a2_2, a2_3);
    std::cout << "a2: " << a2 << std::endl;
    std::cout << "rho: " << rho << std::endl;

    // calculate cx = rho^2 *(a1 * a3)
    cx = (rho * rho) * dot(a1, a3);
    std::cout << "cx: " << cx << std::endl;

    // calculate cy = rho^2 *(a2 * a3)
    cy = (rho * rho) * dot(a2, a3);
    std::cout << "cy: " << cy << std::endl;

    // calculate cos(ø) = ((a1*a3)*(a2*a3))/(||a1*a3||*||a2*a3||)
    Vector3D cross_a1_a3 = cross(a1, a3);
    Vector3D cross_a2_a3 = cross(a2, a3);

    double norm_cross_a1_a3 = norm(cross_a1_a3);
    double norm_cross_a2_a3 = norm(cross_a2_a3);

    std::cout << "cross_a1_a3: " << cross_a1_a3 << std::endl;
    std::cout << "cross_a2_a3: " << cross_a2_a3 << std::endl;

    std::cout << "norm_cross_a1_a3: " << norm_cross_a1_a3 << std::endl;
    std::cout << "norm_cross_a2_a3: " << norm_cross_a2_a3 << std::endl;

    double cos_theta = -1 * dot(cross_a1_a3, cross_a2_a3) / (norm_cross_a1_a3 * norm_cross_a2_a3);

    // Clamp cos_theta to the range [-1, 1]
    if (cos_theta > 1.0) {
        std::cout << "cos_theta > 1.0" << std::endl;
        cos_theta = 1.0;
    } else if (cos_theta < -1.0) {
        std::cout << "cos_theta < 1.0" << std::endl;
        cos_theta = -1.0;
    }

    std::cout << "cos_theta: " << cos_theta << std::endl;

    // calculate alpha= rho^2 * ||a1*a3|| * sin(ø)

    // double sin_theta = sqrt(std::max(0, 1 - cos_theta * cos_theta)); // Ensure numerical stability
    double sin_theta = sqrt(1 - cos_theta * cos_theta);
    std::cout << "sin_theta: " << sin_theta << std::endl;

    double alpha = (rho * rho) * norm_cross_a1_a3 * sin_theta;
    std::cout << "alpha: " << alpha << std::endl;


    // calculate beta = rho^2 * ||a2*a3|| * sin(ø)

    double beta = (rho * rho) * norm_cross_a2_a3 * sin_theta;
    std::cout << "beta: " << beta << std::endl;

    // calculate fx = alpha
    fx = alpha;
    std::cout << "fx: " << fx << std::endl;

    // calculate fy = beta / sin(ø)

    fy = beta / sin_theta;
    std::cout << "fy: " << fy << std::endl;

    // calculate s = -alpha * cot(ø)
    s = -alpha / tan(acos(cos_theta));

    std::cout << "s: " << s << std::endl;

    // ----------------------  extract extrinsic parameters from M.

    // extract b
    double b1 = m_matrix(0, 3);
    double b2 = m_matrix(1, 3);
    double b3 = m_matrix(2, 3);
    Vector3D b = Vector3D(b1, b2, b3);
    std::cout << "b: " << b << std::endl;

    Matrix33 K(fx, s, cx,
               0, fy, cy,
               0, 0, 1);
    std::cout << "K: " << K << std::endl;

    // calculate t= rho*K-1*b
    t = rho * inverse(K) * b;
    std::cout << "t: " << t << std::endl;

    // if t(2) < 0, then rho = -rho
    if (t.z() < 0) {
        rho = -rho;
        t = rho * inverse(K) * b;
    }

    // calculate r1 = (a2*a3)/||a2*a3||
    Vector3D r1 = cross_a2_a3 / norm_cross_a2_a3;
    std::cout << "r1: " << r1 << std::endl;

    // calculate r3= rho * a3
    Vector3D r3 = rho * a3;
    std::cout << "r3: " << r3 << std::endl;

    // calculate r2 = r3 x r1
    Vector3D r2 = cross(r3, r1);
    std::cout << "r2: " << r2 << std::endl;

    // calculate R
    R = Matrix33(r1.x(), r1.y(), r1.z(),
                 r2.x(), r2.y(), r2.z(),
                 r3.x(), r3.y(), r3.z());
    std::cout << "R: " << R << std::endl;


    // ----------------------Optional: test if these values are correct by feeding it a 3D point and checking if the 2D point is correct p=K[R t]P
    std::cout << "Verifying the intermediate results: applying K[R t] on the 3D points" << std::endl;

    Matrix34 Rt_test(R(0, 0), R(0, 1), R(0, 2), t.x(),
             R(1, 0), R(1, 1), R(1, 2), t.y(),
             R(2, 0), R(2, 1), R(2, 2), t.z());


    for (size_t i = 0; i < points_3d.size(); ++i) {
        const Vector3D &test_point = points_3d[i];
        const double x_input = points_2d[i][0];  // Access corresponding 2D point x
        const double y_input = points_2d[i][1];  // Access corresponding 2D point y

        Vector3D s_pi2 = (K * Rt_test) * Vector4D(test_point.x(), test_point.y(), test_point.z(), 1);


        double su2 = s_pi2.x();
        double sv2 = s_pi2.y();
        double sc2 = s_pi2.z();

        double x_rt = su2 / sc2;
        double y_rt = sv2 / sc2;

        std::cout << "Input 2D Point: " << x_input << ", " << y_input <<  std::endl;
        std::cout << "Computed 2D Point: " << x_rt << ", " << y_rt <<  std::endl;
    }

    return true;
}
                                                                                                                                               