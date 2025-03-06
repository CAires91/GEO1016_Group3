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
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx,  /// output: focal length (i.e., K[0][0]).
        double& fy,  /// output: focal length (i.e., K[1][1]).
        double& cx,  /// output: x component of the principal point (i.e., K[0][2]).
        double& cy,  /// output: y component of the principal point (i.e., K[1][2]).
        double& s,   /// output: skew factor (i.e., K[0][1]), which is s = -alpha * cot(theta).
        Matrix33& R, /// output: the 3x3 rotation matrix encoding camera rotation.
        Vector3D& t) /// output：a 3D vector encoding camera translation.
{
    std::cout << "\nTODO: implement the 'calibration()' function in the file 'Calibration/calibration_method.cpp'\n\n";

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
                 "\t    - ensure that your code compiles and can reproduce your results WITHOUT ANY modification.\n\n" << std::flush;

    /// Below are a few examples showing some useful data structures and functions.

    // This is a 1D array of 'double' values. Alternatively, you can use 'double mat[25]' but you cannot change it
    // length. With 'std::vector', you can append/delete/insert elements, and much more. The 'std::vector' can store
    // not only 'double', but also any other types of objects. In case you may want to learn more about 'std::vector'
    // check here: https://en.cppreference.com/w/cpp/container/vector
    std::vector<double> array = {1, 3, 3, 4, 7, 6, 2, 8, 2, 8, 3, 2, 4, 9, 1, 7, 3, 23, 2, 3, 5, 2, 1, 5, 8, 9, 22};
    array.push_back(5); // append 5 to the array (so the size will increase by 1).
    array.insert(array.end(), 10, 3);  // append ten 3 (so the size will grow by 10).

    /// To access the value of an element.
    // double a = array[2];

    /// define a 2D vector/point
    // Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    // Vector3D c(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    // Vector2D p = c.cartesian();

    /// get the Homogeneous coordinates of p
    // Vector3D q = p.homogeneous();

    /// the length of a vector
    // double len = p.length();
    /// the squared length of a vector
    // double sqr_len = p.length2();

    /// the dot product of two vectors
    // double dot_prod = dot(p, q);

    /// the cross product of two vectors
    // Vector cross_prod = cross(c, q);

    /// normalize this vector
    // cross_prod.normalize();

    // Define an m-by-n double valued matrix.
    // Here I use the above array to initialize it. You can also use A(i, j) to initialize/modify/access its elements.
    // const int m = 6, n = 5;
    // Matrix A(m, n, array.data());    // 'array.data()' returns a pointer to the array.
    // std::cout << "M: \n" << A << std::endl;

    /// define a 3 by 4 matrix (and all elements initialized to 0.0)
    // Matrix M(3, 4, 0.0);

    /// set first row by a vector
    // M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    // M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    // Matrix33 B;

    /// define and initialize a 3 by 3 matrix
    // Matrix33 T(1.1, 2.2, 3.3,
    //            0, 2.2, 3.3,
    //            0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    // Matrix34 P(1.1, 2.2, 3.3, 0,
               // 0, 2.2, 3.3, 1,
               // 0, 0, 1, 1);

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    // Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
    // W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    // int num_rows = W.rows();

    /// get the number of columns.
    // int num_cols = W.cols();

    /// get the the element at row 1 and column 2
    // double value = W(1, 2);

    /// get the last column of a matrix
    // Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    // Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    // Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4
    //
    // Matrix U(m, m, 0.0);   // initialized with 0s
    // Matrix S(m, n, 0.0);   // initialized with 0s
    // Matrix V(n, n, 0.0);   // initialized with 0s

    // Compute the SVD decomposition of A
    // svd_decompose(A, U, S, V);

    // Now let's check if the SVD result is correct

    // Check 1: U is orthogonal, so U * U^T must be identity
    // std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;

    // Check 2: V is orthogonal, so V * V^T must be identity
    // std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;

    // Check 3: S must be a diagonal matrix
    // std::cout << "S: \n" << S << std::endl;

    // Check 4: according to the definition, A = U * S * V^T
    // std::cout << "M - U * S * V^T: \n" << A - U * S * transpose(V) << std::endl;

    // Compute the inverse of a matrix
    // Matrix invT;
    // inverse(T, invT);
    // Let's check if the inverse is correct
    // std::cout << "T * invT: \n" << T * invT << std::endl;

    // TODO: the above code just demonstrates some useful data structures and APIs. Please remove all above code in your
    //       final submission.

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

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

    // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    int num_2d_points = points_2d.size();  // Get the number of 2d points
    int num_3d_points = points_3d.size();  // Get the number of 3d points

    if (num_2d_points < 6 || num_3d_points < 6 || num_2d_points != num_3d_points) {
        std::cerr << "Error: Input is invalid. At least 6 correspondences are required" << std::endl;
        return false;
    }

    std::cout << "Input is valid. Proceed with calibration..." << std::endl;

    // TODO: construct the P matrix (so P * m = 0).
    int num_rows_P = 2*num_3d_points;
    int num_cols_P = 12;

    Matrix P(num_rows_P, num_cols_P, 0.0); // Initialize Matrix P as a 2n x 12 matrix filled with zeroes
    for (int i = 0; i < num_3d_points; i++) {
        // Get 3d points (in homogeneous coordinates)
        Vector3D P_i = points_3d[i];    // Use vector3D class
        double X = P_i.x();             // Use x() function
        double Y = P_i.y();
        double Z = P_i.z();

        // Get the 2d points
        Vector2D p_i = points_2d[i];
        double u = p_i.x();
        double v = p_i.y();

        // First row of matrix for the correspondence points
        P.set_row(2*i, {X, Y, Z, 1, 0, 0, 0, 0, -u*X, -u*Y, -u*Z, -u});

        // Second row of matrix for the correspondence points
        P.set_row(2*i+1, {0, 0, 0, 0, X, Y, Z, 1, -v*X, -v*Y, -v*Z, -v});
    }

    // Print the P matrix
    std::cout << "P matrix:\n" << P << std::endl;



    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.
    const int m = 2*num_3d_points, n = 12;  // m is rows of the P matrix, n is the columns of P matrix

    Matrix U(m, m, 0.0);   // initialized with 0s, check slide 31 about SVD for matrix dimension
    Matrix S(m, n, 0.0);   // initialized with 0s
    Matrix V(n, n, 0.0);   // initialized with 0s

    // Compute the SVD decomposition of A
    svd_decompose(P, U, S, V);

    std::cout << "Check 1: U is orthogonal, so U * U^T must be identity" << std::endl;
    std::cout << "U * U^T:\n" << U*transpose(U) << std::endl;

    std::cout << "Check 2: V is orthogonal, so V * V^T must be identity" << std::endl;
    std::cout << "V * V^T:\n" << V*transpose(V) << std::endl;

    std::cout << "Check 3: S is a diagonal matrix" << std::endl;
    std::cout << "S: \n " << S << std::endl;

    std::cout << "Check 4: according to the definition, P = U * S * V^T" << std::endl;
    std::cout << "P - U * S * V^T: \n" << P - U * S * transpose(V) << std::endl;

    // Get the last column of the V matrix to obtain "m" value
    std::cout << "V:\n" << V << std::endl;

    Vector m_vector = V.get_column(V.cols() - 1);

    // build m matrix dimension 3x4
    Matrix m_matrix(3,4,0.0);

    int idx = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            m_matrix[i][j] = m_vector[idx++];
        }
    }

    std::cout << "m_matrix:\n" << m_matrix << std::endl;

    // ---------------------------------------------------------CHECK-------------------------------------------------

    // get va3 from m_matrix:

    double a3_1 = m_matrix(2, 0);
    double a3_2 = m_matrix(2, 1);
    double a3_3 = m_matrix(2, 2);
    Vector3D a3 = Vector3D(a3_1, a3_2, a3_3);
    std::cout << "vector a3: \n" << a3 << std::endl;

    // TODO: Handle rho

    double rho = 1/norm(a3);

    Matrix M = (rho*-1)*m_matrix;

    std::cout << "M: \n" << M << std::endl;

    // // cross-check points between 3d and 2d corresponding points. by multiplying m_matrix

    for (const Vector3D& points: points_3d) {
        std::cout << "\nPoint: \n" << points << std::endl;

        Matrix point_test_matrix(4, 1, 1.0);
        // slide 24, add 1 more row so it's matrix dimension is match with the m_matrix dimension in order to be able for multiplication

        for (int i = 0; i < 3; i++) {
            point_test_matrix(i, 0) = points[i];
        }

        Matrix s_pi = M * point_test_matrix;

        double su = s_pi(0, 0);
        double sv = s_pi(1, 0);
        double sc = s_pi(2, 0);

        double x = su / sc;
        double y = sv / sc;

        std::cout << "Computed image coordinates-x: " << x << std::endl;
        std::cout << "Computed image coordinates-y: " << y << std::endl;
    }

    // TODO: extract intrinsic parameters from M.
    // Rho
    std::cout << "rho: " << rho << std::endl;

    // cx
    double a1_1 = m_matrix(0,0);
    double a1_2 = m_matrix(0,1);
    double a1_3 = m_matrix(0,2);
    Vector3D a1 = Vector3D(a1_1, a1_2, a1_3);
    std::cout << "vector a1:" << a1 << std::endl;

    double a2_1 = m_matrix(1,0);
    double a2_2 = m_matrix(1,1);
    double a2_3 = m_matrix(1,2);
    Vector3D a2 = Vector3D(a2_1, a2_2, a2_3);
    std::cout << "vector a2:" << a2 << std::endl;
    std::cout << "rho: " << rho << std::endl;

    // calculate cx = rho^2 *(a1 * a3)
    cx = (rho * rho) * dot(a1, a3);
    std::cout << "cx: " << cx << std::endl;

    // ---------------------------------------------------------CHECKCHECKCHECK-------------------------------------------------

    // cy
    // calculate cy = rho^2 *(a2 * a3)
    cy = (rho * rho) * dot(a2, a3);
    std::cout << "cy: " << cy << std::endl;

    // cos(theta)
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

    // alpha - fx
    // calculate alpha= rho^2 * ||a1*a3|| * sin(ø)

    // double sin_theta = sqrt(std::max(0, 1 - cos_theta * cos_theta)); // Ensure numerical stability
    double sin_theta = sqrt(1 - cos_theta * cos_theta);
    std::cout << "sin_theta: " << sin_theta << std::endl;

    double alpha = (rho * rho) * norm_cross_a1_a3 * sin_theta;
    std::cout << "alpha: " << alpha << std::endl;

    // beta - fy
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



    // TODO: extract extrinsic parameters from M.
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

    // r1
    // calculate r1 = (a2*a3)/||a2*a3||
    Vector3D r1 = cross_a2_a3 / norm_cross_a2_a3;
    std::cout << "r1: " << r1 << std::endl;

    // r3
    // calculate r3= rho * a3
    Vector3D r3 = rho * a3;
    std::cout << "r3: " << r3 << std::endl;

    // r2
    // calculate r2 = r3 x r1
    Vector3D r2 = cross(r3, r1);
    std::cout << "r2: " << r2 << std::endl;

    // calculate R
    R = Matrix33(r1.x(), r1.y(), r1.z(),
                 r2.x(), r2.y(), r2.z(),
                 r3.x(), r3.y(), r3.z());
    std::cout << "R: " << R << std::endl;

    // t??

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

//     // TODO: make sure the recovered parameters are passed to the corresponding variables (fx, fy, cx, cy, s, R, and t)
//
//     std::cout << "\n\tTODO: After you implement this function, please return 'true' - this will trigger the viewer to\n"
//                  "\t\tupdate the rendering using your recovered camera parameters. This can help you to visually check\n"
//                  "\t\tif your calibration is successful or not.\n\n" << std::flush;
//     return false;
// }

















