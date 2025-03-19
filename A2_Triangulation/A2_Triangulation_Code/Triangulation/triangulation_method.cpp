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

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;

Matrix33 normalize_points(
    const std::vector<Vector2D> &points,
    std::vector<Vector3D> &norm_points
    ) {
    int N = points.size();
    if (N < 8) {
        std::cerr << "At least 8 points are required for normalization." << std::endl;
        return Matrix33();
    }
    // Step 1.1 Normalize the image points
    // Step 1.1.1 Compute the centroid (mean) of the image points (translation)
    double mean_x = 0, mean_y = 0;     // initialize the mean of x,y coordinates of the image points
    for (const auto &p : points) {
        mean_x += p.x();    // sum of x coordinates of the image points
        mean_y += p.y();    // sum of y coordinates of the image points
    }
    mean_x /= N;    // mean of x,y coordinates of the image points
    mean_y /= N;

    // Step 1.1.2 Compute the scale factor of the image points (scaling or average distance from origin)
    double scale = 0;   // initialize the scale factor of the image points
    for (const auto &p : points) {
        scale += sqrt(pow(p.x()-mean_x,2) + pow(p.y()-mean_y,2)); // compute distance from the centroid to the image points (euclidean distance)
    }
    scale = sqrt(2)/(scale/N);    // compute the scale factor of the 1st image points

    // Step 1.1.3 Construct Normalization Matrix T
    Matrix33 T(scale, 0, -scale*mean_x,
                0, scale, -scale*mean_y,
                0, 0, 1);   // normalization matrix T for the 2nd image points

    // Step 1.1.4 Normalize the image points
    norm_points.clear();     // clear the normalized points
    for (const auto &p : points) {
        norm_points.push_back(T * p.homogeneous());    // p' = T*p (convert to homogeneous coordinates)
    }
    return T;
}


// TODO : - estimate the fundamental matrix F; (Step 1)
Matrix33 compute_fundamental_matrix(
    const std::vector<Vector2D> &points_0,
    const std::vector<Vector2D> &points_1
    ) {

    int N = points_0.size();
    if (N < 8) {
        std::cerr << "At least 8 points are required for computing the fundamental matrix." << std::endl;
        return Matrix33();  // return empty matrix
    }

    // Step 1.1 Normalize the image points using function normalize_points
    std::vector<Vector3D> norm_points_0, norm_points_1;
    Matrix33 T0 = normalize_points(points_0, norm_points_0);    // normalize the 1st image points
    Matrix33 T1 = normalize_points(points_1, norm_points_1);    // normalize the 2nd image points
    std::cout << "Normalization matrix T0: \n" << T0 << std::endl;
    std::cout << "Normalization matrix T1: \n" << T1 << std::endl;


    // Step 1.2 Construct Matrix W
    Matrix W(N, 9, 0.0);   // initialize matrix W with N rows, 9 columns, filled with 0s
    for (int i = 0; i < N; i++) {
        double u0 = norm_points_0[i].x(), v0 = norm_points_0[i].y();
        double u1 = norm_points_1[i].x(), v1 = norm_points_1[i].y();

        W.set_row(i, {u0*u1, v0*u1, u1, u0*v1, v0*v1, v1, u0, v1, 1}); // fill the i-th row of matrix W
    }
    std::cout << "Matrix W has number of rows: " << W.rows() << " number of columns: " << W.cols() << std::endl;
    std::cout << "First row of matrix W: \n" << W.get_row(0) << std::endl;
    std::cout << "Last row (" << N << ") of matrix W: \n" << W.get_row(N-1) << std::endl;

    // Step 1.3 Solve for f in Wf = 0 using SVD
    int m = W.rows(), n = W.cols();
    Matrix U(m, m, 0.0); // initialize with 0s -> orthogonal matrix, left singular vectors
    Matrix S(m, n, 0.0); // initialize with 0s -> diagonal matrix, singular values
    Matrix V(n, n, 0.0); // initialize with 0s -> orthogonal matrix, right singular vectors, last column gives F


    svd_decompose(W, U, S, V);
    // std::cout << "Matrix U: \n" << U << std::endl; // do not print, too many lines!
    std::cout << "Matrix S: \n" << S << std::endl;
    std::cout << "Matrix V: \n" << V << std::endl;
    Vector F_vector = V.get_column(V.cols()-1);  // get the last column of V matrix
    std::cout << "Fundamental matrix F in vector form: \n" << F_vector << std::endl;  // F in vector form
    std::cout << "F_vector size: " << F_vector.size() << std::endl;


    // Step 1.4 Enforce the rank-2 constraint on F
    // Step 1.4.1 Reshape vector F to 3x3 matrix
    Matrix33 F_matrix;
    int idx = 0;  // Reset idx to 0
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            F_matrix(i, j) = F_vector[idx++];  // Fill in row-major order
        }
    }
    std::cout << "Fundamental matrix F in Matrix33 form: \n" << F_matrix << std::endl;  // F in Matrix33 form

    // Step 1.4.2 Perform SVD on F_matrix
    Matrix U_2(3,3,0.0), S_2(3,3,0.0), V_2(3,3,0.0);
    svd_decompose(F_matrix, U_2, S_2, V_2);

    // set the last/third singular value to 0
    S_2(2, 2) = 0;
    F_matrix = U_2 * S_2 * V_2.transpose();
    std::cout << "F_matrix33 after rank-2 constraint: \n" << F_matrix << std::endl;

    // Convert 'F_matrix' back to 'F_vector_new' (9x1)
    Vector F_vector_new(9);
    idx = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            F_vector_new[idx++] = F_matrix(i,j);
        }
    }
    std::cout << "F matrix new after rank-2 constraint in vector form \n" << F_vector_new << std::endl;

    // Step 1.5 Denormalize F
    // Step 1.5.1 Convert F_vector_new back to Matrix33 form
    Matrix33 F_denorm;
    idx = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            F_denorm(i,j) = F_vector_new[idx++];
        }
    }
    // Step 1.5.2 Apply denormalization F
    F_denorm = T1.transpose()*F_denorm*T0;
    std::cout << "F denormalized in Matrix33 form: \n" << F_denorm << std::endl;

    return F_denorm;    // return the denormalized fundamental matrix F in Matrix33 form
}





/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        double s,                 /// input: the skew factor (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       'triangulation()'.

    std::cout << "\nTODO: implement the 'triangulation()' function in the file 'Triangulation/triangulation_method.cpp'\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tSimilar to the first assignment, basic linear algebra data structures and functions are provided in\n"
                 "\tthe following files:\n"
                 "\t    - Triangulation/matrix.h: handles matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/vector.h: manages vectors of arbitrary sizes and related functions.\n"
                 "\t    - Triangulation/matrix_algo.h: contains functions for determinant, inverse, SVD, linear least-squares...\n"
                 "\tFor more details about these data structures and a complete list of related functions, please\n"
                 "\trefer to the header files mentioned above.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
                 "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
                 "\tFor your final submission, adhere to the following guidelines:\n"
                 "\t    - submit ONLY the 'Triangulation/triangulation_method.cpp' file.\n"
                 "\t    - remove ALL unrelated test code, debugging code, and comments.\n"
                 "\t    - ensure that your code compiles and can reproduce your results WITHOUT ANY modification.\n\n" << std::flush;

    /// Below are a few examples showing some useful data structures and APIs.

    /// define a 2D vector/point
    // Vector2D b(1.1, 2.2);
    //
    // /// define a 3D vector/point
    // Vector3D a(1.1, 2.2, 3.3);
    //
    // /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    // Vector2D p = a.cartesian();
    //
    // /// get the Homogeneous coordinates of p
    // Vector3D q = p.homogeneous();
    //
    // /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    // Matrix33 A;
    //
    // /// define and initialize a 3 by 3 matrix
    // Matrix33 T(1.1, 2.2, 3.3,
    //            0, 2.2, 3.3,
    //            0, 0, 1);
    //
    // /// define and initialize a 3 by 4 matrix
    // Matrix34 M(1.1, 2.2, 3.3, 0,
    //            0, 2.2, 3.3, 1,
    //            0, 0, 1, 1);
    //
    // /// set first row by a vector
    // M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));
    //
    // /// set second column by a vector
    // M.set_column(1, Vector3D(5.5, 5.5, 5.5));
    //
    // /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    // Matrix W(15, 9, 0.0);
    // /// set the first row by a 9-dimensional vector
    // W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>
    //
    // /// get the number of rows.
    // int num_rows = W.rows();
    //
    // /// get the number of columns.
    // int num_cols = W.cols();
    //
    // /// get the the element at row 1 and column 2
    // double value = W(1, 2);
    //
    // /// get the last column of a matrix
    // Vector last_column = W.get_column(W.cols() - 1);
    //
    // /// define a 3 by 3 identity matrix
    // Matrix33 I = Matrix::identity(3, 3, 1.0);
    //
    // /// matrix-vector product
    // Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).
    std::cout << "points_0 has " << points_0.size() << " points" << std::endl;
    std::cout << "points_1 has " << points_1.size() << " points" << std::endl;

    // Matching number of points check
    if (points_0.size() != points_1.size()) {
        std::cerr <<"Error: Mismatch of the number of 2D points from the inputs\n" << std::endl;
        return false;
    }

    // Minimum number of points check
    if (points_0.size()<8) {
        std::cerr <<"Error: Not enough pairs of 2D points. At least 8 pairs are needed" << std::endl;
        return false;
    }


    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;    -> Step #1.1 Assignment, lines 31
    //      - compute the essential matrix E;       ->
    //      - recover rotation R and t.             -> Step #2 Assignment
    Matrix33 F = compute_fundamental_matrix(points_0, points_1);
    if (F.trace() == 0) {
        std::cerr << "Error: Fundamental matrix computation failed" << std::endl;
        return false;

    }


    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    return points_3d.size() > 0;
}