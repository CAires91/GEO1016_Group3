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
#include <tuple>


using namespace easy3d;

void normalization(
    const std::vector<Vector2D> &points,
    std::vector<Vector2D> &points_q, /// output: points transformed
    Matrix33 &T /// output: 3 by 3 transformation matrix
    ) {

    // calculate size of points
    const int num_points = points.size();

    double sum_x = 0.0;
    double sum_y = 0.0;

    // calculate centroid of points_0

    for (int i = 0; i < num_points; i++) {
        sum_x += points[i][0];
        sum_y += points[i][1];
    }
    Vector2D centroid_points;
    centroid_points = Vector2D(sum_x / num_points, sum_y / num_points);

    // translate by the centroid
    double tx = centroid_points[0];
    double ty = centroid_points[1];

    // scale
    double scale = 1;


    double sum_distance(0.0);

    // Accumulate Euclidean distances
    for (const auto &point : points) {
        double dx = point[0] - centroid_points[0];
        double dy = point[1] - centroid_points[1];
        sum_distance += sqrt(dx * dx  + dy * dy); // Accumulate Euclidean distances
    }

    // Calculate the mean distance from the centroid
    double mean_distance = sum_distance/points.size();

    // Calculate the scale factor to normalize the mean distance to sqrt(2)
    scale = sqrt(2)/mean_distance;

    // T
    T = Matrix33(scale, 0, -scale * tx,
                 0, scale, -scale *ty,
                  0, 0, 1);

    // transformation points
    points_q.clear(); // Clear the output vector before adding new points

    for (const auto& p : points) {
        Vector3D p_hom = p.homogeneous();
        Vector3D transformed_point = T * p_hom;
        Vector2D transformed_point2 = transformed_point.cartesian();
        points_q.push_back(transformed_point2);
    }

    // Check if the number of points is the same
    if (points_q.size() != points.size()) {
        throw std::runtime_error("Normalization failed: mismatch in point count.");
    }
}

void construct_W_matrix(
    const std::vector<Vector2D> &points_q_img1, /// input:
    const std::vector<Vector2D> &points_q_img2, /// input:
    Matrix &W /// output:
    ) {

    for (int i = 0; i < points_q_img1.size(); i++) {

        W.set_row(i, {
            (points_q_img1[i].x() * points_q_img2[i].x()),
            (points_q_img1[i].y() * points_q_img2[i].x()),
            points_q_img2[i].x(),
            (points_q_img1[i].x() * points_q_img2[i].y()),
            (points_q_img1[i].y() * points_q_img2[i].y()),
            points_q_img2[i].y(),
            points_q_img1[i].x(),
            points_q_img1[i].y(),
            1});
    }

}

// Function to reconstruct a 3D point
Vector3D reconstruct3D(const Matrix34& M, //M is the projection matrix for the second camera
                       const Matrix& K,
                       const Vector2D& point_0,
                       const Vector2D& point_1){

    // Identity R | t, t = 0,0,0
    Matrix34 Rt0(1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0);

    //M0 is the projection matrix for the first camera
    Matrix34 M0 = K * Rt0;

    Matrix A(4,4,{(point_0[0]*M0(2, 0))-M0(0, 0), (point_0[0]*M0(2, 1))-M0(0, 1), (point_0[0]*M0(2, 2))-M0(0, 2), (point_0[0]*M0(2, 3))-M0(0, 3),
                                  (point_0[1]*M0(2, 0))-M0(1, 0), (point_0[1]*M0(2, 1))-M0(1, 1), (point_0[1]*M0(2, 2))-M0(1, 2), (point_0[1]*M0(2, 3))-M0(1, 3),
                                   (point_1[0]*M(2, 0))-M(0, 0),  (point_1[0]*M(2, 1))-M(0, 1),   (point_1[0]*M(2, 2))-M(0, 2),   (point_1[0]*M(2, 3))- M(0, 3),
                                   (point_1[1]*M(2, 0))-M(1, 0),  (point_1[1]*M(2, 1))-M(1, 1),   (point_1[1]*M(2, 2))-M(1, 2),   (point_1[1]*M(2, 3))- M(1, 3)});

    // SVD for A

    const int m = 4, n = 4;

    Matrix U(m, m, 0.0); // initialized with 0s
    Matrix S(m, n, 0.0); // initialized with 0s
    Matrix V(n, n, 0.0); // initialized with 0s


    // compute the SVD decomposition of A

    svd_decompose(A, U, S, V);

    Vector4D P_4D = V.get_column(V.cols() - 1);

    Vector3D P = P_4D.cartesian();
    return P;

}


// Function to count points in front of both cameras
int countPointsInFront(const std::vector<Vector2D>& points_0,
                       const std::vector<Vector2D>& points_1,
                       const Matrix& M,
                       const Matrix33& R,
                       const Vector3D& t,
                       std::vector<Vector3D>& temp_points_3d,
                       const Matrix& K) {
    int count = 0;
    temp_points_3d.clear(); // Clear previous 3D points

    for (size_t i = 0; i < points_0.size(); i++) {
        // Reconstruct 3D point from first camera
        Vector3D P = reconstruct3D(M, K, points_0[i], points_1[i]);


        // Check if P is in front of camera 1
        if (P.z() > 0) {
            // Convert to second camera's coordinate system: Q = R * P + t
            Vector3D Q = R * P + t;
            // Check if Q is in front of camera 2
            if (Q.z() > 0) {
                count++;
                temp_points_3d.push_back(P); // Save valid 3D point

                }
        }
    }

    return count;
}

Vector2D project3DPoint(const Vector3D& P, const Matrix33& K, const Matrix34& M) {
    Vector4D P_hom = P.homogeneous();
    // std::cout<< "P_hom: " << P_hom << std::endl;

    // Apply the projection matrix
    Vector3D projected_img = M * P_hom;
    // std::cout<< "projected_img: " << projected_img << std::endl;
    // std::cout<< "original: " << P << std::endl;

    // Convert to 2D Euclidean coordinates
    return Vector2D(projected_img.cartesian());
}


std::tuple<double, double, double, double> computeReprojectionError(const std::vector<Vector2D>& points_0,
                                                    const std::vector<Vector2D>& points_1,
                                                    const std::vector<Vector3D>& points_3d,
                                                    const Matrix& K,
                                                    const Matrix34& M0,
                                                    const Matrix34& M) {
    double total_squared_error_0_x = 0.0;
    double total_squared_error_0_y = 0.0;
    double total_squared_error_1_x = 0.0;
    double total_squared_error_1_y = 0.0;
    size_t num_points = points_3d.size();

    for (size_t i = 0; i < num_points; ++i) {
        // Project 3D points back into image space
        Vector3D P = points_3d[i];
        Vector2D reprojected_0 = project3DPoint(P, K, M0);
        Vector2D reprojected_1 = project3DPoint(P, K, M);

        // Compute squared error for both views using direct multiplication
        Vector2D diff_0 = reprojected_0 - points_0[i];
        Vector2D diff_1 = reprojected_1 - points_1[i];

        double error_0_x = diff_0[0] * diff_0[0];
        double error_0_y = diff_0[1] * diff_0[1];
        double error_1_x = diff_1[0] * diff_1[0];
        double error_1_y= diff_1[1] * diff_1[1];

        total_squared_error_0_x += error_0_x;
        total_squared_error_0_y += error_0_y;
        total_squared_error_1_x += error_1_x;
        total_squared_error_1_y += error_1_y;

    }
    // total squared error by using pythagoras sqrt(x^2+y^2)= total error
    double total_squared_error_0 = total_squared_error_0_x + total_squared_error_0_y;
    double total_squared_error_1 = total_squared_error_1_x + total_squared_error_1_y;



    // Compute RMSE

    double rmse_0 = std::sqrt(total_squared_error_0 / ( num_points));
    double rmse_1 = std::sqrt(total_squared_error_1 / ( num_points));

    return std::make_tuple(rmse_0, total_squared_error_0, rmse_1, total_squared_error_1);
}



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



    // check if the input is valid (always good because you never known how others will call your function).
    int num_points0 = points_0.size();
    int num_points1 = points_1.size();

    if (num_points1 < 8 || num_points0 < 8 || num_points1 != num_points0) {
        std::cerr << "Error: Invalid number of points (less than 8 or mismatched points between images).\n";
        return false;
    }
    // if (points_0 == points_1) {
    //     std::cerr << "Error: Input points of both images are identical.\n";
    //     return false;
    // }

    // Step 1: Estimate relative pose of two views. This can be subdivided into

    //      - step 1.1: Normalize the input points.

    std::vector<Vector2D> points_q_img1;
    std::vector<Vector2D> points_q_img2;

    Matrix33 T_img1;
    Matrix33 T_img2;

    normalization(points_0, points_q_img1, T_img1);
    normalization(points_1, points_q_img2, T_img2);

    //      step 1.2: estimate the fundamental matrix F;

    //  Construct W
    int num_rows_img1 = points_0.size();
    Matrix W(num_rows_img1, 9, 0.0);
    construct_W_matrix(points_q_img1, points_q_img2, W);

    // SVD for W
    const int m = num_rows_img1, n = 9;
    Matrix U(m, m, 0.0);
    Matrix S(m, n, 0.0);
    Matrix V(n, n, 0.0);

    // compute the SVD decomposition of W
    svd_decompose(W, U, S, V);

    /// get the last column of the V matrix

    Vector f = V.get_column(V.cols() - 1);

    Matrix33 F(0.0);

    int idx = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            F[i][j] = f[idx++];
        }
    }

    const int m_2 = 3, n_2 = 3;

    Matrix U_2(m_2, m_2, 0.0);
    Matrix S_2(m_2, n_2, 0.0);
    Matrix V_2(n_2, n_2, 0.0);

    // SVD for ^ F
    svd_decompose(F, U_2, S_2, V_2);

    S_2(2,2) = 0;


    Matrix33 Fq;

    Fq = U_2 * S_2 * transpose(V_2);

    Matrix33 F_denormalized;

    F_denormalized = transpose(T_img2) * Fq * T_img1;

    double det_F = determinant(F_denormalized);
    if (std::abs(det_F) > 1e-6) {
        std::cerr << "Warning: Fundamental matrix determinant is not close to zero. Possible issue with computation.\n";
        std::cout<< "Determinant of F: " << det_F << std::endl;
    }


    Matrix33 K(fx, s, cx,
                0, fy, cy,
                0, 0, 1);

    //      compute the essential matrix E;

    Matrix33 E;

    E = transpose(K) * F_denormalized * K;

    //      recover rotation R and t.

    Matrix33 W_2(0.0, -1.0 , 0.0,
                1.0, 0.0, 0.0,
                0.0, 0.0, 1.0);

    Matrix33 Z(0.0, 1.0 , 0.0,
            -1.0, 0.0, 0.0,
            0.0, 0.0, 0.0);

    //      SVD of E

    constexpr int m_3 = 3, n_3 = 3;

    Matrix U_3(m_3, m_3, 0.0); // initialized with 0s
    Matrix S_3(m_3, n_3, 0.0); // initialized with 0s
    Matrix V_3(n_3, n_3, 0.0); // initialized with 0s

    svd_decompose(E, U_3, S_3, V_3); //@attention V is returned (instead of V^T)
    Vector u_3 = U_3.get_column(U_3.cols() - 1);
    Vector3D t_var1 = u_3;
    Vector3D t_var2 = -u_3;


    // R

    Matrix33 R_var1;
    R_var1 = determinant(U_3 * W_2 * transpose(V_3)) * U_3 * W_2 * transpose(V_3);
    Matrix33 R_var2;
    R_var2 = determinant(U_3 * transpose(W_2) * transpose(V_3)) * U_3 * transpose(W_2) * transpose(V_3);


    // Define the four possible Rt combinations
    Matrix34 Rt1(R_var1(0, 0), R_var1(0, 1), R_var1(0, 2), t_var1.x(),
                 R_var1(1, 0), R_var1(1, 1), R_var1(1, 2), t_var1.y(),
                 R_var1(2, 0), R_var1(2, 1), R_var1(2, 2), t_var1.z());

    Matrix34 Rt2(R_var1(0, 0), R_var1(0, 1), R_var1(0, 2), t_var2.x(),
                 R_var1(1, 0), R_var1(1, 1), R_var1(1, 2), t_var2.y(),
                 R_var1(2, 0), R_var1(2, 1), R_var1(2, 2), t_var2.z());

    Matrix34 Rt3(R_var2(0, 0), R_var2(0, 1), R_var2(0, 2), t_var1.x(),
                 R_var2(1, 0), R_var2(1, 1), R_var2(1, 2), t_var1.y(),
                 R_var2(2, 0), R_var2(2, 1), R_var2(2, 2), t_var1.z());

    Matrix34 Rt4(R_var2(0, 0), R_var2(0, 1), R_var2(0, 2), t_var2.x(),
                 R_var2(1, 0), R_var2(1, 1), R_var2(1, 2), t_var2.y(),
                 R_var2(2, 0), R_var2(2, 1), R_var2(2, 2), t_var2.z());

    // Compute the projection matrices
    Matrix M1 = K * Rt1;
    Matrix M2 = K * Rt2;
    Matrix M3 = K * Rt3;
    Matrix M4 = K * Rt4;

    // Define R and t variants
    std::vector<Matrix33> R_variants = {R_var1, R_var1, R_var2, R_var2};
    std::vector<Vector3D> t_variants = {t_var1, t_var2, t_var1, t_var2};
    std::vector<Matrix> M_matrices = {M1, M2, M3, M4};

    int max_count = 0;
    Matrix33 best_R;
    Vector3D best_t;
    Matrix best_M;
    std::vector<Vector3D> best_points_3d;

    // Iterate through all M matrices and count valid points
    for (int i = 0; i < 4; i++) {
        std::vector<Vector3D> temp_points_3d;
        int count = countPointsInFront(points_0, points_1, M_matrices[i], R_variants[i], t_variants[i], temp_points_3d, K);

        if (count > max_count) {
            max_count = count;
            best_R = R_variants[i];
            best_t = t_variants[i];
            best_M = M_matrices[i];
            best_points_3d = temp_points_3d;

        }
    }
    if (max_count < points_0.size()/2) {
        std::cerr << "Warning: Less than half of the points are in front of both cameras. Possible issue with reconstruction.\n";
        std::cout<< "Number of points in front: " << max_count << std::endl;
    }

    //check if the determinant of the rotation matrix is 1
    double detR = determinant(best_R);
    if (std::abs(detR - 1.0) > 1e-6) {
        std::cerr << "Warning: Rotation matrix determinant is not 1. Possible issue with decomposition.\n";
        std::cout<< "Determinant of R: " << detR << std::endl;
    }

    R = best_R;
    t = best_t;
    points_3d = best_points_3d;

    // Identity R | t, t = 0,0,0
    Matrix34 Rt0(1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0);

    Matrix33 K_real(1000, 0, 320,
                0, 1000, 240,
                0, 0, 1);
    //M0 is the projection matrix for the first camera

    Matrix34 M0 = K_real * Rt0;
    std::cout<< "M0: " << M0 << std::endl;

    // Compute RMSE and total squared error
    auto [rmse_0, total_error_0, rmse_1, total_error_1] = computeReprojectionError(points_0, points_1, points_3d, K_real, M0, best_M);

    // Print results
    std::cout << "Reprojection RMSE camera 1: " << rmse_0 << std::endl;
    std::cout << "Total Squared Error camera 1: " << total_error_0 << std::endl;

    std::cout << "Reprojection RMSE camera 2: " << rmse_1 << std::endl;
    std::cout << "Total Squared Error camera 2: " << total_error_1 << std::endl;


    return points_3d.size() > 0;
}
