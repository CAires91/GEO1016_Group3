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

void normalization(
    const std::vector<Vector2D> &points, 
    std::vector<Vector2D> &points_q,
    Matrix33 &T) {
    
    const int num_points = points.size();
    
    // calculate centroid of the points
    double sum_x = 0.0, sum_y = 0.0;
    for (int i = 0; i < num_points; i++) {
        sum_x += points[i][0];
        sum_y += points[i][1];
    }
    Vector2D centroid(sum_x / num_points, sum_y / num_points);

    // Accumulate Euclidean distances from de centroid_point
    double sum_distance = 0.0;
    for (const auto &point : points) {
        double dx = point[0] - centroid[0];
        double dy = point[1] - centroid[1];
        sum_distance += sqrt(dx * dx  + dy * dy);
    }

    // Calculate the scale factor
    double scale = sqrt(2) / (sum_distance / num_points);

    // Set the transformation matrix T
    T = Matrix33(scale, 0, -scale * centroid[0],
                 0, scale, -scale * centroid[1],
                 0, 0, 1);

    // Apply the transformation to the points
    points_q.clear(); 
    for (const auto& p : points) {
        Vector3D p_hom = p.homogeneous();
        Vector3D transformed_point = T * p_hom;
        points_q.push_back(transformed_point.cartesian());
    }

    // Check if the number of points remains unchanged
    if (points_q.size() != points.size()) {
        throw std::runtime_error("Normalization failed: mismatch in point count.");
    }
}

void construct_W_matrix(
    const std::vector<Vector2D> &points_q_img1, 
    const std::vector<Vector2D> &points_q_img2, 
    Matrix &W) {

for (size_t i = 0; i < points_q_img1.size(); ++i) {
    const auto& p1 = points_q_img1[i];
    const auto& p2 = points_q_img2[i];

    W.set_row(i, {
        p1.x() * p2.x(), p1.y() * p2.x(), p2.x(),
        p1.x() * p2.y(), p1.y() * p2.y(), p2.y(),
        p1.x(), p1.y(), 1
    });
}

}


Vector3D reconstruct3D(const Matrix34& M, 
                       const Matrix& K,   
                       const Vector2D& point_0, 
                       const Vector2D& point_1){

    // Create matrix Rt0 from the Identity matrix of R and t where all is 0
    Matrix34 Rt0(1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0);

    //create the projection matrix for the first camera
    Matrix34 M0 = K * Rt0;
    Matrix A(4, 4, {
        (point_0[0] * M0(2, 0)) - M0(0, 0), (point_0[0] * M0(2, 1)) - M0(0, 1), (point_0[0] * M0(2, 2)) - M0(0, 2), (point_0[0] * M0(2, 3)) - M0(0, 3),
        (point_0[1] * M0(2, 0)) - M0(1, 0), (point_0[1] * M0(2, 1)) - M0(1, 1), (point_0[1] * M0(2, 2)) - M0(1, 2), (point_0[1] * M0(2, 3)) - M0(1, 3),
        (point_1[0] * M(2, 0)) - M(0, 0), (point_1[0] * M(2, 1)) - M(0, 1), (point_1[0] * M(2, 2)) - M(0, 2), (point_1[0] * M(2, 3)) - M(0, 3),
        (point_1[1] * M(2, 0)) - M(1, 0), (point_1[1] * M(2, 1)) - M(1, 1), (point_1[1] * M(2, 2)) - M(1, 2), (point_1[1] * M(2, 3)) - M(1, 3)
    });
    // SVD
    const int m = 4, n = 4;
    Matrix U(m, m, 0.0);
    Matrix S(m, n, 0.0); 
    Matrix V(n, n, 0.0); 
    svd_decompose(A, U, S, V);
    
    //get last column of V
    Vector4D P_4D = V.get_column(V.cols() - 1);
    return P_4D.cartesian();
}

int countPointsInFront(const std::vector<Vector2D>& points_0,
                       const std::vector<Vector2D>& points_1,
                       const Matrix& M,
                       const Matrix33& R,
                       const Vector3D& t,
                       std::vector<Vector3D>& temp_points_3d,
                       const Matrix& K) {
    
    int count = 0;         
    temp_points_3d.clear();
    
    // Reconstruct 3D point from first camera
    for (size_t i = 0; i < points_0.size(); i++) {
        Vector3D P = reconstruct3D(M, K, points_0[i], points_1[i]);


        // Check if P is in front of camera 1
        if (P.z() > 0) {
            Vector3D Q = R * P + t; // Convert to second camera's coordinate system
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
    
    // Apply the projection matrix
    Vector3D projected_img = M * P_hom;
    
    // Convert to 2D Euclidean coordinates
    return Vector2D(projected_img.x() / projected_img.z(), projected_img.y() / projected_img.z());
}


std::pair<double, double> computeReprojectionError(const std::vector<Vector2D>& points_0,
                                                    const std::vector<Vector2D>& points_1,
                                                    const std::vector<Vector3D>& points_3d,
                                                    const Matrix& K,
                                                    const Matrix34& M0,
                                                    const Matrix34& M) {
    double total_squared_error = 0.0;
    size_t num_points = points_3d.size();

    for (size_t i = 0; i < num_points; ++i) {
        // Project 3D point back into image space
        Vector3D P = points_3d[i];
        Vector2D reprojected_0 = project3DPoint(P, K, M0);
        Vector2D reprojected_1 = project3DPoint(P, K, M);

        // Compute squared error for both views using direct multiplication
        Vector2D diff_0 = reprojected_0 - points_0[i];
        Vector2D diff_1 = reprojected_1 - points_1[i];

        double error_0 = diff_0[0] * diff_0[0] + diff_0[1] * diff_0[1];
        double error_1 = diff_1[0] * diff_1[0] + diff_1[1] * diff_1[1];

        total_squared_error += error_0 + error_1;
    }

    // Compute RMSE
    double rmse = std::sqrt(total_squared_error / (2 * num_points));
    return {rmse, total_squared_error};
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
    int num_points0 = points_0.size();
    int num_points1 = points_1.size();

    if (num_points1 < 8 || num_points0 < 8 ) {
        std::cerr << "Error: Not enough input points.\n";
        return false;
    }

    if (num_points1 != num_points0){
        std::cerr << "Error: Number of points of the input layers are not the same.\n";
        return false;

    // Step 1: Estimate relative pose of the two views
        
    // Normalize the input points.
    std::vector<Vector2D> points_q_img1;
    std::vector<Vector2D> points_q_img2;

    Matrix33 T_img1;
    Matrix33 T_img2;

    normalization(points_0, points_q_img1, T_img1);
    normalization(points_1, points_q_img2, T_img2);

    // estimate the fundamental matrix F;

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

    // Constraint enforcement
    const int m_2 = 3, n_2 = 3;
    Matrix U_2(m_2, m_2, 0.0);
    Matrix S_2(m_2, n_2, 0.0);
    Matrix V_2(n_2, n_2, 0.0);
    svd_decompose(F, U_2, S_2, V_2);
    S_2(2,2) = 0;

    Matrix33 Fq;
    Fq = U_2 * S_2 * transpose(V_2);
        
    // denormalization
    Matrix33 F_denormalized;
    F_denormalized = transpose(T_img2) * Fq * T_img1;

    double det_F = determinant(F_denormalized);
    if (std::abs(det_F) > 1e-6) {
        std::cerr << "Warning: Fundamental matrix determinant is not close to zero. Possible issue with computation.\n";
        std::cout<< "Determinant of F: " << det_F << std::endl;
    }

    // compute the essential matrix E;
    Matrix33 K(fx, s, cx,
             0, fy, cy,
             0, 0, 1);     
    Matrix33 E = transpose(K) * F_denormalized * K;

    // recover rotation R and t.
    Matrix33 W_2(0.0, -1.0 , 0.0,
                1.0, 0.0, 0.0,
                0.0, 0.0, 1.0);

    Matrix33 Z(0.0, 1.0 , 0.0,
            -1.0, 0.0, 0.0,
            0.0, 0.0, 0.0);

    // SVD of E
    constexpr int m_3 = 3, n_3 = 3;
    Matrix U_3(m_3, m_3, 0.0); 
    Matrix S_3(m_3, n_3, 0.0); 
    Matrix V_3(n_3, n_3, 0.0);
    svd_decompose(E, U_3, S_3, V_3);
        
    // get the 2 possible t vectors
    Vector u_3 = U_3.get_column(U_3.cols() - 1);
    Vector3D t_var1 = u_3;
    Vector3D t_var2 = -u_3;

    // Get the 2 possible R matrixes
    Matrix33 R_var1 = determinant(U_3 * W_2 * transpose(V_3)) * U_3 * W_2 * transpose(V_3);
    Matrix33 R_var2 = determinant(U_3 * transpose(W_2) * transpose(V_3)) * U_3 * transpose(W_2) * transpose(V_3);
    
    // Define the possible Rt combinations
    std::vector<Matrix34> Rt_combinations(4);
    std::vector<Matrix33> R_variants = {R_var1, R_var1, R_var2, R_var2};
    std::vector<Vector3D> t_variants = {t_var1, t_var2, t_var1, t_var2};

    for (int i = 0; i < 4; ++i) {
        Rt_combinations[i] = Matrix34(
            R_variants[i](0, 0), R_variants[i](0, 1), R_variants[i](0, 2), t_variants[i].x(),
            R_variants[i](1, 0), R_variants[i](1, 1), R_variants[i](1, 2), t_variants[i].y(),
            R_variants[i](2, 0), R_variants[i](2, 1), R_variants[i](2, 2), t_variants[i].z()
        );
    }    

    // Compute the projection matrices
    std::vector<Matrix> M_matrices(4);
    for (int i = 0; i < 4; ++i) {
        M_matrices[i] = K * Rt_combinations[i];
    }
        
    // Initialize variables to track the best results
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
    R = best_R;
    t = best_t;
    points_3d = best_points_3d;

    if (max_count < points_0.size()/2) {
        std::cerr << "Warning: Less than half of the points are in front of both cameras. Possible issue with reconstruction.\n";
        std::cout<< "Number of points in front: " << max_count << std::endl;
    }
        
    double detR = determinant(best_R);
    if (std::abs(detR - 1.0) > 1e-6) {
        std::cerr << "Warning: Rotation matrix determinant is not 1. Possible issue with decomposition.\n";
        std::cout<< "Determinant of R: " << detR << std::endl;
    }

    // Create matrix Rt0 from the Identity matrix of R and t where all is 0
    Matrix34 Rt0(1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0);
        
    //create the projection matrix for the first camera
    Matrix34 M0 = K * Rt0;

    // Compute RMSE and total squared error
    std::pair<double, double> result = computeReprojectionError(points_0, points_1, points_3d, K, M0, best_M);

    std::cout << "Reprojection RMSE: " << result.first << std::endl;
    std::cout << "Total Squared Error: " << result.second << std::endl;

    return points_3d.size() > 0;
}
