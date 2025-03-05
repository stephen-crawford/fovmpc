//
// Created by stephencrawford on 3/3/25.
//

#include <splines/detail/CatmullRomOperations.h>
namespace splines {

    template <typename T>
    math::Row<T> catmullRomBasis(T t) {
        math::Row<T> basis(4);
        T t2 = t * t;
        T t3 = t2 * t;

        basis(0) = -0.5 * t3 + t2 - 0.5 * t;
        basis(1) =  1.5 * t3 - 2.5 * t2 + 1.0;
        basis(2) = -1.5 * t3 + 2.0 * t2 + 0.5 * t;
        basis(3) =  0.5 * t3 - 0.5 * t2;

        return basis;
    }


    template <typename T>
    math::Matrix<T> catmullRomCoefficientMatrix() {
        math::Matrix<T> M(4, 4);

        M << -0.5,  1.5, -1.5,  0.5,
              1.0, -2.5,  2.0, -0.5,
             -0.5,  0.0,  0.5,  0.0,
              0.0,  1.0,  0.0,  0.0;

        return M;
    }

    template <typename T>
    math::Row<T> catmullRomSegment(T t, const std::vector<Eigen::Vector2<T>>& controlPoints) {
        if (controlPoints.size() != 4) {
            throw std::runtime_error("catmullRomSegment: Requires exactly 4 control points.");
        }

        // Compute centripetal parameterization
        std::vector<T> t_values = centripetalParameterization(controlPoints);

        // Ensure t is within the valid range
        if (t < t_values[1] || t > t_values[2]) {
            throw std::runtime_error("catmullRomSegment: t is out of valid range.");
        }

        // Compute local normalized τ
        T tau = (t - t_values[1]) / (t_values[2] - t_values[1]);

        // Compute basis with the normalized τ
        math::Row<T> basis = catmullRomBasis(tau);

        // Construct the 4x2 matrix of control points
        math::Matrix<T> P(4, 2);
        for (size_t i = 0; i < 4; ++i) {
            P.row(i) = controlPoints[i];
        }

        // Compute interpolated position
        return basis * P;
    }

    template <typename T>
    std::vector<T> centripetalParameterization(const std::vector<Eigen::Vector2<T>>& points) {
    std::vector<T> t(points.size(), 0);
    for (size_t i = 1; i < points.size(); ++i) {
        t[i] = t[i - 1] + std::sqrt((points[i] - points[i - 1]).norm());
    }

    // Normalize t to be in [0, 1]
    for (size_t i = 1; i < points.size(); ++i) {
        t[i] /= t.back();
    }

    return t;
}

} // splines
