//
// Created by stephencrawford on 3/3/25.
//
#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "splines/detail/CatmullRomOperations.h"
#include "splines/curves/CatmullRomSpline.h"

int main() {
    using Vector2d = Eigen::Vector2d;

    // Define control points as Vector2d
    std::vector<Vector2d> controlPoints = {
        {0.0, 0.0},
        {1.0, 2.0},
        {2.0, 3.0},
        {4.0, 1.0}
    };

    // Create a Catmull-Rom spline
    splines::CatmullRomSpline<double, 2> CatmullRomSegment(controlPoints);

    // Evaluate the spline at t = 0.5
    double t = 0.5;
    Vector2d result = spline.eval(t);

    // Print result
    std::cout << "Interpolated Point: (" << result.x() << ", " << result.y() << ")" << std::endl;

    return 0;
}
