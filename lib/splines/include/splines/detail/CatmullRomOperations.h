//
// Created by stephencrawford on 3/3/25.
//

#ifndef CATMULLROMOPERATIONS_H
#define CATMULLROMOPERATIONS_H

#include <math/Types.h>

namespace splines {
    template <typename T>
    math::Row<T> catmullRomBasis(T t);

    template <typename T>
    math::Matrix<T> catmullRomCoefficientMatrix();

    template <typename T>
    math::Row<T> catmullRomSegment(T t, const math::Matrix<T>& controlPoints);
}

#endif //CATMULLROMOPERATIONS_H
