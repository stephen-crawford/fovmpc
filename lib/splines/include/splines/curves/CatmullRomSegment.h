//
// Created by stephencrawford on 3/3/25.
//

#ifndef CATMULLROMSEGMENT_H
#define CATMULLROMSEGMENT_H
#include "SingleParameterCurve.h"

/// Catmull-Rom splines work different from Bezier--they require exactly four control points per segment
/// You can make an arbitrary length spline by chaining segments

    template <typename T, unsigned int DIM>
    class CatmullRomSegment : public splines::SingleParameterCurve<T, DIM> {
    public:
        using Base = splines::SingleParameterCurve<T, DIM>;
        using VectorDIM = math::VectorDIM<T, DIM>;
        using MaximumDerivativeMagnitude = typename Base::MaximumDerivativeMagnitude;

        CatmullRomSegment(T max_parameter, const std::vector<VectorDIM>& control_points);

        // append a control point to CatmullRom segment
        void appendControlPoint(const VectorDIM& control_point);

        // get the control point with the given index. return status is not OK if
        // the index is out of range.
        const VectorDIM getControlPoint(std::size_t control_point_idx) const;

        // get all the control points of the CatmullRom segment as as list
        const std::vector<VectorDIM> getControlPoints() const;

        T max_parameter() const override;
        void set_max_parameter(T max_parameter) override;
        void scaleMaxParameter(T scaling_factor) override;

        // Evaluate {derivative_degree}^{th} derivative of the curve at parameter.
        // return status is not ok if parameter is out of range [0, max_parameter_]
        // or curve does not have any control points
        virtual VectorDIM eval(
                T parameter, uint64_t derivative_degree) const override;

        /**
         * @brief Returns the maximum magnitude of the derivative_degree^{th}
         * derivative of the curve as well as the parameter at which the maximum
         * occurs.
         *
         * @param derivative_degree The degree of the derivative.
         * @return absl::StatusOr<MaximumDerivativeMagnitude> The maximum magnitude
         * - parameter pair
         */
        virtual MaximumDerivativeMagnitude
        maximumDerivativeMagnitude(uint64_t derivative_degree) const override;

    private:
        // max parameter of the curve
        T max_parameter_;

        // control points of the bezier curve
        std::vector<VectorDIM> control_points_;
    };

} // splines

#endif //CATMULLROMSEGMENT_H
