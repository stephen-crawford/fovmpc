//
// Created by stephencrawford on 3/3/25.
//

#include <stdexcept>
#include <splines/curves/CatmullRomSegment.h>

namespace splines {
    template<typename T, unsigned int DIM>
    class CatmullRomSegment {
    public:
        using VectorDIM = Eigen::Matrix<T, DIM, 1>;

        CatmullRomSegment(const std::vector<VectorDIM>& control_points)
            : control_points_(control_points) {
            if (control_points_.size() != 4) {
                throw std::runtime_error("CatmullRomSegment requires exactly 4 control points.");
            }
        }

        VectorDIM eval(T t) const {
            math::Row<T> basis = catmullRomBasis(t);
            math::Matrix<T> P(4, DIM);

            for (size_t i = 0; i < 4; ++i) {
                P.row(i) = control_points_[i];
            }

            return basis * P;
        }

    private:
        std::vector<VectorDIM> control_points_;
    };

} // namespace splines
