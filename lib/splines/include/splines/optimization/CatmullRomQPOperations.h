#ifndef CATMULL_ROM_QP_OPERATIONS_H
#define CATMULL_ROM_QP_OPERATIONS_H

#include <vector>
#include <memory>
#include "CatmullRomSpline.h"

namespace splines {

    template <typename T, unsigned int DIM>
    class CatmullRomQPOperations {
    public:
        struct Params {
            std::size_t num_control_points_;
            T max_parameter_;
        };

        explicit CatmullRomQPOperations(Params &p);

        std::size_t numDecisionVariables() const;
        std::unique_ptr<CatmullRomSpline<T, DIM>> generateCurveFromSolution(const std::vector<T> &decision_variables) const;
        
        using Row = Eigen::Matrix<T, 1, Eigen::Dynamic>;
        Row evalBasisRow(unsigned int dimension, T parameter, uint64_t derivative_degree) const;

        struct CostAddition {
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> quadratic_term;
            Eigen::Matrix<T, Eigen::Dynamic, 1> linear_term;
            T constant;
        };
        CostAddition integratedSquaredDerivativeCost(uint64_t derivative_degree, T lambda) const;
        CostAddition evalCost(T parameter, uint64_t derivative_degree, const Eigen::Matrix<T, DIM, 1> &target, T lambda) const;
        
        struct LinearConstraint {
            Row coefficients;
            T lower_bound;
            T upper_bound;
        };
        std::vector<LinearConstraint> evalConstraint(T parameter, uint64_t derivative_degree, const Eigen::Matrix<T, DIM, 1>& target) const;
        std::vector<LinearConstraint> evalBound(T parameter, uint64_t derivative_degree, const Eigen::Matrix<T, DIM, 1>& LB, const Eigen::Matrix<T, DIM, 1>& UB) const;
        
        void set_max_parameter(T new_max_parameter);
        T max_parameter() const;

    private:
        std::size_t num_control_points_;
        T max_parameter_;
    };

    // Explicit template instantiations
    extern template class CatmullRomQPOperations<float, 2U>;
    extern template class CatmullRomQPOperations<double, 2U>;
    extern template class CatmullRomQPOperations<float, 3U>;
    extern template class CatmullRomQPOperations<double, 3U>;

} // namespace splines

#endif // CATMULL_ROM_QP_OPERATIONS_H