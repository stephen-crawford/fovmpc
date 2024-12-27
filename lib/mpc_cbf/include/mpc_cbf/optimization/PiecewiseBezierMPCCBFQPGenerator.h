//
// Created by lishuo on 9/21/24.
//

#ifndef MPC_CBF_PIECEWISEBEZIERMPCCBFQPGENERATOR_H
#define MPC_CBF_PIECEWISEBEZIERMPCCBFQPGENERATOR_H

#include <mpc/optimization/PiecewiseBezierMPCQPGenerator.h>
#include <mpc_cbf/optimization/PiecewiseBezierMPCCBFQPOperations.h>

namespace mpc_cbf {
    template <typename T, unsigned int DIM>
    class PiecewiseBezierMPCCBFQPGenerator {
    public:
        using PiecewiseBezierMPCQPGenerator = mpc::PiecewiseBezierMPCQPGenerator<T, DIM>;
        using PiecewiseBezierMPCCBFQPOperations = mpc_cbf::PiecewiseBezierMPCCBFQPOperations<T, DIM>;
        using PiecewiseBezierMPCQPOperations = mpc::PiecewiseBezierMPCQPOperations<T, DIM>;
        using CostAddition = typename PiecewiseBezierMPCCBFQPOperations::CostAddition;
        using LinearConstraint = typename PiecewiseBezierMPCCBFQPOperations::LinearConstraint;
        using State = typename PiecewiseBezierMPCCBFQPOperations::State;
        using Vector = typename PiecewiseBezierMPCCBFQPOperations::Vector;
        using Row = typename PiecewiseBezierMPCCBFQPOperations::Row;
        void addPiecewise(std::unique_ptr<PiecewiseBezierMPCCBFQPOperations> &&piecewise_mpc_cbf_operations_ptr, int num_neighbors=0, bool slack_mode=false);
        // return a reference to the problem instance generated by this generator
        qpcpp::Problem<T>& problem();
        // return a pointer to the PiecewiseBezierMPCQPGenerator API
        std::shared_ptr<PiecewiseBezierMPCQPGenerator> piecewise_mpc_qp_generator_ptr();

        void addSlackCost(const std::vector<double>& slack_weights);

        void addSafetyCBFConstraint(const State& current_state, const Vector& other_pos, T slack_value=0);
        void addFovLBConstraint(const State& current_state, const Vector& other_pos, T slack_value=0);
        void addFovRBConstraint(const State& current_state, const Vector& other_pos, T slack_value=0);

        void addPredSafetyCBFConstraints(const std::vector<State>& pred_states, const Vector& other_pos, const std::vector<T>& slack_values);
        void addPredFovLBConstraints(const std::vector<State>& pred_states, const Vector& other_pos, const std::vector<T>& slack_values);
        void addPredFovRBConstraints(const std::vector<State>& pred_states, const Vector& other_pos, const std::vector<T>& slack_values);

        void addSafetyCBFConstraintWithSlackVariables(const State& current_state, const Vector& other_pos);
        void addFovLBConstraintWithSlackVariables(const State& current_state, const Vector& other_pos);
        void addFovRBConstraintWithSlackVariables(const State& current_state, const Vector& other_pos);

//        void addPredSafetyCBFConstraintsWithSlackVariables(const std::vector<State>& pred_states, const Vector& other_pos);
        void addPredFovLBConstraintsWithSlackVariables(const std::vector<State>& pred_states, const Vector& other_pos);
        void addPredFovRBConstraintsWithSlackVariables(const std::vector<State>& pred_states, const Vector& other_pos);

        // adds given cost_addition for piecewise to the generated qp.
        void addCostAdditionForSlackVariables(const CostAddition& cost_addition);

        // adds given linear_constraint for piecewise to the generated qp.
        void addLinearConstraintForPiecewiseWithSlackVariables(const LinearConstraint& linear_constraint,
                                                               const Row &slack_coefficients);

    private:
        std::shared_ptr<PiecewiseBezierMPCQPGenerator> piecewise_mpc_qp_generator_ptr_ = std::make_shared<PiecewiseBezierMPCQPGenerator>();
        std::unique_ptr<PiecewiseBezierMPCCBFQPOperations> piecewise_mpc_cbf_operations_ptr_;
        // slack variables, extra to the curve variables
        std::vector<qpcpp::Variable<T>*> slack_variables_;

    };

} // mpc_cbf

#endif //MPC_CBF_PIECEWISEBEZIERMPCCBFQPGENERATOR_H
