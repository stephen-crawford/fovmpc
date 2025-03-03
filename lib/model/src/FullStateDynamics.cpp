//
// Created by stephencrawford on 2/28/25.
//

#include <model/FullStateDynamics.h>

namespace model {
    template <typename T, unsigned int DIM>
    StatePropagator<T> FullStateDynamics<T, DIM>::get_A0(int K) {
        StatePropagator<T> A0;
        A0.pos_ = Matrix::Zero(dim_ * K, 3 * dim_);
        A0.rpy_ = Matrix::Zero(dim_ * K, 3 * dim_);
        A0.vel_ = Matrix::Zero(dim_ * K, 3 * dim_);
        A0.w_ = Matrix::Zero(dim_ * K, 3 * dim_);

        Matrix new_row = Matrix::Zero(3 * dim_, 3 * dim_);
        Matrix prev_row = Matrix::Identity(3 * dim_, 3 * dim_);

        int idx = 0;

        for (int  k = 0; k < K; ++k) {
            new_row = A_ * prev_row;
            A0.pos_.middleRows(idx, 3) = new_row.middleRows(0, 3);
            A0.rpy_.middleRows(idx, 3) = new_row.middleRows(3, 3);
            A0.vel_.middleRows(idx, 3) = new_row.middleRows(6, 3);
            A0.wvel_.middleRows(idx, 3) = new_row.middleRows(9, 3);
            prev_row = new_row;
            idx += 3;
        }

        return A0;
    }

    template <typename T, unsigned int DIM>
    StatePropagator<T> FullStateDynamics<T, DIM>::get_lambda(int K) {
        StatePropagator<T> Lambda;
        Lambda.pos_ = Matrix::Zero(dim_ * K, dim_ * K);
        Lambda.rpy_ = Matrix::Zero(dim_ * K, dim_ * K);
        Lambda.vel_ = Matrix::Zero(dim_ * K, dim_ * K);
        Lambda.wvel_ = Matrix::Zero(dim_ * K, dim_ * K);

        Matrix prev_row = Matrix::Zero(2 * dim_, dim_ * K);
        Matrix new_row = Matrix::Zero(2 * dim_, dim_ * K);
        Matrix add_b = Matrix::Zero(2 * dim_, dim_ * K);

        int idx = 0;

        for (int k = 0; k < K; ++k) {
            add_b << Matrix::Zero(B_.rows(), B_.cols() * (k)), B_,
                    Matrix::Zero(B_.rows(), B_.cols() * (K - k - 1));
            new_row = A_ * prev_row + add_b;
            Lambda.pos_.middleRows(idx, 3) = new_row.middleRows(0, 3);
            Lambda.pos_.middleRows(idx, 3) = new_row.middleRows(3, 3);
            Lambda.vel_.middleRows(idx, 3) = new_row.middleRows(6, 3);
            Lambda.pos_.middleRows(idx, 3) = new_row.middleRows(12, 3);
            prev_row = new_row;
            idx += 3;
        }

        return Lambda;
    }

    template <typename T, unsigned int DIM>
    State<T, DIM> FullStateDynamics<T, DIM>::applyInput(const State &state,
                                                       const Vector &u) {
        Vector x_t0 = Vector::Zero(2*DIM);
        x_t0 << state.pos_, state.vel_;

        assert(A_.cols() == x_t0.rows());
        assert(B_.cols() == u.rows());
        Vector x_t1 = A_ * x_t0 + B_ * u;
        State new_state = {x_t1.segment(0, DIM), x_t1.segment(DIM, DIM)};
        return new_state;
    }

    template class FullStateDynamics<double, 3U>;
    template class FullStateDynamics<float, 3U>;
    template class FullStateDynamics<double, 2U>;
    template class FullStateDynamics<float, 2U>;

} // model