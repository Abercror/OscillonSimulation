#pragma once
#include <functional>
#include <omp.h>

template <typename Traits>
class NumericalMethods {
protected:
    using Scalar    = typename Traits::Scalar;
    using StateType = typename Traits::StateType;
    using Func      = typename Traits::Func;
    using Func2  = typename Traits::Func2;

public:
    void leapfrog2ndOrder(StateType &yVariable, StateType &dyVariable, StateType &d2yVariable, Scalar const &delta, auto secondDifferentialFunc){
        auto &y = yVariable;
        auto &dy = dyVariable;
        auto &ddy = d2yVariable;
        auto halfDelta = delta / Scalar(2);

        secondDifferentialFunc();

        if constexpr (!std::is_same_v<std::decay_t<StateType>, Scalar>) {
            #pragma omp parallel for
            for (int i = 0; i < y.field.size(); ++i){
                dy[i] += ddy[i] * halfDelta;
                y[i] += dy[i] * delta;
            }
        }
        else {
            dy += ddy * halfDelta;
            y += dy * delta;
        }

        secondDifferentialFunc();
        if constexpr (!std::is_same_v<std::decay_t<StateType>, Scalar>) {
            #pragma omp parallel for
            for (int i = 0; i < y.field.size(); ++i){
                dy[i] += ddy[i] * halfDelta;
            }
        }
        else {
            dy += ddy * halfDelta;
        }

        // if constexpr (!std::is_same_v<std::decay_t<StateType>, Scalar>) {
        //     for (int i = 0; i < y.field.size(); ++i){
        //         secondDifferentialFunc();
        //         dy[i] += ddy[i] * halfDelta;
        //         y[i] += dy[i] * delta;
        //         secondDifferentialFunc();
        //         dy[i] += ddy[i] * halfDelta;
        //     }
        // }
        // else {
        //     secondDifferentialFunc();
        //     dy += ddy * halfDelta;
        //     y += dy * delta;
        //     secondDifferentialFunc();
        //     dy += ddy * halfDelta;
        // }
    }
};