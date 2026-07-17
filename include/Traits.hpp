#pragma once
#include "LatticeType.hpp"
#include "InflatonSpacetimeVariablesStruct.hpp"
#include <functional>



template <typename T>
struct VectorTraits {
    using Scalar = T;
    using varStruct = InflatonSpacetimeVariables<Scalar>;
    using StateType = lattice<T>;
    using Func = std::function<void(const StateType&, StateType&)>;
    using Func2 = std::function<void(Func, varStruct)>;
};


template <typename T> 
struct ScalarTraits {
    using Scalar = T;
    using varStruct = InflatonSpacetimeVariables<Scalar>;
    using StateType = T;
    using Func = std::function<void(const T&, T&)>;
    using Func2 = std::function<void(Func, varStruct)>;
};