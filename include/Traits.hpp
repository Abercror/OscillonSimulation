#pragma once
#include "LatticeType.hpp"
#include "OscillonSpacetimeVariablesStruct.hpp"
#include <functional>



template <typename T>
struct VectorTraits {
    using Scalar = T;
    using varStruct = OscillonSpacetimeVariables<Scalar>;
    using StateType = lattice<T>;
    using Func = std::function<void(const Scalar&, Scalar&)>;
    using Func2 = std::function<void(Func, varStruct)>;
};


template <typename T> 
struct ScalarTraits {
    using Scalar = T;
    using varStruct = OscillonSpacetimeVariables<Scalar>;
    using StateType = T;
    using Func = std::function<void(const T&, T&)>;
    using Func2 = std::function<void(Func, varStruct)>;
};