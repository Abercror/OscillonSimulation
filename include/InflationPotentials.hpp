#pragma once
#include <cstddef>

template <typename StateType, typename Scalar>
struct InflationPotentials{

    static void tModel(StateType const &phiValues, StateType &potential){
        #pragma omp parallel for
        for (std::size_t i = 0; i < phiValues.size(); ++i) {
            Scalar const phi = phiValues[i];

            potential[i] = std::tanh(phi) * std::tanh(phi);
        }
    }


    static void eModel(StateType const &phiValues, StateType &potential){
        #pragma omp parallel for
        for (std::size_t i = 0; i < phiValues.size(); ++i) {
            Scalar const phi = phiValues[i];
            Scalar const expPhi = std::exp(-phi);

            potential[i] = (1 - expPhi) * (1 - expPhi);
        }
    }

    static void axionCosine(StateType const &phiValues, StateType &potential){
        #pragma omp parallel for
        for (std::size_t i = 0; i < phiValues.size(); ++i) {
            Scalar const phi = phiValues[i];

            potential[i] = 1 - std::cos(phi);
        }
    }
};


template <typename StateType, typename Scalar>
struct DifferentialInflationPotentials{

    static void tModel(StateType const &phiValues, StateType &potential){
        Scalar const c1 = Scalar(17) / Scalar(45);
        Scalar const c2 = Scalar(45) / Scalar(315);

        #pragma omp parallel for
        for (std::size_t i = 0; i < phiValues.size(); ++i) {
            Scalar const phi = phiValues[i];
            Scalar const sechPhi = 1 / std::cosh(phi);

            potential[i] = std::tanh(phi) * sechPhi * sechPhi;
        }
    }


    static void eModel(StateType const &phiValues, StateType &potential){
        #pragma omp parallel for
        for (std::size_t i = 0; i < phiValues.size(); ++i) {
            Scalar const phi = phiValues[i];
            Scalar const expPhi = std::exp(-phi);

            potential[i] = 2 * expPhi * (1 - expPhi);
        }
    }


    static void axionCosine(StateType const &phiValues, StateType &potential) {
        #pragma omp parallel for
        for (std::size_t i = 0; i < phiValues.size(); ++i) {
            Scalar const phi = phiValues[i];

            potential[i] = std::sin(phi);
        }
    }
};


template <typename StateType, typename Scalar>
struct SecondDifferentialInflationPotentials {
    static void tModel(StateType const &phiValues, StateType &potential){
        Scalar const c1 = Scalar(17) / Scalar(45);
        Scalar const c2 = Scalar(45) / Scalar(315);

        #pragma omp parallel for
        for (std::size_t i = 0; i < phiValues.size(); ++i) {
            Scalar const phi = phiValues[i];
            Scalar const sechPhi = 1 / std::cosh(phi);

            potential[i] = std::tanh(phi) * sechPhi * sechPhi;
        }
    }


    static void eModel(StateType const &phiValues, StateType &potential){
        #pragma omp parallel for
        for (std::size_t i = 0; i < phiValues.size(); ++i) {
            Scalar const phi = phiValues[i];
            Scalar const expPhi = std::exp(-phi);

            potential[i] = 2 * expPhi * (1 - expPhi);
        }
    }


    static void axionCosine(StateType const &phiValues, StateType &potential) {
        #pragma omp parallel for
        for (std::size_t i = 0; i < phiValues.size(); ++i) {
            Scalar const phi = phiValues[i];

            potential[i] = std::cos(phi);
        }
    }
};
    