#pragma once
#include <filesystem>
#include "OscillonClass.hpp"
#include "SpacetimeClass.hpp"
#include "InflationPotentials.hpp"
#include "OscillonSpacetimeVariablesStruct.hpp"
#include <string>
#include <unordered_map>


template <typename VectorTraits, typename ScalarTraits>
class Simulation{
public:
    using VectorT = VectorTraits;
    using ScalarT = ScalarTraits;
    using Scalar = typename VectorT::Scalar;
    using Func = typename VectorT::Func;
    using StateType = typename VectorT::StateType;

    Oscillon<VectorT> m_oscillonField;
    SpacetimeParameters<ScalarT> m_spacetimeParameters;
    Func m_inflationPotential;
    Func m_inflationPotentialDerivative;
    Scalar m_timeDelta;
    std::string m_inflationModel;


    Simulation(Oscillon<VectorT> &oscillonField, SpacetimeParameters<ScalarT> &spacetimeParameters, Scalar timeDelta, std::string const &inflationModel): m_oscillonField(oscillonField), m_spacetimeParameters(spacetimeParameters), m_timeDelta(timeDelta), m_inflationModel(inflationModel) {}

    void inflationaryPotentialSelector(){
        const std::unordered_map<std::string, Func> inflationaryPotentials = {
            {"T-Model", InflationPotentials<Scalar>::tModel},
            {"E-model", InflationPotentials<Scalar>::eModel},
            {"Axion Cosine", InflationPotentials<Scalar>::axionCosine},
        };
        const std::unordered_map<std::string, Func> inflationaryPotentialDerivatives = {
            {"T-Model", DifferentiatedInflationPotentials<Scalar>::tModel},
            {"E-model", DifferentiatedInflationPotentials<Scalar>::eModel},
            {"Axion Cosine", DifferentiatedInflationPotentials<Scalar>::axionCosine},
        };
        this->m_inflationPotential = inflationaryPotentials.at(this->m_inflationModel);
        this->m_inflationPotentialDerivative = inflationaryPotentialDerivatives.at(this->m_inflationModel);
    }


    void run(Scalar const &totalSteps){
        OscillonSpacetimeVariables<Scalar> structVariables;
        auto &oscillon = this->m_oscillonField;
        auto &spacetimeParameters = this->m_spacetimeParameters;
        this->inflationaryPotentialSelector();
        oscillon.setLength(totalSteps);
        spacetimeParameters.setLength(totalSteps);

        for(Scalar stepCount = 0; stepCount < totalSteps; ++stepCount){
            structVariables.getOscillonSpacetimeVariables(oscillon, spacetimeParameters);

            spacetimeParameters.updateSpacetime(this->m_timeDelta, structVariables);

            oscillon.updateOscillon(stepCount, this->m_timeDelta, this->m_inflationPotential, this->m_inflationPotentialDerivative, structVariables);
        }
    }
};
