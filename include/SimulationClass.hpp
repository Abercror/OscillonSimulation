#pragma once
#include <filesystem>
#include "InflatonFieldClass.hpp"
#include "SpacetimeClass.hpp"
#include "InflationPotentials.hpp"
#include "InflatonSpacetimeVariablesStruct.hpp"
#include <string>
#include <unordered_map>
#include <utility>


template <typename VectorTraits, typename ScalarTraits>
class Simulation{
protected:
    using Scalar = typename VectorTraits::Scalar;
    using Func = typename VectorTraits::Func;
    using StateType = typename VectorTraits::StateType;

private:
    InflatonField<VectorTraits> m_inflatonField;
    SpacetimeParameters<ScalarTraits> m_spacetimeParameters;
    Func m_inflationPotential;
    Func m_inflationPotentialDerivative;
    Func m_inflationPotentialSecondDerivative;
    Scalar m_timeDelta;
    std::string m_inflationModel;
    int m_gridSize;

public:
    Simulation(
        InflatonField<VectorTraits> &inflatonField,
        SpacetimeParameters<ScalarTraits> &spacetimeParameters,
        Scalar timeDelta,
        std::string inflationModel,
        int const gridSize
        ):
    m_inflatonField(inflatonField),
    m_spacetimeParameters(spacetimeParameters),
    m_timeDelta(timeDelta),
    m_inflationModel(std::move(inflationModel)),
    m_gridSize(gridSize) {}

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Inflationary Model Selector
    ///
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void inflationaryPotentialSelector(){
        const std::unordered_map<std::string, Func> inflationaryPotentials = {
            {"T-Model", InflationPotentials<StateType, Scalar>::tModel},
            {"E-model", InflationPotentials<StateType, Scalar>::eModel},
            {"Axion Cosine", InflationPotentials<StateType, Scalar>::axionCosine},
        };
        const std::unordered_map<std::string, Func> inflationaryPotentialDerivatives = {
            {"T-Model", DifferentialInflationPotentials<StateType, Scalar>::tModel},
            {"E-model", DifferentialInflationPotentials<StateType, Scalar>::eModel},
            {"Axion Cosine", DifferentialInflationPotentials<StateType, Scalar>::axionCosine},
        };
        const std::unordered_map<std::string, Func> inflationaryPotentialSecondDerivatives = {
            {"T-Model", SecondDifferentialInflationPotentials<StateType, Scalar>::tModel},
            {"E-model", SecondDifferentialInflationPotentials<StateType, Scalar>::eModel},
            {"Axion Cosine", SecondDifferentialInflationPotentials<StateType, Scalar>::axionCosine},
        };
        this->m_inflationPotential = inflationaryPotentials.at(this->m_inflationModel);
        this->m_inflationPotentialDerivative = inflationaryPotentialDerivatives.at(this->m_inflationModel);
        this->m_inflationPotentialSecondDerivative = inflationaryPotentialSecondDerivatives.at(this->m_inflationModel);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Run Function
    ///
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    void run(int const &totalSteps, Scalar const &initialPhiValue, hsize_t const &bufferSize, hsize_t const &writingInterval){
        std::cout << "entered run function" << std::endl;

        InflatonSpacetimeVariables<Scalar> structVariables;

        auto &inflatonField = this->m_inflatonField;
        auto &spacetimeParameters = this->m_spacetimeParameters;

        this->inflationaryPotentialSelector();

        inflatonField.setLength(bufferSize);
        spacetimeParameters.setLength(bufferSize);

        inflatonField.initialConditions(initialPhiValue, this->m_inflationPotential, this->m_inflationPotentialDerivative, this->m_inflationPotentialSecondDerivative);
        Scalar averageInflatonPotential = inflatonField.getAverageInflationPotential();
        spacetimeParameters.initialConditions(averageInflatonPotential);

        hsize_t const arrayLength = m_gridSize * m_gridSize * m_gridSize;
        hsize_t const bufferWritingInterval = bufferSize * writingInterval;
        hsize_t fieldStart[2] = {0, 0};
        hsize_t fieldCount[2] = {1, arrayLength};
        hsize_t doubleStart[1] = {0};
        hsize_t doubleCount[1] = {bufferSize};
        hsize_t index = 0;

        for (int stepCount = 1; stepCount < totalSteps; ++stepCount){
            structVariables.getInflatonSpacetimeVariables(inflatonField, spacetimeParameters);
            spacetimeParameters.updateSpacetime(this->m_timeDelta, structVariables);
            inflatonField.updateInflatonField(stepCount, this->m_timeDelta, this->m_inflationPotential, this->m_inflationPotentialDerivative, structVariables);

            if (stepCount % writingInterval == 0) {
                spacetimeParameters.writeToBuffer(index);
                inflatonField.writeToBuffer(index);
                index = (index + 1) % bufferSize;
            }

            if (stepCount % bufferWritingInterval == 0) {
                inflatonField.writeToFile(fieldStart, fieldCount, doubleStart, doubleCount);
                spacetimeParameters.writeToFile(doubleStart, doubleCount);
                fieldStart[0] += bufferSize;
                doubleStart[0] += bufferSize;
            }
        }
    }
};