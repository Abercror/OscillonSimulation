#pragma once
#include "NumericalMethodsClass.hpp"
#include <vector>
#include <cmath>


template <typename Scalar>
struct SpacetimeParametersData{
    std::vector<Scalar> scaleFactorValues;
    std::vector<Scalar> dScaleFactorValues;
    std::vector<Scalar> d2ScaleFactorValues;
    std::vector<Scalar> hubbleParameterValues;

    void reserveMemory(Scalar const &totalSteps){
        scaleFactorValues.reserve(totalSteps);
        dScaleFactorValues.reserve(totalSteps);
        d2ScaleFactorValues.reserve(totalSteps);
        hubbleParameterValues.reserve(totalSteps);
    }
};

template <typename Traits>
class SpacetimeParameters: public NumericalMethods<Traits>{
public:
    using Scalar = typename Traits::Scalar;
    using StateType = typename Traits::StateType;
    using Func = typename Traits::Func;

private:
    Scalar m_scaleFactor;
    Scalar m_dScaleFactor;
    Scalar m_d2ScaleFactor;
    Scalar m_hubbleParameter;
    SpacetimeParametersData<Scalar> m_history;

public:

    SpacetimeParameters(
        Scalar scaleFactor,
        Scalar dScaleFactor,
        Scalar d2ScaleFactor,
        Scalar hubbleParameter,
        SpacetimeParametersData<Scalar> history): m_scaleFactor(scaleFactor), m_dScaleFactor(dScaleFactor), m_d2ScaleFactor(d2ScaleFactor), m_hubbleParameter(hubbleParameter), m_history(history) {}


    Scalar getScaleFactor() const { return this->m_scaleFactor; }
    Scalar getHubbleParameter() const { return this->m_hubbleParameter; }
    SpacetimeParameters<Scalar> getHistory() const { return this->m_history; }

    void setLength(Scalar const &totalSteps){
        this->m_history.reserveMemory(totalSteps);
    }

    auto accelerationEquationSecondDerivative(auto const &structVariables){
        Scalar const &scaleFactor = this->m_scaleFactor;
        auto &averagePhiDerivative = structVariables.m_averagePhiDerivative;
        auto &averagePressure = structVariables.m_averagePressure;

        this->m_d2ScaleFactor = - (scaleFactor / Scalar(3)) * (averagePhiDerivative + Scalar(3) * averagePressure);
    }


    void hubbleParameter(Scalar const &energyDensity){
        this->m_hubbleParameter = std::sqrt(energyDensity/Scalar(3));
        // this->m_hubbleParameter = this->m_dScaleFactor / this->m_scaleFactor;
    }

    void writeToHistory(){
        this->m_history.scaleFactorValues.push_back(this->m_scaleFactor);
        this->m_history.dScaleFactorValues.push_back(this->m_dScaleFactor);
        this->m_history.d2ScaleFactorValues.push_back(this->m_d2ScaleFactor);
        this->m_history.hubbleParameterValues.push_back(this->m_hubbleParameter);
    }

    void updateSpacetime(Scalar const &timeDelta, auto const &structVariables){
        this->hubbleParameter(structVariables.m_energyDensity);
        this->leapfrog2ndOrder(this->m_scaleFactor, this->m_dScaleFactor, this->m_d2ScaleFactor, timeDelta, [this, &structVariables](){ return this->accelerationEquationSecondDerivative(structVariables);});
        this->writeToHistory();
    }
};
