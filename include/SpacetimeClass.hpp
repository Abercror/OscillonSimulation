#pragma once
#include "NumericalMethodsClass.hpp"
#include "SimulationDataSets.hpp"
#include <vector>
#include <cmath>
#include <H5Cpp.h>
#include <iostream>


template <typename Scalar>
struct SpacetimeParametersData{
    std::vector<Scalar> scaleFactorValues;
    std::vector<Scalar> dScaleFactorValues;
    std::vector<Scalar> d2ScaleFactorValues;
    std::vector<Scalar> hubbleParameterValues;

    void resizeBuffer(hsize_t const &bufferSize){
        scaleFactorValues.resize(bufferSize);
        dScaleFactorValues.resize(bufferSize);
        d2ScaleFactorValues.resize(bufferSize);
        hubbleParameterValues.resize(bufferSize);
    }

    void wipeBuffer() {
        scaleFactorValues.clear();
        dScaleFactorValues.clear();
        d2ScaleFactorValues.clear();
        hubbleParameterValues.clear();
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
    SpacetimeParametersData<Scalar> m_buffer;
    SimulationDataSets m_dataSets;

public:
    SpacetimeParameters(SimulationDataSets &dataSets) {
        StateType scaleFactor = 0.0;
        StateType dScaleFactor = 0.0;
        StateType d2ScaleFactor = 0.0;
        StateType hubbleParameter = 0.0;
        SpacetimeParametersData<Scalar> buffer;

        this->m_scaleFactor = scaleFactor;
        this->m_dScaleFactor = dScaleFactor;
        this->m_d2ScaleFactor = d2ScaleFactor;
        this->m_hubbleParameter = hubbleParameter;
        this->m_buffer = buffer;

        this->m_dataSets = dataSets;

    }


    Scalar getScaleFactor() const { return this->m_scaleFactor; }
    Scalar getHubbleParameter() const { return this->m_hubbleParameter; }
    SpacetimeParametersData<Scalar> getBuffer() const { return this->m_buffer; }

    void setLength(hsize_t const &bufferSize){
        this->m_buffer.resizeBuffer(bufferSize);
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

    void writeToBuffer(auto const &index){
        this->m_buffer.scaleFactorValues[index] = this->m_scaleFactor;
        this->m_buffer.dScaleFactorValues[index] = this->m_dScaleFactor;
        this->m_buffer.d2ScaleFactorValues[index] = this->m_d2ScaleFactor;
        this->m_buffer.hubbleParameterValues[index] = this->m_hubbleParameter;
    }
    
    void writeDoubleToFile(DataSet &dataSet, auto const &buffer, auto const &start, auto const &count) {
        DataSpace fileSpace = dataSet.getSpace();
        fileSpace.selectHyperslab(H5S_SELECT_SET, count, start);
        DataSpace memorySpace(1, count);
        dataSet.write(buffer.data(), PredType::NATIVE_DOUBLE, memorySpace, fileSpace);
    }

    void writeToFile(auto const &doubleStart, auto const &doubleCount) {
        auto data = this->m_dataSets;
        auto buffer = this->m_buffer;
        writeDoubleToFile(data.scaleFactor, buffer.scaleFactorValues, doubleStart, doubleCount);
        writeDoubleToFile(data.scaleFactorDerivative, buffer.dScaleFactorValues, doubleStart, doubleCount);
        writeDoubleToFile(data.scaleFactorSecondDerivative, buffer.d2ScaleFactorValues, doubleStart, doubleCount);
        writeDoubleToFile(data.hubbleParameter, buffer.hubbleParameterValues, doubleStart, doubleCount);
    }

    void updateSpacetime(Scalar const &timeDelta, auto const &structVariables){
        this->hubbleParameter(structVariables.m_energyDensity);
        this->leapfrog2ndOrder(this->m_scaleFactor, this->m_dScaleFactor, this->m_d2ScaleFactor, timeDelta, [this, &structVariables](){ return this->accelerationEquationSecondDerivative(structVariables);});
    }

};
