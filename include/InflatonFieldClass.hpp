#pragma once
#include "NumericalMethodsClass.hpp"
#include "SimulationDataSets.hpp"
#include <vector>
#include <array>
#include <functional>
#include <cstddef>
#include <H5Cpp.h>


template <typename StateType, typename Scalar>
struct inflatonFieldData{
    std::vector<StateType> phiValues;
    std::vector<StateType> inflationPotentialValues;
    std::vector<StateType> energyDensityValues;
    std::vector<StateType> energyOverdensityValues;
    std::vector<Scalar> averageEnergyDensity;
    
    void resizeBuffer(hsize_t const &bufferSize){
        phiValues.resize(bufferSize);
        inflationPotentialValues.resize(bufferSize);
        energyDensityValues.resize(bufferSize);
        energyOverdensityValues.resize(bufferSize);
        averageEnergyDensity.resize(bufferSize);
    }

    void wipeBuffer() {
        phiValues.clear();
        inflationPotentialValues.clear();
        energyDensityValues.clear();
        energyOverdensityValues.clear();
        averageEnergyDensity.clear();
    }
};


template <typename Traits>
class InflatonField: public NumericalMethods<Traits> {
protected:
    using Scalar = typename Traits::Scalar;
    using StateType = typename Traits::StateType;
    using Func = typename Traits::Func;


private:
    StateType m_phi;
    StateType m_dPhi;
    StateType m_d2Phi;
    StateType m_phiLaplacian;
    StateType m_inflationPotential;
    StateType m_inflationPotentialDerivative;
    StateType m_inflationPotentialSecondDerivative;
    StateType m_energyDensity;
    StateType m_energyOverDensity;
    Scalar m_averageEnergyDensity;
    StateType m_kineticEnergyDensity;
    StateType m_potentialEnergyDensity;
    Scalar m_pressure;
    Scalar m_gridSize;
    inflatonFieldData<StateType, Scalar> m_buffer;
    SimulationDataSets m_dataSets;

public:

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Constructor
    ///
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    InflatonField(Scalar gridSize, SimulationDataSets &dataSets) {
        StateType phi(gridSize, gridSize, gridSize);
        StateType dPhi(gridSize, gridSize, gridSize);
        StateType d2Phi(gridSize, gridSize, gridSize);
        StateType phiLaplacian(gridSize, gridSize, gridSize);
        StateType inflationPotential(gridSize, gridSize, gridSize);
        StateType inflationPotentialDerivative(gridSize, gridSize, gridSize);
        StateType inflationPotentialSecondDerivative(gridSize, gridSize, gridSize);
        StateType energyDensity(gridSize, gridSize, gridSize);
        Scalar averageEnergyDensity = 0.0;
        StateType energyOverDensity(gridSize, gridSize, gridSize);
        StateType kineticEnergyDensity(gridSize, gridSize, gridSize);
        StateType potentialEnergyDensity(gridSize, gridSize, gridSize);
        Scalar pressure = 0.0;
        inflatonFieldData<StateType, Scalar> buffer;

        this->m_phi = phi;
        this->m_dPhi = dPhi;
        this->m_d2Phi = d2Phi;
        this->m_phiLaplacian = phiLaplacian;
        this->m_inflationPotential = inflationPotential;
        this->m_inflationPotentialDerivative = inflationPotentialDerivative;
        this->m_inflationPotentialSecondDerivative = inflationPotentialSecondDerivative;
        this->m_energyDensity = energyDensity;
        this->m_averageEnergyDensity = averageEnergyDensity;
        this->m_energyOverDensity = energyOverDensity;
        this->m_kineticEnergyDensity = kineticEnergyDensity;
        this->m_potentialEnergyDensity = potentialEnergyDensity;
        this->m_pressure = pressure;
        this->m_buffer = buffer;

        this->m_gridSize = gridSize;
        this->m_dataSets = dataSets;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Getter and Setter Functions
    ///
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Scalar getPhiDerivativeMean() const {
        auto averagePhiDerivative = this->m_dPhi.mean();
        return averagePhiDerivative;
    }
    Scalar getPressure() const { return this->m_pressure; }
    Scalar getAverageEnergyDensity() const { return this->m_averageEnergyDensity; }
    Scalar getAverageInflationPotential() const {
        auto averagePotential = this->m_inflationPotential.mean();
        return averagePotential;
    }
    inflatonFieldData<StateType, Scalar> getBuffer() const { return this->m_buffer; }

    void setLength(hsize_t const &bufferSize){
        this->m_buffer.resizeBuffer(bufferSize);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Initial Conditions
    ///
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void initialConditions(Scalar const &initialPhiValue, Func const &inflationPotentialFunc, Func const &inflationPotentialDerivativeFunc, Func const &inflationPotentialSecondDerivativeFunc) {
        for (int i = 0; i < this->m_phi.size(); ++i) {
            this->m_phi[i] = initialPhiValue;

        }

        inflationPotentialFunc(this->m_phi, this->m_inflationPotential);
        inflationPotentialDerivativeFunc(this->m_phi, this->m_inflationPotentialDerivative);
        inflationPotentialSecondDerivativeFunc(this->m_phi, this->m_inflationPotentialSecondDerivative);
        this->m_energyDensity = this->m_inflationPotential;
        this->m_averageEnergyDensity = this->m_energyDensity.mean();
        this->m_pressure = -this->m_averageEnergyDensity;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Equations of Motion Functions
    ///
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void phiSecondDifferential(Func const &potentialDerivative, auto const &structVariables){
        StateType const &phi = this->m_phi;
        StateType const &dPhi = this->m_dPhi;
        StateType const &laplacianPhi = this->m_phiLaplacian;

        potentialDerivative(phi, this->m_inflationPotentialDerivative);

        StateType &inflationPotentialDerivative = this->m_inflationPotentialDerivative;

        Scalar const &hubbleParameter = structVariables.m_hubbleParameter;
        Scalar const &scaleFactor = structVariables.m_scaleFactor;

        Scalar const scaleFactor2 = scaleFactor * scaleFactor;
        Scalar const hubbleTerm = Scalar(3) * hubbleParameter;

        StateType &secondDerivative = this->m_d2Phi;
        #pragma omp parallel for
        for (std::size_t i = 0; i < phi.size(); ++i) {
            secondDerivative[i] = laplacianPhi[i] / scaleFactor2 - hubbleTerm * dPhi[i] - inflationPotentialDerivative[i];
        }
    }

    void determineLaplacian(Scalar const &deltaPosition){
        StateType const &phi = this->m_phi;
        StateType &laplacian = this->m_phiLaplacian;

        Scalar deltaPosition2 = deltaPosition * deltaPosition;

        #pragma omp parallel for collapse(3) schedule(static)
        for(std::size_t z = 1; z < this->m_gridSize-1; ++z) {
            for(std::size_t y = 1; y < this->m_gridSize-1; ++y){
                for (std::size_t x = 1; x < this->m_gridSize-1; ++x){
                    laplacian(x, y, z) = ((phi(x+1, y, z) - Scalar(2) * phi(x, y , z) + phi(x-1, y, z)) / deltaPosition2) + ((phi(x, y+1, z) - Scalar(2) * phi(x, y , z) + phi(x, y-1, z)) / deltaPosition2) +((phi(x, y, z+1) - Scalar(2) * phi(x, y , z) + phi(x, y, z-1)) / deltaPosition2);
                    }
                }
            }
        }

    void potentialKineticEnergyDensity(Scalar const &stepCount, Func const &inflationPotentialFunc, auto const &structVariables) {
        StateType const &phi = this->m_phi;
        StateType &inflationPotential = this->m_inflationPotential;
        Scalar const &scaleFactor = structVariables.m_scaleFactor;
        Scalar const scaleFactorTwo = Scalar(2) * scaleFactor;

        inflationPotentialFunc(phi, inflationPotential);

        #pragma omp parallel for
        for (std::size_t i = 0; i < phi.size(); ++i) {
            this->m_energyDensity[i] = this->m_dPhi[i] * this->m_dPhi[i] / Scalar(2) + this->m_phiLaplacian[i] * this->m_phiLaplacian[i] / scaleFactorTwo + this->m_inflationPotential[i];
        }

        this->m_averageEnergyDensity = this->m_energyDensity.mean();

        #pragma omp parallel for
        for (std::size_t i = 0; i < phi.size(); ++i) {
            this->m_energyOverDensity[i] = this->m_energyDensity[i] / this->m_averageEnergyDensity;
        }

        #pragma omp parallel for
        for (std::size_t i = 0; i < phi.size(); ++i) {
            this->m_kineticEnergyDensity[i] = m_inflationPotential[i] * std::cos(stepCount) * std::cos(stepCount);
            this->m_potentialEnergyDensity[i] = m_energyDensity[i] - m_kineticEnergyDensity[i];
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Writing Data
    ///
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void writeToBuffer(auto const &index) {
        this->m_buffer.phiValues[index] = this->m_phi;
        this->m_buffer.energyDensityValues[index] = this->m_energyDensity;
        this->m_buffer.energyOverdensityValues[index] = this->m_energyOverDensity;
        this->m_buffer.inflationPotentialValues[index] = this->m_inflationPotential;
        this->m_buffer.averageEnergyDensity[index] = this->m_averageEnergyDensity;
    }


    void writeFieldsToFile(H5::DataSet &dataSet, auto &buffer, auto const &fieldStart, auto const &fieldCount) {
        hsize_t start[2] = {fieldStart[0], fieldStart[1]};
        auto size = buffer.size();

        for (std::size_t i = 0; i < size; ++i) {
            H5::DataSpace fileSpace = dataSet.getSpace();
            fileSpace.selectHyperslab(H5S_SELECT_SET, fieldCount, start);
            H5::DataSpace memorySpace(2, fieldCount);
            dataSet.write(buffer[i].data(), H5::PredType::NATIVE_DOUBLE, memorySpace, fileSpace);

            start[0] +=1;
        }
    }

    void testing(H5::DataSet &dataSet, auto &buffer, auto const &fieldStart, auto const &fieldCount, auto const &name) {
        H5::DataSpace fileSpace = dataSet.getSpace();

        hsize_t dims[2];
        fileSpace.getSimpleExtentDims(dims);

        std::cout
            << " name = " << name
            << " dims = [" << dims[0] << ", " << dims[1] << "]"
            << " start = [" << fieldStart[0] << ", " << fieldStart[1] << "]"
            << " count = [" << fieldCount[0] << ", " << fieldCount[1] << "]"
            << std::endl;
    }

    void writeDoubleToFile(H5::DataSet &dataSet, auto &buffer, auto const &start, auto const &count) {
        H5::DataSpace fileSpace = dataSet.getSpace();
        fileSpace.selectHyperslab(H5S_SELECT_SET, count,  start);
        H5::DataSpace memorySpace(1, count);
        dataSet.write(buffer.data(), H5::PredType::NATIVE_DOUBLE, memorySpace, fileSpace);
    }

    void writeToFile(auto const &fieldStart, auto const &fieldCount, auto const &doubleStart, auto const &doubleCount) {
        auto data = this->m_dataSets;
        auto buffer = this->m_buffer;
        writeFieldsToFile(data.inflatonField, buffer.phiValues, fieldStart, fieldCount);
        writeFieldsToFile(data.energyDensity, buffer.energyDensityValues, fieldStart, fieldCount);
        writeFieldsToFile(data.inflatonPotential, buffer.inflationPotentialValues, fieldStart, fieldCount);
        writeDoubleToFile(data.averageEnergyDensity, buffer.averageEnergyDensity, doubleStart, doubleCount);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Update Function
    ///
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void updateInflatonField(Scalar const &stepCount, Scalar const &timeDelta, Func &inflationPotential, Func &inflationPotentialDerivative, auto const &structVariables){
        auto deltaPosition = Scalar(1)/(this->m_gridSize) * structVariables.m_scaleFactor;
        this->determineLaplacian(deltaPosition);
        this->leapfrog2ndOrder(this->m_phi, this->m_dPhi, this->m_d2Phi, timeDelta, [this, &inflationPotentialDerivative, &structVariables](){ return this->phiSecondDifferential(inflationPotentialDerivative, structVariables); });
        this->potentialKineticEnergyDensity(stepCount, inflationPotential, structVariables);
    }
};
