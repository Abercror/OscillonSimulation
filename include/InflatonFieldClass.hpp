#pragma once
#include "NumericalMethodsClass.hpp"
#include "SimulationDataSets.hpp"
#include <vector>
#include <array>
#include <functional>
#include <cstddef>
#include <H5Cpp.h>
using namespace H5;

template <typename StateType, typename Scalar>
struct inflatonFieldData{
    std::vector<StateType> phiValues;
    std::vector<StateType> inflationPotentialValues;
    std::vector<StateType> energyDensityValues;
    std::vector<StateType> energyOverdensityValues;
    std::vector<Scalar> averageEnergyDensity;

    // std::vector<StateType> dPhiValues;
    // std::vector<StateType> d2PhiValues;
    // std::vector<StateType> phiLaplacianValues;
    // std::vector<StateType> inflationPotentialDerivativeValues;
    // std::vector<Scalar> pressureValues;

    
    void resizeBuffer(hsize_t const &bufferSize){
        phiValues.resize(bufferSize);
        inflationPotentialValues.resize(bufferSize);
        energyDensityValues.resize(bufferSize);
        energyOverdensityValues.resize(bufferSize);
        averageEnergyDensity.resize(bufferSize);

        // dPhiValues.resize(bufferSize);
        // phiLaplacianValues.resize(bufferSize);
        // inflationPotentialDerivativeValues.resize(bufferSize);
        // pressureValues.resize(bufferSize);
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
    InflatonField(Scalar gridSize, SimulationDataSets &dataSets) {
        StateType phi(gridSize, gridSize, gridSize);
        StateType dPhi(gridSize, gridSize, gridSize);
        StateType d2Phi(gridSize, gridSize, gridSize);
        StateType phiLaplacian(gridSize, gridSize, gridSize);
        StateType inflationPotential(gridSize, gridSize, gridSize);
        StateType inflationPotentialDerivative(gridSize, gridSize, gridSize);
        StateType energyDensity(gridSize, gridSize, gridSize);
        Scalar averageEnergyDensity = 0;
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


    Scalar getPhiDerivativeMean() const {
        auto averagePhiDerivative = this->m_dPhi.mean();
        return averagePhiDerivative;
    }
    Scalar getPressure() const { return m_pressure; }
    Scalar getAverageEnergyDensity() const { return m_averageEnergyDensity; }
    inflatonFieldData<StateType, Scalar> getBuffer() const { return this->m_buffer; }

    void setLength(hsize_t const &bufferSize){
        this->m_buffer.resizeBuffer(bufferSize);
    }


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


    void writeToBuffer(auto const &index) {
        this->m_buffer.phiValues[index] = this->m_phi;
        this->m_buffer.energyDensityValues[index] = this->m_energyDensity;
        this->m_buffer.energyOverdensityValues[index] = this->m_energyOverDensity;
        this->m_buffer.inflationPotentialValues[index] = this->m_inflationPotential;
        this->m_buffer.averageEnergyDensity[index] = this->m_averageEnergyDensity;

        // this->m_buffer.dPhiValues.push_back(this->m_dPhi);
        // this->m_buffer.phiLaplacianValues.push_back(this->m_phiLaplacian);
        // this->m_buffer.pressureValues.push_back(this->m_pressure);
        // this->m_buffer.inflationPotentialDerivativeValues.push_back(this->m_inflationPotentialDerivative);
        // this->m_buffer.kineticEnergyDensity.push_back(this->m_kineticEnergyDensity);
        // this->m_buffer.potentialEnergyDensity.push_back(this->m_potentialEnergyDensity);
    }


    void writeFieldsToFile(DataSet &dataSet, auto &buffer, auto const &fieldStart, auto const &fieldCount) {
        hsize_t start[2] = {fieldStart[0], fieldStart[1]};
        auto size = buffer.size();

        for (std::size_t i = 0; i < size; ++i) {
            DataSpace fileSpace = dataSet.getSpace();
            fileSpace.selectHyperslab(H5S_SELECT_SET, fieldCount, start);
            DataSpace memorySpace(2, fieldCount);
            dataSet.write(buffer[i].data(), PredType::NATIVE_DOUBLE, memorySpace, fileSpace);

            start[0] +=1;
        }
    }

    void testing(DataSet &dataSet, auto &buffer, auto const &fieldStart, auto const &fieldCount, auto const &name) {
        DataSpace fileSpace = dataSet.getSpace();

        hsize_t dims[2];
        fileSpace.getSimpleExtentDims(dims);

        std::cout
            << " name = " << name
            << " dims = [" << dims[0] << ", " << dims[1] << "]"
            << " start = [" << fieldStart[0] << ", " << fieldStart[1] << "]"
            << " count = [" << fieldCount[0] << ", " << fieldCount[1] << "]"
            << std::endl;
    }

    void writeDoubleToFile(DataSet &dataSet, auto &buffer, auto const &start, auto const &count) {
        DataSpace fileSpace = dataSet.getSpace();
        fileSpace.selectHyperslab(H5S_SELECT_SET, count,  start);
        DataSpace memorySpace(1, count);
        dataSet.write(buffer.data(), PredType::NATIVE_DOUBLE, memorySpace, fileSpace);
    }

    void writeToFile(auto const &fieldStart, auto const &fieldCount, auto const &doubleStart, auto const &doubleCount) {
        auto data = this->m_dataSets;
        auto buffer = this->m_buffer;
        writeFieldsToFile(data.inflatonField, buffer.phiValues, fieldStart, fieldCount);
        // testing(data.inflatonField, buffer.phiValues, fieldStart,  fieldCount, "inflaton field");
        writeFieldsToFile(data.energyDensity, buffer.energyDensityValues, fieldStart, fieldCount);
        // testing(data.energyDensity, buffer.energyDensityValues, fieldStart,  fieldCount, "energy density");
        writeFieldsToFile(data.inflatonPotential, buffer.inflationPotentialValues, fieldStart, fieldCount);
        // testing(data.inflatonPotential, buffer.inflationPotentialValues, fieldStart,  fieldCount, "potential");
        writeDoubleToFile(data.averageEnergyDensity, buffer.averageEnergyDensity, doubleStart, doubleCount);
    }


    void updateInflatonField(Scalar const &stepCount, Scalar const &timeDelta, Func &inflationPotential, Func &inflationPotentialDerivative, auto const &structVariables){
        auto deltaPosition = Scalar(1)/(this->m_gridSize) * structVariables.m_scaleFactor;
        this->determineLaplacian(deltaPosition);
        this->leapfrog2ndOrder(this->m_phi, this->m_dPhi, this->m_d2Phi, timeDelta, [this, &inflationPotentialDerivative, &structVariables](){ return this->phiSecondDifferential(inflationPotentialDerivative, structVariables); });
        this->potentialKineticEnergyDensity(stepCount, inflationPotential, structVariables);
    }
};
