#pragma once
#include "NumericalMethodsClass.hpp"
#include "GlobalVariables.hpp"
#include <vector>
#include <functional>
#include <cstddef>
#include <H5Cpp.h>
using namespace H5;

template <typename StateType, typename Scalar>
struct inflatonFieldData{
    std::vector<StateType> phiValues;
    std::vector<StateType> dPhiValues;
    std::vector<StateType> d2PhiValues;
    std::vector<StateType> phiLaplacianValues;
    std::vector<StateType> inflationPotentialValues;
    std::vector<StateType> inflationPotentialDerivativeValues;
    std::vector<StateType> energyDensityValues;
    std::vector<StateType> energyOverdensityValues;
    std::vector<Scalar> pressureValues;
    
    void reserveMemory(Scalar const &totalSteps){
        phiValues.reserve(totalSteps);
        dPhiValues.reserve(totalSteps);
        phiLaplacianValues.reserve(totalSteps);
        inflationPotentialValues.reserve(totalSteps);
        inflationPotentialDerivativeValues.reserve(totalSteps);
        energyDensityValues.reserve(totalSteps);
        energyOverdensityValues.reserve(totalSteps);
        pressureValues.reserve(totalSteps);
    }
};


template <typename Traits>
class InflatonField: public NumericalMethods<Traits> {
protected:
    using Scalar = typename Traits::Scalar;
    using Func = typename Traits::Func;
    using StateType = typename Traits::StateType;

private:
    StateType m_phi;
    StateType m_dPhi;
    StateType m_d2Phi;
    StateType m_phiLaplacian;
    StateType m_inflationPotential;
    StateType m_inflationPotentialDerivative;
    StateType m_energyDensity;
    StateType m_energyOverdensities;
    Scalar m_averageEnergyDensity;
    StateType m_kineticEnergyDensity;
    StateType m_potentialEnergyDensity;
    Scalar m_pressure;
    Scalar m_gridSize;
    inflatonFieldData<StateType, Scalar> m_history;

public:
    InflatonField(
        StateType phi,
        StateType dPhi,
        StateType d2Phi,
        StateType phiLaplacian,
        StateType inflationPotential,
        StateType inflationPotentialDerivative,
        StateType energyDensity,
        StateType energyOverdensity,
        StateType kineticEnergyDensity,
        StateType potentialEnergyDensity,
        Scalar gridSize,
        inflatonFieldData<StateType, Scalar> history):
            m_phi(phi),
            m_dPhi(dPhi),
            m_d2Phi(d2Phi),
            m_phiLaplacian(phiLaplacian),
            m_inflationPotential(inflationPotential),
            m_inflationPotentialDerivative(inflationPotentialDerivative),
            m_energyDensity(energyDensity),
            m_energyOverdensities(energyOverdensity),
            m_kineticEnergyDensity(kineticEnergyDensity),
            m_potentialEnergyDensity(potentialEnergyDensity),
            m_gridSize(gridSize),
            m_history(history) {}


    Scalar getPhiDerivativeMean() const {
        auto averagePhiDerivative = this->m_dPhi.mean();
        return averagePhiDerivative;
    }
    Scalar getPressure() const { return m_pressure; }
    Scalar getAverageEnergyDensity() const { return m_averageEnergyDensity; }

    void setLength(Scalar const &totalSteps){
        this->m_history.reserveMemory(totalSteps);
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
            this->m_energyOverdensities[i] = this->m_energyDensity[i] / this->m_averageEnergyDensity;
        }

        #pragma omp parallel for
        for (std::size_t i = 0; i < phi.size(); ++i) {
            this->m_kineticEnergyDensity[i] = m_inflationPotential[i] * std::cos(stepCount) * std::cos(stepCount);
            this->m_potentialEnergyDensity[i] = m_energyDensity[i] - m_kineticEnergyDensity[i];
        }
    }


    void writeToHistory() {
        // this->m_history.phiValues.push_back(this->m_phi);
        // this->m_history.dPhiValues.push_back(this->m_dPhi);
        // this->m_history.phiLaplacianValues.push_back(this->m_phiLaplacian);
        // this->m_history.pressureValues.push_back(this->m_pressure);
        this->m_history.energyDensityValues.push_back(this->m_energyDensity);
        this->m_history.energyOverdensityValues.push_back(this->m_energyOverdensities);
        // this->m_history.inflationPotentialValues.push_back(this->m_inflationPotential);
        // this->m_history.inflationPotentialDerivativeValues.push_back(this->m_inflationPotentialDerivative);
        // this->m_history.kineticEnergyDensity.push_back(this->m_kineticEnergyDensity);
        // this->m_history.potentialEnergyDensity.push_back(this->m_potentialEnergyDensity);
    }


    int writeFieldsToFile(Scalar const &bufferSize) {

        int data[bufferSize];

        for (int i = 0; i < bufferSize+1; i++) {
            data[i] = this->m_history.phiValues[i];
        }

        try {
            Exception::dontPrint();

            H5File file(FILE_NAME, H5F_ACC_RDWR);

            DataSet dataSet = file.openDataSet(DATASET_NAME);

            dataSet.write(data, PredType::NATIVE_DOUBLE);
        }

        catch (FileIException &error) {
            error.printErrorStack();
            return -1;
        }

        // catch failure caused by the DataSet operations
        catch (DataSetIException &error) {
            error.printErrorStack();
            return -1;
        }

        return 0;
    }


    void updateInflatonField(Scalar const &stepCount, Scalar const &timeDelta, Func &inflationPotential, Func &inflationPotentialDerivative, auto const &structVariables, int const bufferSize){
        auto deltaPosition = Scalar(1)/(this->m_gridSize) * structVariables.m_scaleFactor;
        this->determineLaplacian(deltaPosition);
        this->leapfrog2ndOrder(this->m_phi, this->m_dPhi, this->m_d2Phi, timeDelta, [this, &inflationPotentialDerivative, &structVariables](){ return this->phiSecondDifferential(inflationPotentialDerivative, structVariables); });
        this->potentialKineticEnergyDensity(stepCount, inflationPotential, structVariables);
        this->writeToHistory();
        if (stepCount % bufferSize == 0) {
            int readWrite = this->writeFieldsToFile(bufferSize);
            if (readWrite != 0) {
                std::cout << "Failed to write Data" << std::endl;
            }
        }
    }
};
