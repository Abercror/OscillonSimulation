#pragma once
#include "NumericalMethodsClass.hpp"
#include <vector>
#include <functional>
#include <iostream>


template <typename StateType, typename Scalar>
struct OscillonData{
    std::vector<StateType> phiValues;
    std::vector<StateType> dPhiValues;
    std::vector<StateType> d2PhiValues;
    std::vector<StateType> phiLaplacianValues;
    std::vector<StateType> inflationPotentialValues;
    std::vector<StateType> inflationPotentialDerivativeValues;
    std::vector<Scalar> energyDensityValues;
    std::vector<Scalar> pressureValues;
    
    void reserveMemory(Scalar const &totalSteps){
        phiValues.reserve(totalSteps);
        dPhiValues.reserve(totalSteps);
        phiLaplacianValues.reserve(totalSteps);
        inflationPotentialValues.reserve(totalSteps);
        inflationPotentialDerivativeValues.reserve(totalSteps);
        energyDensityValues.reserve(totalSteps);
        pressureValues.reserve(totalSteps);
    }
};


template <typename Traits>
class Oscillon: public NumericalMethods<Traits> {
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
    Scalar m_energyDensity;
    // Scalar m_kineticEnergyDensity;
    // Scalar m_potentialEnergyDensity;
    Scalar m_pressure;
    Scalar m_gridSize;
    OscillonData<StateType, Scalar> m_history;

public:
    Oscillon(
        StateType phi,
        StateType dPhi,
        StateType d2Phi,
        StateType phiLaplacian,
        StateType inflationPotential,
        StateType inflationPotentialDerivative,
        Scalar gridSize,
        OscillonData<StateType, Scalar> history): m_phi(phi), m_dPhi(dPhi), m_d2Phi(d2Phi), m_phiLaplacian(phiLaplacian), m_inflationPotential(inflationPotential), m_inflationPotentialDerivative(inflationPotential), m_gridSize(gridSize), m_history(history) {}


    Scalar getPhiDerivativeMean() const {
        auto averagePhiDerivative = this->m_dPhi.mean();
        return averagePhiDerivative;
    }
    Scalar getPressure() const { return m_pressure; }
    Scalar getEnergyDensity() const { return m_energyDensity; }

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
        for (int i = 0; i < phi.field.size(); ++i) {
            secondDerivative[i] = laplacianPhi[i] / scaleFactor2 - hubbleTerm * dPhi[i] - inflationPotentialDerivative[i];

        }
    }

    void determineLaplacian(Scalar const &deltaPosition){
        StateType const &phi = this->m_phi;
        StateType &laplacian = this->m_phiLaplacian;

        #pragma omp parallel for
        for(int z = 1; z < this->m_gridSize-1; ++z){
            for(int y = 1; y < this->m_gridSize-1; ++y){
                for (int x = 1; x < this->m_gridSize-1; ++x)
                laplacian(x, y, z) = ((phi(x+1, y, z) - Scalar(2) * phi(x, y , z) + phi(x-1, y, z)) / deltaPosition) + (Scalar(2) / (x * deltaPosition)) * ((phi(x+1, y, z) - phi(x-1, y, z))/(deltaPosition));
            }
        }
    }

    void potentialKineticEnergyDensity(Scalar const &stepCount, Func const &inflationPotentialFunc){
        StateType const &phi = this->m_phi;
        StateType &inflationPotential = this->m_inflationPotential;

        inflationPotentialFunc(phi, inflationPotential);

        this->m_energyDensity = this->m_inflationPotential.mean();

        // kineticEnergyDensity = averagePotential * std::pow(std::cos(frequency * time + phaseDifference), Scalar(2));
        // potentialEnergyDensity = averagePotential * std::pow(std::sin(frequency * time + phaseDifference), Scalar(2));
        // energyDensity = kineticEnergyDensity + potentialEnergyDensity;

        // this->m_kineticEnergyDensity = kineticEnergyDensity;
        // this->m_potentialEnergyDensity = potentialEnergyDensity;
    }

    void writeToHistory(){
        this->m_history.phiValues.push_back(this->m_phi);
        this->m_history.dPhiValues.push_back(this->m_dPhi);
        this->m_history.phiLaplacianValues.push_back(this->m_phiLaplacian);
        this->m_history.pressureValues.push_back(this->m_pressure);
        this->m_history.energyDensityValues.push_back(this->m_energyDensity);
        this->m_history.inflationPotentialValues.push_back(this->m_inflationPotential);
        this->m_history.inflationPotentialDerivativeValues.push_back(this->m_inflationPotentialDerivative);
        // this->m_history.kineticEnergyDensity.push_back(this->m_kineticEnergyDensity);
        // this->m_history.potentialEnergyDensity.push_back(this->m_potentialEnergyDensity);
    }

    void updateOscillon(Scalar const &stepCount, Scalar const &timeDelta, Func &inflationPotential, Func &inflationPotentialDerivative, auto const &structVariables){
        auto deltaPosition = Scalar(1)/(this->m_gridSize);
        this->determineLaplacian(deltaPosition);
        this->leapfrog2ndOrder(this->m_phi, this->m_dPhi, this->m_d2Phi, timeDelta, [this, &inflationPotentialDerivative, &structVariables](){ return this->phiSecondDifferential(inflationPotentialDerivative, structVariables); });
        this->potentialKineticEnergyDensity(stepCount, inflationPotential);
        this->writeToHistory();
    }
};
