#pragma once


template <typename Scalar>
struct OscillonSpacetimeVariables {
    Scalar m_averagePhiDerivative, m_averagePressure, m_energyDensity, m_scaleFactor, m_hubbleParameter;
    
    void getOscillonSpacetimeVariables(auto const &oscillon, auto const &spacetimeParameters){
        m_averagePhiDerivative = oscillon.getPhiDerivativeMean();
        m_averagePressure = oscillon.getPressure();
        m_energyDensity = oscillon.getEnergyDensity();
        m_scaleFactor = spacetimeParameters.getScaleFactor();
        m_hubbleParameter = spacetimeParameters.getHubbleParameter();
    }
};