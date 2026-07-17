#pragma once

template <typename Scalar>
struct InflatonSpacetimeVariables {
    Scalar m_averagePhiDerivative, m_averagePressure, m_energyDensity, m_scaleFactor, m_hubbleParameter;
    
    void getInflatonSpacetimeVariables(auto const &inflatonField, auto const &spacetimeParameters){
        m_averagePhiDerivative = inflatonField.getPhiDerivativeMean();
        m_averagePressure = inflatonField.getPressure();
        m_energyDensity = inflatonField.getAverageEnergyDensity();
        m_scaleFactor = spacetimeParameters.getScaleFactor();
        m_hubbleParameter = spacetimeParameters.getHubbleParameter();
    }
};