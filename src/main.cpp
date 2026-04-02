#include "SpacetimeClass.hpp"
#include "OscillonClass.hpp"
#include "SimulationClass.hpp"
#include "Traits.hpp"
#include <iostream>
#include <fstream>
#include <chrono>


std::string inflationModel = "T-Model";
constexpr int gridSize = 100;
constexpr double deltaTime = 0.001;
constexpr int totalSteps = 1000;

using oscillonTypes = VectorTraits<double>;
using spacetimeTypes = ScalarTraits<double>;


Oscillon<oscillonTypes>initialisingOscillon(){
    oscillonTypes::StateType phi(gridSize, gridSize, gridSize);
    oscillonTypes::StateType dPhi(gridSize, gridSize, gridSize);
    oscillonTypes::StateType d2Phi(gridSize, gridSize, gridSize);
    oscillonTypes::StateType phiLaplacian(gridSize, gridSize, gridSize);
    oscillonTypes::StateType inflationPotential(gridSize, gridSize, gridSize);
    oscillonTypes::StateType inflationPotentialDerivative(gridSize, gridSize, gridSize);
    OscillonData<oscillonTypes::StateType, oscillonTypes::Scalar> history;

    Oscillon<oscillonTypes> oscillonField(phi, dPhi, d2Phi, phiLaplacian, inflationPotential, inflationPotentialDerivative, gridSize, history);

    return oscillonField;
}


SpacetimeParameters<spacetimeTypes> initialisingSpacetime(){
    spacetimeTypes::StateType scaleFactor = 0.0;
    spacetimeTypes::StateType dScaleFactor = 0.0;
    spacetimeTypes::StateType d2ScaleFactor = 0.0;
    spacetimeTypes::StateType hubbleParameter = 0.0;

    SpacetimeParametersData<spacetimeTypes::Scalar> history;

    SpacetimeParameters<spacetimeTypes> spacetime(scaleFactor, dScaleFactor, d2ScaleFactor, hubbleParameter, history);

    return spacetime;
}


int main(){

    std::cout << "Starting" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    auto oscillon = initialisingOscillon();
    auto spacetime = initialisingSpacetime();

    Simulation<oscillonTypes, spacetimeTypes> Sim(oscillon, spacetime, deltaTime, inflationModel);

    Sim.run(totalSteps);

    std::cout << "Simulation Complete" << std::endl;

    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = duration_cast<std::chrono::seconds>(stop - start);

    std::cout << "Simulation Time: " << duration.count() << "seconds" << std::endl;

    return 1;
}