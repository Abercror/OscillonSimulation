#include "SpacetimeClass.hpp"
#include "InflatonFieldClass.hpp"
#include "SimulationClass.hpp"
#include "Traits.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include "GlobalVariables.hpp"

using namespace H5;

std::string inflationModel = "T-Model";
constexpr int gridSize = 100;
constexpr double deltaTime = 0.001;
constexpr int totalSteps = 1000;

using inflatonFieldTypes = VectorTraits<double>;
using spacetimeTypes = ScalarTraits<double>;


InflatonField<inflatonFieldTypes>initialisingInflatonField(){
    inflatonFieldTypes::StateType phi(gridSize, gridSize, gridSize);
    inflatonFieldTypes::StateType dPhi(gridSize, gridSize, gridSize);
    inflatonFieldTypes::StateType d2Phi(gridSize, gridSize, gridSize);
    inflatonFieldTypes::StateType phiLaplacian(gridSize, gridSize, gridSize);
    inflatonFieldTypes::StateType inflationPotential(gridSize, gridSize, gridSize);
    inflatonFieldTypes::StateType inflationPotentialDerivative(gridSize, gridSize, gridSize);
    inflatonFieldTypes::StateType energyDensity(gridSize, gridSize, gridSize);
    inflatonFieldTypes::StateType energyOverdensity(gridSize, gridSize, gridSize);
    inflatonFieldTypes::StateType kineticEnergyDensity(gridSize, gridSize, gridSize);
    inflatonFieldTypes::StateType potentialEnergyDensity(gridSize, gridSize, gridSize);
    inflatonFieldData<inflatonFieldTypes::StateType, inflatonFieldTypes::Scalar> history;

    InflatonField<inflatonFieldTypes> inflatonField(phi, dPhi, d2Phi, phiLaplacian, inflationPotential, inflationPotentialDerivative, energyDensity, energyOverdensity, kineticEnergyDensity, potentialEnergyDensity, gridSize, history);

    return inflatonField;
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


int fileCreation() {

    try {
        Exception::dontPrint();

        H5File file(FILE_NAME, H5F_ACC_TRUNC);

        Group group(file.createGroup(GROUP_NAME));
    }

    catch (FileIException &error) {
        error.printErrorStack();
        return -1;
    }

    catch (GroupIException &error) {
        error.printErrorStack();
        return -1;
    }

    return 0;
}

// int dataDetCreation() {
//     try {
//         Exception::dontPrint();
//
//         hsize_t arrayLength = gridSize*gridSize*gridSize;
//
//         hsize_t dims[1] = {
//             arrayLength
//         };
//
//         DataSpace dataspace(1, arrayLength);
//
//         H5File file(FILE_NAME, H5F_ACC_RDWR);
//         DataSet dataset = file.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace);
//
//     }
// }


int main(){

    std::cout << "Starting" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    auto inflatonField = initialisingInflatonField();
    auto spacetime = initialisingSpacetime();


    Simulation<inflatonFieldTypes, spacetimeTypes> Sim(inflatonField, spacetime, deltaTime, inflationModel);

    Sim.run(totalSteps);


    std::cout << "Simulation Complete" << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Simulation Time: " << duration.count() << "seconds" << std::endl;

    return 1;
}