#include "SpacetimeClass.hpp"
#include "InflatonFieldClass.hpp"
#include "SimulationClass.hpp"
#include "Traits.hpp"
#include "SimulationDataSets.hpp"
#include <H5Cpp.h>
#include <iostream>
#include <fstream>
#include <chrono>


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Simulation Configuration
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string fileName = "simulationData.h5";
std::string inflationModel = "T-Model";
int constexpr gridSize = 100;
double constexpr deltaTime = 0.001;
int constexpr totalSteps = 1000;
hsize_t constexpr bufferSize = 10;
hsize_t constexpr writingInterval = 1;
double initialPhiValue = 10;

using inflatonFieldTypes = VectorTraits<double>;
using spacetimeTypes = ScalarTraits<double>;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// HDF5 File Setup
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


SimulationDataSets fileCreation() {

    H5::H5File const file(fileName, H5F_ACC_TRUNC);

    auto constexpr arraySize = gridSize*gridSize*gridSize;
    hsize_t constexpr arrayLength = arraySize;
    hsize_t constexpr numberOfTimeSteps = totalSteps/writingInterval;

    hsize_t fieldDimensions[2] = {
        numberOfTimeSteps, arrayLength
    };

    H5::DataSpace const fieldSpace(2, fieldDimensions);

    H5::DataSet inflatonField = file.createDataSet("/inflatonField", H5::PredType::NATIVE_DOUBLE, fieldSpace);
    H5::DataSet energyDensity = file.createDataSet("/energyDensity", H5::PredType::NATIVE_DOUBLE, fieldSpace);
    H5::DataSet inflatonPotential = file.createDataSet("/inflatonPotential", H5::PredType::NATIVE_DOUBLE, fieldSpace);

    hsize_t doubleDimensions[1] {
        numberOfTimeSteps
    };

    H5::DataSpace const doubleSpace(1, doubleDimensions);

    H5::DataSet averageEnergyDensity = file.createDataSet("/averageEnergyDensity", H5::PredType::NATIVE_DOUBLE, doubleSpace);
    H5::DataSet scaleFactor = file.createDataSet("/scaleFactor", H5::PredType::NATIVE_DOUBLE, doubleSpace);
    H5::DataSet scaleFactorDerivative = file.createDataSet("/scaleFactorDerivative", H5::PredType::NATIVE_DOUBLE, doubleSpace);
    H5::DataSet scaleFactorSecondDerivative = file.createDataSet("/scaleFactorSecondDerivative", H5::PredType::NATIVE_DOUBLE, doubleSpace);
    H5::DataSet hubbleParameter = file.createDataSet("/hubbleParameter", H5::PredType::NATIVE_DOUBLE, doubleSpace);

    SimulationDataSets data;

    data.inflatonField = inflatonField;
    data.energyDensity = energyDensity;
    data.inflatonPotential = inflatonPotential;
    data.averageEnergyDensity = averageEnergyDensity;
    data.scaleFactor = scaleFactor;
    data.scaleFactorDerivative = scaleFactorDerivative;
    data.scaleFactorSecondDerivative = scaleFactorSecondDerivative;
    data.hubbleParameter = hubbleParameter;

    return data;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Main Function
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(){

    std::cout << "Starting" << std::endl;
    auto const start = std::chrono::high_resolution_clock::now();

    SimulationDataSets dataSets = fileCreation();

    std::cout << "created files" << std::endl;

    InflatonField<inflatonFieldTypes> inflatonField(gridSize, dataSets);
    SpacetimeParameters<spacetimeTypes> spacetime(dataSets);

    std::cout << "initialised objects" << std::endl;

    Simulation<inflatonFieldTypes, spacetimeTypes> Sim(inflatonField, spacetime, deltaTime, inflationModel, gridSize);

    std::cout << "initialised simulation" << std::endl;

    Sim.run(totalSteps, initialPhiValue, bufferSize, writingInterval);

    std::cout << "Simulation Complete" << std::endl;
    const auto stop = std::chrono::high_resolution_clock::now();
    const auto duration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Simulation Time: " << duration.count() << "seconds" << std::endl;

    return 1;
}