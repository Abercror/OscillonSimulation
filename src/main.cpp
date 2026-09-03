#include "SpacetimeClass.hpp"
#include "InflatonFieldClass.hpp"
#include "SimulationClass.hpp"
#include "Traits.hpp"
#include "SimulationDataSets.hpp"
#include <H5Cpp.h>
#include <iostream>
#include <fstream>
#include <chrono>


using namespace H5;

std::string fileName = "simulationData.h5";
std::string inflationModel = "T-Model";
int constexpr gridSize = 100;
double constexpr deltaTime = 0.001;
int constexpr totalSteps = 1000;
hsize_t constexpr bufferSize = 100;

using inflatonFieldTypes = VectorTraits<double>;
using spacetimeTypes = ScalarTraits<double>;


SimulationDataSets fileCreation() {

    H5File const file(fileName, H5F_ACC_TRUNC);

    auto constexpr arraySize = gridSize*gridSize*gridSize;
    hsize_t constexpr arrayLength = arraySize;
    hsize_t constexpr numberOfSteps = totalSteps;

    hsize_t fieldDimensions[2] = {
        numberOfSteps, arrayLength
    };

    DataSpace const fieldSpace(2, fieldDimensions);

    DataSet inflatonField = file.createDataSet("/inflatonField", PredType::NATIVE_DOUBLE, fieldSpace);
    DataSet energyDensity = file.createDataSet("/energyDensity", PredType::NATIVE_DOUBLE, fieldSpace);
    DataSet energyOverDensity = file.createDataSet("/energyOverDensity", PredType::NATIVE_DOUBLE, fieldSpace);
    DataSet inflatonPotential = file.createDataSet("/inflatonPotential", PredType::NATIVE_DOUBLE, fieldSpace);

    hsize_t doubleDimensions[1] {
        numberOfSteps
    };

    DataSpace const doubleSpace(1, doubleDimensions);

    DataSet averageEnergyDensity = file.createDataSet("/averageEnergyDensity", PredType::NATIVE_DOUBLE, doubleSpace);
    DataSet scaleFactor = file.createDataSet("/scaleFactor", PredType::NATIVE_DOUBLE, doubleSpace);
    DataSet scaleFactorDerivative = file.createDataSet("/scaleFactorDerivative", PredType::NATIVE_DOUBLE, doubleSpace);
    DataSet scaleFactorSecondDerivative = file.createDataSet("/scaleFactorSecondDerivative", PredType::NATIVE_DOUBLE, doubleSpace);
    DataSet hubbleParameter = file.createDataSet("/hubbleParameter", PredType::NATIVE_DOUBLE, doubleSpace);

    SimulationDataSets data;

    data.inflatonField = inflatonField;
    data.energyDensity = energyDensity;
    data.energyOverDensity = energyOverDensity;
    data.inflatonPotential = inflatonPotential;
    data.averageEnergyDensity = averageEnergyDensity;
    data.scaleFactor = scaleFactor;
    data.scaleFactorDerivative = scaleFactorDerivative;
    data.scaleFactorSecondDerivative = scaleFactorSecondDerivative;
    data.hubbleParameter = hubbleParameter;

    return data;

}



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

    Sim.run(totalSteps, bufferSize);


    std::cout << "Simulation Complete" << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Simulation Time: " << duration.count() << "seconds" << std::endl;

    return 1;
}