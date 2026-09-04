#pragma once
#include "H5Cpp.h"



struct SimulationDataSets {
    H5::DataSet inflatonField;
    H5::DataSet energyDensity;
    H5::DataSet inflatonPotential;
    H5::DataSet averageEnergyDensity;
    H5::DataSet scaleFactor;
    H5::DataSet scaleFactorDerivative;
    H5::DataSet scaleFactorSecondDerivative;
    H5::DataSet hubbleParameter;
};