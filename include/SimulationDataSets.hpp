#pragma once
#include "H5Cpp.h"

using namespace H5;

struct SimulationDataSets {
    DataSet inflatonField;
    DataSet energyDensity;
    DataSet inflatonPotential;
    DataSet averageEnergyDensity;
    DataSet scaleFactor;
    DataSet scaleFactorDerivative;
    DataSet scaleFactorSecondDerivative;
    DataSet hubbleParameter;
};