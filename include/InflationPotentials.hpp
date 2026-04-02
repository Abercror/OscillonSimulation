#pragma once

template <typename StateType, typename Scalar>
struct InflationPotentials{
    static void tModel(StateType const &phiValues, StateType &potential){
        Scalar const c1 = Scalar(17) / Scalar(45);
        Scalar const c2 = Scalar(45) / Scalar(315);

        #pragma omp parallel for
        for (int i = 0; i < phiValues.field.size(); ++i) {
            Scalar const phi = phiValues[i];
            Scalar const phi2 = phi * phi;
            Scalar const phi4 = phi2 * phi2;
            Scalar const phi6 = phi4 * phi2;

            potential[i] = phi2/Scalar(2) * (Scalar(1) - phi2 * Scalar(2) + phi4 * c1  - phi6 * c2);
        }
    }

    
    static void eModel(StateType const &phiValues, StateType &potential){
        Scalar const c1 = Scalar(7) / Scalar(3);
        Scalar const c2 = Scalar(62) / Scalar(45);

        #pragma omp parallel for
        for (int i = 0; i < phiValues.field.size(); ++i) {
            Scalar const phi = phiValues[i];
            Scalar const phi2 = phi * phi;
            Scalar const phi3 = phi * phi2;
            Scalar const phi4 = phi3 * phi;

            potential[i] = phi2 / Scalar(2) * (1 - phi / Scalar(2) + phi2 * c1 - phi3 / Scalar(2) + phi4 * c2);
        }
    }

    static void axionCosine(StateType const &phiValues, StateType &potential){
        Scalar const c1 = Scalar(1) / Scalar(12);
        Scalar const c2 = Scalar(1) / Scalar(27);

        #pragma omp parallel for
        for (int i = 0; i < phiValues.field.size(); ++i) {
            Scalar const phi = phiValues[i];
            Scalar const phi2 = phi * phi;
            Scalar const phi4 = phi2 * phi2;

            potential[i] = phi2 / Scalar(2) * (1 - phi2 * c1 + phi4 * c2);
        }
    }
};


template <typename StateType, typename Scalar>
struct DifferentiatedInflationPotentials{
    static void tModel(StateType const &phiValues, StateType &potentialDerivative){
        Scalar const c1 = Scalar(4) / Scalar(3);
        Scalar const c2 = Scalar(17) / Scalar(15);
        Scalar const c3 = Scalar(284) / Scalar(315);

        #pragma omp parallel for
        for (int i = 0; i < phiValues.field.size(); ++i) {
            Scalar const phi = phiValues[i];
            Scalar const phi2 = phi * phi;
            Scalar const phi3 = phi2 * phi;
            Scalar const phi5 = phi3 * phi2;
            Scalar const phi7 = phi5 * phi2;

            potentialDerivative[i] = phi - phi3 * c1 + phi5 * c2 - phi7 * c3;
        }
    }

    static void eModel(StateType const &phiValues, StateType &potentialDerivative) {
        Scalar const c1 = Scalar(14) / Scalar(3);
        Scalar const c2 = Scalar(62) / Scalar(15);

        #pragma omp parallel for
        for (int i = 0; i < phiValues.field.size(); ++i) {
            Scalar const phi = phiValues[i];
            Scalar const phi2 = phi * phi;
            Scalar const phi3 = phi2 * phi;
            Scalar const phi4 = phi3 * phi;
            Scalar const phi5 = phi3 * phi2;

            potentialDerivative[i] = phi - phi2 * Scalar(3) + phi3 * c1 - phi4 * Scalar(5) + phi5 * c2;
        }
    }

    static void axionCosine(StateType const &phiValues, StateType &potentialDerivative) {
        Scalar const c1 = Scalar(14) / Scalar(3);
        Scalar const c2 = Scalar(62) / Scalar(15);

        #pragma omp parallel for
        for (int i = 0; i < phiValues.field.size(); ++i) {
            Scalar const phi = phiValues[i];
            Scalar const phi2 = phi * phi;
            Scalar const phi3 = phi2 * phi;
            Scalar const phi5 = phi3 * phi2;

            potentialDerivative[i] = phi - phi3 / Scalar(6) + phi5 / Scalar(120);
        }
    }
};
    