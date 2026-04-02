#pragma once

template <typename Scalar>
struct InflationPotentials{
    static void tModel(Scalar const &phi, Scalar &potential){
        Scalar const c1 = Scalar(17) / Scalar(45);
        Scalar const c2 = Scalar(45) / Scalar(315);

        Scalar const phi2 = phi * phi;
        Scalar const phi4 = phi2 * phi2;
        Scalar const phi6 = phi4 * phi2;

        potential = phi2/Scalar(2) * (Scalar(1) - phi2 * Scalar(2) + phi4 * c1  - phi6 * c2);

    }

    
    static void eModel(Scalar const &phi, Scalar &potential){
        Scalar const c1 = Scalar(7) / Scalar(3);
        Scalar const c2 = Scalar(62) / Scalar(45);

        Scalar const phi2 = phi * phi;
        Scalar const phi3 = phi * phi2;
        Scalar const phi4 = phi3 * phi;

        potential = phi2 / Scalar(2) * (1 - phi / Scalar(2) + phi2 * c1 - phi3 / Scalar(2) + phi4 * c2);

    }

    static void axionCosine(Scalar const &phi, Scalar &potential){
        Scalar const c1 = Scalar(1) / Scalar(12);
        Scalar const c2 = Scalar(1) / Scalar(27);

        Scalar const phi2 = phi * phi;
        Scalar const phi4 = phi2 * phi2;

        potential = phi2 / Scalar(2) * (1 - phi2 * c1 + phi4 * c2);

    }
};


template <typename Scalar>
struct DifferentiatedInflationPotentials {
    static void tModel(Scalar const &phi, Scalar &potentialDerivative){
        Scalar const c1 = Scalar(4) / Scalar(3);
        Scalar const c2 = Scalar(17) / Scalar(15);
        Scalar const c3 = Scalar(284) / Scalar(315);

        Scalar const phi2 = phi * phi;
        Scalar const phi3 = phi2 * phi;
        Scalar const phi5 = phi3 * phi2;
        Scalar const phi7 = phi5 * phi2;

        potentialDerivative = phi - phi3 * c1 + phi5 * c2 - phi7 * c3;
    }

    static void eModel(Scalar const &phi, Scalar &potentialDerivative) {
        Scalar const c1 = Scalar(14) / Scalar(3);
        Scalar const c2 = Scalar(62) / Scalar(15);

        Scalar const phi2 = phi * phi;
        Scalar const phi3 = phi2 * phi;
        Scalar const phi4 = phi3 * phi;
        Scalar const phi5 = phi3 * phi2;

        potentialDerivative = phi - phi2 * Scalar(3) + phi3 * c1 - phi4 * Scalar(5) + phi5 * c2;

}

    static void axionCosine(Scalar const &phi, Scalar &potentialDerivative) {
        Scalar const c1 = Scalar(14) / Scalar(3);
        Scalar const c2 = Scalar(62) / Scalar(15);

        Scalar const phi2 = phi * phi;
        Scalar const phi3 = phi2 * phi;
        Scalar const phi5 = phi3 * phi2;

        potentialDerivative = phi - phi3 / Scalar(6) + phi5 / Scalar(120);
    }

};
    