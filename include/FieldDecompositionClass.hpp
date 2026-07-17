#pragma once
#include "NumericalMethodsClass.hpp"
#include <cstddef>

template <typename StateType, typename Scalar>
struct DecomposedFieldData {


};

template <typename Traits>
class DecomposedField : public NumericalMethods<Traits> {
protected:
    using Scalar = typename Traits::Scalar;
    using Func = typename Traits::Func;
    using StateType = typename Traits::StateType;

private:
    StateType m_deltaPhi;
    StateType m_dDeltaPhi;
    StateType m_d2DeltaPhi;
    DecomposedFieldData<StateType, Scalar> m_history;

public:
    DecomposedField(StateType const deltaPhi, StateType const dDeltaPhi, StateType const d2DeltaPhi): m_deltaPhi(deltaPhi), m_dDeltaPhi(dDeltaPhi), m_d2DeltaPhi(d2DeltaPhi) {}


    void fieldDecomposition(StateType const &phi, StateType const &dPhi, StateType const &d2Phi) {
        Scalar const averagePhi = phi.mean();
        Scalar const averageDPhi = dPhi.mean();
        Scalar const averageD2Phi = d2Phi.mean();

        #pragma omp parallel for
        for (std::size_t i = 0; i < phi.size(); ++i) {
            this->m_deltaPhi[i] = phi[i] - averagePhi;
            this->m_dDeltaPhi[i] = dPhi[i] - averageDPhi;
            this->m_d2DeltaPhi[i] = d2Phi[i] - averageD2Phi;
        }
    }

    void energyOverdensity(Scalar const &energyDensity) {
        
    }

};