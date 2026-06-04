#pragma once
#include "Labs/FinalProject/FluidSimulator.h"

namespace VCX::Labs::Final {
    struct CGSimulator : public Simulator {
        virtual void solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift) override;
    };
}