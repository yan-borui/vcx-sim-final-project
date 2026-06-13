#pragma once

#include "Labs/FinalProject/VariationalCoupledSimulator.h"

namespace VCX::Labs::Final {
    // Keep the standalone Variational demo on the same Batty pressure
    // projection used by the coupled, sub-grid, and wall-separation cases.
    struct VariationalSimulator : public VariationalCoupledSimulator {
    };
} // namespace VCX::Labs::Final
