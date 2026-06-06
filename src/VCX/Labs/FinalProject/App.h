#pragma once

#include <vector>
#include "Engine/app.h"
#include "Labs/FinalProject/FluidSimulator.h"
#include "Labs/FinalProject/CaseCoupled.h"
#include "Labs/FinalProject/CaseSubgrid.h"
#include "Labs/Common/UI.h"

namespace VCX::Labs::FluidSimulation {
    class App : public Engine::IApp {
    private:
        Common::UI _ui;
        CaseFluid _caseFluid;
        CaseSubgrid _caseSubgrid;
        std::size_t _caseId = 0;
        std::vector<std::reference_wrapper<Common::ICase>> _cases = { _caseFluid, _caseSubgrid };

    public:
        App();
        void OnFrame() override;
    };
}
