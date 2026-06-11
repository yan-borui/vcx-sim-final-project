#pragma once

#include <vector>
#include "Engine/app.h"
#include "Labs/FinalProject/FluidSimulator.h"
#include "Labs/FinalProject/CaseCoupled.h"
#include "Labs/FinalProject/CaseFreeSurfaceSeparation.h"
#include "Labs/FinalProject/CaseSubgrid.h"
#include "Labs/FinalProject/CaseVariation.h"
#include "Labs/Common/UI.h"

namespace VCX::Labs::FluidSimulation {
    class App : public Engine::IApp {
    private:
        Common::UI _ui;
        CaseFluid _caseFluid;
        CaseSubgrid _caseSubgrid;
        CaseFreeSurfaceSeparation _caseFreeSurfaceSeparation;
        CaseVariation _caseVariation;
        std::size_t _caseId = 0;
        std::vector<std::reference_wrapper<Common::ICase>> _cases = {
            _caseVariation,
            _caseFluid,
            _caseSubgrid,
            _caseFreeSurfaceSeparation,
        };

    public:
        App();
        void OnFrame() override;
    };
}
