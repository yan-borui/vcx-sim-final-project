#pragma once

#include <vector>
#include "Engine/app.h"
#include "Labs/FinalProject/FluidSimulator.h"
#include "Labs/FinalProject/CaseCoupled.h"
#include "Labs/FinalProject/CaseFluidRigid_zly.h"
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
        Final::CaseFluidRigid _caseFluidRigid;
        std::size_t _caseId = 0;
        std::vector<std::reference_wrapper<Common::ICase>> _cases = {
            _caseFluid,
            _caseVariation,
            _caseSubgrid,
            _caseFreeSurfaceSeparation,
            _caseFluidRigid,
        };

    public:
        App();
        void OnFrame() override;
    };
}
