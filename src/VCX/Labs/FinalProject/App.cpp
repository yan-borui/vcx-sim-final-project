#include "Labs/FinalProject/App.h"

namespace VCX::Labs::FluidSimulation {
    App::App():
        _ui(Labs::Common::UIOptions {}),
        _caseFluid({ Assets::ExampleScene::Fluid }),
        _caseSubgrid({ Assets::ExampleScene::Fluid }),
        _caseFreeSurfaceSeparation({ Assets::ExampleScene::Fluid }),
        _caseVariation({ Assets::ExampleScene::Fluid })
        {
    }

    void App::OnFrame() {
        _ui.Setup(_cases, _caseId);
    }
}
