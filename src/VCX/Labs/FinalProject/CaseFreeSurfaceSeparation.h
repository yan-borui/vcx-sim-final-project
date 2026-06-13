#pragma once

#include "Engine/GL/Frame.hpp"
#include "Engine/GL/Program.h"
#include "Engine/Sphere.h"
#include "Labs/Common/ICase.h"
#include "Labs/Common/OrbitCameraManager.h"
#include "Labs/FinalProject/FreeSurfaceSeparationSimulator.h"
#include "Labs/Scene/Content.h"
#include "Labs/Scene/SceneObject.h"

namespace VCX::Labs::FluidSimulation {
    class CaseFreeSurfaceSeparation : public Common::ICase {
    public:
        CaseFreeSurfaceSeparation(
            std::initializer_list<Assets::ExampleScene> && scenes);

        std::string_view const GetName() override {
            return "Free-Surface Wall Separation";
        }
        void OnSetupPropsUI() override;
        Common::CaseRenderResult OnRender(
            std::pair<std::uint32_t, std::uint32_t> const desiredSize) override;
        void OnProcessInput(ImVec2 const & pos) override;

    private:
        void ResetSystem();

        std::vector<Assets::ExampleScene> const _scenes;
        Engine::GL::UniqueProgram              _program;
        Engine::GL::UniqueProgram              _lineProgram;
        Engine::GL::UniqueRenderFrame          _frame;
        Rendering::SceneObject                 _sceneObject;
        Engine::GL::UniqueIndexedRenderItem    _boundaryItem;
        Common::OrbitCameraManager             _cameraManager;

        std::size_t _sceneIdx { 0 };
        bool        _recompute { true };
        bool        _uniformDirty { true };
        bool        _stopped { false };
        float       _timeStep { 0.008f };
        float       _boundaryWidth { 2.0f };
        int         _resolution { 24 };

        Engine::Model                         _sphere;
        std::optional<Rendering::ModelObject>  _fluidParticles;
        Final::FreeSurfaceSeparationSimulator _simulation;
    };
} // namespace VCX::Labs::FluidSimulation
