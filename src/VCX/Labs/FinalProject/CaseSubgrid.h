#pragma once

#include "Engine/GL/Frame.hpp"
#include "Engine/GL/Program.h"
#include "Engine/Sphere.h"
#include "Labs/Common/ICase.h"
#include "Labs/Common/OrbitCameraManager.h"
#include "Labs/FinalProject/RigidBody.h"
#include "Labs/FinalProject/SubgridSimulator.h"
#include "Labs/Scene/Content.h"
#include "Labs/Scene/SceneObject.h"

namespace VCX::Labs::FluidSimulation {
    class CaseSubgrid : public Common::ICase {
    public:
        CaseSubgrid(std::initializer_list<Assets::ExampleScene> && scenes);

        std::string_view const GetName() override { return "Sub-grid Accurate Flow"; }
        void OnSetupPropsUI() override;
        Common::CaseRenderResult OnRender(
            std::pair<std::uint32_t, std::uint32_t> const desiredSize) override;
        void OnProcessInput(ImVec2 const & pos) override;

    private:
        void ResetSystem();

        std::vector<Assets::ExampleScene> const _scenes;

        Engine::GL::UniqueProgram           _program;
        Engine::GL::UniqueProgram           _lineProgram;
        Engine::GL::UniqueProgram           _flatProgram;
        Engine::GL::UniqueRenderFrame       _frame;
        Rendering::SceneObject              _sceneObject;
        Engine::GL::UniqueIndexedRenderItem _boundaryItem;
        Engine::GL::UniqueIndexedRenderItem _rigidBodyItem;
        Common::OrbitCameraManager          _cameraManager;

        std::size_t _sceneIdx { 0 };
        bool        _recompute { true };
        bool        _uniformDirty { true };
        bool        _stopped { false };
        float       _timeStep { 0.016f };
        float       _boundaryWidth { 2.0f };
        int         _resolution { 18 };

        Engine::Model           _sphere;
        std::optional<Rendering::ModelObject> _fluidParticles;
        Final::SubgridSimulator _simulation;
        Final::RigidBody        _body;
    };
} // namespace VCX::Labs::FluidSimulation
