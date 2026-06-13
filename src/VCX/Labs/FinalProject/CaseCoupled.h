#pragma once

#include <memory>

#include "Engine/GL/Frame.hpp"
#include "Engine/GL/Program.h"
#include "Engine/GL/UniformBlock.hpp"
#include "Engine/loader.h"
#include "Engine/Sphere.h"
#include "Labs/Common/ICase.h"
#include "Labs/Common/ImageRGB.h"
#include "Labs/Common/OrbitCameraManager.h"
#include "Labs/FinalProject/APICSimulator.h"
#include "Labs/FinalProject/CGSimulator.h"
#include "Labs/FinalProject/FluidSimulator.h"
#include "Labs/FinalProject/MeshSDF.h"
#include "Labs/FinalProject/VariationalCoupledSimulator.h"
#include "Labs/Scene/Content.h"
#include "Labs/Scene/SceneObject.h"
#include "RigidBody.h"

namespace VCX::Labs::FluidSimulation {

    class CaseFluid : public Common::ICase {
    public:
        CaseFluid(std::initializer_list<Assets::ExampleScene> && scenes);

        virtual std::string_view const GetName() override { return "Coupled Simulation"; }

        virtual void                     OnSetupPropsUI() override;
        virtual Common::CaseRenderResult OnRender(std::pair<std::uint32_t, std::uint32_t> const desiredSize) override;
        virtual void                     OnProcessInput(ImVec2 const & pos) override;
        void                             OnProcessMouseControl(glm::vec3 mouseDelta);

    private:
        enum class DemoShape {
            Box    = 0,
            Sphere = 1,
            Mesh   = 2
        };

        std::vector<Assets::ExampleScene> const _scenes;

        Engine::GL::UniqueProgram         _program;
        Engine::GL::UniqueProgram         _lineprogram;
        Engine::GL::UniqueRenderFrame     _frame;
        VCX::Labs::Rendering::SceneObject _sceneObject;
        std::size_t                       _sceneIdx { 0 };
        bool                              _recompute { true };
        bool                              _uniformDirty { true };
        int                               _msaa { 2 };
        int                               _useBlinn { 0 };
        float                             _shininess { 32 };
        float                             _ambientScale { 1 };
        bool                              _useGammaCorrection { true };
        int                               _attenuationOrder { 2 };
        int                               _bumpMappingPercent { 20 };

        Engine::GL::UniqueIndexedRenderItem _BoundaryItem;
        Common::OrbitCameraManager          _cameraManager;
        float                               _BndWidth { 2.0f };
        bool                                _stopped { false };
        float                               _dt { 0.016f };
        Engine::Model                       _sphere;
        std::optional<Rendering::ModelObject> _fluidParticles;
        Engine::Model                       _rigidSphere;
        int                                 _res { 16 };
        float                               _r           = 0.0f;
        int                                 numofSpheres = 0;
        Final::Simulator *                  _sim         = nullptr;
        Final::VariationalCoupledSimulator  _simulation;
        Final::APICSimulator                _apicSimulation;
        Final::CGSimulator                  _cgsimulation;
        Final::RigidBody                    _body;
        Engine::GL::UniqueIndexedRenderItem _RigidBodyItem;
        Engine::GL::UniqueIndexedRenderItem _MeshBodyItem;
        Engine::GL::UniqueProgram           _flatProgram;
        bool                                _useAPIC { false };
        bool                                _useCG { false };
        DemoShape                           _demoShape { DemoShape::Box };
        std::shared_ptr<Final::MeshSDF>     _meshSDF;
        bool                                _meshReady { false };

        char const *          GetSceneName(std::size_t const i) const { return VCX::Labs::Rendering::Content::SceneNames[std::size_t(_scenes[i])].c_str(); }
        Engine::Scene const & GetScene(std::size_t const i) const { return VCX::Labs::Rendering::Content::Scenes[std::size_t(_scenes[i])]; }
        void                  ResetSystem();
        void                  LoadMeshIfNeeded();
    };
} // namespace VCX::Labs::FluidSimulation
