#pragma once
#include "Engine/GL/Frame.hpp"
#include "Engine/GL/Program.h"
#include "Engine/GL/RenderItem.h"
#include "Labs/FinalProject/FluidSimulator_zly.h"
#include "Labs/Common/ICase.h"
#include "Labs/Common/ImageRGB.h"
#include "Labs/Common/OrbitCameraManager.h"

namespace VCX::Labs::Final {
    class CaseFluidRigid : public Common::ICase {
    public:
        CaseFluidRigid();
        std::string_view const   GetName() override { return "Fluid-Rigid Coupling 3D"; }
        void                     OnSetupPropsUI() override;
        Common::CaseRenderResult OnRender(std::pair<std::uint32_t, std::uint32_t> const desiredSize) override;
        void                     OnProcessInput(ImVec2 const & pos) override;

    private:
        Engine::GL::UniqueProgram           _program;
        Engine::GL::UniqueRenderFrame       _frame;
        Engine::Camera                      _camera { .Eye = glm::vec3(3, 4, 5) };
        Common::OrbitCameraManager          _cameraManager;
        Engine::GL::UniqueIndexedRenderItem _triItem, _lineItem;

        FluidSimulator3D                        _sim;
        bool                                    _running = false, _initialized = false;
        int                                     _substeps = 2, _gridRes = 20;
        float                                   _domainW = 2.0f;
        std::pair<std::uint32_t, std::uint32_t> _frameSize { 800, 600 };

        std::vector<glm::vec3> _verts, _lineVerts;
        std::vector<uint32_t>  _triIdx, _lineIdx;

        void ResetSimulation();
        void BuildFluidSurface();
        void AddSphere(glm::vec3 c, float r, int seg = 10);
        void AddWireBox(glm::vec3 lo, glm::vec3 hi);
        void AddQuadFace(glm::vec3 c, float hs, int axis, int side);
    };
} // namespace VCX::Labs::GettingStarted
