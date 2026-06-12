#include "Labs/FinalProject/CaseFreeSurfaceSeparation.h"

#include "Engine/app.h"

namespace VCX::Labs::FluidSimulation {
    namespace {
        constexpr std::uint32_t SeparationPassConstantsBinding = 3;

        const std::vector<glm::vec3> BoundaryVertices = {
            { -0.5f, -0.5f, -0.5f },
            {  0.5f, -0.5f, -0.5f },
            {  0.5f,  0.5f, -0.5f },
            { -0.5f,  0.5f, -0.5f },
            { -0.5f, -0.5f,  0.5f },
            {  0.5f, -0.5f,  0.5f },
            {  0.5f,  0.5f,  0.5f },
            { -0.5f,  0.5f,  0.5f },
        };

        const std::vector<std::uint32_t> BoundaryIndices = {
            0, 1, 1, 2, 2, 3, 3, 0,
            4, 5, 5, 6, 6, 7, 7, 4,
            0, 4, 1, 5, 2, 6, 3, 7,
        };
    }

    CaseFreeSurfaceSeparation::CaseFreeSurfaceSeparation(
        std::initializer_list<Assets::ExampleScene> && scenes):
        _scenes(scenes),
        _program(Engine::GL::UniqueProgram({
            Engine::GL::SharedShader("assets/shaders/fluid.vert"),
            Engine::GL::SharedShader("assets/shaders/fluid.frag"),
        })),
        _lineProgram(Engine::GL::UniqueProgram({
            Engine::GL::SharedShader("assets/shaders/flat.vert"),
            Engine::GL::SharedShader("assets/shaders/flat.frag"),
        })),
        _sceneObject(SeparationPassConstantsBinding),
        _boundaryItem(
            Engine::GL::VertexLayout().Add<glm::vec3>(
                "position", Engine::GL::DrawFrequency::Stream, 0),
            Engine::GL::PrimitiveType::Lines) {
        _cameraManager.AutoRotate = false;
        _program.BindUniformBlock("PassConstants", SeparationPassConstantsBinding);

        _boundaryItem.UpdateElementBuffer(BoundaryIndices);
        _boundaryItem.UpdateVertexBuffer(
            "position", Engine::make_span_bytes<glm::vec3>(BoundaryVertices));
        ResetSystem();
    }

    void CaseFreeSurfaceSeparation::OnSetupPropsUI() {
        if (ImGui::Button("Reset System"))
            ResetSystem();
        ImGui::SameLine();
        if (ImGui::Button(_stopped ? "Start Simulation" : "Stop Simulation"))
            _stopped = ! _stopped;

        ImGui::Spacing();
        ImGui::TextWrapped(
            "Batty et al. Section 4: a water blob hits the left wall. "
            "The active set switches between no-penetration contact and "
            "a zero-pressure separating surface.");
        if (ImGui::Checkbox("Wall Separation", &_simulation.enableWallSeparation))
            ResetSystem();

        ImGui::Text("Candidate cells: %d", _simulation.wallSeparationCandidates);
        ImGui::Text("Separating wall cells: %d", _simulation.separatingWallCells);
        ImGui::Text("Separating wall faces: %d", _simulation.separatingWallFaces);
        ImGui::Text("Wall-contact particles: %d", _simulation.wallContactParticles);
        ImGui::Text("Average wall distance: %.3f", _simulation.averageLeftWallDistance);
        ImGui::Text(
            "Active-set iterations: %d",
            _simulation.wallSeparationActiveSetIterations);
        ImGui::Text("Min contact pressure: %.5f", _simulation.minimumContactPressure);
        ImGui::Text("Pressure residual: %.3e", _simulation.pressureResidual);
        ImGui::Checkbox("Highlight Wall Contact", &_simulation.highlightWallContact);

        ImGui::Spacing();
        ImGui::SliderFloat("Time Step", &_timeStep, 0.001f, 0.016f, "%.3f");
        ImGui::SliderFloat("FLIP Ratio", &_simulation.m_fRatio, 0.0f, 1.0f, "%.2f");
        ImGui::SliderInt("Pressure Iters", &_simulation.numPressureIters, 20, 300);

        const char * colorModes[] = { "Default", "Velocity", "Density", "Pressure" };
        int          colorMode    = static_cast<int>(_simulation.m_colorMode);
        if (ImGui::Combo("Color", &colorMode, colorModes, IM_ARRAYSIZE(colorModes)))
            _simulation.m_colorMode =
                static_cast<Final::Simulator::ColorMode>(colorMode);
    }

    Common::CaseRenderResult CaseFreeSurfaceSeparation::OnRender(
        std::pair<std::uint32_t, std::uint32_t> const desiredSize) {
        if (_recompute) {
            _recompute = false;
            _sceneObject.ReplaceScene(
                Rendering::Content::Scenes[std::size_t(_scenes[_sceneIdx])]);
            _cameraManager.Save(_sceneObject.Camera);
        }

        if (! _stopped)
            _simulation.SimulateTimestep(_timeStep);

        _frame.Resize(desiredSize);
        _cameraManager.Update(_sceneObject.Camera);

        float const     aspect = float(desiredSize.first) / float(desiredSize.second);
        glm::mat4 const projection = _sceneObject.Camera.GetProjectionMatrix(aspect);
        glm::mat4 const view       = _sceneObject.Camera.GetViewMatrix();
        _sceneObject.PassConstantsBlock.Update(
            &Rendering::SceneObject::PassConstants::Projection, projection);
        _sceneObject.PassConstantsBlock.Update(
            &Rendering::SceneObject::PassConstants::View, view);
        _sceneObject.PassConstantsBlock.Update(
            &Rendering::SceneObject::PassConstants::ViewPosition,
            _sceneObject.Camera.Eye);

        if (_uniformDirty) {
            _uniformDirty = false;
            _program.GetUniforms().SetByName("u_AmbientScale", 1.0f);
            _program.GetUniforms().SetByName("u_UseBlinn", 0);
            _program.GetUniforms().SetByName("u_Shininess", 32.0f);
            _program.GetUniforms().SetByName("u_UseGammaCorrection", 1);
            _program.GetUniforms().SetByName("u_AttenuationOrder", 2);
        }

        gl_using(_frame);
        glEnable(GL_DEPTH_TEST);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        _lineProgram.GetUniforms().SetByName("u_Projection", projection);
        _lineProgram.GetUniforms().SetByName("u_View", view);
        _lineProgram.GetUniforms().SetByName("u_Color", glm::vec3(1.0f));
        glLineWidth(_boundaryWidth);
        _boundaryItem.Draw({ _lineProgram.Use() });
        glLineWidth(1.0f);

        Rendering::ModelObject particles(
            _sphere, _simulation.m_particlePos, _simulation.m_particleColor);
        particles.Mesh.Draw(
            { _program.Use() },
            _sphere.Mesh.Indices.size(),
            0,
            _simulation.m_iNumSpheres);

        glDisable(GL_DEPTH_TEST);
        return Common::CaseRenderResult {
            .Fixed     = false,
            .Flipped   = true,
            .Image     = _frame.GetColorAttachment(),
            .ImageSize = desiredSize,
        };
    }

    void CaseFreeSurfaceSeparation::OnProcessInput(ImVec2 const & pos) {
        _cameraManager.ProcessInput(_sceneObject.Camera, pos);
    }

    void CaseFreeSurfaceSeparation::ResetSystem() {
        bool const enableWallSeparation = _simulation.enableWallSeparation;
        _simulation.setupScene(_resolution);
        _simulation.enableWallSeparation = enableWallSeparation;
        _sphere = Engine::Model(
            Engine::Sphere(6, _simulation.m_particleRadius));
    }
} // namespace VCX::Labs::FluidSimulation
