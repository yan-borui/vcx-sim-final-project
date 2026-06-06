#include "Labs/FinalProject/CaseSubgrid.h"

#include "Engine/app.h"

namespace VCX::Labs::FluidSimulation {
    namespace {
        constexpr std::uint32_t SubgridPassConstantsBinding = 2;

        const std::vector<glm::vec3> BoundaryVertices = {
            { -0.5f, -0.5f, -0.5f }, { 0.5f, -0.5f, -0.5f },
            {  0.5f,  0.5f, -0.5f }, {-0.5f,  0.5f, -0.5f },
            { -0.5f, -0.5f,  0.5f }, { 0.5f, -0.5f,  0.5f },
            {  0.5f,  0.5f,  0.5f }, {-0.5f,  0.5f,  0.5f },
        };

        const std::vector<std::uint32_t> BoundaryIndices = {
            0, 1, 1, 2, 2, 3, 3, 0,
            4, 5, 5, 6, 6, 7, 7, 4,
            0, 4, 1, 5, 2, 6, 3, 7,
        };

        const std::vector<glm::vec3> BoxVertices = {
            { -0.5f, -0.5f,  0.5f }, { 0.5f, -0.5f,  0.5f },
            {  0.5f,  0.5f,  0.5f }, {-0.5f,  0.5f,  0.5f },
            { -0.5f, -0.5f, -0.5f }, { 0.5f, -0.5f, -0.5f },
            {  0.5f,  0.5f, -0.5f }, {-0.5f,  0.5f, -0.5f },
        };

        const std::vector<std::uint32_t> BoxIndices = {
            0, 1, 2, 0, 2, 3, 1, 5, 6, 1, 6, 2,
            5, 4, 7, 5, 7, 6, 4, 0, 3, 4, 3, 7,
            3, 2, 6, 3, 6, 7, 4, 5, 1, 4, 1, 0,
        };
    }

    CaseSubgrid::CaseSubgrid(std::initializer_list<Assets::ExampleScene> && scenes) :
        _scenes(scenes),
        _program(Engine::GL::UniqueProgram({
            Engine::GL::SharedShader("assets/shaders/fluid.vert"),
            Engine::GL::SharedShader("assets/shaders/fluid.frag"),
        })),
        _lineProgram(Engine::GL::UniqueProgram({
            Engine::GL::SharedShader("assets/shaders/flat.vert"),
            Engine::GL::SharedShader("assets/shaders/flat.frag"),
        })),
        _flatProgram(Engine::GL::UniqueProgram({
            Engine::GL::SharedShader("assets/shaders/flat.vert"),
            Engine::GL::SharedShader("assets/shaders/flat.frag"),
        })),
        _sceneObject(SubgridPassConstantsBinding),
        _boundaryItem(
            Engine::GL::VertexLayout().Add<glm::vec3>("position", Engine::GL::DrawFrequency::Stream, 0),
            Engine::GL::PrimitiveType::Lines),
        _rigidBodyItem(
            Engine::GL::VertexLayout().Add<glm::vec3>("position", Engine::GL::DrawFrequency::Stream, 0),
            Engine::GL::PrimitiveType::Triangles) {
        _cameraManager.AutoRotate = false;
        _program.BindUniformBlock("PassConstants", SubgridPassConstantsBinding);

        _boundaryItem.UpdateElementBuffer(BoundaryIndices);
        _boundaryItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(BoundaryVertices));
        _rigidBodyItem.UpdateElementBuffer(BoxIndices);
        _rigidBodyItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(BoxVertices));

        ResetSystem();
    }

    void CaseSubgrid::OnSetupPropsUI() {
        if (ImGui::Button("Reset System"))
            ResetSystem();
        ImGui::SameLine();
        if (ImGui::Button(_stopped ? "Start Simulation" : "Stop Simulation"))
            _stopped = ! _stopped;

        ImGui::Spacing();
        ImGui::TextWrapped("Volume-weighted MAC pressure projection from Batty et al.");
        if (ImGui::Checkbox("Sub-grid weights", &_simulation.useSubgridWeights))
            ResetSystem();

        ImGui::Text("Cut faces: %d", _simulation.partiallyOpenFaceCount);
        ImGui::Text("Min fraction: %.3f", _simulation.minimumOpenFaceFraction);

        ImGui::Spacing();
        ImGui::SliderFloat("Time Step", &_timeStep, 0.001f, 0.03f, "%.3f");
        ImGui::SliderFloat("FLIP Ratio", &_simulation.m_fRatio, 0.0f, 1.0f, "%.2f");
        ImGui::SliderInt("Pressure Iters", &_simulation.numPressureIters, 20, 300);

        const char * colorModes[] = { "Default", "Velocity", "Density", "Pressure" };
        int colorMode = static_cast<int>(_simulation.m_colorMode);
        if (ImGui::Combo("Color", &colorMode, colorModes, IM_ARRAYSIZE(colorModes)))
            _simulation.m_colorMode = static_cast<Final::Simulator::ColorMode>(colorMode);
    }

    Common::CaseRenderResult CaseSubgrid::OnRender(
        std::pair<std::uint32_t, std::uint32_t> const desiredSize) {
        if (_recompute) {
            _recompute = false;
            _sceneObject.ReplaceScene(Rendering::Content::Scenes[std::size_t(_scenes[_sceneIdx])]);
            _cameraManager.Save(_sceneObject.Camera);
        }

        if (! _stopped)
            _simulation.SimulateTimestep(_timeStep);

        _frame.Resize(desiredSize);
        _cameraManager.Update(_sceneObject.Camera);

        float const aspect = float(desiredSize.first) / float(desiredSize.second);
        glm::mat4 const projection = _sceneObject.Camera.GetProjectionMatrix(aspect);
        glm::mat4 const view = _sceneObject.Camera.GetViewMatrix();
        _sceneObject.PassConstantsBlock.Update(&Rendering::SceneObject::PassConstants::Projection, projection);
        _sceneObject.PassConstantsBlock.Update(&Rendering::SceneObject::PassConstants::View, view);
        _sceneObject.PassConstantsBlock.Update(&Rendering::SceneObject::PassConstants::ViewPosition, _sceneObject.Camera.Eye);

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

        Rendering::ModelObject particles(_sphere, _simulation.m_particlePos, _simulation.m_particleColor);
        particles.Mesh.Draw(
            { _program.Use() },
            _sphere.Mesh.Indices.size(),
            0,
            _simulation.m_iNumSpheres);

        glm::mat4 const bodyModel =
            glm::translate(glm::mat4(1.0f), _body.position)
            * glm::mat4_cast(_body.orientation)
            * glm::scale(glm::mat4(1.0f), _body.dim);
        _flatProgram.GetUniforms().SetByName("u_Projection", projection);
        _flatProgram.GetUniforms().SetByName("u_View", view * bodyModel);
        _flatProgram.GetUniforms().SetByName("u_Color", _body.color);
        _rigidBodyItem.Draw({ _flatProgram.Use() });

        glDisable(GL_DEPTH_TEST);

        return Common::CaseRenderResult {
            .Fixed = false,
            .Flipped = true,
            .Image = _frame.GetColorAttachment(),
            .ImageSize = desiredSize,
        };
    }

    void CaseSubgrid::OnProcessInput(ImVec2 const & pos) {
        _cameraManager.ProcessInput(_sceneObject.Camera, pos);
    }

    void CaseSubgrid::ResetSystem() {
        float const h = 1.0f / float(_resolution);
        _body.Reset(
            { 0.03f, 0.0f, 0.0f },
            { 0.0f, 0.0f, 0.0f },
            { 0.18f, 1.0f - 3.0f * h, 1.0f - 2.0f * h },
            1.0f,
            { 0.85f, 0.55f, 0.15f });
        _body.isStatic = true;
        _body.ComputeInertia();

        _simulation.m_body = &_body;
        _simulation.setupScene(_resolution);
        _sphere = Engine::Model(Engine::Sphere(6, _simulation.m_particleRadius));
    }
}
