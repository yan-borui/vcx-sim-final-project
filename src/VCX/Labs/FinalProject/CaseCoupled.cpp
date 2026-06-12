#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <spdlog/spdlog.h>

#include "Engine/app.h"
#include "Labs/Common/ImGuiHelper.h"
#include "Labs/FinalProject/CaseCoupled.h"
#include "Labs/FinalProject/RenderBindings.h"

namespace VCX::Labs::FluidSimulation {
    const std::vector<glm::vec3> vertex_pos = {
        glm::vec3(-0.5f, -0.5f, -0.5f),
        glm::vec3(0.5f, -0.5f, -0.5f),
        glm::vec3(0.5f, 0.5f, -0.5f),
        glm::vec3(-0.5f, 0.5f, -0.5f),
        glm::vec3(-0.5f, -0.5f, 0.5f),
        glm::vec3(0.5f, -0.5f, 0.5f),
        glm::vec3(0.5f, 0.5f, 0.5f),
        glm::vec3(-0.5f, 0.5f, 0.5f)
    };

    const std::vector<std::uint32_t> line_index = {
        0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4, 0, 4, 1, 5, 2, 6, 3, 7
    };

    CaseFluid::CaseFluid(std::initializer_list<Assets::ExampleScene> && scenes):
        _scenes(scenes),
        _program(
            Engine::GL::UniqueProgram({ Engine::GL::SharedShader("assets/shaders/fluid.vert"),
                                        Engine::GL::SharedShader("assets/shaders/fluid.frag") })),
        _lineprogram(
            Engine::GL::UniqueProgram({ Engine::GL::SharedShader("assets/shaders/flat.vert"),
                                        Engine::GL::SharedShader("assets/shaders/flat.frag") })),
        _flatProgram(
            Engine::GL::UniqueProgram({ Engine::GL::SharedShader("assets/shaders/flat.vert"),
                                        Engine::GL::SharedShader("assets/shaders/flat.frag") })),
        _RigidBodyItem(
            Engine::GL::VertexLayout().Add<glm::vec3>(
                "position",
                Engine::GL::DrawFrequency::Stream,
                0),
            Engine::GL::PrimitiveType::Triangles),
        _MeshBodyItem(
            Engine::GL::VertexLayout().Add<glm::vec3>(
                "position",
                Engine::GL::DrawFrequency::Static,
                0),
            Engine::GL::PrimitiveType::Triangles),
        _sceneObject(CoupledPassConstantsBinding),
        _BoundaryItem(
            Engine::GL::VertexLayout().Add<glm::vec3>(
                "position",
                Engine::GL::DrawFrequency::Stream,
                0),
            Engine::GL::PrimitiveType::Lines) {
        _cameraManager.AutoRotate = false;

        _program.BindUniformBlock(
            "PassConstants",
            CoupledPassConstantsBinding);
        _lineprogram.BindUniformBlock(
            "PassConstants",
            CoupledPassConstantsBinding);

        _program.GetUniforms().SetByName("u_DiffuseMap", 0);
        _program.GetUniforms().SetByName("u_SpecularMap", 1);
        _program.GetUniforms().SetByName("u_HeightMap", 2);

        _BoundaryItem.UpdateElementBuffer(line_index);
        _BoundaryItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(vertex_pos));

        const std::vector<glm::vec3> box_v = {
            { -.5, -.5,  .5 },
            {  .5, -.5,  .5 },
            {  .5,  .5,  .5 },
            { -.5,  .5,  .5 },
            { -.5, -.5, -.5 },
            {  .5, -.5, -.5 },
            {  .5,  .5, -.5 },
            { -.5,  .5, -.5 }
        };
        const std::vector<std::uint32_t> box_idx = {
            0, 1, 2, 0, 2, 3, 1, 5, 6, 1, 6, 2, 5, 4, 7, 5, 7, 6, 4, 0, 3, 4, 3, 7, 3, 2, 6, 3, 6, 7, 4, 5, 1, 4, 1, 0
        };

        _RigidBodyItem.UpdateElementBuffer(box_idx);
        _RigidBodyItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(box_v));

        ResetSystem();
    }

    void CaseFluid::LoadMeshIfNeeded() {
        if (_meshReady)
            return;

        auto mesh = Engine::LoadSurfaceMesh("assets/models/suzanne.obj", true);
        mesh.NormalizePositions();
        _MeshBodyItem.UpdateElementBuffer(mesh.Indices);
        _MeshBodyItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(mesh.Positions));

        _meshSDF = std::make_shared<Final::MeshSDF>();
        _meshReady = _meshSDF->LoadOBJ("assets/models/suzanne.obj");
        if (! _meshReady)
            spdlog::warn("Failed to load assets/models/suzanne.obj for irregular solid SDF.");
    }

    void CaseFluid::OnSetupPropsUI() {
        if (ImGui::Button("Reset System"))
            ResetSystem();
        ImGui::SameLine();
        if (ImGui::Button(_stopped ? "Start Simulation" : "Stop Simulation"))
            _stopped = ! _stopped;

        ImGui::Spacing();

        bool solverChanged = false;
        if (ImGui::Checkbox("Use APIC", &_useAPIC)) {
            _useCG        = false;
            solverChanged = true;
        }
        ImGui::SameLine();
        if (ImGui::Checkbox("Use CG", &_useCG)) {
            _useAPIC      = false;
            solverChanged = true;
        }
        if (solverChanged)
            ResetSystem();

        int          shapeMode    = static_cast<int>(_demoShape);
        const char * shapeModes[] = { "Box", "Sphere", "suzanne" };
        if (ImGui::Combo("Rigid Shape", &shapeMode, shapeModes, IM_ARRAYSIZE(shapeModes))) {
            _demoShape = static_cast<DemoShape>(shapeMode);
            ResetSystem();
        }

        if (! _useAPIC && ! _useCG) {
            ImGui::Separator();
            ImGui::Text("Variational boundary projection");
            ImGui::Checkbox("Sub-grid weights", &_simulation.useSubgridWeights);
            ImGui::Checkbox("Wall separation", &_simulation.enableWallSeparation);
            ImGui::Text("Pressure residual: %.3e", _simulation.pressureResidual);
            ImGui::Text(
                "Wall KKT: %.3e (%d QP steps)",
                _simulation.wallSeparationKktResidual,
                _simulation.wallSeparationIterations);
        } else if (_useAPIC) {
            ImGui::Text("APIC transfer with current rigid boundary");
        } else if (_useCG) {
            ImGui::Text("CG pressure projection with current rigid boundary");
        }

        if (_demoShape == DemoShape::Mesh) {
            ImGui::Text("Irregular solid: suzanne.obj %s", _meshReady ? "loaded" : "pending");
        }

        if (ImGui::SliderInt("Grid Resolution", &_res, 10, 24))
            ResetSystem();

        ImGui::SliderFloat("Time Step", &_dt, 0.001f, 0.03f);

        ImGui::Separator();
        ImGui::Text("Visualization Enhancement");
        const char * colorModes[]     = { "Default", "Velocity", "Density", "Pressure" };
        int          currentColorMode = static_cast<int>(_sim->m_colorMode);
        if (ImGui::Combo("Color Mode", &currentColorMode, colorModes, IM_ARRAYSIZE(colorModes))) {
            _sim->m_colorMode = static_cast<Final::Simulator::ColorMode>(currentColorMode);
        }

        ImGui::Separator();
        ImGui::Checkbox("separateParticles", &_sim->separateParticles);
        ImGui::Checkbox("compensateDrift", &_sim->compensateDrift);
        ImGui::SliderFloat("Flip Ratio", &_sim->m_fRatio, 0.0f, 1.0f);
        ImGui::SliderFloat("compensateDriftWeight", &_sim->compensateDriftWeight, 0.0f, 1.0f);
        ImGui::SliderFloat("overRelaxation", &_sim->overRelaxation, 0.3f, 2.0f);
        ImGui::SliderInt("numPressureIters", &_sim->numPressureIters, 5, 1000);
        ImGui::SliderInt("numParticleIters", &_sim->numParticleIters, 3, 10);
    }

    Common::CaseRenderResult CaseFluid::OnRender(std::pair<std::uint32_t, std::uint32_t> const desiredSize) {
        if (_recompute) {
            _recompute = false;
            _sceneObject.ReplaceScene(GetScene(_sceneIdx));
            _cameraManager.Save(_sceneObject.Camera);
        }

        if (! _stopped) {
            _body.position += _body.velocity * _dt;
            if (glm::length(_body.angularVelocity) > 0.001f) {
                glm::vec3 const axis = glm::normalize(_body.angularVelocity);
                float const angle = glm::length(_body.angularVelocity) * _dt;
                _body.orientation =
                    glm::normalize(glm::angleAxis(angle, axis) * _body.orientation);
            }

            float const boundRadius = _body.BoundingRadius();
            for (int direction = 0; direction < 3; ++direction) {
                if (_body.position[direction] - boundRadius < -0.5f) {
                    _body.position[direction] = -0.5f + boundRadius;
                    _body.velocity[direction] *= -0.5f;
                }
                if (_body.position[direction] + boundRadius > 0.5f) {
                    _body.position[direction] = 0.5f - boundRadius;
                    _body.velocity[direction] *= -0.5f;
                }
            }

            if (_body.position.y < 0.1f) {
                _body.velocity *= 0.96f;
                _body.angularVelocity *= 0.95f;
            }

            _body.velocity += glm::vec3(0.0f, -9.81f, 0.0f) * _dt;

            _sim->m_boxPos = _body.position;
            _sim->m_boxDim = _body.dim;
            _sim->SimulateTimestep(_dt);

            _body.velocity +=
                _sim->m_feedbackForce / std::max(_body.mass, 1e-6f) * _dt;
            _body.angularVelocity +=
                _body.GetInertiaWorldInv() * _sim->m_feedbackTorque * _dt;

            if (_useCG && _body.shape == Final::RigidBody::ShapeType::Mesh) {
                if (! std::isfinite(_body.velocity.x)
                    || ! std::isfinite(_body.velocity.y)
                    || ! std::isfinite(_body.velocity.z))
                    _body.velocity = glm::vec3(0.0f);
                if (! std::isfinite(_body.angularVelocity.x)
                    || ! std::isfinite(_body.angularVelocity.y)
                    || ! std::isfinite(_body.angularVelocity.z))
                    _body.angularVelocity = glm::vec3(0.0f);

                float const maxSpeed = 4.0f;
                float const speed    = glm::length(_body.velocity);
                if (speed > maxSpeed)
                    _body.velocity *= maxSpeed / speed;

                float const maxAngularSpeed = 20.0f;
                float const angularSpeed    = glm::length(_body.angularVelocity);
                if (angularSpeed > maxAngularSpeed)
                    _body.angularVelocity *= maxAngularSpeed / angularSpeed;
            }
        }

        _frame.Resize(desiredSize);
        _cameraManager.Update(_sceneObject.Camera);

        auto proj = _sceneObject.Camera.GetProjectionMatrix((float(desiredSize.first) / desiredSize.second));
        auto view = _sceneObject.Camera.GetViewMatrix();
        _sceneObject.PassConstantsBlock.Update(&VCX::Labs::Rendering::SceneObject::PassConstants::Projection, proj);
        _sceneObject.PassConstantsBlock.Update(&VCX::Labs::Rendering::SceneObject::PassConstants::View, view);
        _sceneObject.PassConstantsBlock.Update(&VCX::Labs::Rendering::SceneObject::PassConstants::ViewPosition, _sceneObject.Camera.Eye);

        gl_using(_frame);
        glEnable(GL_DEPTH_TEST);
        glClearColor(0.1f, 0.1f, 0.1f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        _lineprogram.GetUniforms().SetByName("u_Projection", proj);
        _lineprogram.GetUniforms().SetByName("u_View", view);
        _lineprogram.GetUniforms().SetByName("u_Model", glm::mat4(1.0f));
        _lineprogram.GetUniforms().SetByName("u_Color", glm::vec3(1.0f, 1.0f, 1.0f));

        _BoundaryItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(vertex_pos));
        glLineWidth(_BndWidth);
        _BoundaryItem.Draw({ _lineprogram.Use() });
        glLineWidth(1.f);

        if (_uniformDirty) {
            _uniformDirty = false;
            _program.GetUniforms().SetByName("u_AmbientScale", _ambientScale);
            _program.GetUniforms().SetByName("u_UseBlinn", _useBlinn);
            _program.GetUniforms().SetByName("u_Shininess", _shininess);
            _program.GetUniforms().SetByName("u_UseGammaCorrection", int(_useGammaCorrection));
            _program.GetUniforms().SetByName("u_AttenuationOrder", _attenuationOrder);
            _program.GetUniforms().SetByName("u_BumpMappingBlend", _bumpMappingPercent * .01f);
        }

        Rendering::ModelObject fluidParticles = Rendering::ModelObject(_sphere, _sim->m_particlePos, _sim->m_particleColor);
        auto const &           material       = _sceneObject.Materials[0];
        _program.GetUniforms().SetByName("u_FluidProjection", proj);
        _program.GetUniforms().SetByName("u_FluidView", view);
        fluidParticles.Mesh.Draw(
            { material.Albedo.Use(), material.MetaSpec.Use(), material.Height.Use(), _program.Use() },
            _sphere.Mesh.Indices.size(),
            0,
            _sim->m_iNumSpheres);

        glm::mat4 bodyModel = glm::translate(glm::mat4(1.0f), _body.position)
            * glm::mat4_cast(_body.orientation)
            * glm::scale(glm::mat4(1.0f), _body.dim);

        _flatProgram.GetUniforms().SetByName("u_Projection", proj);
        _flatProgram.GetUniforms().SetByName("u_Color", _body.color);
        _flatProgram.GetUniforms().SetByName("u_View", view);
        _flatProgram.GetUniforms().SetByName("u_Model", bodyModel);

        if (_body.shape == Final::RigidBody::ShapeType::Sphere) {
            Rendering::ModelObject rigidSphere = Rendering::ModelObject(_rigidSphere, std::vector<glm::vec3> { glm::vec3(0.0f) });
            rigidSphere.Mesh.Draw(
                { _flatProgram.Use() },
                _rigidSphere.Mesh.Indices.size(),
                0,
                1);
        } else if (_body.shape == Final::RigidBody::ShapeType::Mesh && _meshReady) {
            _MeshBodyItem.Draw({ _flatProgram.Use() });
        } else {
            _RigidBodyItem.Draw({ _flatProgram.Use() });
        }

        glDepthFunc(GL_LEQUAL);
        glDepthFunc(GL_LESS);
        glDisable(GL_DEPTH_TEST);

        return Common::CaseRenderResult {
            .Fixed     = false,
            .Flipped   = true,
            .Image     = _frame.GetColorAttachment(),
            .ImageSize = desiredSize,
        };
    }

    void CaseFluid::OnProcessInput(ImVec2 const & pos) {
        auto & io  = ImGui::GetIO();
        auto & cam = _sceneObject.Camera;

        if (! io.KeyCtrl) {
            _cameraManager.ProcessInput(cam, pos);
            return;
        }

        if (io.KeyCtrl && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
            ImVec2 delta = io.MouseDelta;

            if (std::abs(delta.x) > 0.1f || std::abs(delta.y) > 0.1f) {
                glm::vec3 fwd   = glm::normalize(cam.Target - cam.Eye);
                glm::vec3 right = glm::normalize(glm::cross(fwd, cam.Up));
                glm::vec3 up    = cam.Up;

                float     strength = 0.002f;
                glm::vec3 move     = (right * delta.x - up * delta.y) * strength;

                _body.velocity += move / std::max(Engine::GetDeltaTime(), 0.001f);
            } else {
                _body.velocity = glm::vec3(0.0f);
            }
        }
    }

    void CaseFluid::ResetSystem() {
        if (_useCG) {
            _sim = (Final::Simulator *) &_cgsimulation;
        } else if (_useAPIC) {
            _sim = (Final::Simulator *) &_apicSimulation;
        } else {
            _sim = (Final::Simulator *) &_simulation;
        }

        _sim->setupScene(_res);

        if (_demoShape == DemoShape::Mesh) {
            LoadMeshIfNeeded();
            _simulation.volumeSamplesPerAxis = 2;
            _body.Reset(
                { 0.0f, 0.32f, 0.0f },
                { 0.0f, 0.0f, 0.0f },
                { 0.36f, 0.36f, 0.36f },
                0.35f,
                { 0.90f, 0.82f, 0.62f },
                Final::RigidBody::ShapeType::Mesh);
            if (_meshReady)
                _body.SetMeshSDF(_meshSDF);
        } else if (_demoShape == DemoShape::Sphere) {
            _simulation.volumeSamplesPerAxis = 4;
            _body.Reset(
                { 0.0f, 0.26f, 0.0f },
                { 0.0f, 0.0f, 0.0f },
                { 0.26f, 0.26f, 0.26f },
                0.10f,
                { 0.92f, 0.58f, 0.20f },
                Final::RigidBody::ShapeType::Sphere);
        } else {
            _simulation.volumeSamplesPerAxis = 4;
            _body.Reset(
                { 0.0f, 0.30f, 0.0f },
                { 0.0f, 0.0f, 0.0f },
                { 0.30f, 0.30f, 0.30f },
                0.30f,
                { 0.8f, 0.2f, 0.2f },
                Final::RigidBody::ShapeType::Box);
        }

        _sim->m_body  = &_body;
        numofSpheres  = _sim->m_iNumSpheres;
        _r            = _sim->m_particleRadius;
        _sphere       = Engine::Model(Engine::Sphere(6, _r));
        _rigidSphere  = Engine::Model(Engine::Sphere(20, 0.5f));
        _uniformDirty = true;
    }
} // namespace VCX::Labs::FluidSimulation
