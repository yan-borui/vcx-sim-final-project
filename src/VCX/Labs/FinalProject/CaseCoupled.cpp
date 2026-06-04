#include <spdlog/spdlog.h>
#include "Engine/app.h"
#include "Labs/FinalProject/CaseCoupled.h"
#include "Labs/Common/ImGuiHelper.h"
#include <iostream>

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
    const std::vector<std::uint32_t> line_index = { 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4, 0, 4, 1, 5, 2, 6, 3, 7 }; // line index

    CaseFluid::CaseFluid(std::initializer_list<Assets::ExampleScene> && scenes) :
        _scenes(scenes),
        _program(
            Engine::GL::UniqueProgram({
                Engine::GL::SharedShader("assets/shaders/fluid.vert"),
                Engine::GL::SharedShader("assets/shaders/fluid.frag") })),
        _lineprogram(
            Engine::GL::UniqueProgram({
                Engine::GL::SharedShader("assets/shaders/flat.vert"),
                Engine::GL::SharedShader("assets/shaders/flat.frag") })),
        _flatProgram(
            Engine::GL::UniqueProgram({ 
                Engine::GL::SharedShader("assets/shaders/flat.vert"), 
                Engine::GL::SharedShader("assets/shaders/flat.frag") })),
        _RigidBodyItem(
            Engine::GL::VertexLayout().Add<glm::vec3>("position", 
                Engine::GL::DrawFrequency::Stream , 0), 
                Engine::GL::PrimitiveType::Triangles),
        _sceneObject(1),
        _BoundaryItem(Engine::GL::VertexLayout()
            .Add<glm::vec3>("position", Engine::GL::DrawFrequency::Stream , 0), Engine::GL::PrimitiveType::Lines){ 
        _cameraManager.AutoRotate = false;

        _program.BindUniformBlock("PassConstants", 1);
        _lineprogram.BindUniformBlock("PassConstants", 1);

        _program.GetUniforms().SetByName("u_DiffuseMap" , 0);
        _program.GetUniforms().SetByName("u_SpecularMap", 1);
        _program.GetUniforms().SetByName("u_HeightMap"  , 2);

        _BoundaryItem.UpdateElementBuffer(line_index);
        _BoundaryItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(vertex_pos));
        
        const std::vector<glm::vec3> box_v = {
            {-.5,-.5,.5},{.5,-.5,.5},{.5,.5,.5},{-.5,.5,.5},
            {-.5,-.5,-.5},{.5,-.5,-.5},{.5,0.5,-.5},{-.5,0.5,-.5}
        };
        const std::vector<std::uint32_t> box_idx = { 
            0,1,2, 0,2,3, 1,5,6, 1,6,2, 5,4,7, 5,7,6, 4,0,3, 4,3,7, 3,2,6, 3,6,7, 4,5,1, 4,1,0 
        };

        _RigidBodyItem.UpdateElementBuffer(box_idx);
        _RigidBodyItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(box_v));

        ResetSystem();
    }

    void CaseFluid::OnSetupPropsUI() {
        if(ImGui::Button("Reset System")) 
            ResetSystem();
        ImGui::SameLine();
        if(ImGui::Button(_stopped ? "Start Simulation":"Stop Simulation"))
            _stopped = ! _stopped;
        ImGui::Spacing();
        if (ImGui::Checkbox("Use APIC", &_useAPIC)) {
            _useCG = false;
            ResetSystem();
        }
        ImGui::SameLine();
        if (ImGui::Checkbox("Use CG", &_useCG)) {
            _useAPIC = false;
            ResetSystem();
        }
        ImGui::SliderFloat("Time Step", &_dt, 0.001f, 0.03f);
        ImGui::Separator();
        ImGui::Text("Visualization Enhancement");
        const char* colorModes[] = { "Default", "Velocity", "Density", "Pressure" };
        int currentColorMode = static_cast<int>(_sim->m_colorMode);
        if (ImGui::Combo("Color Mode", &currentColorMode, colorModes, IM_ARRAYSIZE(colorModes))) {
            _sim->m_colorMode = static_cast<Final::Simulator::ColorMode>(currentColorMode);
        }
        ImGui::Separator();
        ImGui::Checkbox("separateParticles", &_sim->separateParticles);
        ImGui::Checkbox("compensateDrift", &_sim->compensateDrift);
        ImGui::SliderFloat("Flip Ratio", &_sim->m_fRatio, 0.0f, 1.0f);
        ImGui::SliderFloat("compensateDriftWeight", &_sim->compensateDriftWeight, 0.0f, 1.0f);
        ImGui::SliderFloat("overRelaxation", &_sim->overRelaxation, 0.3f, 2.0f);
        ImGui::SliderInt("numPressureIters", &_sim->numPressureIters, 5,1000);
        ImGui::SliderInt("numParticleIters", &_sim->numParticleIters, 3,10);
    }


    Common::CaseRenderResult CaseFluid::OnRender(std::pair<std::uint32_t, std::uint32_t> const desiredSize) {
        if (_recompute) {
            _recompute = false;
            _sceneObject.ReplaceScene(GetScene(_sceneIdx));
            _cameraManager.Save(_sceneObject.Camera);
        }
        
        if (! _stopped) {
            // A. 先更新流体位置
            _sim->m_boxPos = _body.position;
            _sim->m_boxDim = _body.dim;
            _sim->SimulateTimestep(_dt);

            // B. 计算受力
            glm::vec3 gravityForce = glm::vec3(0.0f, -9.81f, 0.0f) * _body.mass;
            glm::vec3 totalForce = gravityForce + _sim->m_feedbackForce; // 来自流体的反馈力

            // C. 更新速度 (v = v + a*dt)
            glm::vec3 acceleration = totalForce / _body.mass;
            _body.velocity += acceleration * _dt;

            // D. 应用阻尼 (水中阻力)
            if (_body.position.y < 0.1f) {
                _body.velocity *= 0.96f; // 水中运动变慢
                _body.angularVelocity *= 0.95f;
            }

            // E. 更新位置 (p = p + v*dt)
            _body.position += _body.velocity * _dt;

            // F. 刚体旋转
            glm::vec3 angularAccel = _body.GetInertiaWorldInv() * _sim->m_feedbackTorque;
            _body.angularVelocity += angularAccel * _dt;
            if (glm::length(_body.angularVelocity) > 0.001f) {
                glm::vec3 axis = glm::normalize(_body.angularVelocity);
                float angle = glm::length(_body.angularVelocity) * _dt;
                _body.orientation = glm::normalize(glm::angleAxis(angle, axis) * _body.orientation);
            }

            // H. 边界碰撞检测
            glm::vec3 hD = _body.dim * 0.5f;
            float b = 0.5f;
            if (_body.position.x - hD.x < -b) { _body.position.x = -b + hD.x; _body.velocity.x *= -0.5f; } 
            if (_body.position.x + hD.x >  b) { _body.position.x =  b - hD.x; _body.velocity.x *= -0.5f; }
            if (_body.position.y - hD.y < -b) { _body.position.y = -b + hD.y; _body.velocity.y *= -0.5f; }
            if (_body.position.y + hD.y >  b) { _body.position.y =  b - hD.y; _body.velocity.y *= -0.5f; }
            if (_body.position.z - hD.z < -b) { _body.position.z = -b + hD.z; _body.velocity.z *= -0.5f; }
            if (_body.position.z + hD.z >  b) { _body.position.z =  b - hD.z; _body.velocity.z *= -0.5f; }
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
        
        _lineprogram.GetUniforms().SetByName("u_Projection", _sceneObject.Camera.GetProjectionMatrix((float(desiredSize.first) / desiredSize.second)));
        _lineprogram.GetUniforms().SetByName("u_View"      , _sceneObject.Camera.GetViewMatrix());
        _lineprogram.GetUniforms().SetByName("u_Model"     , glm::mat4(1.0f));
        _lineprogram.GetUniforms().SetByName("u_Color"     , glm::vec3(1.0f, 1.0f, 1.0f));

        _BoundaryItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(vertex_pos));
            glLineWidth(_BndWidth);
            _BoundaryItem.Draw({ _lineprogram.Use() });
            glLineWidth(1.f);
        
        if (_uniformDirty) {
            _uniformDirty = false;
            _program.GetUniforms().SetByName("u_AmbientScale"      , _ambientScale);
            _program.GetUniforms().SetByName("u_UseBlinn"          , _useBlinn);
            _program.GetUniforms().SetByName("u_Shininess"         , _shininess);
            _program.GetUniforms().SetByName("u_UseGammaCorrection", int(_useGammaCorrection));
            _program.GetUniforms().SetByName("u_AttenuationOrder"  , _attenuationOrder);            
            _program.GetUniforms().SetByName("u_BumpMappingBlend"  , _bumpMappingPercent * .01f);            

        }
        
        // Rendering::ModelObject m = Rendering::ModelObject(_sphere,_simulation.Positions);
        Rendering::ModelObject m = Rendering::ModelObject(_sphere, _sim->m_particlePos, _sim->m_particleColor);
        auto const & material    = _sceneObject.Materials[0];
        m.Mesh.Draw({ material.Albedo.Use(),  material.MetaSpec.Use(), material.Height.Use(),_program.Use() },
            _sphere.Mesh.Indices.size(), 0, _sim->m_iNumSpheres);
        
        _flatProgram.GetUniforms().SetByName("u_Projection", _sceneObject.Camera.GetProjectionMatrix((float(desiredSize.first) / desiredSize.second)));
        _flatProgram.GetUniforms().SetByName("u_View"      , _sceneObject.Camera.GetViewMatrix());
        
        glm::mat4 bodyModel = glm::translate(glm::mat4(1.0f), _body.position) * 
                              glm::mat4_cast(_body.orientation) * 
                              glm::scale(glm::mat4(1.0f), _body.dim);

        glm::mat4 View = view * bodyModel;

        _flatProgram.GetUniforms().SetByName("u_Projection", proj);
        _flatProgram.GetUniforms().SetByName("u_View", View);
        _flatProgram.GetUniforms().SetByName("u_Color", _body.color);
        _RigidBodyItem.Draw({ _flatProgram.Use() });

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

    void CaseFluid::OnProcessInput(ImVec2 const& pos) {
        // _cameraManager.ProcessInput(_sceneObject.Camera, pos);
        auto & io = ImGui::GetIO();
        auto & cam = _sceneObject.Camera; 

        if (!io.KeyCtrl) {
            _cameraManager.ProcessInput(cam, pos);
            return; 
        }

        // 如果按住 Ctrl 且点击了左键
        if (io.KeyCtrl && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
            ImVec2 delta = io.MouseDelta;

            if (std::abs(delta.x) > 0.1f || std::abs(delta.y) > 0.1f) {
                glm::vec3 fwd   = glm::normalize(cam.Target - cam.Eye);
                glm::vec3 right = glm::normalize(glm::cross(fwd, cam.Up));
                glm::vec3 up    = cam.Up; 
                
                float strength = 0.002f; 
                glm::vec3 move = (right * delta.x - up * delta.y) * strength;

                _body.position += move;
                _body.velocity = move / std::max(Engine::GetDeltaTime(), 0.001f);
                
            } else {
                _body.velocity = glm::vec3(0.0f);
            }
        } 
    }

    void CaseFluid::ResetSystem(){
        if (_useCG) {
            _sim = (Final::Simulator*)&_cgsimulation;
        } else if (_useAPIC) {
            _sim = (Final::Simulator*)&_apicSimulation;
        } else {
            _sim = (Final::Simulator*)&_simulation;
        }
        
        _sim->setupScene(_res);
        _body.Reset({0, 0.3f, 0}, {0,0,0}, {0.3f, 0.3f, 0.3f}, 0.3f, {0.8f, 0.2f, 0.2f});
        _sim->m_body = &_body;
        numofSpheres = _sim->m_iNumSpheres;
        _r = _sim->m_particleRadius; //cell size
        _sphere = Engine::Model(Engine::Sphere(6, _r));
    }
}