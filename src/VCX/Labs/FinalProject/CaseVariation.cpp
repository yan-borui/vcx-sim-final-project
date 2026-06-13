#include <spdlog/spdlog.h>
#include "Engine/app.h"
#include "Labs/FinalProject/CaseVariation.h"
#include "Labs/FinalProject/RenderBindings.h"
#include "Labs/Common/ImGuiHelper.h"
#include <iostream>

namespace VCX::Labs::FluidSimulation {
    static const std::vector<glm::vec3> vertex_pos = {
            glm::vec3(-0.5f, -0.5f, -0.5f),
            glm::vec3(0.5f, -0.5f, -0.5f),  
            glm::vec3(0.5f, 0.5f, -0.5f),  
            glm::vec3(-0.5f, 0.5f, -0.5f), 
            glm::vec3(-0.5f, -0.5f, 0.5f),  
            glm::vec3(0.5f, -0.5f, 0.5f),   
            glm::vec3(0.5f, 0.5f, 0.5f),   
            glm::vec3(-0.5f, 0.5f, 0.5f)
    };
    static const std::vector<std::uint32_t> line_index = { 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4, 0, 4, 1, 5, 2, 6, 3, 7 }; // line index

    CaseVariation::CaseVariation(std::initializer_list<Assets::ExampleScene> && scenes) :
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
        _sceneObject(VariationPassConstantsBinding),
        _BoundaryItem(Engine::GL::VertexLayout()
            .Add<glm::vec3>("position", Engine::GL::DrawFrequency::Stream , 0), Engine::GL::PrimitiveType::Lines){ 
        _cameraManager.AutoRotate = false;

        _program.BindUniformBlock(
            "PassConstants",
            VariationPassConstantsBinding);
        _lineprogram.BindUniformBlock(
            "PassConstants",
            VariationPassConstantsBinding);

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

    void CaseVariation::OnSetupPropsUI() {
        if(ImGui::Button("Reset System")) 
            ResetSystem();
        ImGui::SameLine();
        if(ImGui::Button(_stopped ? "Start Simulation":"Stop Simulation"))
            _stopped = ! _stopped;
        ImGui::Spacing();
        ImGui::SliderFloat("Time Step", &_dt, 0.001f, 0.03f);
        ImGui::Separator();
        ImGui::Text("Visualization Enhancement");
        const char* colorModes[] = { "Default", "Velocity", "Density", "Pressure" };
        int currentColorMode = static_cast<int>(_sim->m_colorMode);
        if (ImGui::Combo("Color Mode", &currentColorMode, colorModes, IM_ARRAYSIZE(colorModes))) {
            _sim->m_colorMode = static_cast<Final::VariationalSimulator::ColorMode>(currentColorMode);
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


    Common::CaseRenderResult CaseVariation::OnRender(std::pair<std::uint32_t, std::uint32_t> const desiredSize) {
        if (_recompute) {
            _recompute = false;
            _sceneObject.ReplaceScene(GetScene(_sceneIdx));
            _cameraManager.Save(_sceneObject.Camera);
        }
        
        if (! _stopped) {
            _sim->m_boxPos = _body.position;
            _sim->m_boxDim = _body.dim;
            _sim->SimulateTimestep(_dt);
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
        _program.GetUniforms().SetByName("u_FluidProjection", proj);
        _program.GetUniforms().SetByName("u_FluidView", view);
        m.Mesh.Draw({ material.Albedo.Use(),  material.MetaSpec.Use(), material.Height.Use(),_program.Use() },
            _sphere.Mesh.Indices.size(), 0, _sim->m_iNumSpheres);
        
        _flatProgram.GetUniforms().SetByName("u_Projection", _sceneObject.Camera.GetProjectionMatrix((float(desiredSize.first) / desiredSize.second)));
        _flatProgram.GetUniforms().SetByName("u_View"      , _sceneObject.Camera.GetViewMatrix());
        
        glm::mat4 bodyModel = glm::translate(glm::mat4(1.0f), _body.position) * 
                              glm::mat4_cast(_body.orientation) * 
                              glm::scale(glm::mat4(1.0f), _body.dim);

        _flatProgram.GetUniforms().SetByName("u_Projection", proj);
        _flatProgram.GetUniforms().SetByName("u_View", view);
        _flatProgram.GetUniforms().SetByName("u_Model", bodyModel);
        _flatProgram.GetUniforms().SetByName("u_Color", _body.color);
        _RigidBodyItem.Draw({ _flatProgram.Use() });

        glDepthFunc(GL_LEQUAL);
        glDepthFunc(GL_LESS);
        glDisable(GL_DEPTH_TEST);

        glBindVertexArray(0);
        glUseProgram(0);

        return Common::CaseRenderResult {
            .Fixed     = false,
            .Flipped   = true,
            .Image     = _frame.GetColorAttachment(),
            .ImageSize = desiredSize,
        };
    }

    void CaseVariation::OnProcessInput(ImVec2 const& pos) {
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

    void CaseVariation::ResetSystem(){
        _sim = (Final::VariationalSimulator*)&_simulation;
        
        
        _sim->setupScene(_res);
        std::fill(_sim->m_vel.begin(), _sim->m_vel.end(), glm::vec3(0.0f));
    std::fill(_sim->m_pre_vel.begin(), _sim->m_pre_vel.end(), glm::vec3(0.0f));
        _body.Reset({0, 0.3f, 0}, {0,0,0}, {0.3f, 0.3f, 0.3f}, 1.0f, {0.8f, 0.2f, 0.2f});
        _sim->m_body = &_body;
        numofSpheres = _sim->m_iNumSpheres;
        _r = _sim->m_particleRadius; //cell size
        _sphere = Engine::Model(Engine::Sphere(6, _r));
    }
}
