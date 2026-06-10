#include "Labs/0-GettingStarted/CaseFluidRigid_zly.h"
#include "Labs/Common/ImGuiHelper.h"
#include <glm/gtc/matrix_transform.hpp>
#include <imgui.h>

namespace VCX::Labs::GettingStarted {

    CaseFluidRigid::CaseFluidRigid():
        _program(Engine::GL::UniqueProgram({ Engine::GL::SharedShader("assets/shaders/flat.vert"),
                                             Engine::GL::SharedShader("assets/shaders/flat.frag") })),
        _triItem(Engine::GL::VertexLayout().Add<glm::vec3>("position", Engine::GL::DrawFrequency::Stream, 0), Engine::GL::PrimitiveType::Triangles),
        _lineItem(Engine::GL::VertexLayout().Add<glm::vec3>("position", Engine::GL::DrawFrequency::Stream, 0), Engine::GL::PrimitiveType::Lines) {
        _cameraManager.AutoRotate = false;
        _cameraManager.Save(_camera);
    }

    void CaseFluidRigid::ResetSimulation() {
        int ny = int(_gridRes * 1.5f);
        _sim.Init(_gridRes, ny, _gridRes, _domainW);
        _sim.dt = 0.004f;
        _sim.AddDefaultBodies();
        _running     = false;
        _initialized = true;
    }

    void CaseFluidRigid::OnSetupPropsUI() {
        if (ImGui::CollapsingHeader("Simulation", ImGuiTreeNodeFlags_DefaultOpen)) {
            if (! _initialized) {
                if (ImGui::Button("Initialize")) ResetSimulation();
            } else {
                if (ImGui::Button(_running ? "Pause" : "Run")) _running = ! _running;
                ImGui::SameLine();
                if (ImGui::Button("Reset")) ResetSimulation();
            }
            ImGui::SliderInt("Grid NX", &_gridRes, 10, 28);
            ImGui::SliderFloat("Domain W", &_domainW, 1, 3);
            ImGui::SliderInt("Sub-steps", &_substeps, 1, 4);
            if (_initialized) {
                ImGui::SliderFloat("Gravity", &_sim.gravity, -20, 0);
                ImGui::Text("Grid: %dx%dx%d  dx=%.3f", _sim.nx, _sim.ny, _sim.nz, _sim.dx);
                ImGui::Text("Particles: %d", (int) _sim.particles.size());
            }
        }
        if (_initialized && ImGui::CollapsingHeader("Rigid Bodies", ImGuiTreeNodeFlags_DefaultOpen)) {
            for (size_t i = 0; i < _sim.bodies.size(); i++) {
                auto & b = _sim.bodies[i];
                ImGui::PushID((int) i);
                if (ImGui::TreeNode(("Body " + std::to_string(i)).c_str())) {
                    ImGui::ColorEdit3("Color", &b.color.x);
                    if (ImGui::SliderFloat("Density", &b.density, 50, 5000)) b.ComputeMass();
                    ImGui::SliderFloat("Radius", &b.radius, .02f, .25f);
                    ImGui::Text("Vel: (%.2f,%.2f,%.2f)", b.velocity.x, b.velocity.y, b.velocity.z);
                    ImGui::TreePop();
                }
                ImGui::PopID();
            }
        }
        ImGui::Spacing();
        ImGui::TextWrapped(
            "3D fluid-rigid coupling (Batty07). Particles rendered as small spheres. "
            "Orange=floats, Green=neutral, Red=sinks, Yellow=very light.");
    }

    // ©¤©¤ build particle octahedra (small spheres) ©¤©¤
    void CaseFluidRigid::BuildFluidSurface() {
        float  r = _sim.dx * 0.18f; // particle visual radius
        size_t n = _sim.particles.size();
        _verts.reserve(n * 6);
        _triIdx.reserve(n * 24);
        for (auto & p : _sim.particles) {
            uint32_t b = (uint32_t) _verts.size();
            _verts.push_back({ p.x, p.y + r, p.z }); // 0 top
            _verts.push_back({ p.x, p.y - r, p.z }); // 1 bottom
            _verts.push_back({ p.x + r, p.y, p.z }); // 2 right
            _verts.push_back({ p.x - r, p.y, p.z }); // 3 left
            _verts.push_back({ p.x, p.y, p.z + r }); // 4 front
            _verts.push_back({ p.x, p.y, p.z - r }); // 5 back
            _triIdx.insert(_triIdx.end(), { b, b + 2, b + 4, b, b + 4, b + 3, b, b + 3, b + 5, b, b + 5, b + 2, b + 1, b + 4, b + 2, b + 1, b + 3, b + 4, b + 1, b + 5, b + 3, b + 1, b + 2, b + 5 });
        }
    }

    void CaseFluidRigid::AddSphere(glm::vec3 c, float r, int seg) {
        uint32_t base = (uint32_t) _verts.size();
        int      sl   = seg * 2;
        _verts.push_back({ c.x, c.y + r, c.z });
        for (int j = 1; j < seg; j++) {
            float phi = glm::pi<float>() * j / seg, sp = std::sin(phi), cp = std::cos(phi);
            for (int i = 0; i < sl; i++) {
                float th = glm::pi<float>() * 2 * i / sl;
                _verts.push_back({ c.x + r * sp * std::cos(th), c.y + r * cp, c.z + r * sp * std::sin(th) });
            }
        }
        _verts.push_back({ c.x, c.y - r, c.z });
        for (int i = 0; i < sl; i++) _triIdx.insert(_triIdx.end(), { base, base + 1 + i, base + 1 + (i + 1) % sl });
        for (int j = 0; j < seg - 2; j++)
            for (int i = 0; i < sl; i++) {
                uint32_t a = base + 1 + j * sl + i, b_ = base + 1 + j * sl + (i + 1) % sl, cc = base + 1 + (j + 1) * sl + i, d = base + 1 + (j + 1) * sl + (i + 1) % sl;
                _triIdx.insert(_triIdx.end(), { a, cc, b_, b_, cc, d });
            }
        uint32_t south = base + 1 + (seg - 1) * sl;
        for (int i = 0; i < sl; i++) _triIdx.insert(_triIdx.end(), { south, base + 1 + (seg - 2) * sl + (i + 1) % sl, base + 1 + (seg - 2) * sl + i });
    }

    void CaseFluidRigid::AddWireBox(glm::vec3 lo, glm::vec3 hi) {
        uint32_t b = (uint32_t) _lineVerts.size();
        _lineVerts.push_back({ lo.x, lo.y, lo.z });
        _lineVerts.push_back({ hi.x, lo.y, lo.z });
        _lineVerts.push_back({ hi.x, hi.y, lo.z });
        _lineVerts.push_back({ lo.x, hi.y, lo.z });
        _lineVerts.push_back({ lo.x, lo.y, hi.z });
        _lineVerts.push_back({ hi.x, lo.y, hi.z });
        _lineVerts.push_back({ hi.x, hi.y, hi.z });
        _lineVerts.push_back({ lo.x, hi.y, hi.z });
        _lineIdx.insert(_lineIdx.end(), { b, b + 1, b + 1, b + 2, b + 2, b + 3, b + 3, b, b + 4, b + 5, b + 5, b + 6, b + 6, b + 7, b + 7, b + 4, b, b + 4, b + 1, b + 5, b + 2, b + 6, b + 3, b + 7 });
    }

    void CaseFluidRigid::AddQuadFace(glm::vec3, float, int, int) {} // unused now

    Common::CaseRenderResult CaseFluidRigid::OnRender(
        std::pair<std::uint32_t, std::uint32_t> const desiredSize) {
        _frameSize = desiredSize;
        if (! _initialized) ResetSimulation();
        if (_running)
            for (int s = 0; s < _substeps; s++) _sim.Step();

        _frame.Resize(desiredSize);
        _cameraManager.Update(_camera);
        _program.GetUniforms().SetByName("u_Projection", _camera.GetProjectionMatrix(float(desiredSize.first) / desiredSize.second));
        _program.GetUniforms().SetByName("u_View", _camera.GetViewMatrix());

        gl_using(_frame);
        glEnable(GL_DEPTH_TEST);

        // ©¤©¤ fluid particles as small octahedra ©¤©¤
        _verts.clear();
        _triIdx.clear();
        BuildFluidSurface();
        if (! _triIdx.empty()) {
            _program.GetUniforms().SetByName("u_Color", glm::vec3(.15f, .45f, .95f));
            _triItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(_verts));
            _triItem.UpdateElementBuffer(_triIdx);
            _triItem.Draw({ _program.Use() });
        }

        // ©¤©¤ rigid bodies ©¤©¤
        for (auto & b : _sim.bodies) {
            _verts.clear();
            _triIdx.clear();
            AddSphere(b.position, b.radius, 8);
            if (! _triIdx.empty()) {
                _program.GetUniforms().SetByName("u_Color", b.color);
                _triItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(_verts));
                _triItem.UpdateElementBuffer(_triIdx);
                _triItem.Draw({ _program.Use() });
            }
        }

        // ©¤©¤ wireframe ©¤©¤
        _lineVerts.clear();
        _lineIdx.clear();
        AddWireBox({ 0, 0, 0 }, { _sim.domainW, _sim.domainH, _sim.domainD });
        _program.GetUniforms().SetByName("u_Color", glm::vec3(.5f, .5f, .6f));
        _lineItem.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec3>(_lineVerts));
        _lineItem.UpdateElementBuffer(_lineIdx);
        glLineWidth(1.5f);
        _lineItem.Draw({ _program.Use() });
        glLineWidth(1.f);

        return Common::CaseRenderResult { .Fixed = false, .Flipped = true, .Image = _frame.GetColorAttachment(), .ImageSize = desiredSize };
    }

    void CaseFluidRigid::OnProcessInput(ImVec2 const & pos) {
        _cameraManager.ProcessInput(_camera, pos);
    }

} // namespace VCX::Labs::GettingStarted
