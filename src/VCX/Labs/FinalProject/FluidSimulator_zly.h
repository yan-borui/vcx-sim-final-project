#pragma once
#include <algorithm>
#include <cmath>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <random>
#include <vector>

namespace VCX::Labs::Final {

    struct RigidBody3D {
        glm::vec3 position { 0 }, velocity { 0 }, dragForce { 0 }, color { 1 };
        float     density = 1000, radius = 0.15f, mass = 0;
        bool      dragged = false;
        void      ComputeMass() { mass = density * 4.f / 3.f * glm::pi<float>() * radius * radius * radius; }
        float     SDF(glm::vec3 p) const { return glm::length(p - position) - radius; }
    };

    class FluidSimulator3D {
    public:
        int                             nx = 20, ny = 30, nz = 20;
        float                           dx, domainW = 2, domainH = 3, domainD = 2;
        float                           fluidRho = 1000, gravity = -9.8f, dt = 0.004f;
        std::vector<float>              u, v, w, pressure;
        std::vector<float>              wFU, wFV, wFW; // fluid fractions
        std::vector<float>              aU, aV, aW;    // pre-computed face coefficients
        std::vector<std::vector<float>> bSU, bSV, bSW; // per-body solid fractions
        enum Cell { FLUID,
                    AIR };
        std::vector<Cell>        cell;
        std::vector<glm::vec3>   particles;
        std::vector<RigidBody3D> bodies;
        int                      seedCounter = 0;

        int       ui(int i, int j, int k) const { return (i * ny + j) * nz + k; }
        int       vi(int i, int j, int k) const { return (i * (ny + 1) + j) * nz + k; }
        int       wi(int i, int j, int k) const { return (i * ny + j) * (nz + 1) + k; }
        int       ci(int i, int j, int k) const { return (i * ny + j) * nz + k; }
        glm::vec3 uPos(int i, int j, int k) const { return { i * dx, (j + .5f) * dx, (k + .5f) * dx }; }
        glm::vec3 vPos(int i, int j, int k) const { return { (i + .5f) * dx, j * dx, (k + .5f) * dx }; }
        glm::vec3 wPos(int i, int j, int k) const { return { (i + .5f) * dx, (j + .5f) * dx, k * dx }; }

        void Init(int _nx, int _ny, int _nz, float _dw) {
            nx      = _nx;
            ny      = _ny;
            nz      = _nz;
            domainW = _dw;
            dx      = domainW / nx;
            domainH = ny * dx;
            domainD = nz * dx;
            int uS = (nx + 1) * ny * nz, vS = nx * (ny + 1) * nz, wS = nx * ny * (nz + 1), cS = nx * ny * nz;
            u.assign(uS, 0);
            v.assign(vS, 0);
            w.assign(wS, 0);
            pressure.assign(cS, 0);
            cell.assign(cS, AIR);
            wFU.assign(uS, 0);
            wFV.assign(vS, 0);
            wFW.assign(wS, 0);
            aU.assign(uS, 0);
            aV.assign(vS, 0);
            aW.assign(wS, 0);
            SeedFluid();
            MarkCells();
            seedCounter = 0;
        }
        void SeedFluid() {
            particles.clear();
            std::mt19937                          rng(42);
            std::uniform_real_distribution<float> jit(.15f, .85f);
            float                                 fluidTop = domainH * .7f;
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++) {
                        if ((j + .5f) * dx < fluidTop) {
                            for (int p = 0; p < 2; p++)
                                particles.push_back({ (i + jit(rng)) * dx, (j + jit(rng)) * dx, (k + jit(rng)) * dx });
                        }
                    }
        }
        void AddDefaultBodies() {
            bodies.clear();
            auto  add = [&](glm::vec3 p, float d, float r, glm::vec3 c) {
            RigidBody3D b;b.position=p;b.density=d;b.radius=r;b.color=c;b.ComputeMass();bodies.push_back(b); };
            float cx = domainW * .5f, cz = domainD * .5f;
            add({ cx * .5f, domainH * .35f, cz }, 300, domainW * .065f, { 1, .45f, .15f });
            add({ cx, domainH * .45f, cz }, 1000, domainW * .055f, { .2f, .9f, .4f });
            add({ cx * 1.5f, domainH * .5f, cz }, 3000, domainW * .06f, { .9f, .15f, .2f });
            //add({ cx * .7f, domainH * .25f, cz * 1.3f }, 100, domainW * .05f, { 1, .9f, .1f });
        }

        void Step() {
            MarkCells();
            Advect();
            AddBodyForces();
            ComputeWeights();
            PrecomputeCoeffs();
            SolvePressure();
            ApplyPressureGradient();
            UpdateBodyVelocities();
            AdvectParticles();
            MarkCells();
            AdvanceBodies();
            if (++seedCounter % 5 == 0) ReseedParticles(); // reseed every 5 steps
        }

        void MarkCells() {
            std::fill(cell.begin(), cell.end(), AIR);
            for (auto & p : particles) {
                int i             = glm::clamp(int(p.x / dx), 0, nx - 1);
                int j             = glm::clamp(int(p.y / dx), 0, ny - 1);
                int k             = glm::clamp(int(p.z / dx), 0, nz - 1);
                cell[ci(i, j, k)] = FLUID;
            }
        }

        // ── reseed particles to prevent fluid loss ──
        void ReseedParticles() {
            // count per cell
            std::vector<uint8_t> cnt(nx * ny * nz, 0);
            for (auto & p : particles) {
                int i = glm::clamp(int(p.x / dx), 0, nx - 1);
                int j = glm::clamp(int(p.y / dx), 0, ny - 1);
                int k = glm::clamp(int(p.z / dx), 0, nz - 1);
                if (cnt[ci(i, j, k)] < 255) cnt[ci(i, j, k)]++;
            }
            std::mt19937                          rng(seedCounter * 13 + 7);
            std::uniform_real_distribution<float> jit(.2f, .8f);
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++) {
                        int idx = ci(i, j, k);
                        if (cell[idx] == FLUID && cnt[idx] < 1) {
                            particles.push_back({ (i + jit(rng)) * dx, (j + jit(rng)) * dx, (k + jit(rng)) * dx });
                        }
                    }
            // remove excess (cap at 4 per cell)
            if (particles.size() > size_t(nx * ny * nz * 4)) {
                std::fill(cnt.begin(), cnt.end(), 0);
                std::vector<glm::vec3> kept;
                kept.reserve(particles.size());
                for (auto & p : particles) {
                    int i = glm::clamp(int(p.x / dx), 0, nx - 1);
                    int j = glm::clamp(int(p.y / dx), 0, ny - 1);
                    int k = glm::clamp(int(p.z / dx), 0, nz - 1);
                    if (cnt[ci(i, j, k)] < 3) {
                        cnt[ci(i, j, k)]++;
                        kept.push_back(p);
                    }
                }
                particles = std::move(kept);
            }
        }

        // ── sampling ──
        float SampleF(const std::vector<float> & f, int W, int H, int D, float fi, float fj, float fk) const {
            fi       = glm::clamp(fi, 0.f, float(W - 1) - 1e-4f);
            fj       = glm::clamp(fj, 0.f, float(H - 1) - 1e-4f);
            fk       = glm::clamp(fk, 0.f, float(D - 1) - 1e-4f);
            int   i0 = int(fi), j0 = int(fj), k0 = int(fk), i1 = std::min(i0 + 1, W - 1), j1 = std::min(j0 + 1, H - 1), k1 = std::min(k0 + 1, D - 1);
            float s = fi - i0, t = fj - j0, r = fk - k0;
            auto  I = [&](int a, int b, int c) -> int { return (a * H + b) * D + c; };
            return (1 - s) * ((1 - t) * ((1 - r) * f[I(i0, j0, k0)] + r * f[I(i0, j0, k1)]) + t * ((1 - r) * f[I(i0, j1, k0)] + r * f[I(i0, j1, k1)]))
                + s * ((1 - t) * ((1 - r) * f[I(i1, j0, k0)] + r * f[I(i1, j0, k1)]) + t * ((1 - r) * f[I(i1, j1, k0)] + r * f[I(i1, j1, k1)]));
        }
        float SU(float x, float y, float z) const { return SampleF(u, nx + 1, ny, nz, x / dx, y / dx - .5f, z / dx - .5f); }
        float SV(float x, float y, float z) const { return SampleF(v, nx, ny + 1, nz, x / dx - .5f, y / dx, z / dx - .5f); }
        float SW(float x, float y, float z) const { return SampleF(w, nx, ny, nz + 1, x / dx - .5f, y / dx - .5f, z / dx); }

        glm::vec3 TraceBack(glm::vec3 p, float _dt) const {
            glm::vec3 v0  = { SU(p.x, p.y, p.z), SV(p.x, p.y, p.z), SW(p.x, p.y, p.z) };
            glm::vec3 mid = glm::clamp(p - .5f * _dt * v0, glm::vec3(0), glm::vec3(domainW, domainH, domainD));
            glm::vec3 v1  = { SU(mid.x, mid.y, mid.z), SV(mid.x, mid.y, mid.z), SW(mid.x, mid.y, mid.z) };
            return glm::clamp(p - _dt * v1, glm::vec3(0), glm::vec3(domainW, domainH, domainD));
        }

        void Advect() {
            auto uN = u, vN = v, wN = w;
            for (int i = 0; i <= nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++) {
                        auto b          = TraceBack(uPos(i, j, k), dt);
                        uN[ui(i, j, k)] = SU(b.x, b.y, b.z);
                    }
            for (int i = 0; i < nx; i++)
                for (int j = 0; j <= ny; j++)
                    for (int k = 0; k < nz; k++) {
                        auto b          = TraceBack(vPos(i, j, k), dt);
                        vN[vi(i, j, k)] = SV(b.x, b.y, b.z);
                    }
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k <= nz; k++) {
                        auto b          = TraceBack(wPos(i, j, k), dt);
                        wN[wi(i, j, k)] = SW(b.x, b.y, b.z);
                    }
            u = uN;
            v = vN;
            w = wN;
        }

        void AddBodyForces() {
            for (int i = 0; i < nx; i++)
                for (int j = 0; j <= ny; j++)
                    for (int k = 0; k < nz; k++) v[vi(i, j, k)] += dt * gravity;
            for (auto & b : bodies)
                if (! b.dragged) b.velocity.y += dt * gravity;
        }

        // ── weights with bounding-box acceleration ──
        float SolidFrac(glm::vec3 fp, int ax, const RigidBody3D & bd) const {
            float d = glm::length(fp - bd.position);
            if (d > bd.radius + dx * 1.2f) return 0;
            if (d < bd.radius - dx * 1.2f) return 1;
            int N = 3, ins = 0;
            for (int a = 0; a < N; a++)
                for (int b = 0; b < N; b++) {
                    float     ta = ((a + .5f) / N - .5f) * dx, tb = ((b + .5f) / N - .5f) * dx;
                    glm::vec3 sp = fp;
                    if (ax == 0) {
                        sp.y += ta;
                        sp.z += tb;
                    } else if (ax == 1) {
                        sp.x += ta;
                        sp.z += tb;
                    } else {
                        sp.x += ta;
                        sp.y += tb;
                    }
                    if (bd.SDF(sp) < 0) ins++;
                }
            return float(ins) / (N * N);
        }

        void ComputeWeights() {
            int uS = (nx + 1) * ny * nz, vS = nx * (ny + 1) * nz, wS = nx * ny * (nz + 1);
            wFU.assign(uS, 1);
            wFV.assign(vS, 1);
            wFW.assign(wS, 1);
            bSU.resize(bodies.size());
            bSV.resize(bodies.size());
            bSW.resize(bodies.size());
            for (size_t bi = 0; bi < bodies.size(); bi++) {
                bSU[bi].assign(uS, 0);
                bSV[bi].assign(vS, 0);
                bSW[bi].assign(wS, 0);
            }

            for (int i = 0; i <= nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++) {
                        int       idx = ui(i, j, k);
                        glm::vec3 fp  = uPos(i, j, k);
                        float     ts  = 0;
                        for (size_t bi = 0; bi < bodies.size(); bi++) {
                            float sf     = SolidFrac(fp, 0, bodies[bi]);
                            bSU[bi][idx] = sf;
                            ts += sf;
                        }
                        wFU[idx] = std::max(1 - glm::clamp(ts, 0.f, 1.f), 0.f);
                    }
            for (int i = 0; i < nx; i++)
                for (int j = 0; j <= ny; j++)
                    for (int k = 0; k < nz; k++) {
                        int       idx = vi(i, j, k);
                        glm::vec3 fp  = vPos(i, j, k);
                        float     ts  = 0;
                        for (size_t bi = 0; bi < bodies.size(); bi++) {
                            float sf     = SolidFrac(fp, 1, bodies[bi]);
                            bSV[bi][idx] = sf;
                            ts += sf;
                        }
                        wFV[idx] = std::max(1 - glm::clamp(ts, 0.f, 1.f), 0.f);
                    }
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k <= nz; k++) {
                        int       idx = wi(i, j, k);
                        glm::vec3 fp  = wPos(i, j, k);
                        float     ts  = 0;
                        for (size_t bi = 0; bi < bodies.size(); bi++) {
                            float sf     = SolidFrac(fp, 2, bodies[bi]);
                            bSW[bi][idx] = sf;
                            ts += sf;
                        }
                        wFW[idx] = std::max(1 - glm::clamp(ts, 0.f, 1.f), 0.f);
                    }
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++) {
                    wFU[ui(0, j, k)]  = 0;
                    wFU[ui(nx, j, k)] = 0;
                }
            for (int i = 0; i < nx; i++)
                for (int k = 0; k < nz; k++) {
                    wFV[vi(i, 0, k)]  = 0;
                    wFV[vi(i, ny, k)] = 0;
                }
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++) {
                    wFW[wi(i, j, 0)]  = 0;
                    wFW[wi(i, j, nz)] = 0;
                }
        }

        void PrecomputeCoeffs() {
            float dx4 = dx * dx * dx * dx;
            for (int i = 0; i <= nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++) {
                        int   idx = ui(i, j, k);
                        float a   = dt * wFU[idx] * dx / fluidRho;
                        for (size_t bi = 0; bi < bodies.size(); bi++) {
                            float sf = bSU[bi][idx];
                            a += dt * sf * sf * dx4 / bodies[bi].mass;
                        }
                        aU[idx] = a;
                    }
            for (int i = 0; i < nx; i++)
                for (int j = 0; j <= ny; j++)
                    for (int k = 0; k < nz; k++) {
                        int   idx = vi(i, j, k);
                        float a   = dt * wFV[idx] * dx / fluidRho;
                        for (size_t bi = 0; bi < bodies.size(); bi++) {
                            float sf = bSV[bi][idx];
                            a += dt * sf * sf * dx4 / bodies[bi].mass;
                        }
                        aV[idx] = a;
                    }
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k <= nz; k++) {
                        int   idx = wi(i, j, k);
                        float a   = dt * wFW[idx] * dx / fluidRho;
                        for (size_t bi = 0; bi < bodies.size(); bi++) {
                            float sf = bSW[bi][idx];
                            a += dt * sf * sf * dx4 / bodies[bi].mass;
                        }
                        aW[idx] = a;
                    }
        }

        void SolvePressure() {
            int   cS  = nx * ny * nz;
            float dx2 = dx * dx;
            pressure.assign(cS, 0);
            std::vector<float> rhs(cS, 0);
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++) {
                        if (cell[ci(i, j, k)] != FLUID) continue;
                        float r = 0;
                        r += wFU[ui(i, j, k)] * dx2 * u[ui(i, j, k)];
                        r -= wFU[ui(i + 1, j, k)] * dx2 * u[ui(i + 1, j, k)];
                        r += wFV[vi(i, j, k)] * dx2 * v[vi(i, j, k)];
                        r -= wFV[vi(i, j + 1, k)] * dx2 * v[vi(i, j + 1, k)];
                        r += wFW[wi(i, j, k)] * dx2 * w[wi(i, j, k)];
                        r -= wFW[wi(i, j, k + 1)] * dx2 * w[wi(i, j, k + 1)];
                        for (size_t bi = 0; bi < bodies.size(); bi++) {
                            auto & b = bodies[bi];
                            r += bSU[bi][ui(i, j, k)] * dx2 * b.velocity.x;
                            r -= bSU[bi][ui(i + 1, j, k)] * dx2 * b.velocity.x;
                            r += bSV[bi][vi(i, j, k)] * dx2 * b.velocity.y;
                            r -= bSV[bi][vi(i, j + 1, k)] * dx2 * b.velocity.y;
                            r += bSW[bi][wi(i, j, k)] * dx2 * b.velocity.z;
                            r -= bSW[bi][wi(i, j, k + 1)] * dx2 * b.velocity.z;
                        }
                        rhs[ci(i, j, k)] = r;
                    }
            float omega = 1.75f;
            for (int it = 0; it < 100; it++) {
                float maxR = 0;
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        for (int k = 0; k < nz; k++) {
                            if (cell[ci(i, j, k)] != FLUID) continue;
                            float d = aU[ui(i + 1, j, k)] + aU[ui(i, j, k)] + aV[vi(i, j + 1, k)] + aV[vi(i, j, k)] + aW[wi(i, j, k + 1)] + aW[wi(i, j, k)];
                            if (d < 1e-12f) continue;
                            float off = 0;
                            if (i + 1 < nx && cell[ci(i + 1, j, k)] == FLUID) off += aU[ui(i + 1, j, k)] * pressure[ci(i + 1, j, k)];
                            if (i > 0 && cell[ci(i - 1, j, k)] == FLUID) off += aU[ui(i, j, k)] * pressure[ci(i - 1, j, k)];
                            if (j + 1 < ny && cell[ci(i, j + 1, k)] == FLUID) off += aV[vi(i, j + 1, k)] * pressure[ci(i, j + 1, k)];
                            if (j > 0 && cell[ci(i, j - 1, k)] == FLUID) off += aV[vi(i, j, k)] * pressure[ci(i, j - 1, k)];
                            if (k + 1 < nz && cell[ci(i, j, k + 1)] == FLUID) off += aW[wi(i, j, k + 1)] * pressure[ci(i, j, k + 1)];
                            if (k > 0 && cell[ci(i, j, k - 1)] == FLUID) off += aW[wi(i, j, k)] * pressure[ci(i, j, k - 1)];
                            float np              = (rhs[ci(i, j, k)] + off) / d;
                            float old             = pressure[ci(i, j, k)];
                            pressure[ci(i, j, k)] = (1 - omega) * old + omega * np;
                            maxR                  = std::max(maxR, std::abs(pressure[ci(i, j, k)] - old));
                        }
                if (maxR < 1e-4f) break;
            }
        }

        void ApplyPressureGradient() {
            for (int i = 1; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++) {
                        bool l = cell[ci(i - 1, j, k)] == FLUID, r = cell[ci(i, j, k)] == FLUID;
                        if (l || r) u[ui(i, j, k)] -= dt / fluidRho * ((r ? pressure[ci(i, j, k)] : 0) - (l ? pressure[ci(i - 1, j, k)] : 0)) / dx;
                    }
            for (int i = 0; i < nx; i++)
                for (int j = 1; j < ny; j++)
                    for (int k = 0; k < nz; k++) {
                        bool b = cell[ci(i, j - 1, k)] == FLUID, t = cell[ci(i, j, k)] == FLUID;
                        if (b || t) v[vi(i, j, k)] -= dt / fluidRho * ((t ? pressure[ci(i, j, k)] : 0) - (b ? pressure[ci(i, j - 1, k)] : 0)) / dx;
                    }
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 1; k < nz; k++) {
                        bool b = cell[ci(i, j, k - 1)] == FLUID, f = cell[ci(i, j, k)] == FLUID;
                        if (b || f) w[wi(i, j, k)] -= dt / fluidRho * ((f ? pressure[ci(i, j, k)] : 0) - (b ? pressure[ci(i, j, k - 1)] : 0)) / dx;
                    }
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++) {
                    u[ui(0, j, k)]  = 0;
                    u[ui(nx, j, k)] = 0;
                }
            for (int i = 0; i < nx; i++)
                for (int k = 0; k < nz; k++) {
                    v[vi(i, 0, k)]  = 0;
                    v[vi(i, ny, k)] = 0;
                }
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++) {
                    w[wi(i, j, 0)]  = 0;
                    w[wi(i, j, nz)] = 0;
                }
        }

        void UpdateBodyVelocities() {
            float dx2 = dx * dx;
            for (size_t bi = 0; bi < bodies.size(); bi++) {
                auto & b = bodies[bi];
                if (b.dragged) continue;
                glm::vec3 f(0);
                for (int i = 0; i <= nx; i++)
                    for (int j = 0; j < ny; j++)
                        for (int k = 0; k < nz; k++) {
                            float sf = bSU[bi][ui(i, j, k)];
                            if (sf < 1e-6f) continue;
                            float pL = (i > 0 && cell[ci(i - 1, j, k)] == FLUID) ? pressure[ci(i - 1, j, k)] : 0;
                            float pR = (i < nx && cell[ci(i, j, k)] == FLUID) ? pressure[ci(i, j, k)] : 0;
                            f.x -= sf * dx2 * (pR - pL);
                        }
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j <= ny; j++)
                        for (int k = 0; k < nz; k++) {
                            float sf = bSV[bi][vi(i, j, k)];
                            if (sf < 1e-6f) continue;
                            float pB = (j > 0 && cell[ci(i, j - 1, k)] == FLUID) ? pressure[ci(i, j - 1, k)] : 0;
                            float pT = (j < ny && cell[ci(i, j, k)] == FLUID) ? pressure[ci(i, j, k)] : 0;
                            f.y -= sf * dx2 * (pT - pB);
                        }
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        for (int k = 0; k <= nz; k++) {
                            float sf = bSW[bi][wi(i, j, k)];
                            if (sf < 1e-6f) continue;
                            float pB = (k > 0 && cell[ci(i, j, k - 1)] == FLUID) ? pressure[ci(i, j, k - 1)] : 0;
                            float pF = (k < nz && cell[ci(i, j, k)] == FLUID) ? pressure[ci(i, j, k)] : 0;
                            f.z -= sf * dx2 * (pF - pB);
                        }
                b.velocity += dt / b.mass * f + dt / b.mass * b.dragForce;
                b.dragForce = glm::vec3(0);
            }
        }

        void AdvectParticles() {
            for (auto & p : particles) {
                for (auto & b : bodies)
                    if (b.SDF(p) < 0) {
                        glm::vec3 n = glm::normalize(p - b.position);
                        p           = b.position + n * (b.radius + dx * .1f);
                    }
                glm::vec3 v0  = { SU(p.x, p.y, p.z), SV(p.x, p.y, p.z), SW(p.x, p.y, p.z) };
                glm::vec3 mid = glm::clamp(p + .5f * dt * v0, glm::vec3(dx * .01f), glm::vec3(domainW, domainH, domainD) - dx * .01f);
                glm::vec3 v1  = { SU(mid.x, mid.y, mid.z), SV(mid.x, mid.y, mid.z), SW(mid.x, mid.y, mid.z) };
                p             = glm::clamp(p + dt * v1, glm::vec3(dx * .01f), glm::vec3(domainW, domainH, domainD) - dx * .01f);
            }
        }

        void AdvanceBodies() {
            for (auto & b : bodies) {
                b.position += dt * b.velocity;
                float m  = b.radius + dx * .05f;
                auto  cl = [&](float & p, float & v, float lo, float hi) {
                    if (p < lo + m) {
                        p = lo + m;
                        if (v < 0.0f)
                            v = 0.0f;
                    }
                    if (p > hi - m) {
                        p = hi - m;
                        if (v > 0.0f)
                            v = 0.0f;
                    }
                };
                cl(b.position.x, b.velocity.x, 0, domainW);
                cl(b.position.y, b.velocity.y, 0, domainH);
                cl(b.position.z, b.velocity.z, 0, domainD);
            }
        }

        int PickBody(glm::vec3 ro, glm::vec3 rd) const {
            float best = 1e30f;
            int   idx  = -1;
            for (int i = (int) bodies.size() - 1; i >= 0; i--) {
                glm::vec3 oc = ro - bodies[i].position;
                float     b = glm::dot(oc, rd), c = glm::dot(oc, oc) - bodies[i].radius * bodies[i].radius;
                float     d = b * b - c;
                if (d >= 0) {
                    float t = -b - std::sqrt(d);
                    if (t > 0 && t < best) {
                        best = t;
                        idx  = i;
                    }
                }
            }
            return idx;
        }
        void DragBody(int i, glm::vec3 t) {
            if (i < 0 || i >= (int) bodies.size()) return;
            bodies[i].dragForce = 600.f * (t - bodies[i].position) - 15.f * bodies[i].velocity;
            bodies[i].dragged   = true;
        }
        void ReleaseDrag(int i) {
            if (i >= 0 && i < (int) bodies.size()) {
                bodies[i].dragged   = false;
                bodies[i].dragForce = glm::vec3(0);
            }
        }
    };
} // namespace VCX::Labs::GettingStarted
