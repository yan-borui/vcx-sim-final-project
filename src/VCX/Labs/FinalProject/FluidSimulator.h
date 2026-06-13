#pragma once

#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <glm/glm.hpp>
#include <iostream>
#include <utility>
#include <vector>
#include "Labs/FinalProject/RigidBody.h"


namespace VCX::Labs::Final {
    struct Simulator {
        enum CellType {
            EMPTY_CELL = 0,
            FLUID_CELL = 1,
            SOLID_CELL = 2
        };

        enum ColorMode {
            Default = 0, Velocity, Density, Pressure
        };
        ColorMode              m_colorMode = Default;
        
        std::vector<glm::vec3> m_particlePos; // Particle m_particlePos
        std::vector<glm::vec3> m_particleVel; // Particle Velocity
        std::vector<glm::vec3> m_particleColor;

        float m_fRatio = 0.95f;
        int   m_iCellX;
        int   m_iCellY;
        int   m_iCellZ;
        float m_h;
        float m_fInvSpacing;
        int   m_iNumCells;

        int   m_iNumSpheres;
        float m_particleRadius;

        std::vector<glm::vec3> m_vel;
        std::vector<glm::vec3> m_pre_vel;
        std::vector<float>     m_near_num[3];

        std::vector<int> m_hashtable;
        std::vector<int> m_hashtableindex;

        std::vector<float> m_p;               // Pressure array
        std::vector<float> m_s;               // 0.0 for solid cells, 1.0 for fluid cells, used to update m_type
        std::vector<float> m_baseS;           // Static tank boundary copied before dynamic rigid voxelization
        std::vector<int>   m_type;            // Flags array (const int EMPTY_CELL = 0; const int FLUID_CELL = 1; const int SOLID_CELL = 2;)
                                              // m_type = SOLID_CELL if m_s == 0.0;
                                              // m_type = FLUID_CELL if has particle and m_s == 1;
                                              // m_type = EMPTY_CELL if has No particle and m_s == 1;
        std::vector<float> m_particleDensity; // Particle Density per cell, saved in the grid cell
        float              m_particleRestDensity;

        glm::vec3 gravity { 0, -9.81f, 0 };
        RigidBody* m_body = nullptr;
        glm::vec3 m_boxPos;
        glm::vec3 m_boxDim;
        glm::vec3 m_boxVel;
        glm::vec3 m_feedbackForce = {0.0f, 0.0f, 0.0f};
        glm::vec3 m_feedbackTorque = {0.0f, 0.0f, 0.0f};

        int   numSubSteps       = 1;
        int   numParticleIters  = 5;
        int   numPressureIters  = 30;
        bool  separateParticles = true;
        float overRelaxation    = 0.5f;
        bool  compensateDrift   = true;
        float compensateDriftWeight = 0.1f;
        bool  voxelizeDynamicBody = true;


        void integrateParticles(float timeStep);
        void pushParticlesApart(int numIters);
        void handleParticleCollisions();
        void updateParticleDensity();
        void extrapolateGridVelocities(int layers);
        std::vector<float> buildParticleLevelSet() const;
        float sampleCellCenteredField(
            std::vector<float> const & field,
            glm::vec3 const &         position) const;
        float ghostFluidPressureScale(
            std::vector<float> const & levelSet,
            int                        fluidCell,
            int                        airCell) const;
        bool      isValidCell(glm::ivec3 const & index) const;
        glm::vec3 cellCenter(glm::ivec3 const & index) const;
        void      rebuildSolidCellsFromBody();

        virtual void        transferVelocities(bool toGrid, float flipRatio);
        virtual void        solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift);
        void        updateParticleColors();
        inline bool isValidVelocity(int i, int j, int k, int dir) {
            if (dir == 0) { // X 方向：检查左侧格子 (i-1) 和当前格子 (i)
                bool leftFluid = (i > 0 && m_type[index2GridOffset({ i - 1, j, k })] == FLUID_CELL);
                bool currentFluid = (i < m_iCellX) && (m_type[index2GridOffset({ i, j, k })] == FLUID_CELL);
                return leftFluid || currentFluid;
            }
            if (dir == 1) { // Y 方向：检查下方格子 (j-1) 和当前格子 (j)
                bool downFluid = (j > 0 && m_type[index2GridOffset({ i, j - 1, k })] == FLUID_CELL);
                bool currentFluid = (j < m_iCellY) && (m_type[index2GridOffset({ i, j, k })] == FLUID_CELL);
                return downFluid || currentFluid;
            }
            if (dir == 2) { // Z 方向：检查后方格子 (k-1) 和当前格子 (k)
                bool backFluid = (k > 0 && m_type[index2GridOffset({ i, j, k - 1 })] == FLUID_CELL);
                bool currentFluid = (k < m_iCellZ) && (m_type[index2GridOffset({ i, j, k })] == FLUID_CELL);
                return backFluid || currentFluid;
            }
            return false;
        }

        inline int index2GridOffset(glm::ivec3 index) {
            return index.x + m_iCellX * (index.y + m_iCellY * index.z);
        }

        void SimulateTimestep(float const dt) {
            float sdt = dt / numSubSteps;

            for (int step = 0; step < numSubSteps; step++) {
                integrateParticles(sdt);
                handleParticleCollisions();
                if (separateParticles)
                    pushParticlesApart(numParticleIters);
                handleParticleCollisions();
                transferVelocities(true, m_fRatio);
                updateParticleDensity();
                solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
                transferVelocities(false, m_fRatio);
            }
            updateParticleColors();
        }

        virtual void setupScene(int res) {
            glm::vec3 tank(1.0f);
            glm::vec3 relWater = { 0.6f, 0.7f, 0.6f };

            float _h      = tank.y / res;
            float point_r = 0.3 * _h;
            float dx      = 1.7 * point_r;
            float dy      = sqrt(3.0) / 2.0 * dx;
            float dz      = dx;

            int numX = floor((relWater.x * tank.x - 2.0 * _h - 2.0 * point_r) / dx);
            int numY = floor((relWater.y * tank.y - 2.0 * _h - 2.0 * point_r) / dy);
            int numZ = floor((relWater.z * tank.z - 2.0 * _h - 2.0 * point_r) / dz);

            // update object member attributes
            m_iNumSpheres    = numX * numY * numZ;
            m_iCellX         = res + 1;
            m_iCellY         = res + 1;
            m_iCellZ         = res + 1;
            m_h              = 1.0 / float(res);
            m_fInvSpacing    = float(res);
            m_iNumCells      = m_iCellX * m_iCellY * m_iCellZ;
            m_particleRadius = 0.7 * point_r; // modified

            // update particle array
            m_particlePos.clear();
            m_particlePos.resize(m_iNumSpheres, glm::vec3(0.0f));
            m_particleVel.clear();
            m_particleVel.resize(m_iNumSpheres, glm::vec3(0.0f));
            m_particleColor.clear();
            m_particleColor.resize(m_iNumSpheres, glm::vec3(1.0f));
            m_hashtable.clear();
            m_hashtable.resize(m_iNumSpheres, 0);
            m_hashtableindex.clear();
            m_hashtableindex.resize(m_iNumCells + 1, 0);

            // update grid array
            m_vel.clear();
            m_vel.resize(m_iNumCells, glm::vec3(0.0f));
            m_pre_vel.clear();
            m_pre_vel.resize(m_iNumCells, glm::vec3(0.0f));
            for (int i = 0; i < 3; ++i) {
                m_near_num[i].clear();
                m_near_num[i].resize(m_iNumCells, 0.0f);
            }

            m_p.clear();
            m_p.resize(m_iNumCells, 0.0);
            m_s.clear();
            m_s.resize(m_iNumCells, 0.0);
            m_baseS.clear();
            m_baseS.resize(m_iNumCells, 0.0);
            m_type.clear();
            m_type.resize(m_iNumCells, 0);
            m_particleDensity.clear();
            m_particleDensity.resize(m_iNumCells, 0.0f);

            // the rest density can be assigned after scene initialization
            m_particleRestDensity = 0.0;

            // create particles
            int p = 0;
            for (int i = 0; i < numX; i++) {
                for (int j = 0; j < numY; j++) {
                    for (int k = 0; k < numZ; k++) {
                        m_particlePos[p++] = glm::vec3(m_h + point_r + dx * i + (j % 2 == 0 ? 0.0 : point_r), m_h + point_r + dy * j, m_h + point_r + dz * k + (j % 2 == 0 ? 0.0 : point_r)) + glm::vec3(-0.5f);
                    }
                }
            }
            // setup grid cells for tank
            int n = m_iCellY * m_iCellZ;
            int m = m_iCellZ;

            for (int i = 0; i < m_iCellX; i++) {
                for (int j = 0; j < m_iCellY; j++) {
                    for (int k = 0; k < m_iCellZ; k++) {
                        float s = 1.0; // fluid
                        if (i == 0 || i >= m_iCellX - 2 || j == 0 || j >= m_iCellY - 2 || k == 0 || k >= m_iCellZ - 2)
                            s = 0.0f; // solid
                        m_s[i * n + j * m + k] = s;
                        m_baseS[i * n + j * m + k] = s;
                    }
                }
            }
        }
    };
} // namespace VCX::Labs::Final
