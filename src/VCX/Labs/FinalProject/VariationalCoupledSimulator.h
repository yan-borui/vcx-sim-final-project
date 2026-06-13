#pragma once

#include <array>
#include <cstdint>
#include <unordered_map>

#include "Labs/FinalProject/FluidSimulator.h"

namespace VCX::Labs::Final {
    struct VariationalCoupledSimulator : public Simulator {
        bool  useSubgridWeights           = true;
        bool  enableWallSeparation        = true;
        int   volumeSamplesPerAxis        = 4;
        float pressureResidual            = 0.0f;
        float wallSeparationKktResidual   = 0.0f;
        float wallSeparationQpResidual    = 0.0f;
        float wallSeparationSpeedResidual = 0.0f;
        float maximumBoundaryPressure     = 0.0f;
        int   wallSeparationIterations    = 0;
        int   wallSeparationCandidateCount = 0;
        bool  pressureSolveSucceeded      = true;

        void setupScene(int res) override;
        void solveIncompressibility(
            int   numIters,
            float dt,
            float overRelaxation,
            bool  compensateDrift) override;

    private:
        struct BoundaryPressureSample {
            int   FaceIndex;
            int   Direction;
            float DivergenceSign;
            float BoundaryWeight;
            float BoundaryPressure;
            bool  Contact;
        };

        std::array<std::vector<float>, 3> _faceOpenFraction;
        std::array<std::vector<float>, 3> _faceFluidFraction;
        std::unordered_map<std::uint64_t, double> _wallPressureWarmStart;

        bool      isValidCell(glm::ivec3 const & cell) const;
        int       gridOffset(glm::ivec3 const & cell) const;
        glm::vec3 pressureCellCenter(glm::ivec3 const & cell) const;
        glm::vec3 faceCenter(glm::ivec3 const & face, int dir) const;
        bool      isTankFaceOpen(glm::ivec3 const & face, int dir) const;
        bool      isSolidPressureCell(glm::ivec3 const & cell) const;
        bool      hasFluidSupportAcrossOpenFace(glm::ivec3 const & cell) const;
        bool      isPressureUnknownCell(glm::ivec3 const & cell) const;
        bool      isWallSeparationCandidate(glm::ivec3 const & cell) const;
        glm::vec2 estimateFaceFractions(
            glm::ivec3 const &         face,
            int                        dir,
            std::vector<float> const & liquidLevelSet) const;
        float     faceWeight(int dir, int faceIdx) const;
        float     wallGhostPressureScale(
            glm::ivec3 const & cell,
            glm::ivec3 const & face,
            int                dir,
            bool               bodyBoundary) const;
        float     solidVelocity(glm::ivec3 const & face, int dir, glm::ivec3 const & solidCell) const;

        void updateFaceFluidFractions(std::vector<float> const & liquidLevelSet);
        void applyRigidBodyFeedback(
            std::vector<BoundaryPressureSample> const & boundarySamples,
            float                                      dt);
    };
} // namespace VCX::Labs::Final
