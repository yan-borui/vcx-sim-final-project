#pragma once

#include <array>

#include "Labs/FinalProject/FluidSimulator.h"

namespace VCX::Labs::Final {
    struct VariationalCoupledSimulator : public Simulator {
        bool  useSubgridWeights      = true;
        bool  enableWallSeparation   = true;
        int   volumeSamplesPerAxis   = 4;
        float pressureResidual       = 0.0f;
        bool  pressureSolveSucceeded = true;

        void setupScene(int res) override;
        void solveIncompressibility(
            int   numIters,
            float dt,
            float overRelaxation,
            bool  compensateDrift) override;

    private:
        std::array<std::vector<float>, 3> _faceFluidFraction;

        bool      isValidCell(glm::ivec3 const & cell) const;
        int       gridOffset(glm::ivec3 const & cell) const;
        glm::vec3 pressureCellCenter(glm::ivec3 const & cell) const;
        glm::vec3 faceCenter(glm::ivec3 const & face, int dir) const;
        bool      isTankFaceOpen(glm::ivec3 const & face, int dir) const;
        bool      isSolidPressureCell(glm::ivec3 const & cell) const;
        bool      isWallSeparationCandidate(glm::ivec3 const & cell) const;
        float     estimateFaceFluidFraction(glm::ivec3 const & face, int dir) const;
        float     faceWeight(int dir, int faceIdx) const;
        float     solidVelocity(glm::ivec3 const & face, int dir, glm::ivec3 const & solidCell) const;

        void updateFaceFluidFractions();
        void applyRigidBodyFeedback(std::vector<int> const & pressureCells);
    };
} // namespace VCX::Labs::Final
