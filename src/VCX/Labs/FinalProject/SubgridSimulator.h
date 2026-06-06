#pragma once

#include <array>

#include "Labs/FinalProject/FluidSimulator.h"

namespace VCX::Labs::Final {
    struct SubgridSimulator : public Simulator {
        bool  useSubgridWeights       = true;
        int   volumeSamplesPerAxis    = 4;
        int   partiallyOpenFaceCount  = 0;
        float minimumOpenFaceFraction = 1.0f;
        float pressureResidual        = 0.0f;
        bool  pressureSolveSucceeded  = true;

        void setupScene(int res) override;
        void solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift) override;

    private:
        std::array<std::vector<float>, 3> _faceFluidFraction;

        void      setupChannelParticles();
        void      updateFaceFluidFractions();
        float     estimateFaceFluidFraction(glm::ivec3 const & face, int dir);
        float     faceWeight(int dir, int idx) const;
        bool      isTankFaceOpen(glm::ivec3 const & face, int dir);
        bool      isSolidPressureCell(glm::ivec3 const & cell);
        bool      isValidCell(glm::ivec3 const & cell) const;
        glm::vec3 pressureCellCenter(glm::ivec3 const & cell) const;
    };
}
