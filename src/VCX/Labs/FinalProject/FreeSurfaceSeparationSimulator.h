#pragma once

#include "Labs/FinalProject/FluidSimulator.h"

namespace VCX::Labs::Final {
    struct FreeSurfaceSeparationSimulator : public Simulator {
        bool  enableWallSeparation              = true;
        int   wallSeparationCandidates          = 0;
        int   separatingWallCells               = 0;
        int   wallSeparationActiveSetIterations = 0;
        int   separatingWallFaces               = 0;
        int   wallContactParticles              = 0;
        float minimumContactPressure            = 0.0f;
        float pressureResidual                  = 0.0f;
        float wallSeparationKktResidual         = 0.0f;
        float averageLeftWallDistance           = 0.0f;
        bool  highlightWallContact              = true;

        void setupScene(int res) override;
        void solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift) override;
        void SimulateTimestep(float dt);

    private:
        void setupSplashParticles();
        void constrainParticlesToTank();
        void updateWallContactDiagnostics();

        bool isValidCell(glm::ivec3 const & cell) const;
        bool isSolidCell(glm::ivec3 const & cell) const;
        bool isWallSeparationCandidate(glm::ivec3 const & cell) const;
        int  gridOffset(glm::ivec3 const & cell) const;
    };
} // namespace VCX::Labs::Final
