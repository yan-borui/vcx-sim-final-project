#pragma once
#include "FluidSimulator.h"

namespace VCX::Labs::Final {
    struct APICSimulator : public Simulator {
        std::vector<glm::mat3> m_particleC; // APIC C array

        virtual void transferVelocities(bool toGrid, float flipRatio) override;
        virtual void setupScene(int res) override;
    };
}