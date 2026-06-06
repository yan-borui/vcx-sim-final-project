#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>

#include "Labs/FinalProject/FreeSurfaceSeparationSimulator.h"
#include "Labs/FinalProject/SubgridSimulator.h"

namespace {
    using VCX::Labs::Final::FreeSurfaceSeparationSimulator;
    using VCX::Labs::Final::Simulator;
    using VCX::Labs::Final::SubgridSimulator;

    void require(bool condition, std::string_view message) {
        if (! condition)
            throw std::runtime_error(std::string(message));
    }

    int offset(Simulator const & simulation, glm::ivec3 const cell) {
        return cell.x
            + simulation.m_iCellX * (cell.y + simulation.m_iCellY * cell.z);
    }

    float divergence(Simulator const & simulation, glm::ivec3 const cell) {
        int const center = offset(simulation, cell);
        return simulation.m_vel[offset(simulation, cell + glm::ivec3(1, 0, 0))].x
            - simulation.m_vel[center].x
            + simulation.m_vel[offset(simulation, cell + glm::ivec3(0, 1, 0))].y
            - simulation.m_vel[center].y
            + simulation.m_vel[offset(simulation, cell + glm::ivec3(0, 0, 1))].z
            - simulation.m_vel[center].z;
    }

    bool isFinite(glm::vec3 const value) {
        return std::isfinite(value.x)
            && std::isfinite(value.y)
            && std::isfinite(value.z);
    }

    void prepareSingleWallCell(
        FreeSurfaceSeparationSimulator & simulation,
        float                            wallVelocity) {
        simulation.setupScene(8);
        std::fill(simulation.m_type.begin(), simulation.m_type.end(), Simulator::EMPTY_CELL);
        std::fill(simulation.m_vel.begin(), simulation.m_vel.end(), glm::vec3(0.0f));

        glm::ivec3 const fluidCell { 1, 3, 3 };
        simulation.m_type[offset(simulation, fluidCell)]  = Simulator::FLUID_CELL;
        simulation.m_vel[offset(simulation, fluidCell)].x = wallVelocity;
    }

    void testStandardWallProjection() {
        FreeSurfaceSeparationSimulator simulation;
        prepareSingleWallCell(simulation, 1.0f);
        simulation.enableWallSeparation = false;

        simulation.solveIncompressibility(200, 0.01f, 1.0f, false);

        glm::ivec3 const fluidCell { 1, 3, 3 };
        require(
            std::abs(simulation.m_vel[offset(simulation, fluidCell)].x) < 1e-5f,
            "standard wall must enforce zero normal velocity");
        require(
            std::abs(divergence(simulation, fluidCell)) < 1e-4f,
            "standard wall projection must remain divergence-free after enforcement");
    }

    void testOutwardWallSeparation() {
        FreeSurfaceSeparationSimulator simulation;
        prepareSingleWallCell(simulation, 1.0f);
        simulation.enableWallSeparation = true;

        simulation.solveIncompressibility(200, 0.01f, 1.0f, false);

        glm::ivec3 const fluidCell { 1, 3, 3 };
        float const      normalVelocity = simulation.m_vel[offset(simulation, fluidCell)].x;
        require(normalVelocity > 1e-4f, "outward wall velocity must be allowed to separate");
        require(simulation.separatingWallFaces == 1, "outward wall face must be classified as separating");
    }

    void testInwardWallContact() {
        FreeSurfaceSeparationSimulator simulation;
        prepareSingleWallCell(simulation, -1.0f);
        simulation.enableWallSeparation = true;

        simulation.solveIncompressibility(200, 0.01f, 1.0f, false);

        glm::ivec3 const fluidCell { 1, 3, 3 };
        float const      normalVelocity = simulation.m_vel[offset(simulation, fluidCell)].x;
        require(normalVelocity >= -1e-5f, "wall separation must prevent penetration");
        require(std::abs(normalVelocity) < 1e-5f, "inward wall velocity must remain in contact");
        require(simulation.separatingWallFaces == 0, "contacting wall face must not be classified as separating");
    }

    void testSeparationScenarioRemainsFiniteAndNonpenetrating() {
        FreeSurfaceSeparationSimulator simulation;
        simulation.setupScene(16);
        simulation.enableWallSeparation = true;

        for (int frame = 0; frame < 40; ++frame)
            simulation.SimulateTimestep(0.008f);

        require(
            std::all_of(
                simulation.m_particlePos.begin(),
                simulation.m_particlePos.end(),
                isFinite),
            "wall separation scenario produced a non-finite particle position");
        require(
            std::all_of(
                simulation.m_particleVel.begin(),
                simulation.m_particleVel.end(),
                isFinite),
            "wall separation scenario produced a non-finite particle velocity");
        require(std::isfinite(simulation.pressureResidual), "wall separation pressure solve failed");

        glm::ivec3 const cellOffsets[] = {
            { -1,  0,  0 },
            {  1,  0,  0 },
            {  0, -1,  0 },
            {  0,  1,  0 },
            {  0,  0, -1 },
            {  0,  0,  1 },
        };
        glm::ivec3 const faceOffsets[] = {
            { 0, 0, 0 },
            { 1, 0, 0 },
            { 0, 0, 0 },
            { 0, 1, 0 },
            { 0, 0, 0 },
            { 0, 0, 1 },
        };
        int const   directions[]      = { 0, 0, 1, 1, 2, 2 };
        float const divergenceSigns[] = { -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f };

        for (int idx = 0; idx < simulation.m_iNumCells; ++idx) {
            if (simulation.m_type[idx] != Simulator::FLUID_CELL)
                continue;
            glm::ivec3 const cell {
                idx % simulation.m_iCellX,
                (idx / simulation.m_iCellX) % simulation.m_iCellY,
                idx / (simulation.m_iCellX * simulation.m_iCellY),
            };
            for (int side = 0; side < 6; ++side) {
                glm::ivec3 const neighbor = cell + cellOffsets[side];
                bool const       solid =
                    neighbor.x < 0 || neighbor.x >= simulation.m_iCellX
                    || neighbor.y < 0 || neighbor.y >= simulation.m_iCellY
                    || neighbor.z < 0 || neighbor.z >= simulation.m_iCellZ
                    || simulation.m_s[offset(simulation, neighbor)] <= 0.0f;
                if (! solid)
                    continue;

                glm::ivec3 const face = cell + faceOffsets[side];
                float const      separationSpeed =
                    -divergenceSigns[side]
                    * simulation.m_vel[offset(simulation, face)][directions[side]];
                require(
                    separationSpeed >= -2e-4f,
                    "wall separation scenario allowed fluid to penetrate a wall");
            }
        }
    }

    float runSubgridChannel(bool useSubgridWeights) {
        int constexpr resolution = 16;
        float const h            = 1.0f / float(resolution);

        VCX::Labs::Final::RigidBody body;
        body.Reset(
            { 0.03f, 0.0f, 0.0f },
            { 0.0f, 0.0f, 0.0f },
            { 0.18f, 1.0f - 3.0f * h, 1.0f - 2.0f * h },
            1.0f,
            { 1.0f, 1.0f, 1.0f });
        body.isStatic = true;
        body.ComputeInertia();

        SubgridSimulator simulation;
        simulation.m_body            = &body;
        simulation.useSubgridWeights = useSubgridWeights;
        simulation.setupScene(resolution);
        for (int frame = 0; frame < 80; ++frame)
            simulation.SimulateTimestep(0.01f);

        require(
            std::all_of(
                simulation.m_particlePos.begin(),
                simulation.m_particlePos.end(),
                isFinite),
            "sub-grid scenario produced a non-finite particle position");
        require(simulation.partiallyOpenFaceCount > 0, "sub-grid scenario did not detect cut MAC volumes");
        require(simulation.pressureSolveSucceeded, "sub-grid pressure solve failed");
        require(std::isfinite(simulation.pressureResidual), "sub-grid pressure residual is not finite");

        float maxX = -std::numeric_limits<float>::infinity();
        for (glm::vec3 const position : simulation.m_particlePos)
            maxX = std::max(maxX, position.x);
        return maxX;
    }

    void testSubgridChannelPreservesHalfCellFlow() {
        float const weightedMaxX = runSubgridChannel(true);
        float const binaryMaxX   = runSubgridChannel(false);

        require(weightedMaxX > 0.12f, "sub-grid weights did not carry fluid through the half-cell gap");
        require(
            weightedMaxX > binaryMaxX + 0.08f,
            "sub-grid and voxelized channel flow are not meaningfully different");
    }
} // namespace

int main() {
    try {
        testStandardWallProjection();
        testOutwardWallSeparation();
        testInwardWallContact();
        testSeparationScenarioRemainsFiniteAndNonpenetrating();
        testSubgridChannelPreservesHalfCellFlow();
    } catch (std::exception const & error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }

    std::cout << "All simulator tests passed.\n";
    return 0;
}
