#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

#include "Labs/FinalProject/FreeSurfaceSeparationSimulator.h"
#include "Labs/FinalProject/SubgridSimulator.h"
#include "Labs/FinalProject/VariationalCoupledSimulator.h"

namespace {
    using VCX::Labs::Final::FreeSurfaceSeparationSimulator;
    using VCX::Labs::Final::Simulator;
    using VCX::Labs::Final::SubgridSimulator;
    using VCX::Labs::Final::VariationalCoupledSimulator;

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

    void testDisconnectedClosedPressureComponent() {
        FreeSurfaceSeparationSimulator simulation;
        simulation.setupScene(8);
        simulation.enableWallSeparation = false;
        std::fill(simulation.m_type.begin(), simulation.m_type.end(), Simulator::EMPTY_CELL);
        std::fill(simulation.m_vel.begin(), simulation.m_vel.end(), glm::vec3(0.0f));

        glm::ivec3 const closedLeft { 4, 4, 4 };
        glm::ivec3 const closedRight { 5, 4, 4 };
        glm::ivec3 const openCell { 2, 2, 2 };
        simulation.m_type[offset(simulation, closedLeft)]  = Simulator::FLUID_CELL;
        simulation.m_type[offset(simulation, closedRight)] = Simulator::FLUID_CELL;
        simulation.m_type[offset(simulation, openCell)]    = Simulator::FLUID_CELL;

        glm::ivec3 const neighbors[] = {
            { -1,  0,  0 },
            {  1,  0,  0 },
            {  0, -1,  0 },
            {  0,  1,  0 },
            {  0,  0, -1 },
            {  0,  0,  1 },
        };
        for (glm::ivec3 const cell : { closedLeft, closedRight }) {
            simulation.m_s[offset(simulation, cell)] = 1.0f;
            for (glm::ivec3 const neighborOffset : neighbors) {
                glm::ivec3 const neighbor = cell + neighborOffset;
                if (neighbor != closedLeft && neighbor != closedRight)
                    simulation.m_s[offset(simulation, neighbor)] = 0.0f;
            }
        }
        simulation.m_vel[offset(simulation, closedRight)].x = 1.0f;

        simulation.solveIncompressibility(200, 0.01f, 1.0f, false);

        require(std::isfinite(simulation.pressureResidual), "disconnected pressure component solve failed");
        require(
            std::abs(divergence(simulation, closedLeft)) < 1e-4f,
            "left cell in closed pressure component is not divergence-free");
        require(
            std::abs(divergence(simulation, closedRight)) < 1e-4f,
            "right cell in closed pressure component is not divergence-free");

        VariationalCoupledSimulator coupled;
        coupled.setupScene(8);
        coupled.enableWallSeparation = false;
        std::fill(coupled.m_type.begin(), coupled.m_type.end(), Simulator::EMPTY_CELL);
        std::fill(coupled.m_vel.begin(), coupled.m_vel.end(), glm::vec3(0.0f));
        coupled.m_type[offset(coupled, closedLeft)]  = Simulator::FLUID_CELL;
        coupled.m_type[offset(coupled, closedRight)] = Simulator::FLUID_CELL;
        coupled.m_type[offset(coupled, openCell)]    = Simulator::FLUID_CELL;
        for (glm::ivec3 const cell : { closedLeft, closedRight }) {
            coupled.m_s[offset(coupled, cell)] = 1.0f;
            for (glm::ivec3 const neighborOffset : neighbors) {
                glm::ivec3 const neighbor = cell + neighborOffset;
                if (neighbor != closedLeft && neighbor != closedRight)
                    coupled.m_s[offset(coupled, neighbor)] = 0.0f;
            }
        }
        coupled.m_vel[offset(coupled, closedRight)].x          = 1.0f;
        coupled.m_particleRestDensity                          = 1.0f;
        coupled.m_particleDensity[offset(coupled, closedLeft)] = 2.0f;
        coupled.solveIncompressibility(200, 0.01f, 1.0f, true);
        require(
            coupled.pressureSolveSucceeded && std::isfinite(coupled.pressureResidual),
            "closed component with drift compensation has an unresolved pressure null-space");
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

    void prepareSingleCoupledWallCell(
        VariationalCoupledSimulator & simulation,
        float                         wallVelocity) {
        simulation.setupScene(8);
        std::fill(simulation.m_type.begin(), simulation.m_type.end(), Simulator::EMPTY_CELL);
        std::fill(simulation.m_vel.begin(), simulation.m_vel.end(), glm::vec3(0.0f));

        glm::ivec3 const fluidCell { 1, 3, 3 };
        simulation.m_type[offset(simulation, fluidCell)]  = Simulator::FLUID_CELL;
        simulation.m_vel[offset(simulation, fluidCell)].x = wallVelocity;
    }

    void testCoupledStandardWallProjection() {
        VariationalCoupledSimulator simulation;
        prepareSingleCoupledWallCell(simulation, 1.0f);
        simulation.enableWallSeparation = false;

        simulation.solveIncompressibility(200, 0.01f, 1.0f, false);

        glm::ivec3 const fluidCell { 1, 3, 3 };
        require(
            std::abs(simulation.m_vel[offset(simulation, fluidCell)].x) < 1e-5f,
            "coupled standard wall must enforce zero normal velocity");
        require(
            std::abs(divergence(simulation, fluidCell)) < 1e-4f,
            "coupled standard wall projection must remain divergence-free");
    }

    void testCoupledOutwardWallSeparation() {
        VariationalCoupledSimulator simulation;
        prepareSingleCoupledWallCell(simulation, 1.0f);
        simulation.enableWallSeparation = true;

        simulation.solveIncompressibility(200, 0.01f, 1.0f, false);

        glm::ivec3 const fluidCell { 1, 3, 3 };
        require(
            simulation.m_vel[offset(simulation, fluidCell)].x > 1e-4f,
            "coupled solver must allow outward wall separation");
    }

    void prepareSingleBodyWallCell(
        VariationalCoupledSimulator & simulation,
        VCX::Labs::Final::RigidBody & body,
        float                         bodyVelocity) {
        simulation.setupScene(8);
        body.Reset(
            { 0.1f, -0.0625f, -0.0625f },
            { bodyVelocity, 0.0f, 0.0f },
            { 0.3f, 0.2f, 0.2f },
            1.0f,
            { 1.0f, 1.0f, 1.0f });
        simulation.m_body = &body;

        std::fill(simulation.m_type.begin(), simulation.m_type.end(), Simulator::EMPTY_CELL);
        std::fill(simulation.m_vel.begin(), simulation.m_vel.end(), glm::vec3(0.0f));
        simulation.m_type[offset(simulation, { 3, 3, 3 })] = Simulator::FLUID_CELL;
    }

    void testCoupledMovingBodyContactAndSeparation() {
        VariationalCoupledSimulator contactSimulation;
        VCX::Labs::Final::RigidBody contactBody;
        prepareSingleBodyWallCell(contactSimulation, contactBody, -0.5f);
        contactSimulation.enableWallSeparation = true;
        contactSimulation.solveIncompressibility(200, 0.01f, 1.0f, false);

        glm::ivec3 const bodyFace { 4, 3, 3 };
        require(
            std::abs(contactSimulation.m_vel[offset(contactSimulation, bodyFace)].x + 0.5f) < 1e-4f,
            "body moving into fluid must keep the normal velocity in contact");
        require(contactSimulation.pressureSolveSucceeded, "moving-body contact pressure solve failed");

        VariationalCoupledSimulator separationSimulation;
        VCX::Labs::Final::RigidBody separationBody;
        prepareSingleBodyWallCell(separationSimulation, separationBody, 0.5f);
        separationSimulation.enableWallSeparation = true;
        separationSimulation.solveIncompressibility(200, 0.01f, 1.0f, false);

        float const fluidVelocity =
            separationSimulation.m_vel[offset(separationSimulation, bodyFace)].x;
        require(
            fluidVelocity - separationBody.velocity.x <= 1e-4f,
            "body moving away from fluid produced an interpenetrating relative velocity");
        require(
            std::abs(fluidVelocity - separationBody.velocity.x) > 1e-3f,
            "body moving away from fluid incorrectly kept a sticking contact");
        require(separationSimulation.pressureSolveSucceeded, "moving-body separation pressure solve failed");
    }

    void testCoupledScenarioRemainsFinite() {
        VariationalCoupledSimulator simulation;
        simulation.setupScene(24);

        VCX::Labs::Final::RigidBody body;
        body.Reset(
            { 0.0f, 0.3f, 0.0f },
            { 0.0f, 0.0f, 0.0f },
            { 0.3f, 0.3f, 0.3f },
            0.3f,
            { 1.0f, 1.0f, 1.0f });
        simulation.m_body = &body;

        float constexpr dt = 0.016f;
        for (int frame = 0; frame < 80; ++frame) {
            simulation.SimulateTimestep(dt);
            if (! simulation.pressureSolveSucceeded) {
                throw std::runtime_error(
                    "coupled scenario pressure solve failed at frame "
                    + std::to_string(frame));
            }

            glm::vec3 const force =
                glm::vec3(0.0f, -9.81f, 0.0f) * body.mass
                + simulation.m_feedbackForce;
            body.velocity += force / body.mass * dt;
            if (body.position.y < 0.1f) {
                body.velocity *= 0.96f;
                body.angularVelocity *= 0.95f;
            }
            body.position += body.velocity * dt;
            body.angularVelocity +=
                body.GetInertiaWorldInv() * simulation.m_feedbackTorque * dt;
            if (glm::length(body.angularVelocity) > 0.001f) {
                glm::vec3 const axis  = glm::normalize(body.angularVelocity);
                float const     angle = glm::length(body.angularVelocity) * dt;
                body.orientation =
                    glm::normalize(glm::angleAxis(angle, axis) * body.orientation);
            }

            glm::vec3 const halfDimension = body.dim * 0.5f;
            for (int direction = 0; direction < 3; ++direction) {
                if (body.position[direction] - halfDimension[direction] < -0.5f) {
                    body.position[direction] = -0.5f + halfDimension[direction];
                    body.velocity[direction] *= -0.5f;
                }
                if (body.position[direction] + halfDimension[direction] > 0.5f) {
                    body.position[direction] = 0.5f - halfDimension[direction];
                    body.velocity[direction] *= -0.5f;
                }
            }
        }

        require(std::isfinite(simulation.pressureResidual), "coupled pressure residual is not finite");
        require(isFinite(body.position), "coupled scenario produced a non-finite body position");
        require(isFinite(body.velocity), "coupled scenario produced a non-finite body velocity");
        require(isFinite(simulation.m_feedbackForce), "coupled scenario produced a non-finite feedback force");
        require(
            std::all_of(
                simulation.m_particlePos.begin(),
                simulation.m_particlePos.end(),
                isFinite),
            "coupled scenario produced a non-finite particle position");
    }
} // namespace

int main(int argc, char ** argv) {
    bool const quick =
        argc > 1 && std::string_view(argv[1]) == "--quick";

    try {
        testStandardWallProjection();
        testOutwardWallSeparation();
        testInwardWallContact();
        testDisconnectedClosedPressureComponent();
        testCoupledStandardWallProjection();
        testCoupledOutwardWallSeparation();
        testCoupledMovingBodyContactAndSeparation();
        if (! quick) {
            testSeparationScenarioRemainsFiniteAndNonpenetrating();
            testSubgridChannelPreservesHalfCellFlow();
            testCoupledScenarioRemainsFinite();
        }
    } catch (std::exception const & error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }

    std::cout << "All simulator tests passed.\n";
    return 0;
}
