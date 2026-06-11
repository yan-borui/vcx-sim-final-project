#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>

#include "Labs/FinalProject/FreeSurfaceSeparationSimulator.h"
#include "Labs/FinalProject/SubgridSimulator.h"
#include "Labs/FinalProject/VariationalCoupledSimulator.h"

namespace {
    using VCX::Labs::Final::FreeSurfaceSeparationSimulator;
    using VCX::Labs::Final::RigidBody;
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

    struct CoupledScenarioMetrics {
        float MaxBodySpeed       = 0.0f;
        float MaxAngularSpeed    = 0.0f;
        float MaxParticleSpeed   = 0.0f;
        float MaxFeedbackImpulse = 0.0f;
    };

    CoupledScenarioMetrics runDefaultCoupledScenario(int frameCount) {
        VariationalCoupledSimulator simulation;
        simulation.setupScene(16);

        RigidBody body;
        body.Reset(
            { 0.0f, 0.3f, 0.0f },
            { 0.0f, 0.0f, 0.0f },
            { 0.3f, 0.3f, 0.3f },
            0.3f,
            { 0.8f, 0.2f, 0.2f });
        simulation.m_body = &body;

        CoupledScenarioMetrics metrics;
        float constexpr       dt = 0.016f;
        for (int frame = 0; frame < frameCount; ++frame) {
            body.position += body.velocity * dt;
            if (glm::length(body.angularVelocity) > 0.001f) {
                glm::vec3 const axis  = glm::normalize(body.angularVelocity);
                float const     angle = glm::length(body.angularVelocity) * dt;
                body.orientation =
                    glm::normalize(glm::angleAxis(angle, axis) * body.orientation);
            }

            float const boundRadius = body.BoundingRadius();
            for (int direction = 0; direction < 3; ++direction) {
                if (body.position[direction] - boundRadius < -0.5f) {
                    body.position[direction] = -0.5f + boundRadius;
                    body.velocity[direction] *= -0.5f;
                }
                if (body.position[direction] + boundRadius > 0.5f) {
                    body.position[direction] = 0.5f - boundRadius;
                    body.velocity[direction] *= -0.5f;
                }
            }

            if (body.position.y < 0.1f) {
                body.velocity *= 0.96f;
                body.angularVelocity *= 0.95f;
            }

            body.velocity += glm::vec3(0.0f, -9.81f, 0.0f) * dt;

            simulation.SimulateTimestep(dt);
            require(simulation.pressureSolveSucceeded, "default coupled pressure solve failed");

            body.velocity += simulation.m_feedbackForce / body.mass * dt;
            body.angularVelocity +=
                body.GetInertiaWorldInv() * simulation.m_feedbackTorque * dt;

            metrics.MaxBodySpeed =
                std::max(metrics.MaxBodySpeed, glm::length(body.velocity));
            metrics.MaxAngularSpeed =
                std::max(metrics.MaxAngularSpeed, glm::length(body.angularVelocity));
            metrics.MaxFeedbackImpulse = std::max(
                metrics.MaxFeedbackImpulse,
                glm::length(simulation.m_feedbackForce) * dt);
            for (glm::vec3 const velocity : simulation.m_particleVel)
                metrics.MaxParticleSpeed =
                    std::max(metrics.MaxParticleSpeed, glm::length(velocity));

            require(isFinite(body.position), "default coupled body position is not finite");
            require(isFinite(body.velocity), "default coupled body velocity is not finite");
            require(
                std::all_of(
                    simulation.m_particlePos.begin(),
                    simulation.m_particlePos.end(),
                    isFinite),
                "default coupled particle position is not finite");
        }
        return metrics;
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
            "standard wall projection must remain divergence-free");
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
        require(
            simulation.minimumContactPressure >= -1e-5f,
            "contacting wall face must not keep negative pressure");
        require(simulation.separatingWallFaces == 0, "contacting wall face must not separate");
    }

    void testNegativeWallPressureSeparates() {
        FreeSurfaceSeparationSimulator simulation;
        prepareSingleWallCell(simulation, 0.0f);
        simulation.enableWallSeparation = true;

        glm::ivec3 const fluidCell { 1, 3, 3 };
        simulation.m_vel[offset(simulation, fluidCell + glm::ivec3(1, 0, 0))].x = 1.0f;
        simulation.solveIncompressibility(200, 0.01f, 1.0f, false);

        require(
            simulation.m_vel[offset(simulation, fluidCell)].x > 1e-4f,
            "negative wall pressure must release suction instead of sticking");
        require(
            simulation.separatingWallFaces >= 1,
            "negative wall pressure must classify at least one wall face as separating");
    }

    void testGhostFluidUsesInterfaceDistance() {
        FreeSurfaceSeparationSimulator simulation;
        simulation.setupScene(8);
        simulation.enableWallSeparation = false;
        std::fill(simulation.m_type.begin(), simulation.m_type.end(), Simulator::EMPTY_CELL);
        std::fill(simulation.m_vel.begin(), simulation.m_vel.end(), glm::vec3(0.0f));

        glm::ivec3 const fluidCell { 3, 3, 3 };
        simulation.m_type[offset(simulation, fluidCell)] = Simulator::FLUID_CELL;
        glm::ivec3 const solidNeighbors[] = {
            { 2, 3, 3 },
            { 3, 2, 3 },
            { 3, 4, 3 },
            { 3, 3, 2 },
            { 3, 3, 4 },
        };
        for (glm::ivec3 const neighbor : solidNeighbors)
            simulation.m_s[offset(simulation, neighbor)] = 0.0f;

        glm::vec3 const fluidCenter =
            (glm::vec3(fluidCell) + glm::vec3(0.5f)) * simulation.m_h
            - glm::vec3(0.5f);
        simulation.m_iNumSpheres = 1;
        simulation.m_particlePos.assign(1, fluidCenter);
        simulation.m_particleVel.assign(1, glm::vec3(0.0f));
        simulation.m_vel[offset(simulation, fluidCell + glm::ivec3(1, 0, 0))].x = 1.0f;

        simulation.solveIncompressibility(200, 0.01f, 1.0f, false);

        require(
            std::abs(simulation.m_p[offset(simulation, fluidCell)] + 0.75f) < 1e-4f,
            "ghost-fluid pressure did not use the reconstructed interface distance");
    }

    void testWallSeparationActiveSetConverges() {
        FreeSurfaceSeparationSimulator simulation;
        simulation.setupScene(24);
        simulation.enableWallSeparation = true;

        for (int frame = 0; frame < 20; ++frame) {
            simulation.SimulateTimestep(0.008f);
            require(
                simulation.wallSeparationActiveSetIterations < 128,
                "wall separation active set hit the iteration cap");
            require(
                simulation.minimumContactPressure >= -1e-4f,
                "wall separation kept negative contact pressure");
            require(
                std::all_of(
                    simulation.m_particlePos.begin(),
                    simulation.m_particlePos.end(),
                    isFinite),
                "wall separation produced a non-finite particle position");
            require(std::isfinite(simulation.pressureResidual), "wall separation solve failed");
        }
    }

    float runSubgridChannel(bool useSubgridWeights) {
        int constexpr resolution = 16;
        float const h            = 1.0f / float(resolution);

        RigidBody body;
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
        require(
            simulation.partiallyOpenFaceCount > 0,
            "sub-grid scenario did not detect cut MAC volumes");
        require(simulation.pressureSolveSucceeded, "sub-grid pressure solve failed");
        require(std::isfinite(simulation.pressureResidual), "sub-grid residual is not finite");

        float maxX = -std::numeric_limits<float>::infinity();
        for (glm::vec3 const position : simulation.m_particlePos)
            maxX = std::max(maxX, position.x);
        return maxX;
    }

    void testSubgridChannelPreservesHalfCellFlow() {
        float const weightedMaxX = runSubgridChannel(true);
        float const binaryMaxX   = runSubgridChannel(false);

        require(weightedMaxX > 0.12f, "sub-grid weights did not carry fluid through the gap");
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

    void testCoupledWallSeparation() {
        VariationalCoupledSimulator outward;
        prepareSingleCoupledWallCell(outward, 1.0f);
        outward.enableWallSeparation = true;
        outward.solveIncompressibility(200, 0.01f, 1.0f, false);
        require(
            outward.m_vel[offset(outward, { 1, 3, 3 })].x > 1e-4f,
            "coupled solver must allow outward wall separation");

        VariationalCoupledSimulator inward;
        prepareSingleCoupledWallCell(inward, -1.0f);
        inward.enableWallSeparation = true;
        inward.solveIncompressibility(200, 0.01f, 1.0f, false);
        require(
            std::abs(inward.m_vel[offset(inward, { 1, 3, 3 })].x) < 1e-5f,
            "coupled solver must prevent wall penetration");
    }

    void testCoupledNegativeWallPressureSeparates() {
        VariationalCoupledSimulator simulation;
        prepareSingleCoupledWallCell(simulation, 0.0f);
        simulation.enableWallSeparation = true;

        glm::ivec3 const fluidCell { 1, 3, 3 };
        simulation.m_vel[offset(simulation, fluidCell + glm::ivec3(1, 0, 0))].x = 1.0f;
        simulation.solveIncompressibility(200, 0.01f, 1.0f, false);

        require(
            simulation.m_vel[offset(simulation, fluidCell)].x > 1e-4f,
            "coupled solver must release negative-pressure wall suction");
    }

    void prepareSingleBodyWallCell(
        VariationalCoupledSimulator & simulation,
        RigidBody &                   body,
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
        glm::ivec3 const bodyFace { 4, 3, 3 };

        VariationalCoupledSimulator contactSimulation;
        RigidBody                   contactBody;
        prepareSingleBodyWallCell(contactSimulation, contactBody, -0.5f);
        contactSimulation.enableWallSeparation = true;
        contactSimulation.solveIncompressibility(200, 0.01f, 1.0f, false);
        require(
            contactSimulation.m_feedbackForce.x > 0.0f,
            "contact pressure did not oppose a body moving into fluid");

        VariationalCoupledSimulator separationSimulation;
        RigidBody                   separationBody;
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
    }

    void testCoupledProjectionUsesUpdatedBodyVelocity() {
        VariationalCoupledSimulator simulation;
        RigidBody                   body;
        prepareSingleBodyWallCell(simulation, body, -0.5f);
        simulation.enableWallSeparation = false;

        float constexpr dt = 0.01f;
        simulation.solveIncompressibility(200, dt, 1.0f, false);

        RigidBody projectedBody = body;
        projectedBody.velocity +=
            simulation.m_feedbackForce / projectedBody.mass * dt;
        projectedBody.angularVelocity +=
            projectedBody.GetInertiaWorldInv() * simulation.m_feedbackTorque * dt;

        glm::ivec3 const bodyFace { 4, 3, 3 };
        glm::vec3 const faceCenter =
            (glm::vec3(bodyFace) + glm::vec3(0.0f, 0.5f, 0.5f))
                * simulation.m_h
            - glm::vec3(0.5f);
        float const fluidVelocity =
            simulation.m_vel[offset(simulation, bodyFace)].x;
        float const projectedBodyVelocity =
            projectedBody.GetVelocityAtPoint(faceCenter - projectedBody.position).x;
        require(
            std::abs(fluidVelocity - projectedBodyVelocity) < 1e-4f,
            "coupled projection used the pre-pressure body velocity");
    }

    void testTankPressureDoesNotPushRigidBody() {
        VariationalCoupledSimulator simulation;
        simulation.setupScene(8);
        simulation.enableWallSeparation = false;

        RigidBody body;
        body.Reset(
            { 0.3f, 0.3f, 0.3f },
            { 0.0f, 0.0f, 0.0f },
            { 0.1f, 0.1f, 0.1f },
            0.3f,
            { 0.8f, 0.2f, 0.2f });
        simulation.m_body = &body;

        std::fill(simulation.m_type.begin(), simulation.m_type.end(), Simulator::EMPTY_CELL);
        std::fill(simulation.m_vel.begin(), simulation.m_vel.end(), glm::vec3(0.0f));
        glm::ivec3 const fluidCell { 1, 3, 3 };
        simulation.m_type[offset(simulation, fluidCell)] = Simulator::FLUID_CELL;
        simulation.m_vel[offset(simulation, fluidCell + glm::ivec3(1, 0, 0))].x = -1.0f;

        simulation.solveIncompressibility(200, 0.01f, 1.0f, false);

        require(
            glm::length(simulation.m_feedbackForce) < 1e-6f,
            "tank wall pressure incorrectly pushed a remote rigid body");
        require(
            glm::length(simulation.m_feedbackTorque) < 1e-6f,
            "tank wall pressure incorrectly torqued a remote rigid body");
    }

    void testDefaultCoupledScenarioRemainsStable() {
        CoupledScenarioMetrics const metrics = runDefaultCoupledScenario(240);
        require(
            metrics.MaxBodySpeed < 8.0f,
            "default coupled body gained excessive linear speed");
        require(
            metrics.MaxAngularSpeed < 120.0f,
            "default coupled body gained excessive angular speed");
        require(
            metrics.MaxParticleSpeed < 30.0f,
            "default coupled particles gained excessive speed");
        require(
            metrics.MaxFeedbackImpulse < 0.5f,
            "default coupled pressure feedback produced an excessive impulse");
    }
} // namespace

int main(int argc, char ** argv) {
    bool const quick =
        argc > 1 && std::string_view(argv[1]) == "--quick";

    try {
        testStandardWallProjection();
        testOutwardWallSeparation();
        testInwardWallContact();
        testNegativeWallPressureSeparates();
        testGhostFluidUsesInterfaceDistance();
        testWallSeparationActiveSetConverges();
        testCoupledWallSeparation();
        testCoupledNegativeWallPressureSeparates();
        testCoupledMovingBodyContactAndSeparation();
        testCoupledProjectionUsesUpdatedBodyVelocity();
        testTankPressureDoesNotPushRigidBody();
        if (! quick) {
            testSubgridChannelPreservesHalfCellFlow();
            testDefaultCoupledScenarioRemainsStable();
        }
    } catch (std::exception const & error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }

    std::cout << "All simulator tests passed.\n";
    return 0;
}
