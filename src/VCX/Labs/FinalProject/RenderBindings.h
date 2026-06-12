#pragma once

#include <array>
#include <cstdint>

namespace VCX::Labs::FluidSimulation {
    inline constexpr std::uint32_t CoupledPassConstantsBinding    = 1;
    inline constexpr std::uint32_t SubgridPassConstantsBinding    = 2;
    inline constexpr std::uint32_t SeparationPassConstantsBinding = 3;
    inline constexpr std::uint32_t VariationPassConstantsBinding  = 4;

    constexpr bool PassConstantsBindingsAreUnique() {
        constexpr std::array Bindings {
            CoupledPassConstantsBinding,
            SubgridPassConstantsBinding,
            SeparationPassConstantsBinding,
            VariationPassConstantsBinding,
        };
        for (std::size_t left = 0; left < Bindings.size(); ++left) {
            for (std::size_t right = left + 1; right < Bindings.size(); ++right) {
                if (Bindings[left] == Bindings[right])
                    return false;
            }
        }
        return true;
    }

    static_assert(PassConstantsBindingsAreUnique());
} // namespace VCX::Labs::FluidSimulation
