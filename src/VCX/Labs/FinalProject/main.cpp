#include "Assets/bundled.h"
#include "Labs/FinalProject/App.h"

int main() {
    using namespace VCX;
    return Engine::RunApp<Labs::FluidSimulation::App>(Engine::AppContextOptions {
        .Title         = "VCX-sim Lab4: Final Project - Coupled Fluid-Rigid Simulation",
        .WindowSize    = {1024, 768},
        .FontSize      = 16,
        .IconFileNames = Assets::DefaultIcons,
        .FontFileNames = Assets::DefaultFonts,
    });
}
