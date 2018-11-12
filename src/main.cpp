float rsltion = .1;
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

#include "GrANA/continuous.hpp"
#include "GrANA/grid.hpp"
#include "GrANA/utils.hpp"

int main(int argc, char **argv) {
    // Get positions and associated variables ready.
    if (argc != 4) {
        std::cerr << "Usage: GrANA resolution in_pdb out_pdb" << '\n';
        return 0;
    }
    try {
        rsltion = std::stof(argv[1]);
    } catch (...) {
        std::cerr
            << "Bad input resolution. Please specify a number between .01 and 1"
            << '\n';
    }
    std::vector<int> indices = {
        300, 600, 900, 1200, 1500, 1800, 1240, 400, 500, 700, 800, 1000, 1100};

    GrANA::Molecule prote(argv[2]);

    GrANA::Triangulation incl_area(prote, indices);

    incl_area.draw("aux/ia.pdb");

    return 0;
}
