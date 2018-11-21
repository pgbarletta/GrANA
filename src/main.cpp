float resolution = .1;
#include <cmath>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

#include "GrANA/continuous.hpp"
#include "GrANA/grid.hpp"
#include "GrANA/utils.hpp"

int main(int argc, char **argv) {

    auto [in_pdb, out_pdb, resolution] = GrANA::get_input(argc, argv);

    std::vector<int> indices = {
        300, 600, 900, 1200, 1500, 1800, 1240, 400, 500, 700, 800, 1000, 1100};

    GrANA::Molecule prote(argv[2]);

    GrANA::Triangulation incl_area(prote, indices);
    incl_area.draw("aux/ia.pdb");

    GrANA::GridMolecule Gprote(prote);
    Gprote.draw("aux/g1mtn.pdb");

    GrANA::fill_grid_tetrahedron(Gprote);

    return 0;
}
