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

    GrANA::Molecule prote(in_pdb);

    GrANA::Triangulation incl_area(prote, indices);
    incl_area.draw(out_pdb);

    GrANA::GridMolecule Gprote(prote, resolution);
    Gprote.draw("aux/g1mtn.pdb", resolution);

    printf("_x_min: %i\n", Gprote._x_min);
    printf("_y_min: %i\n", Gprote._y_min);
    printf("_z_min: %i\n", Gprote._z_min);

    printf("_x_max: %i\n", Gprote._x_max);
    printf("_y_max: %i\n", Gprote._y_max);
    printf("_z_max: %i\n", Gprote._z_max);

    auto mtx = GrANA::fill_grid_tetrahedron(Gprote, resolution);
    draw(out_pdb, mtx, Gprote._orig_vtor, resolution);

    return 0;
}
