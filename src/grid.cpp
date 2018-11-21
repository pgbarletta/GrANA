#include "GrANA/grid.hpp"
namespace GrANA {

void GridMolecule::draw(std::string const &ou_fil) {

    FILE *file = std::fopen(ou_fil.c_str(), "w");
    if (file) {
        for (int i = 0; i <= _natoms - 1; ++i) {
            GridPoint const atm(_x[i], _y[i], _z[i]);
            atm.draw(file, i + 1, i + 1, _orig_vtor);
        }
    } else {
        std::cerr << "Could not open " << ou_fil << ". " << '\n';
    }
    std::fclose(file);
    return;
}

void fill_grid_tetrahedron(GridMolecule const &in_mol) {

    float constexpr sq = 1.4142f;

    grid_t const x0 = in_mol._x[in_mol._idx_x[0]];
    grid_t const y0 = in_mol._y[in_mol._idx_x[0]];
    grid_t const z0 = in_mol._z[in_mol._idx_x[0]];

    grid_t const x1 = in_mol._x[in_mol._idx_x[1]];
    grid_t const y1 = in_mol._y[in_mol._idx_x[1]];
    grid_t const z1 = in_mol._z[in_mol._idx_x[1]];

    grid_t const x2 = in_mol._x[in_mol._idx_x[2]];
    grid_t const y2 = in_mol._y[in_mol._idx_x[2]];
    grid_t const z2 = in_mol._z[in_mol._idx_x[2]];

    grid_t const x3 = in_mol._x[in_mol._idx_x[3]];
    grid_t const y3 = in_mol._y[in_mol._idx_x[3]];
    grid_t const z3 = in_mol._z[in_mol._idx_x[3]];

    grid_t const x01 = std::fabs(x1 - x0);
    grid_t const y01 = std::fabs(y1 - y0);
    grid_t const z01 = std::fabs(z1 - z0);

    grid_t const x02 = std::fabs(x2 - x0);
    grid_t const y02 = std::fabs(y2 - y0);
    grid_t const z02 = std::fabs(z2 - z0);

    grid_t const x03 = std::fabs(x3 - x0);
    grid_t const y03 = std::fabs(y3 - y0);
    grid_t const z03 = std::fabs(z3 - z0);

    grid_t const x12 = std::fabs(x2 - x1);
    grid_t const y12 = std::fabs(y2 - y1);
    grid_t const z12 = std::fabs(z2 - z1);

    grid_t const x13 = std::fabs(x3 - x1);
    grid_t const y13 = std::fabs(y3 - y1);
    grid_t const z13 = std::fabs(z3 - z1);

    grid_t const x23 = std::fabs(x3 - x2);
    grid_t const y23 = std::fabs(y3 - y2);
    grid_t const z23 = std::fabs(z3 - z2);

    return;
}

}
