#include "GrANA/grid.hpp"
namespace GrANA {

void GridMolecule::draw(std::string const &ou_fil, float const resolution) {

    FILE *file = std::fopen(ou_fil.c_str(), "w");
    if (file) {
        for (int i = 0; i <= _natoms - 1; ++i) {
            GridPoint const atm(_x[i], _y[i], _z[i]);
            atm.draw(file, i + 1, i + 1, _orig_vtor, resolution);
        }
    } else {
        std::cerr << "Could not open " << ou_fil << ". " << '\n';
    }
    std::fclose(file);
    return;
}

auto fill_grid_tetrahedron(GridMolecule const &in_mol, float const resolution)
    -> std::vector<std::vector<std::vector<grid_t>>> {
    // Inverse square root
    float constexpr isqrt = 0.7071f;
    std::vector<std::vector<std::vector<grid_t>>> mtx(in_mol._z_max,
        std::vector<std::vector<grid_t>>(
            in_mol._y_max, std::vector<grid_t>(in_mol._x_max)));

    for (auto ii = 0; ii < in_mol._natoms; ++ii) {
        grid_t const x = in_mol._x[ii];
        grid_t const y = in_mol._y[ii];
        grid_t const z = in_mol._z[ii];

        // Atom VdW radius, inner square dimension and their difference.
        grid_t const radius =
            static_cast<grid_t>(in_mol._radii[ii] / resolution);
        grid_t const is = static_cast<grid_t>(radius * isqrt);
        [[maybe_unused]] grid_t const leftover = radius - is;

        for (grid_t k = -is; k <= is; ++k) {
            grid_t const iz = z - k;
            for (grid_t j = -is; j <= is; ++j) {
                grid_t const iy = y - j;
                for (grid_t i = -is; i <= is; ++i) {
                    grid_t const ix = x - i;
                    if (ix >= 0 and iy >= 0 and iz >= 0 and
                        ix < in_mol._x_max and iy < in_mol._y_max and
                        iz < in_mol._z_max) {
                        mtx[iz][iy][ix] = 1;

                        // printf("ix: %i iy: %i iz: %i \n", ix, iy, iz);
                    }
                }
            }
        }
    }

    return mtx;
    // for (auto i = 0; i < in_mol._natoms - 3; ++i) {
    //     grid_t const x0 = in_mol._x[in_mol._idx_x[i]];
    //     grid_t const y0 = in_mol._y[in_mol._idx_x[i]];
    //     grid_t const z0 = in_mol._z[in_mol._idx_x[i]];

    //     grid_t const x1 = in_mol._x[in_mol._idx_x[i + 1]];
    //     grid_t const y1 = in_mol._y[in_mol._idx_x[i + 1]];
    //     grid_t const z1 = in_mol._z[in_mol._idx_x[i + 1]];

    //     grid_t const x2 = in_mol._x[in_mol._idx_x[i + 2]];
    //     grid_t const y2 = in_mol._y[in_mol._idx_x[i + 2]];
    //     grid_t const z2 = in_mol._z[in_mol._idx_x[i + 2]];

    //     grid_t const x3 = in_mol._x[in_mol._idx_x[i + 3]];
    //     grid_t const y3 = in_mol._y[in_mol._idx_x[i + 3]];
    //     grid_t const z3 = in_mol._z[in_mol._idx_x[i + 3]];

    //     grid_t const x01 = std::fabs(x1 - x0);
    //     grid_t const y01 = std::fabs(y1 - y0);
    //     grid_t const z01 = std::fabs(z1 - z0);

    //     grid_t const x02 = std::fabs(x2 - x0);
    //     grid_t const y02 = std::fabs(y2 - y0);
    //     grid_t const z02 = std::fabs(z2 - z0);

    //     grid_t const x03 = std::fabs(x3 - x0);
    //     grid_t const y03 = std::fabs(y3 - y0);
    //     grid_t const z03 = std::fabs(z3 - z0);

    //     grid_t const x12 = std::fabs(x2 - x1);
    //     grid_t const y12 = std::fabs(y2 - y1);
    //     grid_t const z12 = std::fabs(z2 - z1);

    //     grid_t const x13 = std::fabs(x3 - x1);
    //     grid_t const y13 = std::fabs(y3 - y1);
    //     grid_t const z13 = std::fabs(z3 - z1);

    //     grid_t const x23 = std::fabs(x3 - x2);
    //     grid_t const y23 = std::fabs(y3 - y2);
    //     grid_t const z23 = std::fabs(z3 - z2);

    //     grid_t const d01 = x01 * x01 + y01 * y01 + z01 * z01;
    //     grid_t const d02 = x02 * x02 + y02 * y02 + z02 * z02;
    //     grid_t const d03 = x03 * x03 + y03 * y03 + z03 * z03;

    //     float const mn = (4 / (resolution * resolution));

    //     if (d01 > mn || d02 > mn || d03 > mn)
    //         printf("%i+%i+%i+%i\n", in_mol._idx_x[i], in_mol._idx_x[i + 1],
    //             in_mol._idx_x[i + 2], in_mol._idx_x[i + 3]);
    // }
}

}
