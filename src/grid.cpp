#include "GrANA/grid.hpp"
namespace GrANA {

GridMolecule::GridMolecule(Molecule const &in_mol, float const resolution) :
    _idx_x(sort_indices(in_mol._x)), _idx_y(sort_indices(in_mol._y)),
    _idx_z(sort_indices(in_mol._z)), _natoms(in_mol._natoms),
    _resolution(resolution) {

    _orig_vtor = Vector(std::floor(in_mol._x[_idx_x[0]]),
        std::floor(in_mol._y[_idx_y[0]]), std::floor(in_mol._z[_idx_z[0]]));

    _x.reserve(_natoms * 3);
    _y.reserve(_natoms * 3);
    _z.reserve(_natoms * 3);
    _radii.reserve(_natoms);

    for (int i = 0; i < _natoms; ++i) {
        grid_t const x =
            cont_to_grid(in_mol._x[i] - _orig_vtor[0], _resolution);
        grid_t const y =
            cont_to_grid(in_mol._y[i] - _orig_vtor[1], _resolution);
        grid_t const z =
            cont_to_grid(in_mol._z[i] - _orig_vtor[2], _resolution);

        if (x < _start[0])
            _start[0] = x;
        if (y < _start[1])
            _start[1] = y;
        if (z < _start[2])
            _start[2] = z;
        if (x > _end[0])
            _end[0] = x;
        if (y > _end[1])
            _end[1] = y;
        if (z > _end[2])
            _end[2] = z;

        _x.push_back(x);
        _y.push_back(y);
        _z.push_back(z);
        _radii.push_back(in_mol._radii[i]);
    }

    _sizes = {_end[0] - _start[0], _end[1] - _start[1], _end[2] - _start[2]};
    _size = std::max({_sizes[0], _sizes[1], _sizes[2]});
    // don't mind the rounding.
    _center = {(_sizes[0]) / 2, (_sizes[1]) / 2, (_sizes[2]) / 2};
}

void GridMolecule::draw(std::string const &ou_fil) {

    FILE *file = std::fopen(ou_fil.c_str(), "w");
    if (file) {
        for (int i = 0; i <= _natoms - 1; ++i) {
            GridPoint const atm(_x[i], _y[i], _z[i]);
            atm.draw(file, i + 1, i + 1, _orig_vtor, _resolution);
        }
    } else {
        std::cerr << "Could not open " << ou_fil << ". " << '\n';
    }
    std::fclose(file);
    return;
}

auto fill_grid_tetrahedron(GridMolecule const &in_mol)
    -> std::vector<std::vector<std::vector<grid_t>>> {
    // Inverse square root
    float constexpr isqrt = 0.7071f;
    std::vector<std::vector<std::vector<grid_t>>> mtx(in_mol._end[2],
        std::vector<std::vector<grid_t>>(
            in_mol._end[1], std::vector<grid_t>(in_mol._end[0])));

    for (auto ii = 0; ii < in_mol._natoms; ++ii) {
        grid_t const x = in_mol._x[ii];
        grid_t const y = in_mol._y[ii];
        grid_t const z = in_mol._z[ii];

        // Atom VdW radius, inner square dimension and their difference.
        grid_t const radius = static_cast<grid_t>(
            std::floor(in_mol._radii[ii] / in_mol._resolution));
        grid_t const is = static_cast<grid_t>(std::floor(radius * isqrt));
        [[maybe_unused]] grid_t const leftover = radius - is;

        for (grid_t k = -is; k <= is; ++k) {
            grid_t const iz = z - k;
            for (grid_t j = -is; j <= is; ++j) {
                grid_t const iy = y - j;
                for (grid_t i = -is; i <= is; ++i) {
                    grid_t const ix = x - i;
                    if (ix >= 0 and iy >= 0 and iz >= 0 and
                        ix < in_mol._end[0] and iy < in_mol._end[1] and
                        iz < in_mol._end[2]) {
                        mtx[iz][iy][ix] = 1;
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
