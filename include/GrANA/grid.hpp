#ifndef GrANA_GRID
#define GrANA_GRID

#include "GrANA/grid_base.hpp"
#include "utils.hpp"
#include <limits>
#include <tuple>

namespace GrANA {

class GridMolecule {
public:
    GridMolecule() = default;
    GridMolecule(Molecule const &in_mol, float const resolution) :
        _idx_x(sort_indices(in_mol._x)), _idx_y(sort_indices(in_mol._y)),
        _idx_z(sort_indices(in_mol._z)), _natoms(in_mol._natoms) {

        _orig_vtor = Vector(std::floor(in_mol._x[_idx_x[0]]),
            std::floor(in_mol._y[_idx_y[0]]), std::floor(in_mol._z[_idx_z[0]]));

        _x.reserve(_natoms * 3);
        _y.reserve(_natoms * 3);
        _z.reserve(_natoms * 3);
        _radii.reserve(_natoms);

        for (int i = 0; i < _natoms; ++i) {
            auto const x =
                cont_to_grid(in_mol._x[i] - _orig_vtor[0], resolution);
            auto const y =
                cont_to_grid(in_mol._y[i] - _orig_vtor[1], resolution);
            auto const z =
                cont_to_grid(in_mol._z[i] - _orig_vtor[2], resolution);

            if (x < _x_min)
                _x_min = x;
            if (y < _y_min)
                _y_min = y;
            if (z < _z_min)
                _z_min = z;
            if (x > _x_max)
                _x_max = x;
            if (y > _y_max)
                _y_max = y;
            if (z > _z_max)
                _z_max = z;

            _x.push_back(x);
            _y.push_back(y);
            _z.push_back(z);
            _radii.push_back(in_mol._radii[i]);
        }
    }

    void draw(std::string const &ou_fil, float const resolution);

    // Indices that sort the atoms along the 3 axes.
    std::vector<int> const _idx_x, _idx_y, _idx_z;
    int const _natoms;
    Vector _orig_vtor{0.0f, 0.0f, 0.0f};
    std::vector<grid_t> _x, _y, _z;
    std::vector<float> _radii;
    grid_t _x_min = std::numeric_limits<int>::max(),
           _y_min = std::numeric_limits<int>::max(),
           _z_min = std::numeric_limits<int>::max(), _x_max = 0, _y_max = 0,
           _z_max = 0;
};

class GridConvexHull {
public:
    GridConvexHull() = default;
    GridConvexHull(ConvexHull const &CH);

    void draw(std::string const &ou_fil);

    std::vector<float> _dots;
};

class GridBool {
public:
    GridBool() = default;
    GridBool(ConvexHull const &CH);
};

auto fill_grid_tetrahedron(GridMolecule const &in_mol, float const resolution)
    -> std::vector<std::vector<std::vector<grid_t>>>;

} // namespace GrANA

#endif // GrANA_GRID