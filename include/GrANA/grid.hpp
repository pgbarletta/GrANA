#ifndef GrANA_GRID
#define GrANA_GRID

#include "GrANA/grid_base.hpp"
#include "utils.hpp"
#include <algorithm>
#include <array>
#include <limits>
#include <tuple>

namespace GrANA {

class GridMolecule {
public:
    GridMolecule() = default;
    GridMolecule(Molecule const &in_mol, float const resolution);

    void draw(std::string const &ou_fil, float const resolution);

    // Indices that sort the atoms along the 3 axes.
    std::vector<int> const _idx_x, _idx_y, _idx_z;
    int const _natoms;
    Vector _orig_vtor{0.0f, 0.0f, 0.0f};
    std::vector<grid_t> _x, _y, _z;
    std::vector<float> _radii;
    grid_t _x_min = std::numeric_limits<int>::max(),
           _y_min = std::numeric_limits<int>::max(),
           _z_min = std::numeric_limits<int>::max(),
           _x_max = std::numeric_limits<int>::min(),
           _y_max = std::numeric_limits<int>::min(),
           _z_max = std::numeric_limits<int>::min();
    std::array<grid_t, 3> _center;
    grid_t _size_x, _size_y, _size_z, _size;
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