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

    void draw(std::string const &ou_fil);

    // Indices that sort the atoms along the 3 axes.
    std::vector<int> const _idx_x, _idx_y, _idx_z;

    // number of atoms.
    int const _natoms;

    // Origin coordinates.
    Vector _orig_vtor{0.0f, 0.0f, 0.0f};

    // Atoms coordinates. Using SoA.
    std::vector<grid_t> _x, _y, _z;

    // Atoms VdW radii
    std::vector<float> _radii;

    // Min and max coordinates for each axis. Can be used to get main cube's
    // vertices coordinates.
    std::array<grid_t, 3> _start{std::numeric_limits<grid_t>::max(),
        std::numeric_limits<grid_t>::max(), std::numeric_limits<grid_t>::max()},
        _end{std::numeric_limits<grid_t>::min(),
            std::numeric_limits<grid_t>::min(),
            std::numeric_limits<grid_t>::min()};

    // Cube's center coordinates.
    std::array<grid_t, 3> _center;

    // _size holds the cube's size and is equal to the max element of _sizes,
    // which holds the sizes in xyz coordinates.
    grid_t _size;
    std::array<grid_t, 3> _sizes;

    // Resolution used to build the GridMolecule.
    float const _resolution = 1.0;
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

auto fill_grid_tetrahedron(GridMolecule const &in_mol)
    -> std::vector<std::vector<std::vector<grid_t>>>;

} // namespace GrANA

#endif // GrANA_GRID