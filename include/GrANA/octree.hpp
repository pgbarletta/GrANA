#ifndef GrANA_OCTREE_H
#define GrANA_OCTREE_H

#include "GrANA/grid.hpp"
#include "GrANA/utils.hpp"
#include <algorithm>
#include <array>

namespace GrANA {

class Octree {
public:
    Octree(Molecule const &in_mol, float const resolution,
        grid_t const bucket_size);

    class Octant {
    public:
        Octant() = default;

        std::array<grid_t, 3> _start, _end;
        int _natoms;
        Octant *_childs[8];
        bool _leaf = true;
    };

    template <typename Iterator>
    auto newOctant(grid_t const x_min, grid_t const y_min, grid_t const z_min,
        grid_t const x_max, grid_t const y_max, grid_t const z_max,
        std::array<grid_t, 3> center, Iterator x_begin, Iterator x_end,
        Iterator y_begin, Iterator y_end, Iterator z_begin, Iterator z_end,
        int const natoms) -> Octant * {

        Octant *octant = new Octant;
        octant->_start = {x_min, y_min, z_min};
        octant->_end = {x_max, y_max, z_max};
        octant->_natoms = natoms;

        auto const x_medio =
            std::lower_bound(x_sorted.begin(), x_sorted.end(), center[0]);
        auto const y_medio =
            std::lower_bound(y_sorted.begin(), y_sorted.end(), center[1]);
        auto const z_medio =
            std::lower_bound(z_sorted.begin(), z_sorted.end(), center[2]);

        return octant;
    }

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

    // Minimum octant size.
    int const _bucket_size = 5;

    // Main octant
    Octant *_root;
};
}

#endif // GrANA_OCTREE_H