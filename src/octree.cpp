#include <GrANA/octree.hpp>

namespace GrANA {

Octree::Octree(GridMolecule const &molecula, grid_t const bucket_size) :
    _x_sorted(reorder(molecula._x, molecula._idx_x)),
    _y_sorted(reorder(molecula._y, molecula._idx_y)),
    _z_sorted(reorder(molecula._z, molecula._idx_z)),
    _bucket_size(bucket_size) {

    _root = newOctant(molecula._start[0], molecula._start[1],
        molecula._start[2], molecula._end[0], molecula._end[1],
        molecula._end[2], molecula._center, _x_sorted.begin(), _x_sorted.end(),
        _y_sorted.begin(), _y_sorted.end(), _z_sorted.begin(), _z_sorted.end());
    1;
    return;
}

// Octree::Octant *Octree::newOctant(grid_t const x_min, grid_t const y_min,
//     grid_t const z_min, grid_t const x_max, grid_t const y_max,
//     grid_t const z_max, std::array<grid_t, 3> center,
//     std::vector<grid_t> x_sorted, std::vector<grid_t> y_sorted,
//     std::vector<grid_t> z_sorted) {

// }
}