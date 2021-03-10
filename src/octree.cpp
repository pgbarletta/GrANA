#include <GrANA/octree.hpp>

namespace GrANA {

// Octree::Octree(
//     Molecule const &in_mol, float const resolution, grid_t const bucket_size)
//     : _idx_x(sort_indices_uint32(in_mol._x)),
//     _idx_y(sort_indices_uint32(in_mol._y)),
//     _idx_z(sort_indices_uint32(in_mol._z)), _natoms(in_mol._natoms),
//     _resolution(resolution), _bucket_size(bucket_size) {

//     _origin = Vector(std::floor(in_mol._x[_idx_x[0]]),
//         std::floor(in_mol._y[_idx_y[0]]), std::floor(in_mol._z[_idx_z[0]]));

//     _x.reserve(_natoms * 3);
//     _y.reserve(_natoms * 3);
//     _z.reserve(_natoms * 3);
//     _radii.reserve(_natoms);

//     for (int i = 0; i < _natoms; ++i) {
//         grid_t const x = cont_to_grid(in_mol._x[i] - _origin[0],
//         _resolution); grid_t const y = cont_to_grid(in_mol._y[i] -
//         _origin[1], _resolution); grid_t const z = cont_to_grid(in_mol._z[i]
//         - _origin[2], _resolution);

//         if (x < _start[0])
//             _start[0] = x;
//         if (y < _start[1])
//             _start[1] = y;
//         if (z < _start[2])
//             _start[2] = z;
//         if (x > _end[0])
//             _end[0] = x;
//         if (y > _end[1])
//             _end[1] = y;
//         if (z > _end[2])
//             _end[2] = z;

//         _x.push_back(x);
//         _y.push_back(y);
//         _z.push_back(z);
//         _radii.push_back(in_mol._radii[i]);
//     }

//     _sizes = {_end[0] - _start[0], _end[1] - _start[1], _end[2] - _start[2]};
//     _size = std::max({_sizes[0], _sizes[1], _sizes[2]});
//     // don't mind the rounding.
//     _center = {(_sizes[0]) / 2, (_sizes[1]) / 2, (_sizes[2]) / 2};

//     // _root = newOctant(molecula._start[0], molecula._start[1],
//     //     molecula._start[2], molecula._end[0], molecula._end[1],
//     //     molecula._end[2], molecula._center, _x_sorted.begin(),
//     //     _x_sorted.end(), _y_sorted.begin(), _y_sorted.end(),
//     //     _z_sorted.begin(), _z_sorted.end());

//     return;
// }

// Octree::Octant *Octree::newOctant(grid_t const x_min, grid_t const y_min,
//     grid_t const z_min, grid_t const x_max, grid_t const y_max,
//     grid_t const z_max, std::array<grid_t, 3> center,
//     std::vector<grid_t> x_sorted, std::vector<grid_t> y_sorted,
//     std::vector<grid_t> z_sorted) {

// }
}