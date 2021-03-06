#ifndef GrANA_GRID_BASE_H
#define GrANA_GRID_BASE_H

#include "GrANA/grid_primitives.hpp"

namespace GrANA {

// class GridPrism {
// public:
//     Prism() = default;

//     // From GrANA::GridPoint
//     Prism(GridPoint const &p0, GridPoint const &p1, GridPoint const &p2,
//         GridPoint const &p3, GridPoint const &p4, GridPoint const &p5,
//         GridPoint const &p6, GridPoint const &p7);

//     // Using minimum and maximum coordinates.
//     Prism(float const x_min, float const y_min, float const z_min,
//         float const x_max, float const y_max, float const z_max);

//     GridPoint &operator[](int const idx) { return _p[idx]; }

//     std::array<GridPoint, 8> _p;
// };

// std::ostream &operator<<(std::ostream &stream, GridPrism const &t);

} // namespace GrANA

#endif // H
