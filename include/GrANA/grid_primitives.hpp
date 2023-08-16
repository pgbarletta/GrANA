#ifndef GrANA_GRID_PRIMITIVES_H
#define GrANA_GRID_PRIMITIVES_H

#include "GrANA/continuous_primitives.hpp"
#include "GrANA/utils.hpp"
#include <array>
#include <cstdint>
#include <fmt/format.h>
#include <fstream>
#include <iostream>
#include <string>

namespace GrANA {
using grid_t = uint32_t;

// Turn a coordinate in the matching matrix grid index.
inline grid_t cont_to_grid(float const x, float const resolution) {
    return static_cast<grid_t>(fabs(x - fmod(x, resolution)) / resolution);
}

// Turn a grid index into a xyz coordinate.
inline float grid_to_cont(grid_t const idx, float const resolution) {
    return static_cast<float>(idx * resolution);
}

class GridPoint {
public:
    GridPoint() = default;

    GridPoint(grid_t const x, grid_t const y, grid_t const z) noexcept :
        _xyz {x, y, z} { }

    GridPoint(Point const &p, float const resolution) noexcept :
        _xyz {cont_to_grid(p[0], resolution), cont_to_grid(p[1], resolution),
            cont_to_grid(p[2], resolution)} { }

    GridPoint(
        Point const &p, float const resolution, Point const &offset) noexcept :
        _xyz {cont_to_grid(p[0] - offset[0], resolution),
            cont_to_grid(p[1] - offset[1], resolution),
            cont_to_grid(p[2] - offset[2], resolution)} { }

    grid_t &operator[](int idx) { return _xyz[idx]; }
    grid_t operator[](int idx) const { return _xyz[idx]; }

    std::array<grid_t, 3> _xyz;
};

inline std::ostream &operator<<(std::ostream &stream, const GridPoint &p) {
    stream << p[0] << " " << p[1] << " " << p[2];
    return stream;
}

// Turn a grid point into a continuous point, given the resolution.
// This is not a constructor in the `Point` class to avoid circular
// dependencies.
[[nodiscard]] inline Point GridPoint_to_point(
    GridPoint const &in_point, float const resolution) {
    return Point(grid_to_cont(in_point[0], resolution),
        grid_to_cont(in_point[1], resolution),
        grid_to_cont(in_point[2], resolution));
}

// Turn a grid point into a continuous point, given the resolution and an
// origin.
[[nodiscard]] inline Point GridPoint_to_point(
    GridPoint const &in_point, float const resolution, Point const &origin) {
    return Point(grid_to_cont(in_point[0], resolution) + origin[0],
        grid_to_cont(in_point[1], resolution) + origin[1],
        grid_to_cont(in_point[2], resolution) + origin[2]);
}

} // namespace GrANA

#endif // GrANA_GRID
