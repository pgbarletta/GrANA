#ifndef GrANA_GRID_PRIMITIVES_H
#define GrANA_GRID_PRIMITIVES_H

#include "GrANA/continuous.hpp"
#include "GrANA/utils.hpp"
#include <array>
#include <fmt/format.h>
#include <fstream>
#include <iostream>
#include <string>

namespace GrANA {
using grid_t = int;

// Turn a coordinate in the matching matrix grid index.
inline grid_t cont_to_grid(float const x, float const resolution) {
    return static_cast<grid_t>(fabs(x - fmod(x, resolution)) / resolution);
}

// Turn a grid index into a xyz coordinate.
inline float grid_to_cont(int const idx, float const resolution) {
    return static_cast<float>(idx * resolution);
}

class GridPoint {
public:
    GridPoint() = default;

    GridPoint(const grid_t x, const grid_t y, const grid_t z) noexcept :
        _xyz{x, y, z} {}

    GridPoint(const Point &p, float const resolution) noexcept :
        _xyz{cont_to_grid(p[0], resolution), cont_to_grid(p[1], resolution),
            cont_to_grid(p[2], resolution)} {}

    // Draw GridPoint as atom.
    void draw(FILE *ou_fil, int idx, int resid, Vector const &orig_vtor,
        float const resolution) const {
        float const fx = grid_to_cont(_xyz[0], resolution) + orig_vtor[0];
        float const fy = grid_to_cont(_xyz[1], resolution) + orig_vtor[1];
        float const fz = grid_to_cont(_xyz[2], resolution) + orig_vtor[2];
        fmt::print(ou_fil,
            "{: <6}{: >5} {: <4s} {:3} {:1}{: >4}    "
            "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {: >2s}\n",
            "HETATM", idx, "H", "GPU", "A", resid, fx, fy, fz, 1.0, 0.0, "H");
        return;
    }

    grid_t &operator[](int idx) { return _xyz[idx]; }
    grid_t operator[](int idx) const { return _xyz[idx]; }

private:
    std::array<grid_t, 3> _xyz;
};

inline std::ostream &operator<<(std::ostream &stream, const GridPoint &p) {
    stream << p[0] << " " << p[1] << " " << p[2];
    return stream;
}

// Turn a grid point into a continuous point, given the resolution.
// This is not a constructor in the `Point` class to avoid circular
// dependencies.
inline Point GridPoint_to_point(GridPoint const &in_point, float const resolution) {
    return Point(grid_to_cont(in_point[0], resolution),
        grid_to_cont(in_point[1], resolution),
        grid_to_cont(in_point[2], resolution));
}

} // namespace GrANA

#endif // GrANA_GRID
