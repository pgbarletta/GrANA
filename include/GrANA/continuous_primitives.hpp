#ifndef GrANA_CONTINUOUS_PRIMITIVES_H
#define GrANA_CONTINUOUS_PRIMITIVES_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
using EPIC = CGAL::Exact_predicates_inexact_constructions_kernel;
using CPoint = EPIC::Point_3;

#include <fmt/format.h>

namespace GrANA {

class Vector {
public:
    Vector() = default;

    Vector(float const x, float const y, float const z) noexcept :
        _vxyz{x, y, z}, _origin{0., 0., 0.} {}

    Vector(float const x, float const y, float const z, float const ox,
        float const oy, float const oz) noexcept :
        _vxyz{x, y, z},
        _origin{ox, oy, oz} {}

    float operator[](int const idx) const { return _vxyz[idx]; }
    float &operator[](int const idx) { return _vxyz[idx]; }

    // Origin access and modification.
    float get_ox() const { return _origin[0]; }
    float get_oy() const { return _origin[1]; }
    float get_oz() const { return _origin[2]; }
    void set_ox(float ox) { _origin[0] = ox; }
    void set_oy(float oy) { _origin[1] = oy; }
    void set_oz(float oz) { _origin[2] = oz; }

private:
    std::array<float, 3> _vxyz, _origin;
};

inline std::ostream &operator<<(std::ostream &stream, Vector const &v) {
    stream << v[0] << " " << v[1] << " " << v[2];
    return stream;
}

// Returns a Vector starting on this same Vector coordinates.
inline Vector operator+(Vector const &lhs, Vector const &rhs) {
    return Vector(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2],
        lhs.get_ox(), lhs.get_oy(), lhs.get_oz());
}

// Returns a Vector starting on this same Vector coordinates.
inline Vector operator-(Vector const &lhs, Vector const &rhs) {
    return Vector(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2],
        lhs.get_ox(), lhs.get_oy(), lhs.get_oz());
}

inline bool operator==(Vector const &lhs, Vector const &rhs) {
    return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2]
        && lhs.get_ox() == rhs.get_ox() && lhs.get_oy() == rhs.get_oy()
        && lhs.get_oz() == rhs.get_oz());
}

// Get the magnitude of the Vector.
inline float norm(Vector const &v) {
    const float dx = v[0] - v.get_ox();
    const float dy = v[1] - v.get_oy();
    const float dz = v[2] - v.get_oz();
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

inline auto determinant(Vector const &v0, Vector const &v1, Vector const &v2)
    -> float {
    // First, compute the det2x2.
    float const m01 = v0[0] * v1[1] - v0[1] * v1[0];
    float const m02 = v0[0] * v2[1] - v0[1] * v2[0];
    float const m12 = v1[0] * v2[1] - v1[1] * v2[0];
    // Now compute the minors of rank 3.
    return m01 * v2[2] - m02 * v1[2] + m12 * v0[2];
}

class Point {
public:
    Point() = default;

    Point(float const x, float const y, float const z) : _xyz{x, y, z} {}

    Point(CPoint const p) :
        _xyz{static_cast<float>(CGAL::to_double(p.x())),
            static_cast<float>(CGAL::to_double(p.y())),
            static_cast<float>(CGAL::to_double(p.z()))} {}

    // Draw atom.
    void draw(FILE *out_file, int const idx, int const resid) {
        fmt::print(out_file,
            "{: <6}{: >5} {: <4s} {:3} {:1}{: >4}    "
            "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {: >2s}\n",
            "HETATM", idx, "H", "GPU", "A", resid, _xyz[0], _xyz[1], _xyz[2],
            1.0, 0.0, "H");
        return;
    }

    float operator[](int const idx) const { return _xyz[idx]; }

private:
    std::array<float, 3> _xyz;
};

inline std::ostream &operator<<(std::ostream &stream, Point const &p) {
    stream << p[0] << " " << p[1] << " " << p[2];
    return stream;
}

// Returns a Vector starting on this Point coordinates.
inline Vector operator-(Point const &lhs, Point const &rhs) {
    return Vector(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2], lhs[0],
        lhs[1], lhs[2]);
}
// Displaces the Point along the Vector.
inline Point operator+(Point const &p, Vector const &v) {
    return Point(p[0] + (v[0] - v.get_ox()), p[1] + (v[1] - v.get_oy()),
        p[2] + (v[2] - v.get_oz()));
}
// Displaces the Point along the Vector.
inline Point operator-(Point const &p, Vector const &v) {
    return Point(p[0] - (v[0] - v.get_ox()), p[1] - (v[1] - v.get_oy()),
        p[2] - (v[2] - v.get_oz()));
}

inline bool operator==(Point const &lhs, Point const &rhs) {
    return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2]);
}

inline Vector point_to_vector(Point const &in_point) {
    return Vector(in_point[0], in_point[1], in_point[2]);
}

// Get the distance between 2 points
inline float distance(Point const &p0, Point const &p1) {
    float const dx = p0[0] - p1[0];
    float const dy = p0[1] - p1[1];
    float const dz = p0[2] - p1[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}
}
#endif // _H