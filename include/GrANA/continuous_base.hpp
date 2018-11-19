#ifndef GrANA_CONTINUOUS_BASE_H
#define GrANA_CONTINUOUS_BASE_H

#include <GrANA/continuous_primitives.hpp>
#include <array>
#include <fmt/format.h>
namespace GrANA {

class Triangle {
public:
    Triangle() = default;

    Triangle(float const x0, float const y0, float const z0, float const x1,
        float const y1, float const z1, float const x2, float const y2,
        float const z2) :
        _p({Point(x0, y0, z0), Point(x1, y1, z1), Point(x2, y2, z2)}) {}

    // From GrANA::Point
    Triangle(Point const &p0, Point const &p1, Point const &p2) :
        _p({p0, p1, p2}) {}

    Point operator[](int const idx) const { return _p[idx]; }

    // Draw triangle.
    void draw(FILE *out_file, int const start_idx, int const resid);

    std::array<Point, 3> _p;
};

inline std::ostream &operator<<(std::ostream &stream, const Triangle &t) {
    stream << t._p[0] << "\t" << t._p[1] << "\t" << t._p[2];
    return stream;
}

class Tetrahedron {
public:
    Tetrahedron() = default;
    Tetrahedron(float const x0, float const y0, float const z0, float const x1,
        float const y1, float const z1, float const x2, float const y2,
        float const z2, float const x3, float const y3, float const z3) :
        _p({Point(x0, y0, z0), Point(x1, y1, z1), Point(x2, y2, z2),
            Point(x3, y3, z3)}) {}

    // From GrANA::Point
    Tetrahedron(
        Point const &p0, Point const &p1, Point const &p2, Point const &p3) :
        _p({p0, p1, p2, p3}) {}

    Point operator[](int const idx) const { return _p[idx]; }

    // Draw tetrahedron.
    void draw(FILE *out_file, int const start_idx, int const resid);

    std::array<Point, 4> _p;
};

inline std::ostream &operator<<(std::ostream &stream, Tetrahedron const &t) {
    stream << t._p[0] << "\t" << t._p[1] << "\t" << t._p[2] << "\t" << t._p[3];
    return stream;
}

inline auto determinant(Tetrahedron const &t) -> float {
    Vector const v10 = t[1] - t[0];
    Vector const v20 = t[2] - t[0];
    Vector const v30 = t[3] - t[0];
    return determinant(v10, v20, v30);
}

class Cube {
public:
    Cube() = default;

    Cube(float const p0x, float const p0y, float const p0z, float const dim);

    // From GrANA::Point.
    Cube(Point const p0, float const dim);

    // Draw cube.
    void draw(FILE *out_file, int const start_idx, int const resid);

    std::array<Point, 8> _p;
    float _dim;
};
std::ostream &operator<<(std::ostream &stream, Cube const &t);

class Prism {
public:
    Prism() = default;

    // From GrANA::Point
    Prism(Point const &p0, Point const &p1, Point const &p2, Point const &p3,
        Point const &p4, Point const &p5, Point const &p6, Point const &p7);

    Point &operator[](int const idx) { return _p[idx]; }

    // Draw prism.
    void draw(FILE *out_file, int const start_idx, int const resid);

    std::array<Point, 8> _p;
};
std::ostream &operator<<(std::ostream &stream, Prism const &t);

}
#endif // _H