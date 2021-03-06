#include "GrANA/continuous_base.hpp"

namespace GrANA {

Cube::Cube(float const p0x, float const p0y, float const p0z, float const dim) {
    _dim = dim;
    _p[0] = Point(p0x, p0y, p0z);
    _p[1] = Point(p0x, p0y, p0z + _dim);
    _p[2] = Point(p0x, p0y + _dim, p0z);
    _p[3] = Point(p0x, p0y + _dim, p0z + _dim);
    _p[4] = Point(p0x + _dim, p0y, p0z);
    _p[5] = Point(p0x + _dim, p0y, p0z + _dim);
    _p[6] = Point(p0x + _dim, p0y + _dim, p0z);
    _p[7] = Point(p0x + _dim, p0y + _dim, p0z + _dim);

    return;
}

// From GrANA::Point.
Cube::Cube(Point const p0, float const dim) {
    _dim = dim;
    _p[0] = p0;
    _p[1] = p0 + Vector(0.f, 0.f, _dim);
    _p[2] = p0 + Vector(0.f, _dim, _dim);
    _p[3] = p0 + Vector(0.f, _dim, 0.f);
    _p[4] = p0 + Vector(_dim, 0.f, 0.f);
    _p[5] = p0 + Vector(_dim, 0.f, _dim);
    _p[6] = p0 + Vector(_dim, _dim, _dim);
    _p[7] = p0 + Vector(_dim, _dim, 0.f);

    return;
}

// From GrANA::Point
Prism::Prism(Point const &p0, Point const &p1, Point const &p2, Point const &p3,
    Point const &p4, Point const &p5, Point const &p6, Point const &p7) {

    _p[0] = p0;
    _p[1] = p1;
    _p[2] = p2;
    _p[3] = p3;
    _p[4] = p4;
    _p[5] = p5;
    _p[6] = p6;
    _p[7] = p7;

    return;
}

// Using minimum and maximum coordinates.
Prism::Prism(float const x_min, float const y_min, float const z_min,
    float const x_max, float const y_max, float const z_max) {

    _p[0] = Point(x_min, y_min, z_min);
    _p[1] = Point(x_max, y_min, z_min);
    _p[2] = Point(x_min, y_max, z_min);
    _p[3] = Point(x_max, y_max, z_min);
    _p[4] = Point(x_min, y_min, z_max);
    _p[5] = Point(x_max, y_min, z_max);
    _p[6] = Point(x_min, y_max, z_max);
    _p[7] = Point(x_max, y_max, z_max);

    return;
}

}
