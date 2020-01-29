#include "GrANA/continuous_base.hpp"

namespace GrANA {
// Draw triangle.
void Triangle::draw(FILE *out_file, int const start_idx, int const resid) {
    const auto i = start_idx;
    const auto j = start_idx + 1;
    const auto k = start_idx + 2;

    _p[0].draw(out_file, i, resid);
    _p[1].draw(out_file, j, resid);
    _p[2].draw(out_file, k, resid);

    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", i, j, k);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", j, k, i);
    return;
}

// Draw tetrahedron.
void Tetrahedron::draw(FILE *out_file, int const start_idx, int const resid) {
    const auto i = start_idx;
    const auto j = start_idx + 1;
    const auto k = start_idx + 2;
    const auto l = start_idx + 3;

    _p[0].draw(out_file, i, resid);
    _p[1].draw(out_file, j, resid);
    _p[2].draw(out_file, k, resid);
    _p[3].draw(out_file, l, resid);

    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", i, j, k, l);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", j, k, l, i);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", k, l, i, j);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", l, i, j, k);
    return;
}

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

// Draw cube.
void Cube::draw(FILE *out_file, int const start_idx, int const resid) {
    const auto i = start_idx;
    const auto j = start_idx + 1;
    const auto k = start_idx + 2;
    const auto l = start_idx + 3;
    const auto ii = start_idx + 4;
    const auto jj = start_idx + 5;
    const auto kk = start_idx + 6;
    const auto ll = start_idx + 7;

    _p[0].draw(out_file, i, resid);
    _p[1].draw(out_file, j, resid);
    _p[2].draw(out_file, k, resid);
    _p[3].draw(out_file, l, resid);
    _p[4].draw(out_file, ii, resid);
    _p[5].draw(out_file, jj, resid);
    _p[6].draw(out_file, kk, resid);
    _p[7].draw(out_file, ll, resid);

    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", i, j, l, ii);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", k, j, l, kk);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", jj, ii, kk, j);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", ll, ii, kk, l);

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

// Draw prism. Can't draw connectivity properly if the prism wasn't constructed
// with proper Point ordering. SO this class is kind of useless.
void Prism::draw(FILE *out_file, int const start_idx, int const resid) {
    const auto i = start_idx;
    const auto j = start_idx + 1;
    const auto k = start_idx + 2;
    const auto l = start_idx + 3;
    const auto ii = start_idx + 4;
    const auto jj = start_idx + 5;
    const auto kk = start_idx + 6;
    const auto ll = start_idx + 7;

    _p[0].draw(out_file, i, resid);
    _p[1].draw(out_file, j, resid);
    _p[2].draw(out_file, k, resid);
    _p[3].draw(out_file, l, resid);
    _p[4].draw(out_file, ii, resid);
    _p[5].draw(out_file, jj, resid);
    _p[6].draw(out_file, kk, resid);
    _p[7].draw(out_file, ll, resid);

    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", i, j, l, ii);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", k, j, l, kk);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", jj, ii, kk, j);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", ll, ii, kk, l);

    return;
}
}
