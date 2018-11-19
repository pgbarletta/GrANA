#ifndef GrANA_GRID
#define GrANA_GRID

#include "GrANA/continuous.hpp"
#include "GrANA/utils.hpp"
#include <array>
#include <fmt/format.h>
#include <fstream>
#include <iostream>
#include <string>

namespace GrANA {

class GridPoint {
public:
    using grid_idx_t = int;

    GridPoint() = default;

    GridPoint(
        const grid_idx_t x, const grid_idx_t y, const grid_idx_t z) noexcept :
        _xyz{x, y, z} {}

    GridPoint(const Point &p) noexcept :
        _xyz{cont_to_grid(p[0]), cont_to_grid(p[1]), cont_to_grid(p[2])} {}

    // Draw GridPoint as atom.
    void draw(FILE *ou_fil, int idx, int resid, Vector const &orig_vtor);

    grid_idx_t &operator[](int idx) { return _xyz[idx]; }
    grid_idx_t operator[](int idx) const { return _xyz[idx]; }

private:
    std::array<grid_idx_t, 3> _xyz;
};

std::ostream &operator<<(std::ostream &stream, const GridPoint &t);

// Turn a grid point into a continuous point, given the resolution.
Point GridPoint_to_point(const GridPoint &in_point);

// Turn a grid point into a continuous point, given the resolution.
GridPoint point_to_GridPoint(const Point &in_point);

class GridMolecule {
public:
    // GridMolecule() noexcept : _natoms(0), _orig_vtor(Vector(0.0f, 0.0f,
    // 0.0f)),
    //     _xyz(nullptr), _in_xyz(nullptr), _radii(nullptr), _in_radii(nullptr)
    //     {};
    GridMolecule() = default;
    GridMolecule(Molecule const &in_mol, Point const &orig_point);

    void draw(std::string const &ou_fil);

    int _natoms{0};
    Vector _orig_vtor{0.0f, 0.0f, 0.0f};
    std::vector<GridPoint> _xyz, _in_xyz;
    std::vector<float> _radii, _in_radii;
};

// class GridTriangulation {
// public:
//     GridTriangulation() = default;

//     GridTriangulation(Triangulation const &T, Point const &orig_point);

//     void draw(const std::string &out_file);

//     int _ntetrahedrons;
//     std::vector<GridTetrahedron> _tetrahedrons;
//     std::vector<GridCube> _bboxes;
// };

class GridConvexHull {
public:
    GridConvexHull() = default;
    GridConvexHull(ConvexHull const &CH);

    void draw(std::string const &ou_fil);

    std::vector<float> _dots;
};

class GridBool {
public:
    GridBool() = default;
    GridBool(ConvexHull const &CH);
};
} // namespace GrANA

#endif // GrANA_GRID
