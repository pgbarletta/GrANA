#ifndef GrANA_GRID
#define GrANA_GRID

#include "GrANA/continuous.hpp"
#include "GrANA/grid_base.hpp"
#include "utils.hpp"
#include <algorithm>
#include <array>
#include <limits>
#include <tuple>

namespace GrANA {

class BoundingBox {
public:
    BoundingBox() = default;

    // BoundingBox(float const xmin, float const ymin, float const zmin,
    //     float const xmax, float const ymax, float const zmax) noexcept :
    //     _xmin(xmin),
    //     _ymin(ymin), _zmin(zmin), _xmax(xmax), _ymax(ymax), _zmax(zmax) {};

    // Box Points and GridPoints
    std::array<Point, 8> _p;
    std::array<GridPoint, 8> _gp;

    // Box center coordinates.
    std::array<grid_t, 3> _center;

    // _size holds the box's size and is equal to the max element of _sizes,
    // which holds the sizes in xyz coordinates.
    grid_t _size;
    std::array<grid_t, 3> _sizes;
};

// Mainly used by GridMolecule to fill its bounding box.
void fill_bounding_box(BoundingBox &bbox, float const xmin, float const ymin,
    float const zmin, float const xmax, float const ymax, float const zmax,
    Point const &origin, float const resolution);

class GridMolecule {
public:
    GridMolecule() = default;
    GridMolecule(Molecule const &in_mol, float const resolution,
        float const bbox_margin);
    // Build a GridMolecule using a given origin.
    GridMolecule(Molecule const &in_mol, float const resolution,
        float const bbox_margin, Point const &origin);

    // Finish building GridMolecule.
    void construct_GridMolecule(Molecule const &in_mol);

    void draw(std::string const &ou_fil);

    // Indices that sort the atoms along the 3 axes.
    std::vector<int> const _idx_x, _idx_y, _idx_z;

    // number of atoms.
    int const _natoms;

    // Resolution used to build the GridMolecule.
    float const _resolution = 1.0;

    // Extra margin for the bounding box.
    float const _bbox_margin = 1.0;

    // Origin coordinates.
    Point _origin {0.0f, 0.0f, 0.0f};

    // Atoms coordinates. Using SoA.
    std::vector<grid_t> _x, _y, _z;

    // Atoms VdW radii
    std::vector<float> _radii;

    // Bouding box
    BoundingBox _bbox;
};

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

auto fill_grid_tetrahedron(GridMolecule const &in_mol)
    -> std::vector<std::vector<std::vector<grid_t>>>;

} // namespace GrANA

#endif // GrANA_GRID