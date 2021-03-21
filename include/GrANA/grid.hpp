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

    // Resolution used to build the BoundingBox.
    float _resolution = 1.0;

    // Origin coordinates.
    Point _origin {0.0f, 0.0f, 0.0f};

    // Extra margin for the BoundingBox.
    unsigned int _bbox_margin = 1;

    // Box Points and GridPoints
    std::array<Point, 8> _p;
    std::array<GridPoint, 8> _gp;

    // Box center coordinates.
    std::array<grid_t, 3> _center;

    // Box dimensions.
    grid_t _dimx, _dimy, _dimz;
};

// Mainly used by GridMolecule to fill its bounding box.
void fill_bounding_box(BoundingBox &bbox, grid_t const xmin, grid_t const ymin,
    grid_t const zmin, grid_t const xmax, grid_t const ymax, grid_t const zmax,
    float const resolution, Point const &origin, unsigned int const margin);

class GridMolecule {
public:
    GridMolecule() = default;
    GridMolecule(Molecule const &in_mol, float const resolution,
        unsigned int const bbox_margin);
    // Build a GridMolecule using a given origin.
    GridMolecule(Molecule const &in_mol, float const resolution,
        unsigned int const bbox_margin, Point const &origin);

    // Finish building GridMolecule.
    void construct_GridMolecule(Molecule const &in_mol);

    void draw(std::string const &ou_fil);

    // Indices that sort the atoms along the 3 axes.
    std::vector<uint32_t> const _idx_x, _idx_y, _idx_z;

    // number of atoms.
    int const _natoms;

    // Resolution used to build the GridMolecule.
    float _resolution = 1.0;

    // Extra margin for the bounding box.
    unsigned int _bbox_margin = 1;

    // Origin coordinates.
    Point _origin {0.0f, 0.0f, 0.0f};

    // Atoms coordinates. Using SoA.
    std::vector<grid_t> _x, _y, _z;

    // Atoms VdW radii
    std::vector<float> _radii;

    // Bouding box
    BoundingBox _bbox;
};

// Doing this in O(N), without any sorting or binary search. We'll see. TODO.
std::vector<uint32_t> get_atoms_in_bbox(
    GridMolecule const &molecule, BoundingBox const &bbox);

class GridConvexHull {
public:
    GridConvexHull() = default;
    GridConvexHull(ConvexHull const &CH);

    void draw(std::string const &ou_fil);

    std::vector<float> _dots;
};

class GridMatrix {
public:
    GridMatrix() = default;

    // Delegated constructor
    GridMatrix(grid_t const x, grid_t const y, grid_t const z,
        float const resolution, Point const origin);

    // Delegating constructor
    GridMatrix(GridMolecule const &gmolecule);

    // GridMatrix dimensions.
    grid_t _dimx, _dimy, _dimz;

    // Number of grid particles.
    grid_t _n;

    // Flattened 3D array that indicates if a grid particle is filled or not.
    std::vector<bool> _bool;

    // Resolution used to build the GridMolecule.
    float _resolution = 1.0;

    // Origin coordinates.
    Point _origin {0.0f, 0.0f, 0.0f};
};

auto fill_grid_tetrahedron(GridMolecule const &in_mol)
    -> std::vector<std::vector<std::vector<grid_t>>>;

} // namespace GrANA

#endif // GrANA_GRID