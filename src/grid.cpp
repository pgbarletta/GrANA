#include "GrANA/grid.hpp"
namespace GrANA {

GridMolecule::GridMolecule(Molecule const &in_mol, float const resolution,
    unsigned int const bbox_margin) :
    _idx_x(sort_indices_uint32(in_mol._x)),
    _idx_y(sort_indices_uint32(in_mol._y)),
    _idx_z(sort_indices_uint32(in_mol._z)), _natoms(in_mol._natoms),
    _resolution(resolution), _bbox_margin(bbox_margin) {

    float const safe_distance = _resolution * _bbox_margin;
    // std::floor for rounding.
    _origin = Point(std::floor(in_mol._x[_idx_x[0]]) - safe_distance,
        std::floor(in_mol._y[_idx_y[0]]) - safe_distance,
        std::floor(in_mol._z[_idx_z[0]]) - safe_distance);

    construct_GridMolecule(in_mol);

    return;
}

// Build a GridMolecule using a given origin.
GridMolecule::GridMolecule(Molecule const &in_mol, float const resolution,
    unsigned int const bbox_margin, Point const &origin) :
    _idx_x(sort_indices_uint32(in_mol._x)),
    _idx_y(sort_indices_uint32(in_mol._y)),
    _idx_z(sort_indices_uint32(in_mol._z)), _natoms(in_mol._natoms),
    _resolution(resolution), _bbox_margin(bbox_margin), _origin(origin) {

    construct_GridMolecule(in_mol);

    return;
}

// Finish building GridMolecule.
void GridMolecule::construct_GridMolecule(Molecule const &in_mol) {
    _x.reserve(_natoms * 3);
    _y.reserve(_natoms * 3);
    _z.reserve(_natoms * 3);
    _radii.reserve(_natoms);

    for (int i = 0; i < _natoms; ++i) {
        grid_t const x = cont_to_grid(in_mol._x[i] - _origin[0], _resolution);
        grid_t const y = cont_to_grid(in_mol._y[i] - _origin[1], _resolution);
        grid_t const z = cont_to_grid(in_mol._z[i] - _origin[2], _resolution);

        _x.push_back(x);
        _y.push_back(y);
        _z.push_back(z);
        _radii.push_back(in_mol._radii[i]);
    }

    // in_mol._(x,y,z)[_idx_(x,y,z)[0]] - _bbox_margin is a bug.
    fill_bounding_box(_bbox, _x[_idx_x[0]] - _bbox_margin,
        _y[_idx_y[0]] - _bbox_margin, _z[_idx_z[0]] - _bbox_margin,
        _x[_idx_x[_idx_x.size() - 1]] + _bbox_margin,
        _y[_idx_y[_idx_x.size() - 1]] + _bbox_margin,
        _z[_idx_z[_idx_x.size() - 1]] + _bbox_margin, _resolution, _origin,
        _bbox_margin);

    return;
}

// Build a GridMolecule from a selection of another GridMolecule. The _bbox
// bounding box will not be its bounding box but its parents.
GridMolecule::GridMolecule(
    GridMolecule const &in_gmol, std::vector<uint32_t> const subset) :
    _natoms(subset.size()),
    _resolution(in_gmol._resolution), _bbox_margin(in_gmol._bbox_margin),
    _origin(in_gmol._origin), _bbox(in_gmol._bbox) {

    _x.reserve(_natoms);
    _y.reserve(_natoms);
    _z.reserve(_natoms);
    _radii.reserve(_natoms);

    for (uint32_t i = 0; static_cast<int>(i) < _natoms; ++i) {
        auto const atm = subset[i];

        _x.push_back(in_gmol._x[atm]);
        _y.push_back(in_gmol._y[atm]);
        _z.push_back(in_gmol._z[atm]);
        _radii.push_back(in_gmol._radii[atm]);
    }

    // TODO. Por ahora no le copio los índices q ordenan _idx_x, _idx_y, _idx_z

    return;
}

// Mainly used by GridMolecule to fill its bounding box.
void fill_bounding_box(BoundingBox &bbox, grid_t const xmin, grid_t const ymin,
    grid_t const zmin, grid_t const xmax, grid_t const ymax, grid_t const zmax,
    float const resolution, Point const &origin, unsigned int const margin) {

    bbox._resolution = resolution;
    bbox._origin = origin;
    bbox._bbox_margin = margin;

    bbox._gp[0] = GridPoint(xmin, ymin, zmin);
    bbox._gp[1] = GridPoint(xmax, ymin, zmin);
    bbox._gp[2] = GridPoint(xmin, ymax, zmin);
    bbox._gp[3] = GridPoint(xmax, ymax, zmin);
    bbox._gp[4] = GridPoint(xmin, ymin, zmax);
    bbox._gp[5] = GridPoint(xmax, ymin, zmax);
    bbox._gp[6] = GridPoint(xmin, ymax, zmax);
    bbox._gp[7] = GridPoint(xmax, ymax, zmax);

    bbox._dimx = bbox._gp[7][0] - bbox._gp[0][0];
    bbox._dimy = bbox._gp[7][1] - bbox._gp[0][1];
    bbox._dimz = bbox._gp[7][2] - bbox._gp[0][2];

    // don't mind the rounding.
    bbox._center = {bbox._dimx / 2, bbox._dimy / 2, bbox._dimz / 2};

    return;
}

// Doing this in O(N), without any sorting or binary search. We'll see. TODO.
GridMolecule get_atoms_in_bbox(
    GridMolecule const &gmolecule, BoundingBox const &bbox) {

    // TODO. Hardcodeo 2 angstroms de margen p/ incluir átomos.
    // grid_t const x_min = bbox._gp[0][0] - (2 / bbox._resolution);
    // grid_t const x_max = bbox._gp[7][0] + (2 / bbox._resolution);

    // grid_t const y_min = bbox._gp[0][1] - (2 / bbox._resolution);
    // grid_t const y_max = bbox._gp[7][1] + (2 / bbox._resolution);

    // grid_t const z_min = bbox._gp[0][2] - (2 / bbox._resolution);
    // grid_t const z_max = bbox._gp[7][2] + (2 / bbox._resolution);

    grid_t const x_min = bbox._gp[0][0];
    grid_t const x_max = bbox._gp[7][0];
    grid_t const y_min = bbox._gp[0][1];
    grid_t const y_max = bbox._gp[7][1];
    grid_t const z_min = bbox._gp[0][2];
    grid_t const z_max = bbox._gp[7][2];

    std::vector<uint32_t> in_bounding_box;

    for (uint32_t i = 0; static_cast<int>(i) < gmolecule._natoms; ++i) {

        if ((gmolecule._x[i] > x_min) && (gmolecule._x[i] < x_max)) {
            if ((gmolecule._y[i] > y_min) && (gmolecule._y[i] < y_max)) {
                if ((gmolecule._z[i] > z_min) && (gmolecule._z[i] < z_max)) {

                    in_bounding_box.push_back(i);
                }
            }
        }
    }

    return {gmolecule, in_bounding_box};
}

// Delegating constructor.
GridMatrix::GridMatrix(
    GridMolecule const &gmolecule, uint32_t const gap_depth) :
    GridMatrix(gmolecule._bbox._dimx + 1, gmolecule._bbox._dimy + 1,
        gmolecule._bbox._dimz + 1, gap_depth, gmolecule._resolution,
        gmolecule._bbox._gp[0], gmolecule._origin) { }

// Delegated constructor.
GridMatrix::GridMatrix(grid_t const x, grid_t const y, grid_t const z,
    uint32_t const gap_depth, float const resolution,
    GridPoint const bbox_vtx_0, Point const molecule_origin) :
    _dimx(x),
    _dimy(y), _dimz(z), _gap_depth(gap_depth), _plane_size(x * y),
    _n(x * y * z), _resolution(resolution), _grid_origin(bbox_vtx_0),
    _molecule_origin(molecule_origin) {

    _origin = GridPoint_to_point(_grid_origin, _resolution, _molecule_origin);

    _bool.reserve(_n);

    // First `_gap_depth` planes are empty.
    _bool.insert(std::end(_bool), _gap_depth * _plane_size, false);
    // Now, the subsequent plantes are empty on the borders, but filled on the
    // inside.
    int const n_planes = _dimz - (2 * _gap_depth);
    for (int i = 0; i < n_planes; i++) {
        // First `_gap_depth` rows are empty.
        _bool.insert(std::end(_bool), _gap_depth * _dimx, false);
        // Now, the next rows begin with `_gap_depth` grid points empty, then
        // the next ones are filled, and the final `_gap_depth` grid points are
        // also empty.
        int const n_rows = _dimy - 2 * _gap_depth;
        int const n_dots = _dimx - 2 * _gap_depth;
        for (int i = 0; i < n_rows; i++) {
            _bool.insert(std::end(_bool), _gap_depth, false);
            _bool.insert(std::end(_bool), n_dots, true);
            _bool.insert(std::end(_bool), _gap_depth, false);
        }
        // Last `_gap_depth` rows are empty.
        _bool.insert(std::end(_bool), _gap_depth * _dimx, false);
    }
    // Last `_gap_depth` planes are empty.
    _bool.insert(std::end(_bool), _gap_depth * _plane_size, false);

    return;
}

// Remove grid points overlapping with atoms. TODO. Parallelize this.
void carve_atoms_in_mtx(
    GrANA::GridMolecule const &in_bounding_box, GridMatrix &mtx) {

    grid_t const _plane_size = mtx._dimx * mtx._dimy;
    std::vector<uint32_t> indices_ocupados;

    for (int a = 0; a < in_bounding_box._natoms; ++a) {
        // The atoms are using a different grid origin that the matrix.
        grid_t const x_off = in_bounding_box._x[a] - mtx._grid_origin[0];
        grid_t const y_off = in_bounding_box._y[a] - mtx._grid_origin[1];
        grid_t const z_off = in_bounding_box._z[a] - mtx._grid_origin[2];

        // What follows is all about turning the x, y, z coordinates of the atom
        // center and all the overlapping grid points into a linear array.
        int const centro = z_off * _plane_size + y_off * mtx._dimx + x_off;
        int const ancho =
            std::ceil(in_bounding_box._radii[a] / in_bounding_box._resolution);

        int const extremo_inferior =
            centro - ancho * _plane_size - ancho * mtx._dimx;
        int const extremo_superior =
            centro + ancho * _plane_size + ancho * mtx._dimx;

        if (extremo_inferior < 0 ||
            extremo_superior >= static_cast<int>(mtx._n)) {
            // The atom is on the outer cells of the box. Its VdW radius
            // protrudes out of the box. Trying to set surrounding cells as
            // empty will cause a segmentation fault on `mtx._bool` vector.

            for (int k = -ancho; k <= ancho; ++k) {
                for (int j = -ancho; j <= ancho; ++j) {

                    int const beg =
                        centro - ancho + k * _plane_size + j * mtx._dimx;
                    int const end = beg + 2 * ancho;

                    for (uint32_t i = beg; static_cast<int>(i) < end; ++i) {
                        if (i < mtx._n) {
                            // if i >= mtx._n, then beg was negative and `i`
                            // overflowed (atome close to the beginnning of the
                            // matrix), or`i` is just too large (atom close to
                            // the end of the matrix).

                            // VdW radii of bonded atoms overlap, hence this is
                            // gonna be a race condition when parallelized. I
                            // don't mind 'cause all I'm doing is setting it to
                            // false.

                            mtx._bool[i] = false;
                        }
                    }
                }
            }
        } else {
            for (int k = -ancho; k <= ancho; ++k) {
                for (int j = -ancho; j <= ancho; ++j) {

                    int const beg = centro + k * _plane_size + j * mtx._dimx;
                    int const end = beg + ancho;

                    for (uint32_t i = beg; static_cast<int>(i) < end; ++i) {
                        // VdW radii of bonded atoms overlap, hence this is
                        // gonna be a race condition when parallelized. I don't
                        // mind 'cause all I'm doing is setting it to false.

                        mtx._bool[i] = false;
                    }
                }
            }
        }
    }

    return;
}

}
