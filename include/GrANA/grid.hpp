#ifndef GrANA_GRID
#define GrANA_GRID

#include "GrANA/grid_base.hpp"
#include "utils.hpp"
#include <tuple>

namespace GrANA {

class GridMolecule {
public:
    GridMolecule() = default;
    GridMolecule(Molecule const &in_mol) :
        _idx_x(sort_indices(in_mol._x)), _idx_y(sort_indices(in_mol._y)),
        _idx_z(sort_indices(in_mol._z)), _natoms(in_mol._natoms) {

        _orig_vtor = Vector(std::floor(in_mol._x[_idx_x[0]]),
            std::floor(in_mol._y[_idx_y[0]]), std::floor(in_mol._z[_idx_z[0]]));

        _x.reserve(_natoms * 3);
        _y.reserve(_natoms * 3);
        _z.reserve(_natoms * 3);
        _radii.reserve(_natoms);

        for (int i = 0; i < _natoms; ++i) {
            _x.push_back(cont_to_grid(in_mol._x[i] - _orig_vtor[0]));
            _y.push_back(cont_to_grid(in_mol._y[i] - _orig_vtor[1]));
            _z.push_back(cont_to_grid(in_mol._z[i] - _orig_vtor[2]));
            _radii.push_back(in_mol._radii[i]);
        }
    }

    void draw(std::string const &ou_fil);

    std::vector<int> const _idx_x;
    std::vector<int> const _idx_y;
    std::vector<int> const _idx_z;
    int const _natoms;
    Vector _orig_vtor{0.0f, 0.0f, 0.0f};
    std::vector<grid_t> _x;
    std::vector<grid_t> _y;
    std::vector<grid_t> _z;
    std::vector<float> _radii;
    // Indices that sort the atoms along the 3 axes.
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

void fill_grid_tetrahedron(std::vector<grid_t> i_x, std::vector<grid_t> g_x,
    std::vector<grid_t> g_y, std::vector<grid_t> g_z, std::vector<float> radii);

} // namespace GrANA

#endif // GrANA_GRID