#include "GrANA/grid.hpp"
namespace GrANA {

// Draw GridPoint as atom.
void GridPoint::draw(FILE *ou_fil, int idx, int resid) {
    const float fx = grid_to_cont(_xyz[0]);
    const float fy = grid_to_cont(_xyz[1]);
    const float fz = grid_to_cont(_xyz[2]);
    fmt::print(ou_fil,
        "{: <6}{: >5} {: <4s} {:3} {:1}{: >4}    "
        "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {: >2s}\n",
        "HETATM", idx, "H", "GPU", "A", resid, fx, fy, fz, 1.0, 0.0, "H");
    return;
}

std::ostream &operator<<(std::ostream &stream, const GridPoint &p) {
    stream << p[0] << " " << p[1] << " " << p[2];
    return stream;
}

// Turn a grid point into a continuous point, given the resolution.
// This is not a constructor in the `Point` class to avoid circular
// dependencies.
Point GridPoint_to_point(const GridPoint &in_point) {
    return Point(grid_to_cont(in_point[0]), grid_to_cont(in_point[1]),
        grid_to_cont(in_point[2]));
}

GridMolecule::GridMolecule(Molecule const &in_mol, Point const &orig_point) {

    _orig_vtor = point_to_vector(orig_point);
    _natoms = in_mol._natoms;

    _xyz.reserve(_natoms * 3);
    _in_xyz.reserve(_natoms * 3);
    _radii.reserve(_natoms);
    _in_radii.reserve(_natoms);

    for (int i = 0; i < _natoms; ++i) {
        _xyz.emplace_back(in_mol._xyz[i] - _orig_vtor);
        _radii.push_back(in_mol._radii[i]);
    }
}

void GridMolecule::draw(const std::string &ou_fil) {

    FILE *file = std::fopen(ou_fil.c_str(), "w");
    if (file) {
        for (int i = 1; i <= _natoms; ++i) {
            _xyz[i - 1].draw(file, i, i);
        }
    } else {
        std::cerr << "Could not open " << ou_fil << ". " << '\n';
    }
    std::fclose(file);
    return;
}

}
