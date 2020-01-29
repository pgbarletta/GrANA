#include "GrANA/continuous.hpp"

#include <iostream>

#include <CGAL/Delaunay_triangulation_3.h>
using Delaunay = CGAL::Delaunay_triangulation_3<EPIC>;
using Finite_cells_iterator = Delaunay::Finite_cells_iterator;

namespace GrANA {

Molecule::Molecule(std::string const &in_filename) {
    chemfiles::Trajectory in_trj(in_filename);
    auto in_frm = in_trj.read();
    auto in_top = in_frm.topology();
    auto in_xyz = in_frm.positions();
    _natoms = in_xyz.size();

    _x.reserve(_natoms * 3);
    _y.reserve(_natoms * 3);
    _z.reserve(_natoms * 3);
    _radii.reserve(_natoms);

    // Get atoms positions and VdW radii.
    int j = 0;
    for (const auto &residuo : in_top.residues()) {
        for (const auto &i : residuo) {
            const auto atom_xyz = in_xyz[i];
            _x.push_back(static_cast<float>(atom_xyz[0]));
            _y.push_back(static_cast<float>(atom_xyz[1]));
            _z.push_back(static_cast<float>(atom_xyz[2]));
            _radii.push_back(in_top[i].vdw_radius().value_or(1.5));
            ++j;
        }
    }
}

// Draw the molecule in the **out_file** path in PDB format.
void Molecule::draw(std::string const &out_file) {

    FILE *file = std::fopen(out_file.c_str(), "w");
    if (file) {
        for (int i = 0; i <= _natoms - 1; ++i) {
            Point const atm(_x[i], _y[i], _z[i]);
            atm.draw(file, i + 1, i + 1);
        }
    } else {
        std::cerr << "Could not open " << out_file << ". " << '\n';
    }
    std::fclose(file);
    return;
}

ConvexHull::ConvexHull(Molecule const &prote, std::vector<int> const &indices) {

    // Get the coordinates of the input indices atoms.
    auto point_vtor = atom_indices_to_points(prote, indices);

    // Get the convex hull.
    Polyhedron con_hul;
    CGAL::convex_hull_3(point_vtor.begin(), point_vtor.end(), con_hul);

    // Turn CGAL's polyhedron holding the convex hull into GrANA triangles.
    P_Facet_const_iterator f_ite = con_hul.facets_begin();
    const P_Facet_const_iterator f_end = con_hul.facets_end();
    _ntriangles = static_cast<int>(std::distance(f_ite, f_end));
    _triangles.reserve(_ntriangles);

    for (int i = 0; f_ite != f_end; ++f_ite) {
        P_Halfedge_around_facet_const_circulator he_ite = f_ite->facet_begin();
        _triangles.emplace_back(he_ite->vertex()->point(),
            (he_ite++)->vertex()->point(), (he_ite++)->vertex()->point());
        ++i;
    }
}

void ConvexHull::draw(std::string const &out_file) {

    FILE *file = std::fopen(out_file.c_str(), "w");
    if (file) {
        int resid = 1;
        for (int i = 0; i < _ntriangles; ++i) {
            const auto start_idx = i * 3 + 1;
            _triangles[i].draw(file, start_idx, resid++);
        }
        for (int i = 0; i < _ntriangles; ++i) {
            const auto j = i * 3 + 1;
            fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n", j, j + 1, j + 2);
            fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n", j + 1, j, j + 2);
        }
    } else {
        std::cerr << "Could not open " << out_file << ". " << '\n';
    }
    std::fclose(file);

    return;
}

Triangulation::Triangulation(
    Molecule const &prote, std::vector<int> const &indices) {

    // Get the coordinates of the input indices atoms.
    auto point_vtor = atom_indices_to_points(prote, indices);

    // Get the convex hull.
    Polyhedron con_hul;
    CGAL::convex_hull_3(point_vtor.begin(), point_vtor.end(), con_hul);

    // Get the triangulation from the convex hull.
    Delaunay T(con_hul.points_begin(), con_hul.points_end());

    // Turn CGAL's Delaunay triangulation into GrANA triangulation.
    auto cell_ite = T.finite_cells_begin();
    auto const cell_end = T.finite_cells_end();
    _ntetrahedrons = static_cast<int>(std::distance(cell_ite, cell_end));
    _tetrahedrons.reserve(_ntetrahedrons);
    _bboxes.reserve(_ntetrahedrons);

    // // Initialize bounding box.
    // constexpr float x = -999.0f;
    // constexpr float y = -999.0f;
    // constexpr float z = -999.0f;
    // bbo_x = cube(point(-x, -y, -z), point(-x, -y, z), point(-x, y, z),
    // 	point(-x, y, -z), point(x, -y, -z), point(x, -y, z),
    // 	point(x, y, z), point(x, y, -z));

    for (; cell_ite != cell_end; ++cell_ite) {
        _tetrahedrons.emplace_back(cell_ite->vertex(0)->point(),
            cell_ite->vertex(1)->point(), cell_ite->vertex(2)->point(),
            cell_ite->vertex(3)->point());
    }

    return;
}

void Triangulation::draw(std::string const &out_file) {

    FILE *file = std::fopen(out_file.c_str(), "w");
    if (file) {
        int resid = 1;
        for (int i = 0; i < _ntetrahedrons; ++i) {
            const auto start_idx = i * 4 + 1;
            _tetrahedrons[i].draw(file, start_idx, resid++);
        }
        for (int i = 0; i < _ntetrahedrons; ++i) {
            const auto j = i * 4 + 1;
            fmt::print(
                file, "CONECT {:>4} {:>4} {:>4}\n", j, j + 1, j + 2, j + 3);
            fmt::print(
                file, "CONECT {:>4} {:>4} {:>4}\n", j + 1, j + 2, j + 3, j);
            fmt::print(
                file, "CONECT {:>4} {:>4} {:>4}\n", j + 2, j + 3, j, j + 1);
            fmt::print(
                file, "CONECT {:>4} {:>4} {:>4}\n", j + 3, j, j + 1, j + 2);
        }
    } else {
        std::cerr << "Could not open " << out_file << ". " << '\n';
    }
    std::fclose(file);

    return;
}

BoundingBox::BoundingBox(Tetrahedron const &in_T) noexcept {
    // in_T._p0.

    return;
}

} // namespace GrANA
