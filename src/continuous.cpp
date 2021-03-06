#include "GrANA/continuous.hpp"

namespace GrANA {

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

} // namespace GrANA
