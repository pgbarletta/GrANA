#ifndef GrANA_CONTINUOUS_H
#define GrANA_CONTINUOUS_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

using EPIC = CGAL::Exact_predicates_inexact_constructions_kernel;
using CPoint = EPIC::Point_3;
using Polyhedron = CGAL::Polyhedron_3<EPIC>;

#include "GrANA/continuous_base.hpp"
#include "GrANA/utils.hpp"
#include "chemfiles.hpp"

using P_Facet_iterator = Polyhedron::Facet_iterator;
using P_Facet_const_iterator = Polyhedron::Facet_const_iterator;
using P_Edge_iterator = Polyhedron::Edge_iterator;
using P_Edge_const_iterator = Polyhedron::Edge_const_iterator;
using P_Halfedge_around_facet_circulator =
    Polyhedron::Halfedge_around_facet_circulator;
using P_Halfedge_around_facet_const_circulator =
    Polyhedron::Halfedge_around_facet_const_circulator;
using P_Vertex_iterator = Polyhedron::Vertex_iterator;
using P_Vertex_const_iterator = Polyhedron::Vertex_const_iterator;
using Delaunay = CGAL::Delaunay_triangulation_3<EPIC>;
using Finite_cells_iterator = Delaunay::Finite_cells_iterator;

namespace GrANA {

class Molecule {
public:
    Molecule() = default;
    Molecule(std::vector<float> &&x, std::vector<float> &&y,
        std::vector<float> &&z, std::vector<float> &&radii) :
        _x(x),
        _y(y), _z(z), _radii(radii), _natoms(radii.size()) {
        _empty = false;
    };

    void draw(std::string const &out_file);

    std::vector<float> _x;
    std::vector<float> _y;
    std::vector<float> _z;
    std::vector<float> _radii;
    int _natoms;
    bool _empty = true;
};

class ConvexHull {
public:
    ConvexHull() = default;
    ConvexHull(Molecule const &prote, std::vector<int> const &indices);

    void draw(std::string const &out_file);

    int _ntriangles;
    std::vector<Triangle> _triangles;
};

class Triangulation {
public:
    Triangulation() = default;
    Triangulation(Molecule const &prote, std::vector<int> const &indices);

    void draw(const std::string &out_file);

    int _ntetrahedrons;
    std::vector<Tetrahedron> _tetrahedrons;
    std::vector<Cube> _bboxes;
};

class PointSoA {
public:
    PointSoA(Triangulation const &T) noexcept;

private:
    float *_x = NULL, *_y = NULL, *_z = NULL;
};

// Get the coordinates of the "indices" atoms as CGAL points.
inline std::vector<CPoint> atom_indices_to_points(
    Molecule const &prote, std::vector<int> const &indices) {

    std::size_t const sz = indices.size();
    std::vector<CPoint> point_set;
    point_set.reserve(sz);

    for (std::size_t i = 0; i < sz; ++i) {
        auto const x = prote._x[indices[i] - 1];
        auto const y = prote._y[indices[i] - 1];
        auto const z = prote._z[indices[i] - 1];
        point_set.emplace_back(x, y, z);
    }
    return point_set;
}
} // namespace GrANA

#endif // _H