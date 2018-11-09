#ifndef GrANA_CONTINUOUS_H
#define GrANA_CONTINUOUS_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
using EPIC = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = EPIC::Point_3;
using Polyhedron = CGAL::Polyhedron_3<EPIC>;

#include "GrANA/primitives.hpp"
#include "GrANA/utils.hpp"
#include "chemfiles.hpp"

using P_Facet_iterator = Polyhedron::Facet_iterator;
using P_Facet_const_iterator = Polyhedron::Facet_const_iterator;
using P_Edge_iterator = Polyhedron::Edge_iterator;
using P_Edge_const_iterator = Polyhedron::Edge_const_iterator;
using P_Halfedge_around_facet_circulator
    = Polyhedron::Halfedge_around_facet_circulator;
using P_Halfedge_around_facet_const_circulator
    = Polyhedron::Halfedge_around_facet_const_circulator;
using P_Vertex_iterator = Polyhedron::Vertex_iterator;
using P_Vertex_const_iterator = Polyhedron::Vertex_const_iterator;

namespace GrANA {
class Molecule {
public:
    Molecule() = default;
    Molecule(std::string const &in_filename);

    ~Molecule() {
        free(_xyz);
        free(_in_xyz);
        free(_radii);
        free(_in_radii);
    }

    void draw(std::string const &out_file);

    int _natoms;
    Point *_xyz = NULL, *_in_xyz = NULL;
    float *_radii = NULL, *_in_radii = NULL;
};

class ConvexHull {
public:
    ConvexHull() = default;
    ConvexHull(Molecule const &prote, std::vector<int> const &indices);

    ~ConvexHull() { free(_triangles); }

    void draw(const std::string &out_file);

    int _ntriangles;
    Triangle *_triangles = NULL;
};

class Triangulation {
public:
    Triangulation() = default;
    Triangulation(Molecule const &prote, std::vector<int> const &indices);

    ~Triangulation() {
        free(_tetrahedrons);
        free(_bboxes);
    }

    void draw(const std::string &out_file);

    int _ntetrahedrons;
    Tetrahedron *_tetrahedrons = NULL;
    Cube *_bboxes = NULL;
};

class BoundingBox {
public:
    BoundingBox() = default;

    BoundingBox(float const xmin, float const ymin, float const zmin,
        float const xmax, float const ymax, float const zmax) noexcept :
        _xmin(xmin),
        _ymin(ymin), _zmin(zmin), _xmax(xmax), _ymax(ymax), _zmax(zmax){};

    BoundingBox(Tetrahedron const &in_tet) noexcept;

    float _xmin, _ymin, _zmin, _xmax, _ymax, _zmax;
};

class PointSoA {
public:
    PointSoA(Triangulation const &T) noexcept;

private:
    float *_x = NULL, *_y = NULL, *_z = NULL;
};

// Get the coordinates of the "indices" atoms as CGAL points.
inline std::vector<Point_3> atom_indices_to_points(
    Molecule const &prote, std::vector<int> const &indices) {

    std::size_t const sz = indices.size();
    std::vector<Point_3> point_set(sz);

    for (std::size_t i = 0; i < sz; ++i) {
        const auto xyz = prote._xyz[indices[i] - 1];
        point_set[sz] = Point_3(xyz[0], xyz[1], xyz[2]);
    }
    return point_set;
}
} // namespace

#endif // _H
