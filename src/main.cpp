#include <cmath>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

#include "GrANA/continuous.hpp"
#include "GrANA/grid.hpp"
#include "GrANA/utils.hpp"
#include <GrANA/PDB.hpp>
#include <GrANA/octree.hpp>

int main(int argc, char **argv) {

    auto [in_pdb, resolution, margin, out_pdb] = GrANA::get_input(argc, argv);

    std::vector<int> indices = {
        300, 600, 900, 1200, 1500, 1800, 1240, 400, 500, 700, 800, 1000, 1100};

    auto const [prote, lig] = GrANA::read_PDB(in_pdb);

    // GrANA::Triangulation incl_area(prote, indices);
    // incl_area.draw(out_pdb);

    GrANA::GridMolecule Gprote(prote, resolution, margin);
    GrANA::GridMolecule Glig(lig, resolution, margin, Gprote._origin);

    GrANA::draw_PDB(Gprote, "gprote.pdb");
    GrANA::draw_PDB(Glig, "glig.pdb");
    GrANA::draw_PDB(Glig._bbox, "bbox.pdb");

    auto const in_bounding_box = get_atoms_in_bbox(Gprote, Glig._bbox);
    GrANA::draw_selection_PDB(Gprote, "in_bbox.pdb", in_bounding_box);

    // auto mtx = GrANA::fill_grid_tetrahedron(Gprote, resolution);
    // draw(out_pdb, mtx, Gprote._origin, resolution);

    return 0;
}
