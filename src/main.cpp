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

    auto [in_pdb, resolution, margin, gap_depth, out_pdb] =
        GrANA::get_input(argc, argv);

    std::vector<int> indices = {
        300, 600, 900, 1200, 1500, 1800, 1240, 400, 500, 700, 800, 1000, 1100};

    auto const [prote, lig] = GrANA::read_PDB(in_pdb);

    GrANA::GridMolecule Gprote(prote, resolution, margin);
    GrANA::GridMolecule Glig(lig, resolution, margin, Gprote._origin);

    GrANA::GridMolecule const in_bounding_box =
        get_atoms_in_bbox(Gprote, Glig._bbox);

    GrANA::GridMatrix mtx(Glig, gap_depth);

    std::cout << "Grid matrix size:  " << mtx._n << '\n';
    GrANA::draw_PDB(mtx, "grilla.pdb");

    GrANA::carve_atoms_in_mtx(in_bounding_box, mtx);

    GrANA::draw_PDB(Gprote, "gprote.pdb");
    GrANA::draw_PDB(Glig, "glig.pdb");
    GrANA::draw_PDB(Gprote._bbox, "pbbox.pdb");
    GrANA::draw_PDB(Glig._bbox, "bbox.pdb");
    GrANA::draw_PDB(in_bounding_box, "in_bbox.pdb");

    GrANA::draw_PDB(mtx, "hueco.pdb");

    for (const auto &each : Gprote._radii) {
        std::cout << each << '\n';
    }

    return 0;
}
