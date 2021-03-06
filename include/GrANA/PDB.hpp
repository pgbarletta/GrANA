#ifndef GrANA_PDB_H
#define GrANA_PDB_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "GrANA/continuous.hpp"
#include "GrANA/grid.hpp"
#include "GrANA/utils.hpp"
#include "chemfiles.hpp"

namespace GrANA {

// Generates a molecule for the protein and another one for the ligand, if
// present.
auto read_PDB(std::string const &in_filename) -> std::tuple<Molecule, Molecule>;

// Write GridPoint as atom to PDB file.
void draw(GridPoint const &gpoint, FILE *ou_fil, int idx, int resid,
    float const resolution, Point const &origin);

// Write GridMolecule to PDB file.
void draw_PDB(GridMolecule const &gmolecula, std::string const &ou_fil);

// Write selected GridPoints from GridMolecule to PDB file.
void draw_selection_PDB(GridMolecule const &gmolecula,
    std::string const &ou_fil, std::vector<uint32_t> const &selection);

// Write range of GridPoints from GridMolecule to PDB file.
void draw_range_PDB(GridMolecule const &gmolecula, std::string const &ou_fil,
    std::tuple<uint32_t, uint32_t> const range);

// Write Point as atom to PDB file.
void draw(Point const &point, FILE *out_file, int const idx, int const resid);

// Write triangle to PDB file.
void draw_PDB(Triangle const &triangle, FILE *out_file, int const start_idx,
    int const resid);

// Write tetrahedron to PDB file.
void draw_PDB(Tetrahedron const &tetrahedro, FILE *out_file,
    int const start_idx, int const resid);

// Write cube to PDB file.
void draw_PDB(
    Cube const &cubo, FILE *out_file, int const start_idx, int const resid);

// Write prism to an existing PDB file. Can't draw connectivity properly if the
// prism wasn't constructed with proper Point ordering. So this class is kind of
// useless.
void draw(
    Prism const &prisma, FILE *out_file, int const start_idx, int const resid);

// Write prism to a new PDB file.
void draw_PDB(Prism const &prisma, std::string const &out_file);

// Write BoundingBox to a new PDB file.
void draw_PDB(BoundingBox const &bbox, std::string const &out_file);

// Write ConvexHull to PDB file.
void draw_PDB(ConvexHull const &ch, std::string const &out_file);

// Write Triangulation to PDB file.
void draw_PDB(Triangulation const &triangulacion, std::string const &out_file);

// Write Molecule to PDB file.
void draw_PDB(Molecule const &molecula, std::string const &out_file);

// Write grid matrix to PDB file.
void draw_PDB(GridMatrix const &mtx, std::string const &out_fil);

}
#endif // _H