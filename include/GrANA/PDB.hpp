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
void draw_PDB(GridPoint const &gpoint, FILE *ou_fil, int idx, int resid,
    Point const &origin, float const resolution);

// Write GridMolecule to PDB file.
void draw_PDB(GridMolecule const &gmolecula, std::string const &ou_fil);

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

// Write prism. Can't draw connectivity properly if the prism wasn't constructed
// with proper Point ordering. So this class is kind of useless.
void draw_PDB(
    Prism const &prisma, FILE *out_file, int const start_idx, int const resid);

// Write prism to a new PDB file.
void draw_PDB(Prism const &prisma, std::string const &out_file);

// Write Triangulation to PDB file.
void draw_PDB(Triangulation const &triangulacion, std::string const &out_file);

// Write ConvexHull to PDB file.
void draw_PDB(ConvexHull const &ch, std::string const &out_file);

// Write Molecule to PDB file.
void draw_PDB(Molecule const &molecula, std::string const &out_file);

// Write grid matrix to PDB file.
void draw_grid_mtx(std::string const &out_fil,
    std::vector<std::vector<std::vector<grid_t>>> const &mtx,
    Point const &origin, float const resolution);

}
#endif // _H