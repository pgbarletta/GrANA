#include "GrANA/PDB.hpp"

namespace GrANA {

auto read_PDB(std::string const &in_filename)
    -> std::tuple<Molecule, Molecule> {
    chemfiles::Trajectory in_trj(in_filename);
    auto in_frm = in_trj.read();
    auto in_top = in_frm.topology();
    auto in_xyz = in_frm.positions();

    std::vector<float> x, x_lig;
    std::vector<float> y, y_lig;
    std::vector<float> z, z_lig;
    std::vector<float> radii, radii_lig;

    int natoms = in_xyz.size();

    x.reserve(natoms * 3);
    y.reserve(natoms * 3);
    z.reserve(natoms * 3);
    radii.reserve(natoms);

    // Get atoms positions and VdW radii.
    for (const auto &residuo : in_top.residues()) {
        bool is_lig = !(
            residuo.get<chemfiles::Property::BOOL>("is_standard_pdb").value());
        if (is_lig) {
            // dealing now with the ligand.
            for (const auto &i : residuo) {
                const auto atom_xyz = in_xyz[i];
                x_lig.push_back(static_cast<float>(atom_xyz[0]));
                y_lig.push_back(static_cast<float>(atom_xyz[1]));
                z_lig.push_back(static_cast<float>(atom_xyz[2]));
                radii_lig.push_back(in_top[i].vdw_radius().value_or(1.5));
            }
        } else {
            for (const auto &i : residuo) {
                const auto atom_xyz = in_xyz[i];
                x.push_back(static_cast<float>(atom_xyz[0]));
                y.push_back(static_cast<float>(atom_xyz[1]));
                z.push_back(static_cast<float>(atom_xyz[2]));
                radii.push_back(in_top[i].vdw_radius().value_or(1.5));
            }
        }
    }
    int natoms_lig = radii_lig.size();

    // Molecule proteina(
    //     std::move(x), std::move(y), std::move(z), std::move(radii));
    if (natoms_lig == 0) {
        // there's no ligand
        // Molecule ligando;
        return {
            {std::move(x), std::move(y), std::move(z), std::move(radii)}, {}};
    } else {
        // Molecule ligando(std::move(x_lig), std::move(y_lig),
        // std::move(z_lig),
        //     std::move(radii_lig));
        return {{std::move(x), std::move(y), std::move(z), std::move(radii)},
            {std::move(x_lig), std::move(y_lig), std::move(z_lig),
                std::move(radii_lig)}};
    }
}

// Write GridPoint as atom to PDB file.
void draw(GridPoint const &gpoint, FILE *ou_fil, int idx, int resid,
    float const resolution, Point const &origin) {
    float const fx = grid_to_cont(gpoint._xyz[0], resolution) + origin[0];
    float const fy = grid_to_cont(gpoint._xyz[1], resolution) + origin[1];
    float const fz = grid_to_cont(gpoint._xyz[2], resolution) + origin[2];
    fmt::print(ou_fil,
        "{: <6}{: >5} {: <4s} {:3} {:1}{: >4}    "
        "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {: >2s}\n",
        "HETATM", idx, "H", "GPU", "A", resid, fx, fy, fz, 1.0, 0.0, "H");
    return;
}

// Write GridMolecule to PDB file.
void draw_PDB(GridMolecule const &gmolecula, std::string const &ou_fil) {

    FILE *file = std::fopen(ou_fil.c_str(), "w");
    if (file) {
        for (int i = 0; i < gmolecula._natoms; ++i) {
            GridPoint const atm(
                gmolecula._x[i], gmolecula._y[i], gmolecula._z[i]);

            draw(atm, file, i + 1, i + 1, gmolecula._resolution,
                gmolecula._origin);
        }
    } else {
        std::cerr << "Could not open " << ou_fil << ". " << '\n';
    }

    fmt::print(file, "END\n");
    std::fclose(file);
    return;
}
// Write selected GridPoints from GridMolecule to PDB file.
void draw_selection_PDB(GridMolecule const &gmolecula,
    std::string const &ou_fil, std::vector<uint32_t> const &selection) {

    FILE *file = std::fopen(ou_fil.c_str(), "w");
    if (file) {
        for (uint32_t const i : selection) {
            GridPoint const atm(
                gmolecula._x[i], gmolecula._y[i], gmolecula._z[i]);

            draw(atm, file, i + 1, i + 1, gmolecula._resolution,
                gmolecula._origin);
        }
    } else {
        std::cerr << "Could not open " << ou_fil << ". " << '\n';
    }

    fmt::print(file, "END\n");
    std::fclose(file);
    return;
}

// Write range of GridPoints from GridMolecule to PDB file.
void draw_range_PDB(GridMolecule const &gmolecula, std::string const &ou_fil,
    std::tuple<uint32_t, uint32_t> const range) {

    FILE *file = std::fopen(ou_fil.c_str(), "w");
    if (file) {
        for (uint32_t i = std::get<0>(range); i < std::get<1>(range); ++i) {
            GridPoint const atm(
                gmolecula._x[i], gmolecula._y[i], gmolecula._z[i]);

            draw(atm, file, i + 1, i + 1, gmolecula._resolution,
                gmolecula._origin);
        }
    } else {
        std::cerr << "Could not open " << ou_fil << ". " << '\n';
    }

    fmt::print(file, "END\n");
    std::fclose(file);
    return;
}

// Write Point as atom to PDB file.
void draw(Point const &point, FILE *out_file, int const idx, int const resid) {
    fmt::print(out_file,
        "{: <6}{: >5} {: <4s} {:3} {:1}{: >4}    "
        "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {: >2s}\n",
        "HETATM", idx, "H", "GPU", "A", resid, point._xyz[0], point._xyz[1],
        point._xyz[2], 1.0, 0.0, "H");
    return;
}

// Write triangle.
void draw(Triangle const &triangulo, FILE *out_file, int const start_idx,
    int const resid) {
    const auto i = start_idx;
    const auto j = start_idx + 1;
    const auto k = start_idx + 2;

    draw(triangulo._p[0], out_file, i, resid);
    draw(triangulo._p[1], out_file, j, resid);
    draw(triangulo._p[2], out_file, k, resid);

    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", i, j, k);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", j, k, i);
    return;
}

// Write tetrahedron.
void draw(Tetrahedron const &tetrahedro, FILE *out_file, int const start_idx,
    int const resid) {
    const auto i = start_idx;
    const auto j = start_idx + 1;
    const auto k = start_idx + 2;
    const auto l = start_idx + 3;

    draw(tetrahedro._p[0], out_file, i, resid);
    draw(tetrahedro._p[1], out_file, j, resid);
    draw(tetrahedro._p[2], out_file, k, resid);
    draw(tetrahedro._p[3], out_file, l, resid);

    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", i, j, k, l);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", j, k, l, i);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", k, l, i, j);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", l, i, j, k);
    return;
}

// Write cube.
void draw_PDB(
    Cube const &cubo, FILE *out_file, int const start_idx, int const resid) {
    const auto i = start_idx;
    const auto j = start_idx + 1;
    const auto k = start_idx + 2;
    const auto l = start_idx + 3;
    const auto ii = start_idx + 4;
    const auto jj = start_idx + 5;
    const auto kk = start_idx + 6;
    const auto ll = start_idx + 7;

    draw(cubo._p[0], out_file, i, resid);
    draw(cubo._p[1], out_file, j, resid);
    draw(cubo._p[2], out_file, k, resid);
    draw(cubo._p[3], out_file, l, resid);
    draw(cubo._p[4], out_file, ii, resid);
    draw(cubo._p[5], out_file, jj, resid);
    draw(cubo._p[6], out_file, kk, resid);
    draw(cubo._p[7], out_file, ll, resid);

    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", i, j, l, ii);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", k, j, l, kk);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", jj, ii, kk, j);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", ll, ii, kk, l);

    return;
}

// Write prism to an existing PDB file. Can't draw connectivity properly if the
// prism wasn't constructed with proper Point ordering. So this class is kind of
// useless.
void draw(
    Prism const &prisma, FILE *out_file, int const start_idx, int const resid) {
    const auto i = start_idx;
    const auto j = start_idx + 1;
    const auto k = start_idx + 2;
    const auto l = start_idx + 3;
    const auto ii = start_idx + 4;
    const auto jj = start_idx + 5;
    const auto kk = start_idx + 6;
    const auto ll = start_idx + 7;

    draw(prisma._p[0], out_file, i, resid);
    draw(prisma._p[1], out_file, j, resid);
    draw(prisma._p[2], out_file, k, resid);
    draw(prisma._p[3], out_file, l, resid);
    draw(prisma._p[4], out_file, ii, resid);
    draw(prisma._p[5], out_file, jj, resid);
    draw(prisma._p[6], out_file, kk, resid);
    draw(prisma._p[7], out_file, ll, resid);

    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", i, j, k, ii);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", j, l, jj);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", k, l, kk);
    fmt::print(out_file, "CONECT {:>4} {:>4}\n", l, ll);

    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", ii, jj, kk);
    fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", ll, jj, kk);

    return;
}

// Write prism to a new PDB file.
void draw_PDB(Prism const &prisma, std::string const &out_file) {

    FILE *file = std::fopen(out_file.c_str(), "w");
    int const n_prismas = 1;
    if (file) {
        int resid = 1;
        for (int i = 0; i < n_prismas; ++i) {
            const auto start_idx = i * 8 + 1;
            draw(prisma, file, start_idx, resid++);
        }
    } else {
        std::cerr << "Could not open " << out_file << ". " << '\n';
    }

    fmt::print(file, "END\n");
    std::fclose(file);
    return;
}

// Write BoundingBox to a new PDB file.
void draw_PDB(BoundingBox const &bbox, std::string const &out_file) {

    FILE *file = std::fopen(out_file.c_str(), "w");
    if (file) {
        for (int i = 0; i < 8; ++i) {
            draw(bbox._gp[i], file, i + 1, 1, bbox._resolution, bbox._origin);
        }
    } else {
        std::cerr << "Could not open " << out_file << ". " << '\n';
    }

    int const i = 1;
    int const j = i + 1;
    int const k = i + 2;
    int const l = i + 3;
    int const ii = i + 4;
    int const jj = i + 5;
    int const kk = i + 6;
    int const ll = i + 7;

    fmt::print(file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", i, j, k, ii);
    fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n", j, l, jj);
    fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n", k, l, kk);
    fmt::print(file, "CONECT {:>4} {:>4}\n", l, ll);

    fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n", ii, jj, kk);
    fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n", ll, jj, kk);

    fmt::print(file, "END\n");
    std::fclose(file);
    return;
}

// Write ConvexHull to PDB file.
void draw_PDB(ConvexHull const &ch, std::string const &out_file) {

    FILE *file = std::fopen(out_file.c_str(), "w");
    if (file) {
        int resid = 1;
        for (int i = 0; i < ch._ntriangles; ++i) {
            const auto start_idx = i * 3 + 1;
            draw(ch._triangles[i], file, start_idx, resid++);
        }
        for (int i = 0; i < ch._ntriangles; ++i) {
            const auto j = i * 3 + 1;
            fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n", j, j + 1, j + 2);
            fmt::print(file, "CONECT {:>4} {:>4} {:>4}\n", j + 1, j, j + 2);
        }
    } else {
        std::cerr << "Could not open " << out_file << ". " << '\n';
    }
    fmt::print(file, "END\n");
    std::fclose(file);

    return;
}

// Write Triangulation to PDB file.
void draw_PDB(Triangulation const &triangulacion, std::string const &out_file) {

    FILE *file = std::fopen(out_file.c_str(), "w");
    if (file) {
        int resid = 1;
        for (int i = 0; i < triangulacion._ntetrahedrons; ++i) {
            const auto start_idx = i * 4 + 1;
            draw(triangulacion._tetrahedrons[i], file, start_idx, resid++);
        }
        for (int i = 0; i < triangulacion._ntetrahedrons; ++i) {
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
    fmt::print(file, "END\n");
    std::fclose(file);

    return;
}

// Write Molecule to PDB file.
void draw_PDB(Molecule const &molecula, std::string const &out_file) {

    FILE *file = std::fopen(out_file.c_str(), "w");
    if (file) {
        for (int i = 0; i <= molecula._natoms - 1; ++i) {
            Point const atm(molecula._x[i], molecula._y[i], molecula._z[i]);
            draw(atm, file, i + 1, i + 1);
        }
    } else {
        std::cerr << "Could not open " << out_file << ". " << '\n';
    }
    fmt::print(file, "END\n");
    std::fclose(file);

    return;
}

// Write grid matrix to PDB file.
void draw_PDB(GridMatrix const &mtx, std::string const &out_fil) {

    grid_t idx = 0;
    grid_t flat_idx = 0;

    FILE *file = std::fopen(out_fil.c_str(), "w");
    if (file) {
        for (grid_t k = 0; k < mtx._dimz; ++k) {
            for (grid_t j = 0; j < mtx._dimy; ++j) {
                for (grid_t i = 0; i < mtx._dimx; ++i) {

                    if (mtx._bool[flat_idx] == true) {
                        ++idx;
                        draw({i, j, k}, file, idx, 1, mtx._resolution,
                            mtx._origin);
                    }
                    ++flat_idx;
                }
            }
        }

    } else {
        std::cerr << "Could not open " << out_fil << ". " << '\n';
    }
    fmt::print(file, "END\n");
    std::fclose(file);

    return;
}

}