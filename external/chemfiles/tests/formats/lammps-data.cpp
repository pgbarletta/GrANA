// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license
#include <fstream>

#include "catch.hpp"
#include "helpers.hpp"
#include "chemfiles.hpp"
using namespace chemfiles;

TEST_CASE("Read files in LAMMPS data format") {
    SECTION("File created by VMD/Topotools") {
        Trajectory file("data/lammps-data/solvated.lmp", 'r', "LAMMPS Data");
        Frame frame = file.read();

        CHECK(frame.size() == 7772);

        CHECK(frame.cell().shape() == UnitCell::ORTHORHOMBIC);
        CHECK(frame.cell().a() == 34.023997999999999);
        CHECK(frame.cell().b() == 34.023998000000006);
        CHECK(frame.cell().c() == 163.03599500000001);

        auto positions = frame.positions();
        CHECK(positions[0] == Vector3D(4.253000, 12.759000, 63.506001));
        CHECK(positions[364] == Vector3D(8.134000, 2.322000, 82.219002));
        CHECK(positions[653] == Vector3D(6.184000, 8.134000, 104.334000));

        auto& topology = frame.topology();
        CHECK(topology.bonds().size() == 6248);
        CHECK(topology[3].name() == "Zn");
        CHECK(topology[3].type() == "Zn");
        CHECK(topology[3].mass() == 65.408997);

        CHECK(topology[12].name() == "C1");
        CHECK(topology[12].type() == "C1");
        CHECK(topology[12].mass() == 42.0);

        // Check the read_step function
        frame = file.read_step(0);
        positions = frame.positions();
        CHECK(positions[0] == Vector3D(4.253000, 12.759000, 63.506001));
        CHECK(positions[364] == Vector3D(8.134000, 2.322000, 82.219002));
        CHECK(positions[653] == Vector3D(6.184000, 8.134000, 104.334000));
    }

    SECTION("File created with LAMMPS") {
        Trajectory file("data/lammps-data/data.body", 'r', "LAMMPS Data");
        Frame frame = file.read();

        CHECK(frame.size() == 100);
        CHECK(frame.cell() == UnitCell(31.064449134, 31.064449134, 1.0));

        auto positions = frame.positions();
        CHECK(positions[0] == Vector3D(-15.5322, -15.5322, 0.0));
        CHECK(positions[22] == Vector3D(-9.31933, -9.31933, 0.0));

        REQUIRE(frame.velocities());
        auto velocities = *frame.velocities();
        CHECK(velocities[5] == Vector3D(1.14438145745, 4.42784814304, 1.75516442452));
        CHECK(velocities[0] == Vector3D(1.02255489961, 2.92322463726, 4.88805110017));
        CHECK(velocities[1] == Vector3D(0.111646059519, 0.474226666855, 0.68604865644));
        CHECK(velocities[42] == Vector3D(4.70147770939, 2.13317266836, 1.29333445263));

        auto& topology = frame.topology();
        CHECK(topology.bonds().size() == 0);
        CHECK(topology[0].mass() == 6);
        CHECK(topology[1].mass() == 4);
        CHECK(topology[2].mass() == 3);

        CHECK(topology[12].name() == "1");
        CHECK(topology[12].type() == "1");
    }
}

TEST_CASE("Write files in LAMMPS data format") {
    const auto EXPECTED_CONTENT =
    "LAMMPS data file -- atom_style full -- generated by chemfiles\n"
    "6 atoms\n"
    "4 bonds\n"
    "4 angles\n"
    "2 dihedrals\n"
    "1 impropers\n"
    "4 atom types\n"
    "3 bond types\n"
    "3 angle types\n"
    "2 dihedral types\n"
    "1 improper types\n"
    "0 5 xlo xhi\n"
    "0 7 ylo yhi\n"
    "0 9 zlo zhi\n"
    "\n"
    "# Pair Coeffs\n"
    "# 1 As\n"
    "# 2 As\n"
    "# 3 B\n"
    "# 4 C\n"
    "\n"
    "# Bond Coeffs\n"
    "# 1 As-B\n"
    "# 2 B-B\n"
    "# 3 B-C\n"
    "\n"
    "# Angle Coeffs\n"
    "# 1 As-B-B\n"
    "# 2 As-B-C\n"
    "# 3 B-B-C\n"
    "\n"
    "# Dihedrals Coeffs\n"
    "# 1 As-B-B-C\n"
    "# 2 C-B-B-C\n"
    "\n"
    "# Impropers Coeffs\n"
    "# 1 As-B-B-C\n"
    "\n"
    "Masses\n"
    "\n"
    "1 25 # As\n"
    "2 74.9216 # As\n"
    "3 10.81 # B\n"
    "4 12.011 # C\n"
    "\n"
    "Atoms # full\n"
    "\n"
    "1 1 1 0 1.1 2.2 3.3 # As \n"
    "2 2 2 0 1.1 2.2 3.3 # As \n"
    "3 3 3 0 1.1 2.2 3.3 # B \n"
    "4 4 4 0 1.1 2.2 3.3 # C \n"
    "5 5 3 0 1.1 2.2 3.3 # B \n"
    "6 6 4 0 1.1 2.2 3.3 # C \n"
    "\n"
    "Velocities\n"
    "\n"
    "1 0.1 0.2 0.3\n"
    "2 0.1 0.2 0.3\n"
    "3 0.1 0.2 0.3\n"
    "4 0.1 0.2 0.3\n"
    "5 0.1 0.2 0.3\n"
    "6 0.1 0.2 0.3\n"
    "\n"
    "Bonds\n"
    "\n"
    "1 1 2 3\n"
    "2 3 3 4\n"
    "3 2 3 5\n"
    "4 3 5 6\n"
    "\n"
    "Angles\n"
    "\n"
    "1 2 2 3 4\n"
    "2 1 2 3 5\n"
    "3 3 3 5 6\n"
    "4 3 4 3 5\n"
    "\n"
    "Dihedrals\n"
    "\n"
    "1 1 2 3 5 6\n"
    "2 2 4 3 5 6\n"
    "\n"
    "Impropers\n"
    "\n"
    "1 1 2 3 4 5\n";



    auto tmpfile = NamedTempPath(".lmp");

    auto topology = Topology();
    topology.add_atom(Atom("As"));
    topology.add_atom(Atom("As"));
    topology.add_atom(Atom("B"));
    topology.add_atom(Atom("C"));
    topology.add_atom(Atom("B"));
    topology.add_atom(Atom("C"));
    topology.add_bond(2, 1);
    topology.add_bond(2, 3);
    topology.add_bond(2, 4);
    topology.add_bond(4, 5);

    topology[0].set_mass(25);

    auto frame = Frame(topology);
    frame.set_cell(UnitCell(5, 7, 9));

    frame.add_velocities();
    auto positions = frame.positions();
    auto velocities = *frame.velocities();
    for(size_t i=0; i<frame.size(); i++) {
        positions[i] = Vector3D(1.1, 2.2, 3.3);
        velocities[i] = Vector3D(0.1, 0.2, 0.3);
    }

    Trajectory(tmpfile, 'w', "LAMMPS Data").write(frame);

    std::ifstream checking(tmpfile);
    std::string content((std::istreambuf_iterator<char>(checking)),
                         std::istreambuf_iterator<char>());
    CHECK(content == EXPECTED_CONTENT);
}