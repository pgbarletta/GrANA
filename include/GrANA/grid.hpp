#ifndef GrANA_GRID
#define GrANA_GRID

#include <string>
#include <fstream>
#include <iostream>
#include <array>
#include <fmt/format.h>
#include "GrANA/continuous.hpp"
#include "GrANA/utils.hpp"

namespace GrANA {

    class GridPoint {
    public:
        using grid_idx_t = int32_t;

        GridPoint() = default;

        GridPoint(const grid_idx_t x, const grid_idx_t y, const grid_idx_t z)
            noexcept : _xyz{x, y, z} {}

        GridPoint(const Point &p) noexcept :
            _xyz{cont_to_grid(p[0]), cont_to_grid(p[1]), cont_to_grid(p[2])} {}

        // Draw GridPoint as atom.
        void draw(FILE *ou_fil, int idx, int resid);

        // // Returns a vector starting on this point coordinates.
        // vector operator-(const point& p) const {
        //     return vector(_x - p[0], _y - p[1], _z - p[2],
        //         _x, _y, _z);
        // }
        // bool operator==(const point &p) const {
        //         return (_x == p[0] && _y == p[1] && _z == p[2]);
        // }
        grid_idx_t &operator[] (int idx) {
        	return _xyz[idx];
        }
        grid_idx_t operator[] (int idx) const {
        	return _xyz[idx];
        }
    private:
        std::array<grid_idx_t, 3> _xyz;
    };
    std::ostream& operator<<(std::ostream &stream, const GridPoint& t);
    // Turn a grid point into a continuous point, given the resolution.
    Point GridPoint_to_point(const GridPoint &in_point);

    // Turn a grid point into a continuous point, given the resolution.
    GridPoint point_to_GridPoint(const Point &in_point);

    class GridMolecule {
    public:
        // GridMolecule() noexcept : _natoms(0), _orig_vtor(Vector(0.0f, 0.0f, 0.0f)),
        //     _xyz(nullptr), _in_xyz(nullptr), _radii(nullptr), _in_radii(nullptr) {};
        GridMolecule() = default;
        GridMolecule(Molecule const &in_mol, Point const &orig_point);

        ~GridMolecule() {
            free(_xyz);
            free(_in_xyz);
            free(_radii);
            free(_in_radii);
        }

        void draw(std::string const &ou_fil);

        int _natoms{0};
        Vector _orig_vtor{0.0f, 0.0f, 0.0f};
        GridPoint *_xyz = nullptr, *_in_xyz = nullptr;
        float *_radii = nullptr, *_in_radii = nullptr;
    };

    class GridConvexHull {
    public:
        GridConvexHull() = default;
        GridConvexHull(ConvexHull const &CH);

        ~GridConvexHull() {
            free(_dots);
        }

        void draw(std::string const &ou_fil);
        
        float *_dots = nullptr;
    };

    class GridBool {
    public:
        GridBool() = default;
        GridBool(ConvexHull const &CH);

    };
} //namespace GrANA

#endif // GrANA_GRID
