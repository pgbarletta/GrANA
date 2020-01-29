#ifndef GrANA_GRID_BASE_H
#define GrANA_GRID_BASE_H

#include "GrANA/grid_primitives.hpp"

namespace GrANA {

void draw(std::string const &out_fil,
    std::vector<std::vector<std::vector<grid_t>>> const &mtx,
    Vector const &orig_vtor, float resolution);

} // namespace GrANA

#endif // GrANA_GRID
