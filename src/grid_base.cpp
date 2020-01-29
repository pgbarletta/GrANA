#include "GrANA/grid_primitives.hpp"
namespace GrANA {

void draw(std::string const &out_fil,
    std::vector<std::vector<std::vector<grid_t>>> const &mtx,
    Vector const &orig_vtor, float const resolution) {
    grid_t const sz_z = mtx.size();
    grid_t const sz_y = mtx[0].size();
    grid_t const sz_x = mtx[0][0].size();

    grid_t idx = 0;
    FILE *file = std::fopen(out_fil.c_str(), "w");
    if (file) {
        for (grid_t k = 0; k < sz_z; ++k) {
            float const fz = grid_to_cont(k, resolution) + orig_vtor[2];
            for (grid_t j = 0; j < sz_y; ++j) {
                float const fy = grid_to_cont(j, resolution) + orig_vtor[1];
                for (grid_t i = 0; i < sz_x; ++i) {
                    if (mtx[k][j][i] != 1) {

                        ++idx;
                        float const fx =
                            grid_to_cont(i, resolution) + orig_vtor[0];
                        fmt::print(file,
                            "{: <6}{: >6} {: <4s} {:3} {:1}{: >3}    "
                            "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {: "
                            ">2s}\n",
                            "ATOM  ", idx, "H", "GPU", "A", 0, fx, fy, fz, 0.0,
                            0.0, "H");
                    }
                }
            }
        }
    } else {
        std::cerr << "Could not open " << out_fil << ". " << '\n';
    }
    std::fclose(file);
    return;
}
}