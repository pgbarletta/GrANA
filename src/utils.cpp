#include "GrANA/utils.hpp"

namespace GrANA {

// Read input.
auto get_input(int argc, char **argv)
    -> std::tuple<std::string, float, uint32_t, uint32_t, std::string> {

    float resolution = 1.f;
    uint32_t margin = 0, gap_depth = 0;
    try {
        if (argc != 6) {
            throw std::invalid_argument(
                "Usage: GrANA in_pdb resolution margin gap_depth out_pdb\n");
        }

        resolution = std::stof(argv[2]);
        margin = std::stoi(argv[3]);
        gap_depth = std::stoi(argv[4]);

    } catch (...) { std::terminate(); }

    return {argv[1], resolution, margin, gap_depth, argv[5]};
}

// Helper function to get the indices of the true elements of a bool array.
// Optimized for large (>500) and sparse bool arrays.
auto get_indices_from_sparse_bool_array(
    bool const in_array[], const int in_size) -> std::vector<int> {

    long const *cin_array = (long *)in_array;
    int const sz = in_size / sizeof(long);
    int const resto = in_size % sizeof(long);

    std::vector<int> indices;
    int const sum = std::accumulate(in_array, in_array + in_size, 0);
    indices.reserve(sum);

    for (int i = 0; i < sz; ++i) {
        if (cin_array[i] != 0) {
            size_t const lo = i * sizeof(long);
            size_t const hi = lo + sizeof(long);
            for (size_t j = lo; j < hi; ++j) {
                if (in_array[j] != 0) {
                    indices.push_back(j);
                }
            }
        }
    }

    for (int i = in_size - resto; i < in_size; ++i) {
        if (in_array[i] != 0) {
            indices.push_back(i);
        }
    }

    return indices;
}
} // namespace GrANA
