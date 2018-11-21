#ifndef GrANA_UTILS
#define GrANA_UTILS

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

extern float resolution;

namespace GrANA {

// Read input.
auto get_input(int argc, char **argv)
    -> std::tuple<std::string, std::string, float>;

// Helper function for getting the indices that sort a vector.
template <typename T>
auto sort_indices(std::vector<T> const &v) -> std::vector<int> {

    // initialize original index locations
    std::vector<int> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indices based on comparing values in v
    sort(
        idx.begin(), idx.end(), [&v](int i1, int i2) { return v[i1] < v[i2]; });
    return idx;
}

// Helper function to get the indices of the true elements of a bool array.
// Optimized for large (>500) and sparse bool arrays.
void get_indices_from_sparse_bool_array(
    bool *in_array, const int n_in, std::vector<int> &indices);

} // namespace GrANA

#endif // GrANA_UTILS
