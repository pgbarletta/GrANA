#ifndef GrANA_UTILS
#define GrANA_UTILS

#include <vector>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <cmath>

extern float rsltion;

namespace GrANA {

// Turn a coordinate in the matching matrix grid index.
inline int32_t cont_to_grid(float x) {
    return static_cast<int32_t>(fabs(x - fmod(x, rsltion)) / rsltion);
}

// Turn a grid index into a xyz coordinate.
inline float grid_to_cont(int32_t idx) {
    return static_cast<float>(idx * rsltion);
}

// Helper function for getting the indices that sort a cont_vector.
template <typename T>
    std::vector<int> sort_indices(const std::vector<T> &v);

// Helper function to get the indices of the true elements of a bool array.
// Optimized for large (>500) and sparse bool arrays.
void get_indices_from_sparse_bool_array(bool *in_array, const int n_in,
    std::vector<int> &indices);

} // namespace GrANA

#endif // GrANA_UTILS
