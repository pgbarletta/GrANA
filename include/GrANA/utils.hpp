#ifndef GrANA_UTILS
#define GrANA_UTILS

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace GrANA {

// Read input.
auto get_input(int argc, char **argv)
    -> std::tuple<std::string, float, float, std::string>;

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

//
template <typename Container_u, typename T_u = typename Container_u::value_type,
    typename Container_v, typename T_v = typename Container_v::value_type>
auto reorder(Container_u const &u, Container_v const &v) -> Container_u {

    auto const sz = u.size();
    if (sz != v.size()) {
        throw std::runtime_error(
            "reorder() fail. u and v don't have same length.\n");
    }
    Container_u w;
    w.reserve(sz);
    for (size_t i = 0; i != sz; ++i) {
        w.push_back(u[v[i]]);
    }

    return w;
}

} // namespace GrANA

#endif // GrANA_UTILS
