#pragma once

#include <filesystem>
#include <fstream>
#include <omp.h>
#include <queue>
#include <random>

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/detail/sketch/hyperloglog.hpp>

namespace hibf::sketch
{

class toolbox
{
public:
    //!\brief type for a node in the clustering tree when for the rearrangement
    struct clustering_node
    {
        // children in the tree
        size_t left;
        size_t right;
        // hll sketch of the union if the node is still a root
        hyperloglog hll;
    };

    //!\brief element of the second priority queue layer of the distance matrix
    struct neighbor
    {
        size_t id;
        double dist;

        bool operator>(neighbor const & other) const
        {
            return dist > other.dist;
        }
    };

    //!\brief type of a min heap based priority queue
    using prio_queue = std::priority_queue<neighbor, std::vector<neighbor>, std::greater<neighbor>>;

    //!\brief entry of the distance matrix that has the id of a cluster with its neighbors in a prio queue
    struct entry
    {
        size_t id;
        prio_queue pq;
    };

    //!\brief type of the distance matrix for the clustering for the rearrangement
    using distance_matrix = std::vector<entry>;

    //!\brief A pointer to kmer counts associated with the above files used to layout user bin into technical bins.
    std::vector<size_t> const & kmer_counts;

    //!\brief HyperLogLog sketches on the k-mer sets of the sequences from the files of filenames.
    std::vector<hyperloglog> const & sketches;

    //!\brief HyperLogLog sketches on the k-mer sets of the sequences from the files of filenames.
    std::vector<size_t> & positions;

public:
    toolbox() = default;                            //!< Defaulted.
    toolbox(toolbox const &) = default;             //!< Defaulted.
    toolbox & operator=(toolbox const &) = default; //!< Defaulted.
    toolbox(toolbox &&) = default;                  //!< Defaulted.
    toolbox & operator=(toolbox &&) = default;      //!< Defaulted.
    ~toolbox() = default;                           //!< Defaulted.

    /*!\brief A sequence of user bins for which filenames and counts are given.
     * \param[in] kmer_counts_ counts of the k-mer sets of the bins corresponding to filenames
     * \param[in] sketches_ sketches of the bins corresponding to filenames
     * \param[in] positions_ The realtive positions of the input information to correctly access sketches and counts.
     */
    toolbox(std::vector<size_t> const & kmer_counts_,
            std::vector<hyperloglog> const & sketches_,
            std::vector<size_t> & positions_) :
        kmer_counts{kmer_counts_},
        sketches{sketches_},
        positions{positions_}
    {}

    //!\brief Sorts filenames and cardinalities by looking only at the cardinalities.
    void sort_by_cardinalities()
    {
        assert(positions.size() <= kmer_counts.size());

        auto cardinality_compare = [this](size_t const index1, size_t const index2)
        {
            return kmer_counts[index1] > kmer_counts[index2];
        };

        std::sort(positions.begin(), positions.end(), cardinality_compare);
    }

    /*!\brief Restore the HLL sketches from the files in hll_dir and target_filenames into target container.
    */
    static void read_hll_files_into(std::filesystem::path const & hll_dir,
                                    std::vector<std::string> const & target_filenames,
                                    std::vector<hyperloglog> & target)
    {
        assert(std::filesystem::exists(hll_dir) && !std::filesystem::is_empty(hll_dir)); // checked in chopper_layout

        target.reserve(target_filenames.size());

        try
        {
            for (auto const & filename : target_filenames)
            {
                std::filesystem::path path = hll_dir / std::filesystem::path(filename).stem();
                path += ".hll";
                std::ifstream hll_fin(path, std::ios::binary);

                if (!hll_fin.good())
                    throw std::runtime_error{"Could not open file " + path.string()};

                // the sketch bits will be automatically read from the files
                target.emplace_back().restore(hll_fin);
            }
        }
        catch (std::runtime_error const & err)
        {
            std::string const chopper_msg{"[CHOPPER LAYOUT ERROR] Something went wrong trying to read the HyperLogLog"
                                          " sketches from files:\n"};
            throw std::runtime_error{chopper_msg + err.what()};
        }
    }

    /*!\brief Estimate the cardinality of the union for a single user bin j with all prior ones j' < j.
     * \param[out] estimates output row
     * \param[in] sketches The hyperloglog sketches of the respective user bins.
     * \param[in] counts The counts/sketch.estimates() of the respective user bins.
     * \param[in] positions The realtive positions of the input information to correctly access sketches and counts.
     * \param[in] j The current user bin (column in the DP matrix)
     *
     * estimates[j_prime] will be the union cardinality estimate of the interval {j_prime, ..., j}.
     */
    static void precompute_union_estimates_for(std::vector<uint64_t> & estimates,
                                               std::vector<hyperloglog> const & sketches,
                                               std::vector<size_t> const & counts,
                                               std::vector<size_t> const & positions,
                                               int64_t const j)
    {
        assert(counts.size() == sketches.size());
        assert(positions.size() <= counts.size());
        assert(estimates.size() == sketches.size()); // Resize happens in precompute_init_interval_union_estimations
        assert(estimates.size() > static_cast<size_t>(j));

        hyperloglog temp_hll = sketches[positions[j]];
        estimates[j] = counts[positions[j]];

        for (int64_t j_prime = j - 1; j_prime >= 0; --j_prime)
            estimates[j_prime] = static_cast<uint64_t>(temp_hll.merge_and_estimate_SIMD(sketches[positions[j_prime]]));
    }

    /*!\brief Estimate the cardinality of the union for each interval [0, j] for all user bins j.
     * \param[out] estimates output row
     * \param[in] sketches The hyperloglog sketches of the respective user bins.
     * \param[in] counts The counts/sketch.estimates() of the respective user bins.
     * \param[in] positions The realtive positions of the input information to correctly access sketches and counts.
     *
     * estimates[j] will be the union cardinality estimate of the interval {0, ..., j}.
     */
    static void precompute_initial_union_estimates(std::vector<uint64_t> & estimates,
                                                   std::vector<hyperloglog> const & sketches,
                                                   std::vector<size_t> const & counts,
                                                   std::vector<size_t> const & positions)
    {
        assert(counts.size() == sketches.size());
        assert(positions.size() <= counts.size());
        assert(sketches.size() > 0u);

        estimates.resize(sketches.size());

        hyperloglog temp_hll = sketches[positions[0]];
        estimates[0] = counts[positions[0]];

        for (size_t j = 1; j < positions.size(); ++j)
            estimates[j] = static_cast<uint64_t>(temp_hll.merge_and_estimate_SIMD(sketches[positions[j]]));
    }

    /*!\brief Estimate the cardinality of the union for a single interval.
     * \param[in] sketches The hyperloglog sketches to be used for estimation.
     * \param[in] positions The realtive positions of the input information to correctly access sketches and counts.
     * \returns The the cardinality of the union for the interval [start, end).
     */
    static uint64_t estimate_interval(std::vector<hyperloglog> const & sketches, std::vector<size_t> const & positions)
    {
        assert(positions.size() <= sketches.size());
        assert(!positions.empty());

        hyperloglog temp_hll = sketches[positions[0]];

        for (size_t i = 1; i < positions.size(); ++i)
            temp_hll.merge(sketches[positions[i]]);

        return temp_hll.estimate();
    }

    /*!\brief Rearrange filenames, sketches and counts such that similar bins are close to each other
     * \param[in] max_ratio the maximal cardinality ratio in the clustering intervals (must be <= 1 and >= 0)
     * \param[in] num_threads the number of threads to use
     */
    void rearrange_bins(double const max_ratio, size_t const num_threads)
    {

        std::vector<size_t> permutation;

        size_t first = 0;
        size_t last = 1;

        while (first < positions.size())
        {
            // size difference is too large or sequence is over -> do the clustering
            if (last == positions.size() || kmer_counts[positions[first]] * max_ratio > kmer_counts[positions[last]])
            {
                // if this is not the first group, we want one bin overlap
                cluster_bins(permutation, first, last, num_threads);
                first = last;
            }
            ++last;
        }

        for (size_t i{0}; i < permutation.size(); ++i)
        {
            size_t swap_index = permutation[i];
            while (swap_index < i)
                swap_index = permutation[swap_index];

            std::swap(positions[i], positions[swap_index]);
        }
    }

    /*!\brief Perform an agglomerative clustering variant on the index range [first:last)
     * \param[in] first id of the first cluster of the interval
     * \param[in] last id of the last cluster of the interval plus one
     * \param[in] num_threads the number of threads to use
     * \param[out] permutation append the new order to this
     */
    void
    cluster_bins(std::vector<size_t> & permutation, size_t const first, size_t const last, size_t const num_threads);

    /*!\brief Randomly swap entries in dist while keeping track of the changes of indices.
     * \param[in] dist the distance matrix (vector of priority queues) to shuffle
     * \param[in] remaining_ids the map with information about which ids remain at which index
     */
    void random_shuffle(distance_matrix & dist, robin_hood::unordered_flat_map<size_t, size_t> & remaining_ids) const
    {
        size_t const n = dist.size();

        std::mt19937_64 gen(0x7E1E5665D46800E5ULL);

        for (size_t i = 0; i < n - 1; ++i)
        {
            std::uniform_int_distribution<size_t> distrib(i, n - 1);
            size_t const swap_i = distrib(gen);

            size_t const id = dist[i].id;
            size_t const swap_id = dist[swap_i].id;

            // swap entries and update the reming ids, because the indices in dist changed
            std::swap(dist[i], dist[swap_i]);
            std::swap(remaining_ids[id], remaining_ids[swap_id]);
        }
    }

    /*!\brief Delete inactive entries out of dist and shrink to fit its size while keeping track of the changes of indices
     * \param[in] dist the distance matrix (vector of priority queues) to prune
     * \param[in] remaining_ids the map with information about which ids remain at which index
     */
    void prune(distance_matrix & dist, robin_hood::unordered_flat_map<size_t, size_t> & remaining_ids) const
    {
        if (dist.empty())
            return;

        // index of the first entry after the valid range
        size_t valid_range_end = 0;
        // index of the first entry before the invalid range
        size_t invalid_range_start = dist.size() - 1;

        while (valid_range_end != invalid_range_start)
        {
            if (remaining_ids.contains(dist[valid_range_end].id))
            {
                ++valid_range_end;
            }
            else if (!remaining_ids.contains(dist[invalid_range_start].id))
            {
                --invalid_range_start;
            }
            else
            {
                // If we arrive here, then valid_range_end has an invalid id
                // and invalid_range_start has a valid id. The correspoding entries should be swapped
                std::swap(dist[valid_range_end], dist[invalid_range_start]);

                // update the index of the valid entry
                remaining_ids.at(dist[valid_range_end].id) = valid_range_end;
            }
        }

        // check the last element between the valid and invalid range
        if (remaining_ids.contains(dist[valid_range_end].id))
        {
            ++valid_range_end;
        }

        // cut off invalid values
        dist.resize(valid_range_end);
    }

    /*!\brief Rotate the previous rightmost bin to the left of the clustering tree
     * \param[in, out] clustering the tree to do the rotation on
     * \param[in] previous_rightmost the id of the node to be rotated to the left
     * \param[in] first the id of the first node in the interval to shift the index
     * \param[in] id the id of the current node
     *
     * If called with the root of the tree, this function recursively calls itself while rotating
     * several subtrees until previous_rightmost is at the very left end of the whole clustering tree.
     *
     * \return whether previous rightmost was in the subtree rooted at id
     */
    bool rotate(std::vector<clustering_node> & clustering,
                size_t const previous_rightmost,
                size_t const first,
                size_t const id) const
    {
        if (id == previous_rightmost) // we are at the leaf that is previous_rightmost (Anchor of the recursion)
            return true;

        clustering_node & curr = clustering[id - first];

        if (curr.left == std::numeric_limits<size_t>::max()) // we are at a leaf that is not previous_rightmost
        {
            return false;
        }
        // nothing to do if previous_rightmost is in the left subtree
        else if (rotate(clustering, previous_rightmost, first, curr.left))
        {
            return true;
        }
        // rotate if previous_rightmost is in the right subtree
        else if (rotate(clustering, previous_rightmost, first, curr.right))
        {
            std::swap(curr.left, curr.right);
            return true;
        }

        // else: previous_rightmost is not in this subtree
        return false;
    }

    /*!\brief Do a recursive traceback to find the order of leaves in the clustering tree
     * \param[in] clustering the tree to do the traceback on
     * \param[out] permutation append the new order to this
     * \param[in] previous_rightmost the id of the node on the left which should be ignored
     * \param[in] first the id of the first node in the interval to shift the index
     * \param[in] id the id of the current node
     *
     * This function traverses the tree in depth-first-search accessing the leaves from left to right.
     * 'Left to right' refers to the order of nodes in `clustering`.
     */
    void trace(std::vector<clustering_node> const & clustering,
               std::vector<size_t> & permutation,
               size_t const previous_rightmost,
               size_t const first,
               size_t const id) const
    {
        clustering_node const & curr = clustering[id - first];

        if (curr.left == std::numeric_limits<size_t>::max()) // I am at a leaf
        {
            if (id != previous_rightmost)
                permutation.push_back(id);
            return;
        }

        trace(clustering, permutation, previous_rightmost, first, curr.left);
        trace(clustering, permutation, previous_rightmost, first, curr.right);
    }
};

} // namespace hibf::sketch
