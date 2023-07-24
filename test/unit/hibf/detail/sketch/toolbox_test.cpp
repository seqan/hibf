#include <gtest/gtest.h> // for Message, TestPartResult, EXPECT_EQ, TestInfo, TEST_F, Test

#include <cinttypes>   // for uint64_t
#include <cstddef>     // for size_t
#include <functional>  // for greater
#include <limits>      // for numeric_limits
#include <string>      // for allocator, basic_string, string
#include <string_view> // for string_view
#include <tuple>       // for tie, make_tuple
#include <vector>      // for vector

#include <hibf/contrib/robin_hood.hpp>        // for unordered_flat_map
#include <hibf/detail/sketch/hyperloglog.hpp> // for hyperloglog
#include <hibf/detail/sketch/toolbox.hpp>     // for clustering_node, entry, precompute_union_estimates_for, cluste...
#include <hibf/test/expect_range_eq.hpp>      // for expect_range_eq, EXPECT_RANGE_EQ

// inherits from toolbox to test private members
struct toolbox_test : public ::testing::Test
{
    std::vector<std::string> test_filenames{"small.fa", "small.fa", "small2.fa", "small2.fa"};
    std::vector<size_t> test_kmer_counts{500, 600, 700, 800};
    std::vector<size_t> test_positions{0, 1, 2, 3};
    std::vector<hibf::sketch::hyperloglog> test_sketches = [this]()
    {
        std::vector<hibf::sketch::hyperloglog> result(test_kmer_counts.size());

        std::vector<std::string> const small_input{
            {"ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"
             "AGCTATTTCAGACCTACACTATCTAGCTTATTCACAAATATTATAACGGCATACGTCTAGTGCTCATCGTGATCTAGCGA"
             "GCTAGCGATCTGATTCACGAGCGTACGTGACGTACGTATCGTACTACGTATCGTACTACATGCATCGATCGACGTAGCTA"
             "TCAGCGTAGCGTACGAGTCAGCTGACTGACGTCGTAGCATCGTACGTAGCGTAGCGATCGAGTCACTTATCGTAGCTAGT"
             "CGACTAGCGTACGTAGTCAGCTATTATGACGAGGCGACTTAGCGACTACGAGCTAGCGAGGAGGCGAGGCGAGCGGACTG"},
            {"ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"
             "AGCTATTTCAGACCTACACTATCTAGCTTATTCACAAATATTATAACGGCATACGTCTAGTGCTCATCGTGATCTAGCGA"
             "ATATCGATCGAGCGAGGCAGGCAGCGATCGAGCGAGCGCATGCAGCGACTAGCTACGACAGCTACTATCAGCAGCGAGCG"
             "GCTAGCGATCTGATTCACGAGCGTACGTGACGTACGTATCGTACTACGTATCGTACTACATGCATCGATCGACGTAGCTA"
             "TCAGCGTAGCGTACGAGTCAGCTGACTGACGTCGTAGCATCGTACGTAGCGTAGCGATCGAGTCACTTATCGTAGCTAGT"
             "CGACTAGCGTACGTAGTCAGCTATTATGACGAGGCGACTTAGCGACTACGAGCTAGCGAGGAGGCGAGGCGAGCGGACTG"},
            {"ACGATCGACTAGGAGCGATTACGACTGACTACATCTAGCTAGCTAGAGATTCTTCAGAGCTTAGCGATCTCGAGCTATCG"
             "AGCTATTTCAGACCTACACTATCTAGCTTATTCACAAATATTATAACGGCATACGTCTAGTGCTCATCGTGATCTAGCGA"
             "GCTAGCGATCTGATTCACGAGCGTACGTGACGTACGTATCGTACTACGTATCGTACTACATGCATCGATCGACGTAGCTA"
             "ATCGATCACGATCAGCGAGCGATATCTTATCGTAGGCATCGAGCATCGAGGAGCGATCTATCTATCTATCATCTATCTAT"
             "TCAGCGTAGCGTACGAGTCAGCTGACTGACGTCGTAGCATCGTACGTAGCGTAGCGATCGAGTCACTTATCGTAGCTAGT"
             "CGACTAGCGTACGTAGTCAGCTATTATGACGAGGCGACTTAGCGACTACGAGCTAGCGAGGAGGCGAGGCGAGCGGACTG"
             "G"}};

        for (std::string_view seq : small_input)
        {
            // we have to go C-style here for the HyperLogLog Interface
            char const * c_seq_it = seq.begin();
            char const * end = c_seq_it + seq.size();

            while (c_seq_it + 21 <= end)
            {
                result[0].add(c_seq_it, 21);
                ++c_seq_it;
            }
        }

        result[1] = result[0];
        result[2] = result[0];
        result[3] = result[0];

        // hibf::sketch::toolbox::read_hll_files_into(data(""), test_filenames, result);
        return result;
    }();
};

TEST_F(toolbox_test, sort_by_cardinalities)
{
    hibf::sketch::toolbox::sort_by_cardinalities(test_sketches, test_kmer_counts, test_positions);

    // filenames do not change
    EXPECT_RANGE_EQ(test_filenames, (std::vector<std::string>{"small.fa", "small.fa", "small2.fa", "small2.fa"}));
    EXPECT_RANGE_EQ(test_kmer_counts, (std::vector<size_t>{500, 600, 700, 800}));
    // only positions change
    EXPECT_RANGE_EQ(test_positions, (std::vector<size_t>{3, 2, 1, 0}));
}

TEST_F(toolbox_test, precompute_union_estimates_for)
{
    std::vector<uint64_t> estimates(4);

    hibf::sketch::toolbox::precompute_union_estimates_for(estimates,
                                                          test_sketches,
                                                          test_kmer_counts,
                                                          test_positions,
                                                          0);
    EXPECT_RANGE_EQ(estimates, (std::vector<uint64_t>{500, 0, 0, 0}));

    hibf::sketch::toolbox::precompute_union_estimates_for(estimates,
                                                          test_sketches,
                                                          test_kmer_counts,
                                                          test_positions,
                                                          1);
    EXPECT_RANGE_EQ(estimates, (std::vector<uint64_t>{658, 600, 0, 0}));

    hibf::sketch::toolbox::precompute_union_estimates_for(estimates,
                                                          test_sketches,
                                                          test_kmer_counts,
                                                          test_positions,
                                                          2);
    EXPECT_RANGE_EQ(estimates, (std::vector<uint64_t>{658, 658, 700, 0}));

    hibf::sketch::toolbox::precompute_union_estimates_for(estimates,
                                                          test_sketches,
                                                          test_kmer_counts,
                                                          test_positions,
                                                          3);
    EXPECT_RANGE_EQ(estimates, (std::vector<uint64_t>{658, 658, 658, 800}));
}

TEST_F(toolbox_test, random_shuffle)
{
    hibf::sketch::toolbox::prio_queue default_pq{};
    hibf::sketch::toolbox::distance_matrix dist{{0, default_pq},
                                                {1, default_pq},
                                                {2, default_pq},
                                                {3, default_pq},
                                                {4, default_pq}};
    robin_hood::unordered_flat_map<size_t, size_t> ids{{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}};

    hibf::sketch::toolbox::random_shuffle(dist, ids);

    // since randomness is seeded, the output is deterministic
    auto [new_pos_0, new_pos_1, new_pos_2, new_pos_3, new_pos_4] = std::make_tuple(3u, 2u, 1u, 0u, 4u);

    EXPECT_EQ(ids[0], new_pos_0);
    EXPECT_EQ(ids[1], new_pos_1);
    EXPECT_EQ(ids[2], new_pos_2);
    EXPECT_EQ(ids[3], new_pos_3);
    EXPECT_EQ(ids[4], new_pos_4);

    EXPECT_EQ(dist[new_pos_0].id, 0u);
    EXPECT_EQ(dist[new_pos_1].id, 1u);
    EXPECT_EQ(dist[new_pos_2].id, 2u);
    EXPECT_EQ(dist[new_pos_3].id, 3u);
    EXPECT_EQ(dist[new_pos_4].id, 4u);
}

TEST_F(toolbox_test, prune)
{
    hibf::sketch::toolbox::prio_queue default_pq{};
    hibf::sketch::toolbox::distance_matrix dist{{0, default_pq},
                                                {1, default_pq},
                                                {2, default_pq},
                                                {3, default_pq},
                                                {4, default_pq}};
    robin_hood::unordered_flat_map<size_t, size_t> remaining_ids{{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}};

    // since remaining_ids contains all_ids, prune shouldn't do anything. All ids are valid.
    hibf::sketch::toolbox::prune(dist, remaining_ids);

    EXPECT_EQ(remaining_ids[0], 0u);
    EXPECT_EQ(remaining_ids[1], 1u);
    EXPECT_EQ(remaining_ids[2], 2u);
    EXPECT_EQ(remaining_ids[3], 3u);
    EXPECT_EQ(remaining_ids[4], 4u);

    EXPECT_EQ(dist.size(), 5u);
    EXPECT_EQ(dist[0].id, 0u);
    EXPECT_EQ(dist[1].id, 1u);
    EXPECT_EQ(dist[2].id, 2u);
    EXPECT_EQ(dist[3].id, 3u);
    EXPECT_EQ(dist[4].id, 4u);

    remaining_ids.erase(1);
    remaining_ids.erase(3);

    // distance entry 1 and 3 are now invalid, since they do not occur in remaining_ids
    // prune() should therefore remove them from dist.
    hibf::sketch::toolbox::prune(dist, remaining_ids);

    EXPECT_EQ(remaining_ids[0], 0u);
    EXPECT_EQ(remaining_ids[2], 2u);
    EXPECT_EQ(remaining_ids[4], 1u);

    EXPECT_EQ(dist.size(), 3u);
    EXPECT_EQ(dist[0].id, 0u);
    EXPECT_EQ(dist[1].id, 4u);
    EXPECT_EQ(dist[2].id, 2u);
}

TEST_F(toolbox_test, rotate)
{
    hibf::sketch::hyperloglog s{5}; // default sketch for every entry in the tree as it is not important for rotate
    auto f = std::numeric_limits<size_t>::max();

    /* test clustering tree
     * The root is at position 0. 'f' means infinity.
     *             (5,6)
     *            /     \
     *        (0,1)     (2,3)
     *       /    \     /    \
     *   (f,f)  (f,f) (f,f) (f,f) the leaves are the UBs to be clustered
     */
    std::vector<hibf::sketch::toolbox::clustering_node> clustering{{f, f, s},
                                                                   {f, f, s},
                                                                   {f, f, s},
                                                                   {f, f, s}, // the leaves come first
                                                                   {5, 6, s},
                                                                   {0, 1, s},
                                                                   {2, 3, s}};

    // previous_rightmost is already at the very left. Nothing has to be rotated.
    hibf::sketch::toolbox::rotate(clustering, 0 /*previous_rightmost*/, 0 /*interval_start*/, 4 /*root_id*/);

    EXPECT_EQ(std::tie(clustering[0].left, clustering[0].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[1].left, clustering[1].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[2].left, clustering[2].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[3].left, clustering[3].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[4].left, clustering[4].right), std::make_tuple(5u, 6u));
    EXPECT_EQ(std::tie(clustering[5].left, clustering[5].right), std::make_tuple(0u, 1u));
    EXPECT_EQ(std::tie(clustering[6].left, clustering[6].right), std::make_tuple(2u, 3u));

    // now the previous_rightmost is within the tree. Rotation should take place
    hibf::sketch::toolbox::rotate(clustering, 2 /*previous_rightmost*/, 0 /*interval_start*/, 4 /*root_id*/);

    EXPECT_EQ(std::tie(clustering[0].left, clustering[0].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[1].left, clustering[1].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[2].left, clustering[2].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[3].left, clustering[3].right), std::tie(f, f));
    EXPECT_EQ(std::tie(clustering[4].left, clustering[4].right), std::make_tuple(6u, 5u));
    EXPECT_EQ(std::tie(clustering[5].left, clustering[5].right), std::make_tuple(0u, 1u));
    EXPECT_EQ(std::tie(clustering[6].left, clustering[6].right), std::make_tuple(2u, 3u));
}

TEST_F(toolbox_test, trace)
{
    hibf::sketch::hyperloglog s{5}; // default sketch for every entry in the tree as it is not important for rotate
    auto f = std::numeric_limits<size_t>::max();

    /* test clustering tree
     * The root is at position 0. 'f' means infinity.
     *             (5,6)
     *            /     \
     *        (1,3)     (2,0)
     *       /    \     /    \
     *   (f,f)  (f,f) (f,f) (f,f) the leaves are the UBs to be clustered
     */
    std::vector<hibf::sketch::toolbox::clustering_node> clustering{{f, f, s},
                                                                   {f, f, s},
                                                                   {f, f, s},
                                                                   {f, f, s}, // the leaves come first
                                                                   {5, 6, s},
                                                                   {1, 3, s},
                                                                   {2, 0, s}};

    std::vector<size_t> permutation{};

    hibf::sketch::toolbox::trace(clustering,
                                 permutation,
                                 2 /*previous_rightmost*/,
                                 0 /*interval_start*/,
                                 4 /*root_id*/);

    EXPECT_RANGE_EQ(permutation, (std::vector<size_t>{1, 3, 0}));
}

TEST_F(toolbox_test, cluster_bins)
{
    { // whole range
        std::vector<size_t> permutation{};
        hibf::sketch::toolbox::cluster_bins(test_sketches,
                                            test_kmer_counts,
                                            test_positions,
                                            permutation,
                                            0 /*interval start*/,
                                            3 /*interval_end*/,
                                            1 /*number of threads*/);
        // index 3 is not part of current permutation so it can participate in "the next interval"
        EXPECT_RANGE_EQ(permutation, (std::vector<size_t>{2, 0, 1}));
    }

    { // intervals
        std::vector<size_t> permutation{};
        hibf::sketch::toolbox::cluster_bins(test_sketches,
                                            test_kmer_counts,
                                            test_positions,
                                            permutation,
                                            0 /*interval start*/,
                                            1 /*interval_end*/,
                                            1 /*number of threads*/);
        EXPECT_RANGE_EQ(permutation, (std::vector<size_t>{0}));
        hibf::sketch::toolbox::cluster_bins(test_sketches,
                                            test_kmer_counts,
                                            test_positions,
                                            permutation,
                                            1 /*interval start*/,
                                            3 /*interval_end*/,
                                            1 /*number of threads*/);
        EXPECT_RANGE_EQ(permutation, (std::vector<size_t>{0, 1, 2}));
    }
}
