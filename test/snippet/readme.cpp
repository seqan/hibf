// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <cinttypes>  // for uint64_t, int64_t
#include <cstddef>    // for size_t
#include <functional> // for function
#include <iostream>   // for operator<<, basic_ostream, cout
#include <ranges>     // for __fn, iota, views
#include <vector>     // for vector

#include <hibf/config.hpp>                                // for insert_iterator, config
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter

int main()
{
    // Let's say we have groups that have data that we find interesting.
    // For example, each file of the RefSeq data set could be such a group.
    // In the context of the HIBF, we call such groups user bins.

    // Given a query, we want to quickly determine which user bins this query is likely to occur in.
    // This is also called Approximate Membership Query (AMQ).

    // In this example, we have three user bins. Each of these user bins is characterized by a range of
    // unsigned integer values. Some popular techniques for obtaining such unsigned integers from
    // biological sequences include k-mers, minimisers, and syncmers.

    // For clarity, we show each user bin individually before copying them to user_bin_data.
    std::vector<uint64_t> user_bin_1{1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u};
    std::vector<uint64_t> user_bin_2{1u, 2u, 3u, 4u, 5u};
    std::vector<uint64_t> user_bin_3{3u, 9u, 11u};
    std::vector<std::vector<uint64_t>> user_bin_data{user_bin_1, user_bin_2, user_bin_3};

    // The HIBF uses a config. There are two required options:
    // 1) The number of user bins: 3 (user_bin_data.size())
    // 2) A function to access the input data.
    //    The signature is (size_t const user_bin_id, seqan::hibf::insert_iterator it). You need to
    //    provide the function body, and the hibf lib will use this function to access the data of each
    //    user bin. When this function is called by the library with a specific user_bin_id, all
    //    unsigned integer values (data) belonging to this user bin have to be assigned to the
    //    seqan::hibf::insert_iterator.
    //    Conveniently, this function can be a lambda, and hence capture data outside the function body.
    auto get_user_bin_data = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        for (auto value : user_bin_data[user_bin_id])
            it = value;
    };

    // Now we can construct a config, any other settings are optional. We have included some interesting
    // settings with their respective default values here.
    seqan::hibf::config config{.input_fn = get_user_bin_data, // required
                               .number_of_user_bins = 3u,     // required
                               .number_of_hash_functions = 2u,
                               .maximum_fpr = 0.05,
                               .threads = 1u};

    // The HIBF constructor will determine a hierarchical layout for the user bins and build the filter.
    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};

    // Now we can search for some query.
    std::vector<uint64_t> query1{3u, 9u, 12u, 14u};

    // For this, we use the membership agent of the HIBF. This agent only needs to be created once and
    // can be reused for multiple subsequent queries.
    // If you are using multiple threads in your app, each thread should have its own membership agent.
    auto agent = hibf.membership_agent();

    // The membership_for function takes the query and a threshold. Here, a threshold of two means that
    // at least (>=) 2 values of the query must be found within a user bin to be a hit.
    // While exact thresholds can be obtained for some approaches such as k-mers, another popular
    // approach is to require at least x% of the values in the query to hit.
    // For example, a threshold of 2 equals 40% of the values in query1 (5 values).
    // This threshold needs to be provided by the user. In general, some care should be taken with the
    // threshold. A low threshold requires a traversal of more parts of the hierarchy and slows down
    // the search.
    // Note that we bind the result with a `&` to avoid copies!
    auto & result1 = agent.membership_for(query1, 2u);

    // query1 hits in user_bin_1 and user_bin_3, which have the IDs 0 and 2, respectively.
    for (int64_t hit_user_bin : result1)
        std::cout << hit_user_bin << ' '; // The results are not sorted: 2 0
    std::cout << '\n';

    // Another query.
    // A query is simply a range of unsigned integer values, e.g., it does not have to be a vector.
    auto query2 = std::views::iota(0u, 15u); // 0,1,2,...,14
    auto & result2 = agent.membership_for(query2, 5u);
    agent.sort_results(); // Sort the results.

    // query2 hits in user_bin_1 and user_bin_2, which have the IDs 0 and 1, respectively.
    for (int64_t hit_user_bin : result2)
        std::cout << hit_user_bin << ' '; // The results are sorted: 0 1
    std::cout << '\n';
}
