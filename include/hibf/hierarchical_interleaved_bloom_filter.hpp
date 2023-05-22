// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <lemon/list_graph.h> /// Must be first include.

#include <ranges>
#include <sstream>

#include <hibf/config.hpp>
#include <hibf/contrib/robin_hood.hpp>
#include <hibf/detail/build/hibf/build_data.hpp>
#include <hibf/detail/build/hibf/compute_kmers.hpp>
#include <hibf/detail/build/hibf/construct_ibf.hpp>
#include <hibf/detail/build/hibf/hierarchical_build.hpp>
#include <hibf/detail/build/hibf/initialise_build_tree.hpp>
#include <hibf/detail/build/hibf/insert_into_ibf.hpp>
#include <hibf/detail/build/hibf/loop_over_children.hpp>
#include <hibf/detail/build/hibf/read_chopper_pack_file.hpp>
#include <hibf/detail/build/hibf/update_parent_kmers.hpp>
#include <hibf/detail/build/hibf/update_user_bins.hpp>
#include <hibf/detail/configuration.hpp>
#include <hibf/detail/data_store.hpp>
#include <hibf/detail/layout/compute_fp_correction.hpp>
#include <hibf/detail/layout/execute.hpp>
#include <hibf/detail/layout/layout.hpp>
#include <hibf/detail/sketch/estimate_kmer_counts.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>
#include <hibf/interleaved_bloom_filter.hpp>

namespace hibf
{

/*!\brief The HIBF binning directory. A data structure that efficiently answers set-membership queries for multiple
 *        bins.
 * \tparam data_layout_mode_ Indicates whether the underlying data type is compressed. See
 *                           [hibf::data_layout](https://docs.seqan.de/seqan/3.0.3/group__submodule__dream__index.html#gae9cb143481c46a1774b3cdf5d9fdb518).
 * \see [hibf::interleaved_bloom_filter][1]
 * \details
 *
 * This class improves the [hibf::interleaved_bloom_filter][1] by adding additional bookkeeping that allows
 * to establish a hierarchical structure. This structure can then be used to split or merge user bins and distribute
 * them over a variable number of technical bins. In the [hibf::interleaved_bloom_filter][1], the number of user bins
 * and technical bins is always the same. This causes performance degradation when there are many user bins or the user
 * bins are unevenly distributed.
 *
 * # Terminology
 *
 * ## Technical Bin
 * A Technical Bin represents an actual bin in the binning directory. In the IBF, it stores its kmers in a single Bloom
 * Filter (which is interleaved with all the other BFs).
 *
 * ## User Bin
 * The user may impose a structure on his sequence data in the form of logical groups (e.g. species). When querying the
 * IBF, the user is interested in an answer that differentiates between these groups.
 *
 * # Hierarchical Interleaved Bloom Filter (HIBF)
 *
 * In constrast to the [hibf::interleaved_bloom_filter][1], the user bins may be split across multiple technical bins
 * , or multiple user bins may be merged into one technical bin. When merging multiple user bins, the HIBF stores
 * another IBF that is built over the user bins constituting the merged bin. This lower-level IBF can then be used
 * to further distinguish between merged bins.
 *
 * In this example, user bin 1 was split into two technical bins. Bins 3, 4, and 5 were merged into a single technical
 * bin, and another IBF was added for the merged bin.
 * \image html hibf.svg width=90%
 *
 * The individual IBFs may have a different number of technical bins and differ in their sizes, allowing an efficient
 * distribution of the user bins.
 *
 * ## Querying
 * To query the Hierarchical Interleaved Bloom Filter for values, call
 * hibf::hierarchical_interleaved_bloom_filter::membership_agent() and use the returned
 * hibf::hierarchical_interleaved_bloom_filter::membership_agent.
 * In contrast to the [hibf::interleaved_bloom_filter][1], the result will consist of indices of user bins.
 *
 * To count the occurrences in each user bin of a range of values in the Hierarchical Interleaved Bloom Filter, call
 * hibf::hierarchical_interleaved_bloom_filter::counting_agent() and use
 * the returned hibf::hierarchical_interleaved_bloom_filter::counting_agent_type.
 *
 * ## Thread safety
 *
 * The Interleaved Bloom Filter promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 * [1]: https://docs.seqan.de/seqan/3.0.3/classseqan3_1_1interleaved__bloom__filter.html
 */
class hierarchical_interleaved_bloom_filter
{
public:
    /*!\brief Bookkeeping for user and technical bins.
    */
    class user_bins
    {
    private:
        //!\brief Contains filenames of all user bins.
        std::vector<std::string> user_bin_filenames;

        /*!\brief Stores for each bin in each IBF of the HIBF the ID of the filename.
        * \details
        * Assume we look up a bin `b` in IBF `i`, i.e. `ibf_bin_to_filename_position[i][b]`.
        * If `-1` is returned, bin `b` is a merged bin, and there is no filename, we need to look into the lower level IBF.
        * Otherwise, the returned value `j` can be used to access the corresponding filename `user_bin_filenames[j]`.
        */
        std::vector<std::vector<int64_t>> ibf_bin_to_filename_position{};

    public:
        //!\brief Returns the number of managed user bins.
        size_t num_user_bins() const noexcept
        {
            return user_bin_filenames.size();
        }

        //!\brief Changes the number of managed IBFs.
        void set_ibf_count(size_t const size)
        {
            ibf_bin_to_filename_position.resize(size);
        }

        //!\brief Changes the number of managed user bins.
        void set_user_bin_count(size_t const size)
        {
            user_bin_filenames.resize(size);
        }

        /*!\brief Returns a vector containing user bin indices for each bin in the `idx`th IBF.
        * \param idx The id of the x-th IBF.
        *
        * \details
        *
        * ### Example
        *
        * \include test/snippet/hibf/bin_indices_of_ibf.cpp
        */
        std::vector<int64_t> & bin_indices_of_ibf(size_t const idx)
        {
            return ibf_bin_to_filename_position[idx];
        }

        //!\brief Returns the filename index of the `ibf_idx`th IBF for bin `bin_idx`.
        int64_t filename_index(size_t const ibf_idx, size_t const bin_idx) const
        {
            return ibf_bin_to_filename_position[ibf_idx][bin_idx];
        }

        /*!\cond DEV
        * \brief Serialisation support function.
        * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
        * \param[in] archive The archive being serialised from/to.
        *
        * \attention These functions are never called directly.
        * \sa https://docs.seqan.de/seqan/3.2.0/group__io.html#serialisation
        */
        template <typename archive_t>
        void serialize(archive_t & archive)
        {
            archive(user_bin_filenames);
            archive(ibf_bin_to_filename_position);
        }
        //!\endcond
    };

    /*!\brief Manages membership queries for the hibf::hierarchical_interleaved_bloom_filter.
    * \see hibf::hierarchical_interleaved_bloom_filter::user_bins::filename_of_user_bin
    * \details
    * In contrast to the [hibf::interleaved_bloom_filter][1], the result will consist of indices of user bins.
    */
    class membership_agent
    {
    private:
        //!\brief The type of the augmented hierarchical_interleaved_bloom_filter.
        using hibf_t = hierarchical_interleaved_bloom_filter;

        //!\brief A pointer to the augmented hierarchical_interleaved_bloom_filter.
        hibf_t const * const hibf_ptr{nullptr};

        //!\brief Helper for recursive membership querying.
        template <std::ranges::forward_range value_range_t>
        void bulk_contains_impl(value_range_t && values, int64_t const ibf_idx, size_t const threshold)
        {
            auto agent = hibf_ptr->ibf_vector[ibf_idx].template counting_agent<uint16_t>();
            auto & result = agent.bulk_count(values);

            uint16_t sum{};

            for (size_t bin{}; bin < result.size(); ++bin)
            {
                sum += result[bin];

                auto const current_filename_index = hibf_ptr->user_bins.filename_index(ibf_idx, bin);

                if (current_filename_index < 0) // merged bin
                {
                    if (sum >= threshold)
                        bulk_contains_impl(values, hibf_ptr->next_ibf_id[ibf_idx][bin], threshold);
                    sum = 0u;
                }
                else if (bin + 1u == result.size() || // last bin
                         current_filename_index
                             != hibf_ptr->user_bins.filename_index(ibf_idx, bin + 1)) // end of split bin
                {
                    if (sum >= threshold)
                        result_buffer.emplace_back(current_filename_index);
                    sum = 0u;
                }
            }
        }

    public:
        /*!\name Constructors, destructor and assignment
        * \{
        */
        membership_agent() = default;                                     //!< Defaulted.
        membership_agent(membership_agent const &) = default;             //!< Defaulted.
        membership_agent & operator=(membership_agent const &) = default; //!< Defaulted.
        membership_agent(membership_agent &&) = default;                  //!< Defaulted.
        membership_agent & operator=(membership_agent &&) = default;      //!< Defaulted.
        ~membership_agent() = default;                                    //!< Defaulted.

        /*!\brief Construct a membership_agent for an existing hierarchical_interleaved_bloom_filter.
        * \private
        * \param hibf The hierarchical_interleaved_bloom_filter.
        */
        explicit membership_agent(hibf_t const & hibf) : hibf_ptr(std::addressof(hibf))
        {}
        //!\}

        //!\brief Stores the result of bulk_contains().
        std::vector<int64_t> result_buffer;

        /*!\name Lookup
        * \{
        */
        /*!\brief Determines set membership of given values, and returns the user bin indices of occurrences.
        * \param[in] values The values to process; must model std::ranges::forward_range.
        * \param[in] threshold Report a user bin if there are at least this many hits.
        *
        * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
        * \attention Sequential calls to this function invalidate the previously returned reference.
        *
        * \details
        *
        * ### Thread safety
        *
        * Concurrent invocations of this function are not thread safe, please create a
        * hibf::hierarchical_interleaved_bloom_filter::membership_agent for each thread.
        */
        template <std::ranges::forward_range value_range_t>
        [[nodiscard]] std::vector<int64_t> const & bulk_contains(value_range_t && values,
                                                                 size_t const threshold) & noexcept
        {
            assert(hibf_ptr != nullptr);

            static_assert(std::ranges::forward_range<value_range_t>, "The values must model forward_range.");
            static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                          "An individual value must be an unsigned integral.");

            result_buffer.clear();

            bulk_contains_impl(values, 0, threshold);

            std::ranges::sort(result_buffer); // TODO: necessary?

            return result_buffer;
        }

        // `bulk_contains` cannot be called on a temporary, since the object the returned reference points to
        // is immediately destroyed.
        template <std::ranges::range value_range_t>
        [[nodiscard]] std::vector<int64_t> const & bulk_contains(value_range_t && values,
                                                                 size_t const threshold) && noexcept = delete;
        //!\}
    };

    //!\brief The type of an individual Bloom filter.
    using ibf_t = hibf::interleaved_bloom_filter<hibf::data_layout::uncompressed>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    hierarchical_interleaved_bloom_filter() = default;                                              //!< Defaulted.
    hierarchical_interleaved_bloom_filter(hierarchical_interleaved_bloom_filter const &) = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter &
    operator=(hierarchical_interleaved_bloom_filter const &) = default;                        //!< Defaulted.
    hierarchical_interleaved_bloom_filter(hierarchical_interleaved_bloom_filter &&) = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter &
    operator=(hierarchical_interleaved_bloom_filter &&) = default; //!< Defaulted.
    ~hierarchical_interleaved_bloom_filter() = default;            //!< Defaulted.

    hierarchical_interleaved_bloom_filter(config const & configuration)
    {
        auto layout = compute_layout(configuration);
        build_index(configuration, layout);
    }
    //!\}

    //!\brief The individual interleaved Bloom filters.
    std::vector<ibf_t> ibf_vector;

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the next IBF.
     * \details
     * Assume we look up a bin `b` in IBF `i`, i.e. `next_ibf_id[i][b]`.
     * If `i` is returned, there is no lower level IBF, bin `b` is hence not a merged bin.
     * If `j != i` is returned, there is a lower level IBF, bin `b` is a merged bin, and `j` is the ID of the lower
     * level IBF in ibf_vector.
     */
    std::vector<std::vector<int64_t>> next_ibf_id;

    //!\brief The underlying user bins.
    user_bins user_bins;

    //!\brief Returns a membership_agent to be used for counting.
    membership_agent membership_agent() const
    {
        return typename hierarchical_interleaved_bloom_filter::membership_agent{*this};
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy hibf::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly.
     * \sa https://docs.seqan.de/seqan/3.2.0/group__io.html#serialisation
     */
    template <hibf::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(ibf_vector);
        archive(next_ibf_id);
        archive(user_bins);
    }
    //!\endcond

protected:
    hibf::layout::layout compute_layout(config const & hibf_config)
    {
        hibf::layout::layout resulting_layout{};

        hibf::configuration chopper_config{.sketch_bits = hibf_config.sketch_bits,
                                           .disable_sketch_output = true,
                                           .tmax = hibf_config.tmax,
                                           .num_hash_functions = hibf_config.number_of_hash_functions,
                                           .false_positive_rate = hibf_config.maximum_false_positive_rate,
                                           .alpha = hibf_config.alpha,
                                           .max_rearrangement_ratio = hibf_config.max_rearrangement_ratio,
                                           .threads = hibf_config.threads,
                                           .disable_estimate_union = hibf_config.disable_estimate_union,
                                           .disable_rearrangement = hibf_config.disable_rearrangement};

        // The output streams facilitate writing the layout file in hierarchical structure.
        // hibf::execute currently writes the filled buffers to the output file.
        std::stringstream output_buffer;
        std::stringstream header_buffer;

        std::vector<std::string> filenames{};
        std::vector<size_t> kmer_counts{};
        std::vector<sketch::hyperloglog> sketches{};

        // dummy init filenames
        filenames.resize(hibf_config.number_of_user_bins);
        for (size_t i = 0; i < hibf_config.number_of_user_bins; ++i)
            filenames[i] = "UB_" + std::to_string(i);

        // compute sketches
        sketches.resize(hibf_config.number_of_user_bins);
        kmer_counts.resize(hibf_config.number_of_user_bins);

        // #pragma omp parallel for schedule(static) num_threads(config.threads)
        robin_hood::unordered_flat_set<uint64_t> kmers;
        for (size_t i = 0; i < hibf_config.number_of_user_bins; ++i)
        {
            hibf::sketch::hyperloglog sketch(hibf_config.sketch_bits);

            kmers.clear();
            hibf_config.input_fn(i, std::inserter(kmers, kmers.begin()));

            for (auto k_hash : kmers)
                sketch.add(reinterpret_cast<char *>(&k_hash), sizeof(k_hash));

            // #pragma omp critical
            sketches[i] = sketch;
        }

        sketch::estimate_kmer_counts(sketches, kmer_counts);

        data_store store{.false_positive_rate = chopper_config.false_positive_rate,
                         .hibf_layout = &resulting_layout,
                         .kmer_counts = kmer_counts,
                         .sketches = sketches};

        size_t const max_hibf_id = hibf::execute(chopper_config, store);
        store.hibf_layout->top_level_max_bin_id = max_hibf_id;

        return *store.hibf_layout; // return layout as string for now, containing the file
    }

    void build_index(config const & hibf_config, hibf::layout::layout & hibf_layout)
    {
        size_t const number_of_ibfs = hibf_layout.max_bins.size() + 1;

        this->ibf_vector.resize(number_of_ibfs);
        this->user_bins.set_ibf_count(number_of_ibfs);
        this->user_bins.set_user_bin_count(hibf_layout.user_bins.size());
        this->next_ibf_id.resize(number_of_ibfs);

        build_data data{.hibf_config = hibf_config};

        initialise_build_tree(hibf_layout, data.ibf_graph, data.node_map);
        lemon::ListDigraph::Node root_node = data.ibf_graph.nodeFromId(0); // root node = top-level IBF node

        size_t const t_max{data.node_map[root_node].number_of_technical_bins};
        data.fp_correction = layout::compute_fp_correction(hibf_config.maximum_false_positive_rate,
                                                           hibf_config.number_of_hash_functions,
                                                           t_max);

        hierarchical_build(root_node, data);
    }

private:
    size_t hierarchical_build(robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                              lemon::ListDigraph::Node const & current_node,
                              build_data & data,
                              bool is_root)
    {
        auto & current_node_data = data.node_map[current_node];

        size_t const ibf_pos{data.request_ibf_idx()};

        std::vector<int64_t> ibf_positions(current_node_data.number_of_technical_bins, ibf_pos);
        std::vector<int64_t> filename_indices(current_node_data.number_of_technical_bins, -1);
        robin_hood::unordered_flat_set<uint64_t> kmers{};

        auto initialise_max_bin_kmers = [this](robin_hood::unordered_flat_set<uint64_t> & kmers,
                                               std::vector<int64_t> & ibf_positions,
                                               std::vector<int64_t> & filename_indices,
                                               lemon::ListDigraph::Node const & node,
                                               build_data & data) -> size_t
        {
            auto & node_data = data.node_map[node];

            if (node_data.favourite_child != lemon::INVALID) // max bin is a merged bin
            {
                // recursively initialize favourite child first
                ibf_positions[node_data.max_bin_index] =
                    hierarchical_build(kmers, node_data.favourite_child, data, false);
                return 1;
            }
            else // max bin is not a merged bin
            {
                // we assume that the max record is at the beginning of the list of remaining records.
                auto const & record = node_data.remaining_records[0];
                compute_kmers(kmers, data, record);
                update_user_bins(filename_indices, record);

                return record.number_of_technical_bins;
            }
        };

        // initialize lower level IBF
        size_t const max_bin_tbs = initialise_max_bin_kmers(kmers, ibf_positions, filename_indices, current_node, data);
        auto && ibf = construct_ibf(parent_kmers, kmers, max_bin_tbs, current_node, data, is_root);
        kmers.clear(); // reduce memory peak

        // parse all other children (merged bins) of the current ibf
        auto loop_over_children = [this](robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                                         hibf::interleaved_bloom_filter<> & ibf,
                                         std::vector<int64_t> & ibf_positions,
                                         lemon::ListDigraph::Node const & current_node,
                                         build_data & data,
                                         bool is_root)
        {
            auto & current_node_data = data.node_map[current_node];
            std::vector<lemon::ListDigraph::Node> children{};

            for (lemon::ListDigraph::OutArcIt arc_it(data.ibf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
                children.emplace_back(data.ibf_graph.target(arc_it));

            if (children.empty())
                return;

            size_t const number_of_mutex = (data.node_map[current_node].number_of_technical_bins + 63) / 64;
            std::vector<std::mutex> local_ibf_mutex(number_of_mutex);

            auto worker = [this,
                           &children,
                           &current_node_data,
                           &parent_kmers,
                           &data,
                           &local_ibf_mutex,
                           &ibf_positions,
                           &ibf,
                           &is_root](auto && index, auto &&)
            {
                auto & child = children[index];

                if (child != current_node_data.favourite_child)
                {
                    robin_hood::unordered_flat_set<uint64_t> kmers{};
                    size_t const ibf_pos = hierarchical_build(kmers, child, data, false);
                    auto parent_bin_index = data.node_map[child].parent_bin_index;
                    {
                        size_t const mutex_id{parent_bin_index / 64};
                        std::lock_guard<std::mutex> guard{local_ibf_mutex[mutex_id]};
                        ibf_positions[parent_bin_index] = ibf_pos;
                        insert_into_ibf(kmers, 1, parent_bin_index, ibf, data.fill_ibf_timer);
                        if (!is_root)
                            update_parent_kmers(parent_kmers, kmers, data.merge_kmers_timer);
                    }
                }
            };

            size_t number_of_threads{};
            std::vector<size_t> indices(children.size());
            std::iota(indices.begin(), indices.end(), size_t{});

            if (is_root)
            {
                // Shuffle indices: More likely to not block each other. Optimal: Interleave
                std::shuffle(indices.begin(), indices.end(), std::mt19937_64{std::random_device{}()});
                number_of_threads = data.hibf_config.threads;
            }
            else
            {
                number_of_threads = 1u;
            }

            seqan3::detail::execution_handler_parallel executioner{number_of_threads};
            executioner.bulk_execute(std::move(worker), std::move(indices), []() {});
        };

        loop_over_children(parent_kmers, ibf, ibf_positions, current_node, data, is_root);

        // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
        size_t const start{(current_node_data.favourite_child != lemon::INVALID) ? 0u : 1u};
        for (size_t i = start; i < current_node_data.remaining_records.size(); ++i)
        {
            auto const & record = current_node_data.remaining_records[i];

            if (is_root && record.number_of_technical_bins == 1) // no splitting needed
            {
                insert_into_ibf(data, record, ibf);
            }
            else
            {
                compute_kmers(kmers, data, record);
                insert_into_ibf(kmers, record.number_of_technical_bins, record.storage_TB_id, ibf, data.fill_ibf_timer);
                if (!is_root)
                    update_parent_kmers(parent_kmers, kmers, data.merge_kmers_timer);
            }

            update_user_bins(filename_indices, record);
            kmers.clear();
        }

        this->ibf_vector[ibf_pos] = std::move(ibf);
        this->next_ibf_id[ibf_pos] = std::move(ibf_positions);
        this->user_bins.bin_indices_of_ibf(ibf_pos) = std::move(filename_indices);

        return ibf_pos;
    }

    size_t hierarchical_build(lemon::ListDigraph::Node const & root_node, build_data & data)
    {
        robin_hood::unordered_flat_set<uint64_t> root_kmers{};
        return hierarchical_build(root_kmers, root_node, data, true);
    }
};

} // namespace hibf
