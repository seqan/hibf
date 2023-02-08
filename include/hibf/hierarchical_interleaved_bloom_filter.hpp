// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <ranges>
#include <sstream>

#include <lemon/list_graph.h> /// Must be first include.

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

// layout
#include <seqan3/search/dream_index/detail/configuration.hpp>
#include <seqan3/search/dream_index/detail/data_store.hpp>
#include <seqan3/search/dream_index/detail/layout/execute.hpp>
#include <seqan3/search/dream_index/detail/layout/layout.hpp>

#include <seqan3/search/dream_index/detail/build/hibf/hierarchical_build.hpp>
#include <seqan3/search/dream_index/detail/build/hibf/read_chopper_pack_file.hpp>


namespace seqan3
{

template <typename input_range_type>
struct hibf_config
{
    static_assert(std::ranges::random_access_range<input_range_type>,
                  "The input data must be a std::random_access_range.");
    static_assert(std::ranges::input_range<std::ranges::range_value_t<std::ranges::range_value_t<input_range_type>>>,
                  "The input data must be a range of a range of a range, whose value type is an std::integral, "
                  "e.g. type `std::vector<std::vector<std::vector<size_t>>>`.");
    static_assert(std::integral<std::ranges::range_value_t<std::ranges::range_value_t<std::ranges::range_value_t<input_range_type>>>>,
                  "The input data must be a range of a range of a range, whose value type is an std::integral, "
                  "e.g. type `std::vector<std::vector<std::vector<size_t>>>`.");

    /*!\name General Configuration
     * \{
     */
    //!\brief A range over inputs. The most inner value type must be an integer (hash value);
    input_range_type input;

    //!\brief The number of hash functions for the IBFs.
    size_t number_of_hash_functions{2};

    //!\brief The desired false positive rate of the IBFs.
    double maximum_false_positive_rate{0.05};

    //!\brief The number of threads to use to compute merged HLL sketches.
    size_t threads{1u};
    //!\}

    /*!\name Layout Configuration
     * \{
     */
    //!\brief The number of bits the HyperLogLog sketch should use to distribute the values into bins.
    uint8_t sketch_bits{12};

    //!\brief The maximum number of technical bins on each IBF in the HIBF.
    uint16_t tmax{64};

    /*\brief A scaling factor to influence the amount of merged bins produced by the algorithm.
     *
     * The higher alpha, the more weight is added artificially to the low level IBFs and thus the optimal
     * solution will contain less merged bins in the end because it costs more to merge bins.
     */
    double alpha{1.2};

    //!\brief The maximal cardinality ratio in the clustering intervals.
    double max_rearrangement_ratio{0.5};

    //!\brief Whether to estimate the union of kmer sets to possibly improve the binning or not.
    bool estimate_union{true};

    //!\brief Whether to do a second sorting of the bins which takes into account similarity or not.
    bool rearrange_user_bins{true};
    //!\}

    /*!\name Build Configuration
     * \{
     */
    // Related to k-mers
    bool disable_cutoffs{false};

    // Related to IBF
    // bool compressed{false};
    //!\}
};

/*!\brief The HIBF binning directory. A data structure that efficiently answers set-membership queries for multiple
 *        bins.
 * \tparam data_layout_mode_ Indicates whether the underlying data type is compressed. See
 *                           [seqan3::data_layout](https://docs.seqan.de/seqan/3.0.3/group__submodule__dream__index.html#gae9cb143481c46a1774b3cdf5d9fdb518).
 * \see [seqan3::interleaved_bloom_filter][1]
 * \details
 *
 * This class improves the [seqan3::interleaved_bloom_filter][1] by adding additional bookkeeping that allows
 * to establish a hierarchical structure. This structure can then be used to split or merge user bins and distribute
 * them over a variable number of technical bins. In the [seqan3::interleaved_bloom_filter][1], the number of user bins
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
 * In constrast to the [seqan3::interleaved_bloom_filter][1], the user bins may be split across multiple technical bins
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
 * In contrast to the [seqan3::interleaved_bloom_filter][1], the result will consist of indices of user bins.
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
template <seqan3::data_layout data_layout_mode_ = seqan3::data_layout::uncompressed>
class hierarchical_interleaved_bloom_filter
{
public:
    // Forward declaration
    class user_bins;

    // Forward declaration
    class membership_agent;

    //!\brief Indicates whether the Interleaved Bloom Filter is compressed.
    static constexpr seqan3::data_layout data_layout_mode = data_layout_mode_;

    //!\brief The type of an individual Bloom filter.
    using ibf_t = seqan3::interleaved_bloom_filter<data_layout_mode_>;

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

    template <typename config_type>
    hierarchical_interleaved_bloom_filter(config_type const & configuration)
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
        return typename hierarchical_interleaved_bloom_filter<data_layout_mode>::membership_agent{*this};
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly.
     * \sa https://docs.seqan.de/seqan/3.2.0/group__io.html#serialisation
     */
    template <seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(ibf_vector);
        archive(next_ibf_id);
        archive(user_bins);
    }
    //!\endcond
protected:

    template <typename config_type>
    seqan3::hibf::layout compute_layout(config_type const & config)
    {
        seqan3::hibf::layout resulting_layout{};

        seqan3::hibf::configuration chopper_config
        {
            .sketch_bits = config.sketch_bits,
            .disable_sketch_output = true,
            .tmax = config.tmax,
            .num_hash_functions = config.number_of_hash_functions,
            .false_positive_rate = config.maximum_false_positive_rate,
            .alpha = config.alpha,
            .max_rearrangement_ratio = config.max_rearrangement_ratio,
            .threads = config.threads,
            .estimate_union = config.estimate_union,
            .rearrange_user_bins = config.rearrange_user_bins
        };

        // The output streams facilitate writing the layout file in hierarchical structure.
        // seqan3::hibf::execute currently writes the filled buffers to the output file.
        std::stringstream output_buffer;
        std::stringstream header_buffer;

        seqan3::hibf::data_store store{.false_positive_rate = chopper_config.false_positive_rate,
                                  .output_buffer = &output_buffer,
                                  .header_buffer = &header_buffer,
                                  .merged_bin_max_ids = &resulting_layout.merged_bin_max_ids};

        size_t const number_of_user_bins = std::ranges::size(config.input);

        // dummy init filenames
        store.filenames.resize(number_of_user_bins);
        for (size_t i = 0; i < number_of_user_bins; ++i)
            store.filenames[i] = "UB_" + std::to_string(i);

        // compute sketches
        store.sketches.resize(number_of_user_bins);
        store.kmer_counts.resize(number_of_user_bins);

// #pragma omp parallel for schedule(static) num_threads(config.threads)
        for (size_t i = 0; i < number_of_user_bins; ++i)
        {
            seqan3::sketch::hyperloglog sketch(config.sketch_bits);

            for (auto && hash_sequence : config.input[i]) // multi range input
                for (auto k_hash : hash_sequence)
                    sketch.add(reinterpret_cast<char *>(&k_hash), sizeof(k_hash));

// #pragma omp critical
            store.sketches[i] = sketch;
// #pragma omp critical
            store.kmer_counts[i] = sketch.estimate();
        }

        size_t const max_hibf_id = seqan3::hibf::execute(chopper_config, store);

        // brief Write the output to the layout file.
        std::stringstream fout{};
        write_layout_header_to(chopper_config, max_hibf_id, store.header_buffer->str(), fout);
        fout << store.output_buffer->str();

        resulting_layout.layout_str = fout.str();

        return resulting_layout; // return layout as string for now, containing the file
    }

    template <typename input_data_type>
    void build_index(hibf_config<input_data_type> const & config, seqan3::hibf::layout layout)
    {
        seqan3::hibf::build_data<data_layout_mode, hibf_config<input_data_type>> data{.hibf_config = config};
        data.hibf = this;

        seqan3::hibf::read_chopper_pack_file(data, layout.layout_str);
        lemon::ListDigraph::Node root = data.ibf_graph.nodeFromId(0); // root node = high level IBF node
        robin_hood::unordered_flat_set<size_t> root_kmers{};

        size_t const t_max{data.node_map[root].number_of_technical_bins};
        data.compute_fp_correction(t_max, config.number_of_hash_functions, config.maximum_false_positive_rate);

        seqan3::hibf::hierarchical_build(root_kmers, root, data, true);
    }
};

/*!\brief Bookkeeping for user and technical bins.
 */
template <seqan3::data_layout data_layout_mode>
class hierarchical_interleaved_bloom_filter<data_layout_mode>::user_bins
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

    //!\brief Returns a vector containing user bin indices for each bin in the `idx`th IBF.
    std::vector<int64_t> & bin_indices_of_ibf(size_t const idx)
    {
        return ibf_bin_to_filename_position[idx];
    }

    //!\brief Returns the filename of the `idx`th user bin.
    std::string & filename_of_user_bin(size_t const idx)
    {
        return user_bin_filenames[idx];
    }

    //!\brief For a pair `(a,b)`, returns a const reference to the filename of the user bin at IBF `a`, bin `b`.
    std::string const & operator[](std::pair<size_t, size_t> const & index_pair) const
    {
        return user_bin_filenames[ibf_bin_to_filename_position[index_pair.first][index_pair.second]];
    }

    /*!\brief Returns a view over the user bin filenames for the `ibf_idx`th IBF.
     *        An empty string is returned for merged bins.
     */
    auto operator[](size_t const ibf_idx) const
    {
        return ibf_bin_to_filename_position[ibf_idx]
             | std::views::transform(
                   [this](int64_t i)
                   {
                       if (i == -1)
                           return std::string{};
                       else
                           return user_bin_filenames[i];
                   });
    }

    //!\brief Returns the filename index of the `ibf_idx`th IBF for bin `bin_idx`.
    int64_t filename_index(size_t const ibf_idx, size_t const bin_idx) const
    {
        return ibf_bin_to_filename_position[ibf_idx][bin_idx];
    }

    /*!\brief Writes all filenames to a stream. Index and filename are tab-separated.
     * \details
     * 0	\<path_to_user_bin_0\>
     * 1	\<path_to_user_bin_1\>
     */
    template <typename stream_t>
    void write_filenames(stream_t & out_stream) const
    {
        size_t position{};
        std::string line{};
        for (auto const & filename : user_bin_filenames)
        {
            line.clear();
            line = '#';
            line += std::to_string(position);
            line += '\t';
            line += filename;
            line += '\n';
            out_stream << line;
            ++position;
        }
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
 * In contrast to the [seqan3::interleaved_bloom_filter][1], the result will consist of indices of user bins.
 */
template <seqan3::data_layout data_layout_mode> // TODO: value_t as template?
class hierarchical_interleaved_bloom_filter<data_layout_mode>::membership_agent
{
private:
    //!\brief The type of the augmented hierarchical_interleaved_bloom_filter.
    using hibf_t = hierarchical_interleaved_bloom_filter<data_layout_mode>;

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
            else if (bin + 1u == result.size() ||                                                    // last bin
                     current_filename_index != hibf_ptr->user_bins.filename_index(ibf_idx, bin + 1)) // end of split bin
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
    [[nodiscard]] std::vector<int64_t> const & bulk_contains(value_range_t && values, size_t const threshold) & noexcept
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

} // namespace seqan3
