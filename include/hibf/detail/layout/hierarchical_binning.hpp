#pragma once

#include <cassert>
#include <cmath>

#include <hibf/detail/configuration.hpp>
#include <hibf/detail/layout/simple_binning.hpp>
#include <hibf/detail/prefixes.hpp>
#include <hibf/next_multiple_of_64.hpp>

namespace hibf::layout
{

class hierarchical_binning
{
private:
    //!\brief The user configuration passed down from the command line.
    configuration const config{};
    //!\brief Stores all data that is needed to compute the layout, e.g. the counts, sketches and the layout::layout.
    data_store * const data{nullptr};

    //!\brief The number of user bins, initialised with the length of user_bin_kmer_counts.
    size_t const num_user_bins{};
    //!\brief The number of technical bins requested by the user.
    size_t const num_technical_bins{};

public:
    hierarchical_binning() = default;                                        //!< Defaulted.
    hierarchical_binning(hierarchical_binning const &) = delete;             //!< Deleted. Would modify same data.
    hierarchical_binning & operator=(hierarchical_binning const &) = delete; //!< Deleted. Would modify same data.
    hierarchical_binning(hierarchical_binning &&) = default;                 //!< Defaulted.
    hierarchical_binning & operator=(hierarchical_binning &&) = delete; //!< Deleted. config has no move assignment.
    ~hierarchical_binning() = default;                                  //!< Defaulted.

    /*!\brief The constructor from user bin names, their kmer counts and a configuration.
     * \param[in, out] data_ Stores all data that is needed to compute the layout.
     * \param[in] config_ A configuration object that holds information from the user that influence the computation.
     *
     * Each entry in the names_ and input vector respectively is considered a user bin (both vectors must have the
     * same length).
     */
    hierarchical_binning(data_store & data_, configuration const & config_) :
        config{config_},
        data{std::addressof(data_)},
        num_user_bins{data->positions.size()},
        num_technical_bins{data->previous.empty() ? config.tmax : needed_technical_bins(num_user_bins)}
    {
        assert(data != nullptr);
    }

    //!\brief Executes the hierarchical binning algorithm and layouts user bins into technical bins.
    size_t execute();

private:
    /*!\brief Returns the number of technical bins given a number of user bins.
     * \param[in] requested_num_ub The number of user bins.
     */
    [[nodiscard]] size_t needed_technical_bins(size_t const requested_num_ub) const
    {
        return std::min<size_t>(next_multiple_of_64(requested_num_ub), config.tmax);
    }

    /*!\brief Returns the maximum number of needed levels when merging `num_ubs_in_merge` many user bins.
     * \param[in] num_ubs_in_merge The number of user bins in the merge.
     */
    [[nodiscard]] size_t max_merge_levels(size_t const num_ubs_in_merge) const
    {
        size_t const lower_lvl_tbs = needed_technical_bins(num_ubs_in_merge);
        double const levels = std::log(num_ubs_in_merge) / std::log(lower_lvl_tbs);
        return static_cast<size_t>(std::ceil(levels));
    }

    /*!\brief Initialize the matrices M (hig_level_ibf), L (low_level_ibfs) and T (trace)
     *
     * \image html hierarchical_dp_init.png
     */
    void initialization(std::vector<std::vector<size_t>> & matrix,
                        std::vector<std::vector<size_t>> & ll_matrix,
                        std::vector<std::vector<std::pair<size_t, size_t>>> & trace);

    /*!\brief Performs the recursion.
     *
     * \image html hierarchical_dp_recursion.png
     *
     * Explanations to the formula:
     *
     * Remember that M (matrix) stores the maximum technical bin size of the high level IBF (HIBF) and
     * L (ll_matrix) stores the sum of all estimate low level ibfs (LIBF) memory footprints
     * (assuming perfect splitting of kmer_content which can obviously not be achieved but `alpha` may be adjusted
     * in order to counteract this underestimation).
     *
     * Now in order to minimize the memory footprint we...
     *
     * 1. ... (\f$v_ij\f$) firstly check if we should split the bin. <br>
     * Therefore we compute for every possible \f$ i' \f$ above \f$ i \f$, the technical bin size if we split
     * \f$ c_j \f$ into \f$ i - i' \f$ bins (\f$ \frac{c_j}{i - i'} \f$). We only take the maximum of the current
     * maximum where I come from (\f$ M_{i',j-1} \f$) and the new technical bin size computed just now.
     * This maximum is the new current *maximum technical bin size* that needs to be multiplied by the current number
     * technical bins \f$ (i + 1) \f$ in order to estimate the memory footprint of the HIBF.
     * Now that we have the memory footprint of the HIBF we also consider the LIBFs memory footprint but nothing
     * changed here, since we split not merge, so we just take \f$ L_{i',j-1} \f$ scaled by alpha.
     *
     * 2. ... (\f$h_ij\f$) secondly check if we should merge the bin with the ones before.<br>
     * Therefore we start by compute for every possible \f$ j' \f$ to the left of \f$ j \f$, the merged bin weight
     * (\f$ \sum_{g = j'}^{j} c_g \f$) of merging together user bins \f$ [j', ...,  j] \f$. This is only possible
     * iff every bin \f$ [j', ...,  j] \f$ was **not splitted**. If we merge those bins starting with user bin
     * \f$ j' \f$ we start the trace at \f$ M_{i-1,j'-1} \f$. Therefore we need to compute the new maximal technical
     * bin size of the HIBF, which is the maximum of the merged bin weight and where we would come from
     * (\f$ \max(M_{i-1,j'-1}, \sum_{g = j'}^{j} c_g) \f$) and multiple it by the number of technical bins so far
     * \f$ (i + 1) \f$ to get the HIBF memory footprint. The LIBFs memory footprint also changes since we introduce a
     * new merged bin. Namely, we add the weight of the new merged bin (\f$ \sum_{g = j'}^{j} c_g \f$) again to the
     * LIBFs memory footprint from where we would come from \f$ L_{i-1,j'-1} \f$ and scale this by alpha. Just adding
     * the combined merged bin weight neglects the fact, that the merged bin weight has to be distributed within the
     * new low level IBF, potentially causing the effective text ratio and thereby the IBF memory footprint to increase.
     * This we cannot know before hand how the data are, we need to accept this as a flaw in the "optimumal result" of
     * this algorithm. It would be too computational intensive to compute the splitting for every possibility.
     *
     */
    void recursion(std::vector<std::vector<size_t>> & matrix,
                   std::vector<std::vector<size_t>> & ll_matrix,
                   std::vector<std::vector<std::pair<size_t, size_t>>> & trace);

    //!\brief Backtracks the trace matrix and writes the resulting binning into the output file.
    size_t backtracking(std::vector<std::vector<std::pair<size_t, size_t>>> const & trace);

    std::string to_string_with_precision(double const value) const
    {
        // TODO std::to_chars after https://github.com/seqan/product_backlog/issues/396
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << value;
        return stream.str();
    };

    data_store initialise_libf_data(size_t const trace_j) const
    {
        data_store libf_data{.false_positive_rate = data->false_positive_rate,
                             .hibf_layout = data->hibf_layout,
                             .kmer_counts = data->kmer_counts,
                             .sketches = data->sketches,
                             .positions = {data->positions[trace_j]},
                             .fp_correction = data->fp_correction};

        return libf_data;
    }

    void process_merged_bin(data_store & libf_data, size_t const bin_id) const
    {
        update_libf_data(libf_data, bin_id);

        // now do the binning for the low-level IBF:
        size_t const lower_max_bin = add_lower_level(libf_data);

        data->hibf_layout->max_bins.emplace_back(libf_data.previous.bin_indices, lower_max_bin);
    }

    void update_libf_data(data_store & libf_data, size_t const bin_id) const
    {
        bool const is_top_level = data->previous.empty();

        libf_data.previous = data->previous;
        libf_data.previous.bin_indices.push_back(bin_id);

#if HIBF_WORKAROUND_GCC_BOGUS_MEMCPY
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wrestrict"
#endif // HIBF_WORKAROUND_GCC_BOGUS_MEMCPY

        libf_data.previous.num_of_bins += (is_top_level ? "" : ";") + std::string{"1"};

#if HIBF_WORKAROUND_GCC_BOGUS_MEMCPY
#    pragma GCC diagnostic pop
#endif // HIBF_WORKAROUND_GCC_BOGUS_MEMCPY
    }

    size_t add_lower_level(data_store & libf_data) const
    {
        // now do the binning for the low-level IBF:
        if (libf_data.positions.size() > config.tmax)
        {
            // recursively call hierarchical binning if there are still too many UBs
            return hierarchical_binning{libf_data, config}.execute(); // return id of maximum technical bin
        }
        else
        {
            // use simple binning to distribute remaining UBs
            return simple_binning{libf_data, 0, config.debug}.execute(); // return id of maximum technical bin
        }
    }

    void update_max_id(size_t & max_id, size_t & max_size, size_t const new_id, size_t const new_size) const
    {
        if (new_size > max_size)
        {
            max_id = new_id;
            max_size = new_size;
        }
    }
};

} // namespace hibf::layout
