// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <iostream>
#include <random>

#include <sharg/parser.hpp>

#include <hibf/build/bin_size_in_bits.hpp>
#include <hibf/contrib/robin_hood.hpp>
#include <hibf/interleaved_bloom_filter.hpp>
#include <hibf/layout/compute_fpr_correction.hpp>
#include <hibf/misc/divide_and_ceil.hpp>

struct config
{
    // Set by user
    size_t kmer_size{12};
    size_t elements{130};
    size_t hash{2};
    size_t repetitions{1};
    double fpr{0.05};

    // Internal
    size_t kmer_max{};
};

class positive_integer_validator
{
public:
    using option_value_type = size_t;

    positive_integer_validator() = default;
    positive_integer_validator(positive_integer_validator const &) = default;
    positive_integer_validator & operator=(positive_integer_validator const &) = default;
    positive_integer_validator(positive_integer_validator &&) = default;
    positive_integer_validator & operator=(positive_integer_validator &&) = default;
    ~positive_integer_validator() = default;

    explicit positive_integer_validator(bool const is_zero_positive_) : is_zero_positive{is_zero_positive_}
    {}

    void operator()(option_value_type const & val) const
    {
        if (!is_zero_positive && !val)
            throw sharg::validation_error{"The value must be a positive integer."};
    }

    std::string get_help_page_message() const
    {
        if (is_zero_positive)
            return "Value must be a positive integer or 0.";
        else
            return "Value must be a positive integer.";
    }

private:
    bool is_zero_positive{false};
};

void init_parser(sharg::parser & parser, config & cfg)
{
    parser.add_option(cfg.kmer_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer",
                                    .description = "The k-mer size.",
                                    .validator = sharg::arithmetic_range_validator{1, 32}});
    parser.add_option(cfg.elements,
                      sharg::config{.short_id = '\0',
                                    .long_id = "elements",
                                    .description = "Number of elements to insert.",
                                    .validator = positive_integer_validator{}});
    parser.add_option(cfg.hash,
                      sharg::config{.short_id = '\0',
                                    .long_id = "hash",
                                    .description = "The number of hash functions to use.",
                                    .validator = sharg::arithmetic_range_validator{1, 5}});
    parser.add_option(cfg.repetitions,
                      sharg::config{.short_id = '\0',
                                    .long_id = "repeats",
                                    .description = "Number of repetitions.",
                                    .validator = positive_integer_validator{}});
    parser.add_option(cfg.fpr,
                      sharg::config{.short_id = '\0',
                                    .long_id = "fpr",
                                    .description = "The desired false positive rate.",
                                    .validator = sharg::arithmetic_range_validator{0.0, 1.0}});
}

double false_positive_rate(size_t const elements, size_t const hash_count, size_t const bin_size)
{
    double const exp_arg = (hash_count * elements) / static_cast<double>(bin_size);
    double const log_arg = 1.0 - std::exp(-exp_arg);
    return std::exp(hash_count * std::log(log_arg));
}

struct stats
{
    size_t min{};
    size_t max{};
    double mean{};
    double variance{};
};

stats get_stats(std::vector<size_t> const & values)
{
    size_t min = std::ranges::min(values);
    size_t max = std::ranges::max(values);

    size_t number_of_values = values.size();
    double mean = std::reduce(values.begin(), values.end()) / static_cast<double>(number_of_values);
    double variance = [&values, mean, number_of_values]()
    {
        if (number_of_values == 1u)
            return std::numeric_limits<double>::quiet_NaN();

        auto helper = [mean, number_of_values](double current, double const value)
        {
            return current + ((value - mean) * (value - mean) / (number_of_values - 1));
        };

        return std::reduce(values.begin(), values.end(), 0.0, helper);
    }();

    return stats{.min = min, .max = max, .mean = mean, .variance = variance};
}

stats get_stats(std::vector<size_t> const & values_, size_t const elements)
{
    std::vector<size_t> values;
    values.reserve(values_.size());
    std::ranges::transform(values_,
                           std::back_inserter(values),
                           [elements](size_t const val)
                           {
                               return elements - val;
                           });

    return get_stats(values);
}

void run(config const & cfg)
{
    size_t const bin_size = seqan::hibf::build::bin_size_in_bits({.fpr = cfg.fpr, //
                                                                  .hash_count = cfg.hash,
                                                                  .elements = cfg.elements});

    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{1u},
                                              seqan::hibf::bin_size{bin_size},
                                              seqan::hibf::hash_function_count{cfg.hash},
                                              true};

    std::vector<size_t> occupancies{};
    occupancies.reserve(cfg.repetitions);

    // Simulation
    {
        auto agent = ibf.containment_agent();

        robin_hood::unordered_set<uint64_t> inserted_values{};
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<uint64_t> distrib(0ULL, cfg.kmer_max);

        for (size_t i{}; i < cfg.repetitions; ++i)
        {
            for (; inserted_values.size() < cfg.elements;)
                inserted_values.emplace(distrib(gen));

            for (uint64_t const value : inserted_values)
                ibf.emplace(value, seqan::hibf::bin_index{0u});

            occupancies.push_back(ibf.occupancy[0]);
            inserted_values.clear();
            ibf.clear(seqan::hibf::bin_index{0u});
        }
    }

    stats const occupancy_stats = get_stats(occupancies);
    stats const collision_stats = get_stats(occupancies, cfg.elements);

    // I don't know why this works reasonably well, but it does.
    // Removing `/ std::sqrt(cfg.hash)` and using `1` hash is a known computation for the expected number of
    // hash collisions.
    // Hence, it is exact for 1 hash, but only a rough estimate for more than 1 hash.
    double const rough_estimate = [&]()
    {
        double const bin_size_1_hash = seqan::hibf::build::bin_size_in_bits({.fpr = cfg.fpr, //
                                                                             .hash_count = 1u,
                                                                             .elements = cfg.elements});

        double estimate_1 = std::exp(cfg.elements * std::log((bin_size_1_hash - 1.0) / bin_size_1_hash));
        double estimate_2 = bin_size_1_hash * (1.0 - estimate_1);
        return (cfg.elements - estimate_2) / std::sqrt(cfg.hash);
    }();

    // This one does make sense, but is rather expensive to compute.
    double const exact_estimate = [&]()
    {
        double estimate = 0.0;
        for (size_t i = 1; i <= cfg.elements; ++i)
            estimate += false_positive_rate(i - 1, cfg.hash, bin_size);
        return estimate;
    }();

    std::stringstream result{};
    result.setf(std::ios::left, std::ios::adjustfield);

    result.width(20);
    result << "";
    result.width(11);
    result << "min";
    result.width(11);
    result << "max";
    result.width(11);
    result << "mean";
    result.width(11);
    result << "variance";
    result << '\n';

    result.width(20);
    result << "occupancy";
    result.width(11);
    result << occupancy_stats.min;
    result.width(11);
    result << occupancy_stats.max;
    result.width(11);
    result << occupancy_stats.mean;
    result.width(11);
    result << occupancy_stats.variance;
    result << '\n';

    result.width(20);
    result << "collisions";
    result.width(11);
    result << collision_stats.min;
    result.width(11);
    result << collision_stats.max;
    result.width(11);
    result << collision_stats.mean;
    result.width(11);
    result << collision_stats.variance;
    result << '\n';

    result.width(20);
    result << "exact estimate";
    result.width(11);
    result << "";
    result.width(11);
    result << "";
    result.width(11);
    result << exact_estimate;
    result.width(11);
    result << "";
    result << '\n';

    result.width(20);
    result << "rough estimate";
    result.width(11);
    result << "";
    result.width(11);
    result << "";
    result.width(11);
    result << rough_estimate;
    result.width(11);
    result << "";
    result << '\n';

    std::cout << result.str();
}

int main(int argc, char ** argv)
{
    sharg::parser parser{"hash_collisions", argc, argv, sharg::update_notifications::off};
    parser.info.author = "Enrico Seiler";
    parser.info.short_copyright = "BSD 3-Clause License";
    parser.info.short_description = "Determines the number of hash collisions in the IBF.";
    config cfg{};
    init_parser(parser, cfg);
    parser.parse();

    cfg.kmer_max = cfg.kmer_size == 32 ? std::numeric_limits<uint64_t>::max() : (1ULL << (2 * cfg.kmer_size)) - 1u;

    if (cfg.elements > cfg.kmer_max && cfg.elements != cfg.kmer_max + 1u)
    {
        std::cout << "[WARNING] Inserting more elements than there are possible k-mers. "
                  << "Setting number of elements to number of possible k-mers.\n\n";
        cfg.elements = cfg.kmer_max + 1u;
    }

    std::cout << "kmer: " << cfg.kmer_size << '\n';
    std::cout << "elements: " << cfg.elements << '\n';
    std::cout << "hash: " << cfg.hash << '\n';
    std::cout << "fpr: " << cfg.fpr << '\n';
    std::cout << "repetitions: " << cfg.repetitions << "\n\n";
    run(cfg);
}
