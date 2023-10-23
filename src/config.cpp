// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cassert>     // for assert
#include <cmath>       // for ceil, sqrt
#include <functional>  // for function
#include <iostream>    // for basic_ostream, operator<<, basic_istream, stringstream, cerr, istream
#include <sstream>     // for basic_stringstream
#include <stdexcept>   // for invalid_argument
#include <string>      // for char_traits, getline, operator<<, string
#include <string_view> // for operator<<, operator==, basic_string_view

#include <hibf/config.hpp>                   // for config
#include <hibf/layout/prefixes.hpp>          // for meta_header, meta_hibf_config_end, meta_hibf_config_start
#include <hibf/misc/next_multiple_of_64.hpp> // for next_multiple_of_64

#include <cereal/archives/json.hpp> // for JSONInputArchive, JSONOutputArchive
#include <cereal/cereal.hpp>        // for make_nvp, InputArchive, OutputArchive

namespace seqan::hibf
{

void config::read_from(std::istream & stream)
{
    std::string line;
    std::stringstream config_str;

    while (std::getline(stream, line) && line != prefix::meta_hibf_config_start)
        ;

    assert(line == prefix::meta_hibf_config_start);

    // TODO ##CONFIG: as prefix
    while (std::getline(stream, line) && line != prefix::meta_hibf_config_end)
    {
        assert(line.size() >= 2);
        assert(std::string_view{line}.substr(0, 1) == seqan::hibf::prefix::meta_header);
        config_str << line.substr(1); // remove seqan::hibf::prefix::meta_header
    }

    assert(line == prefix::meta_hibf_config_end);

    cereal::JSONInputArchive iarchive(config_str);
    iarchive(*this);
}

void config::write_to(std::ostream & stream) const
{
    // write json file to temprorary string stream with cereal
    std::stringstream config_stream{};
    cereal::JSONOutputArchive output(config_stream); // stream to cout
    output(cereal::make_nvp("hibf_config", *this));

    // write config
    stream << prefix::meta_hibf_config_start << '\n';
    std::string line;
    while (std::getline(config_stream, line, '\n'))
        stream << prefix::meta_header << line << '\n';
    stream << prefix::meta_header << "}\n" // last closing bracket isn't written by loop above
           << prefix::meta_hibf_config_end << '\n';
}

void config::validate_and_set_defaults()
{
    if (!input_fn)
        throw std::invalid_argument{"[HIBF CONFIG ERROR] You did not set the required config::input_fn."};

    if (number_of_user_bins == 0u)
        throw std::invalid_argument{"[HIBF CONFIG ERROR] You did not set the required config::number_of_user_bins."};

    if (number_of_hash_functions == 0u || number_of_hash_functions > 5u)
        throw std::invalid_argument{"[HIBF CONFIG ERROR] config::number_of_hash_functions must be in [1,5]."};

    if (maximum_fpr <= 0.0 || maximum_fpr >= 1.0)
        throw std::invalid_argument{"[HIBF CONFIG ERROR] config::maximum_fpr must be in (0.0,1.0)."};

    if (relaxed_fpr <= 0.0 || relaxed_fpr >= 1.0)
        throw std::invalid_argument{"[HIBF CONFIG ERROR] config::relaxed_fpr must be in (0.0,1.0)."};

    if (relaxed_fpr < maximum_fpr)
        throw std::invalid_argument{"[HIBF CONFIG ERROR] config::relaxed_fpr must be "
                                    "greater than or equal to config::maximum_fpr."};

    if (threads == 0u)
        throw std::invalid_argument{"[HIBF CONFIG ERROR] config::threads must be greater than 0."};

    // The following validations are not necessary when construction an interleaved_bloom_filter via config.
    // However, they also shouldn't result in problems.

    if (sketch_bits < 5u || sketch_bits > 32u)
        throw std::invalid_argument{"[HIBF CONFIG ERROR] config::sketch_bits must be in [5,32]."};

    if (tmax == 0) // no tmax was set by the user on the command line
    {
        tmax = seqan::hibf::next_multiple_of_64(std::ceil(std::sqrt(number_of_user_bins)));
    }
    else if (tmax > 18'446'744'073'709'551'552ULL) // next_multiple_of_64 would not fit in size_t. Underflowed by user?
    {
        throw std::invalid_argument{"[HIBF CONFIG ERROR] The maximum possible tmax is 18446744073709551552."};
    }
    else if (tmax % 64 != 0)
    {
        tmax = seqan::hibf::next_multiple_of_64(tmax);
        std::cerr << "[HIBF CONFIG WARNING]: Your requested number of technical bins was not a multiple of 64. "
                  << "Due to the architecture of the HIBF, it will use up space equal to the next multiple of 64 "
                  << "anyway, so we increased your number of technical bins to " << tmax << ".\n";
    }

    if (alpha < 0.0)
        throw std::invalid_argument{"[HIBF CONFIG ERROR] config::alpha must be positive."};

    if (max_rearrangement_ratio < 0.0 || max_rearrangement_ratio > 1.0)
        throw std::invalid_argument{"[HIBF CONFIG ERROR] config::max_rearrangement_ratio must be in [0.0,1.0]."};

    if (disable_estimate_union || max_rearrangement_ratio == 0.0)
        disable_rearrangement = true;
}

} // namespace seqan::hibf
