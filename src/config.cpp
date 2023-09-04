// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/hibf/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <cassert> // for assert
#include <cmath>
#include <iostream>
#include <sstream>     // for basic_istream, basic_ostream, operator<<, basic_stringstream, stringstream
#include <stdexcept>   // for invalid_argument
#include <string>      // for char_traits, getline, operator<<, string
#include <string_view> // for operator<<, operator==, basic_string_view

#include <hibf/config.hpp>              // for config
#include <hibf/layout/prefixes.hpp>     // for meta_header, meta_hibf_config_end, meta_hibf_config_start
#include <hibf/next_multiple_of_64.hpp> // for next_multiple_of_64

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
    if (disable_estimate_union)
        disable_rearrangement = true;

    if (number_of_user_bins == 0u)
        throw std::invalid_argument{"[HIBF CONFIG ERROR] You did not set the required config::number_of_user_bins."};

    if (tmax == 0) // no tmax was set by the user on the command line
    {
        // Set default as sqrt(#samples). Experiments showed that this is a reasonable default.
        if (number_of_user_bins >= 1ULL << 32) // sqrt is bigger than uint16_t
        {
            throw std::invalid_argument{
                "[HIBF CONFIG ERROR] Too many user bins to compute a default tmax. " // GCOVR_EXCL_LINE
                "Please set a tmax manually."};                                      // GCOVR_EXCL_LINE
        }
        else
        {
            tmax = seqan::hibf::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(number_of_user_bins))));
        }
    }
    else if (tmax % 64 != 0)
    {
        tmax = seqan::hibf::next_multiple_of_64(tmax);
        std::cerr << "[HIBF CONFIG WARNING]: Your requested number of technical bins was not a multiple of 64. "
                  << "Due to the architecture of the HIBF, it will use up space equal to the next multiple of 64 "
                  << "anyway, so we increased your number of technical bins to " << tmax << ".\n";
    }
}

} // namespace seqan::hibf
