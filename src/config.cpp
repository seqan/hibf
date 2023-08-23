// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/hibf/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <cassert>
#include <charconv>
#include <iostream>

#include <hibf/config.hpp>
#include <hibf/detail/prefixes.hpp>

#include <cereal/archives/json.hpp>

namespace hibf
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
        assert(std::string_view{line}.substr(0, 1) == hibf::prefix::meta_header);
        config_str << line.substr(1); // remove hibf::prefix::meta_header
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

} // namespace hibf
