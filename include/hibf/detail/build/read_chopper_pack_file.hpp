// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <string> // for string
#include <vector> // for vector

#include <hibf/config.hpp>
#include <hibf/detail/configuration.hpp> // for configuration
#include <hibf/detail/layout/layout.hpp> // for layout

namespace hibf
{

layout::layout read_chopper_pack_file(configuration & config,
                                      std::vector<std::vector<std::string>> & filenames,
                                      std::string const & chopper_pack_filename);

} // namespace hibf
