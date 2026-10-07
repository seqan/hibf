// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <cstddef>   // for size_t
#include <iostream>  // for basic_ostream, operator<<, cout
#include <stdexcept> // for invalid_argument
#include <vector>    // for vector

#include <hibf/config.hpp>        // for config
#include <hibf/layout/layout.hpp> // for layout

using diagnostic = seqan::hibf::layout::layout::diagnostic;

int main()
{
    // Top-level IBF: User bin 0 in technical bin 0, user bin 1 in technical bins 2-62, merged bin in technical bin 63.
    // Technical bin 1 is empty (a note). The IBF below the merged bin only contains user bin 2 (a warning).
    seqan::hibf::layout::layout layout{};
    layout.top_level_max_bin_id = 0u;
    layout.max_bins.push_back({.previous_TB_indices = {63u}, .id = 0u});
    layout.user_bins.push_back(
        {.previous_TB_indices = {}, .storage_TB_id = 0u, .number_of_technical_bins = 1u, .idx = 0u});
    layout.user_bins.push_back(
        {.previous_TB_indices = {}, .storage_TB_id = 2u, .number_of_technical_bins = 61u, .idx = 1u});
    layout.user_bins.push_back(
        {.previous_TB_indices = {63u}, .storage_TB_id = 0u, .number_of_technical_bins = 64u, .idx = 2u});

    seqan::hibf::config const config{.number_of_user_bins = 3u, .tmax = 64u};

    std::cout << "# collect\n";
    //![collect]
    std::vector<diagnostic> diagnostics{};
    bool const valid = layout.validate(config,
                                       [&diagnostics](diagnostic const & finding)
                                       {
                                           diagnostics.push_back(finding);
                                       });
    std::cout << "valid: " << valid << ", findings: " << diagnostics.size() << '\n'; // valid: 1, findings: 2
    //![collect]

    std::cout << "# print\n";
    //![print]
    auto print = [](diagnostic const & finding)
    {
        std::cout << finding << '\n'; // [HIBF LAYOUT <SEVERITY>] <message>
    };
    if (layout.validate(config, print))
        std::cout << "valid\n";
    //![print]

    std::cout << "# ignore\n";
    //![ignore]
    auto ignore_empty_technical_bins = [](diagnostic const & finding)
    {
        if (finding.what != diagnostic::code::empty_technical_bins)
            std::cout << finding << '\n';
    };
    if (layout.validate(config, ignore_empty_technical_bins))
        std::cout << "valid\n";
    //![ignore]

    std::cout << "# strict\n";
    //![strict]
    auto warnings_as_errors = [](diagnostic const & finding)
    {
        if (finding.level >= diagnostic::severity::warning)
            throw std::invalid_argument{finding.message};
    };

    try
    {
        if (layout.validate(config, warnings_as_errors))
            std::cout << "valid\n";
    }
    catch (std::invalid_argument const & exception)
    {
        std::cout << "rejected: " << exception.what() << '\n';
    }
    //![strict]

    std::cout << "# locate\n";
    //![locate]
    auto locate = [](diagnostic const & finding)
    {
        std::cout << "IBF:";
        if (!finding.ibf.has_value())
            std::cout << " none"; // the finding concerns no particular IBF
        else if (finding.ibf->empty())
            std::cout << " Root-IBF";
        else
            for (size_t const technical_bin : *finding.ibf) // the merged bins on the path to the IBF
                std::cout << ' ' << technical_bin;
        if (finding.user_bin.has_value())
            std::cout << ", user bin: " << *finding.user_bin;
        if (finding.technical_bin.has_value())
            std::cout << ", technical bin: " << *finding.technical_bin;
        std::cout << '\n';
    };
    if (layout.validate(config, locate))
        std::cout << "valid\n";
    //![locate]

    std::cout << "# validate\n";
    //![validate]
    try
    {
        layout.validate(seqan::hibf::config{.number_of_user_bins = 4u}); // user bin 3 is missing
    }
    catch (std::invalid_argument const & exception)
    {
        std::cout << exception.what() << '\n';
    }
    //![validate]
}
