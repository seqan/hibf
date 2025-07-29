// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <generator>

#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

#include "fmt.hpp"
#include "inspector.hpp"

namespace seqan::hibf::util
{

// For occupancy, add
// helper("Occupancy", inspector::occupancy); // in this function
// #include <fmt/ranges.h> // at top
// In helper in this function:
// if constexpr (std::integral<decltype(fun(ibf))>)
// {
//     fmt::println(" {}", fmt::group_digits(fun(ibf)));
// }
// else
// {
//     fmt::println(" {}", fun(ibf));
// }
void inspect(seqan::hibf::interleaved_bloom_filter const & ibf,
             size_t const indent = 0u,
             std::optional<size_t> user_bins = std::nullopt)
{
    auto hibf_user_bins = [user_bins](seqan::hibf::interleaved_bloom_filter const &) -> size_t
    {
        return user_bins.value();
    };

    static auto const style = styles::emphasis(fmt::emphasis::bold) | styles::color(fmt::color::medium_purple);
    auto helper = [&](std::string_view const & what, auto && fun)
    {
        fmt::print("{:{}}", "", indent);
        fmt::print(style, fmt::runtime(what));
        fmt::print(style, ":");
        fmt::println(" {}", fmt::group_digits(fun(ibf)));
    };

    using inspector = seqan::hibf::inspector;

    if (user_bins.has_value())
    {
        helper("User bins", hibf_user_bins);
        helper("Active bins", inspector::bin_count);
    }
    else
    {
        helper("User bins", inspector::bin_count);
    }
    helper("Technical bins", inspector::technical_bins);
    helper("Empty bins", inspector::empty_bins);
    helper("Bin words", inspector::bin_words);
    helper("Bin size", inspector::bin_size);
    helper("Total size", inspector::bit_size);
    helper("Hash functions", inspector::hash_function_count);
    helper("Hash shift", inspector::hash_shift);
}

struct tree
{
    uint64_t idx{};
    uint64_t level{};
    seqan::hibf::interleaved_bloom_filter const * ibf{nullptr};
    std::vector<tree> children;

    tree(seqan::hibf::hierarchical_interleaved_bloom_filter const & hibf,
         uint64_t const ibf_idx = 0u,
         uint64_t const ibf_level = 0u) :
        idx{ibf_idx},
        level{ibf_level},
        ibf{std::addressof(hibf.ibf_vector[idx])}
    {
        for (uint64_t const i : hibf.next_ibf_id[idx])
        {
            if (i != idx)
                children.emplace_back(hibf, i, level + 1u);
        }
    }

    std::generator<tree> traverse() const
    {
        co_yield *this;

        auto view = std::views::transform(children,
                                          [](tree const & child)
                                          {
                                              return child.traverse();
                                          });

        co_yield std::ranges::elements_of(std::views::join(view));
    }
};

void inspect(seqan::hibf::hierarchical_interleaved_bloom_filter const & hibf)
{
    static auto const style = styles::emphasis(fmt::emphasis::bold) | styles::color(fmt::color::medium_purple);
    auto helper = [&](std::string_view const & what, auto && fun)
    {
        fmt::print(style, fmt::runtime(what));
        fmt::print(style, ":");
        fmt::println(" {}", fmt::group_digits(fun(hibf)));
    };

    using inspector = seqan::hibf::inspector;

    helper("User bins", inspector::number_of_user_bins);
    helper("Empty bins", inspector::total_empty);
    helper("IBFs", inspector::number_of_ibfs);
    helper("Total size", inspector::total_size);

    tree const tree{hibf};
    size_t const width = std::ceil(std::log10(inspector::number_of_ibfs(hibf) + 1u));
    auto const style2 = styles::emphasis(fmt::emphasis::bold) | styles::color(fmt::color::sienna);
    for (auto const & x : tree.traverse())
    {
        fmt::print("{:{}}", "", x.level * 4u);
        fmt::print(style2, "ID:");
        fmt::print(" {:<{}} ", x.idx, width);
        fmt::print(style2, "Level:");
        fmt::print(" {} ", x.level);
        fmt::print(style2, "Children:");
        fmt::print(" {}\n", x.children.size());
        size_t const user_bins = inspector::user_bins(hibf, x.idx);
        inspect(hibf.ibf_vector[x.idx], x.level * 4u + 2u, user_bins);
    }
}

} // namespace seqan::hibf::util
