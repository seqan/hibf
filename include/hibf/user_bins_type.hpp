// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

#include <hibf/migration/cereal.hpp>

namespace hibf
{

//!\brief Bookkeeping for user and technical bins.
class user_bins_type
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

} // namespace hibf
