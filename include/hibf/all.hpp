// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Meta-header.
 */

/*!\defgroup ibf IBF
 * \brief The Interleaved Bloom Filter: seqan::hibf::interleaved_bloom_filter
 */

/*!\defgroup hibf HIBF
 * \brief The Hierarchical Interleaved Bloom Filter: seqan::hibf::hierarchical_interleaved_bloom_filter
 *
 * The Hierarchical Interleaved Bloom Filter is a data structure that provides fast answers to set-membership queries
 * for multiple samples/set/user-bins.
 *
 * The data structure is implemented in the seqan::hibf::hierarchical_interleaved_bloom_filter and improves its
 * predecessor the seqan::hibf::interleaved_bloom_filter data structure by more efficient storage of unevenly sized
 * samples/sets/user-bins and fast queries of thousands of samples.
 *
 * ### Cite
 *
 * *Mehringer, Svenja, et al. "Hierarchical Interleaved Bloom Filter: enabling ultrafast, approximate sequence queries."
 * Genome Biology 24.1 (2023): 1-25.*
 *
 * ### Example
 *
 * \include test/snippet/hibf/hierarchical_interleaved_bloom_filter.cpp
 *
 * \see seqan::hibf::hierarchical_interleaved_bloom_filter
 * \see seqan::hibf::interleaved_bloom_filter
 * \see Official paper: https://doi.org/10.1186/s13059-023-02971-4
 */

/*!\defgroup hibf_sketch Sketching
 * \brief Includes API to compute HyperLoglog sketches for the cardinality/size estimates for the layout algorithm.
 * \ingroup hibf
 */

/*!\defgroup hibf_sketch_toolbox Toolbox
 * \brief Algorithms for union estimation and similarity rearrangement executed before computing the layout.
 * \ingroup hibf_sketch
 */

/*!\defgroup hibf_layout Layout
 * \brief Includes all API relevant to computing a hierarchical layout for building an HIBF.
 * \ingroup hibf
 */

/*!\defgroup hibf_build Build
 * \brief Includes all API relevant to build an HIBF from a given layout.
 * \ingroup hibf
 */
#pragma once

#include <hibf/interleaved_bloom_filter.hpp>
