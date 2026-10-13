<!--
SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

# Changelog {#about_changelog}

[TOC]

This changelog contains a top-level entry for each release with sections on new features, API changes and notable
bug-fixes (not all bug-fixes will be listed).

<!--
The following API changes should be documented as such:
  * a previously experimental interface now being marked as stable
  * an interface being removed
  * syntactical changes to an interface (e.g. renaming or reordering of files, functions, parameters)
  * semantic changes to an interface (e.g. a function's result is now always one larger) [DANGEROUS!]

If possible, provide tooling that performs the changes, e.g. a shell-script.
-->

# 1.0.0

## New features

## Notable Bug-fixes

* HIBF queries no longer miss a split user bin if the counts of its technical bins sum to more than the counter type
  can hold, e.g., more than 65535 for membership queries. Counting agents report such counts as the counter type's
  maximum instead of wrapping around.
* `config::operator==` now also compares `number_of_hash_functions` and `track_occupancy`.
* `interleaved_bloom_filter::clear` now resets the occupancy of the cleared bins when given a range of bins, as it
  already did for a single bin.
* Adding a `bit_vector` to a `counting_vector` saturates the counts at the maximum of the counter type instead of
  wrapping around. Hence, the counting agents of the IBF and HIBF report such counts as the maximum, and HIBF
  membership queries with more than 65535 values no longer miss user bins.

## API changes
