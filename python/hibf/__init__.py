# SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

"""Python bindings for the Hierarchical Interleaved Bloom Filter (HIBF) library."""

from importlib.metadata import PackageNotFoundError, version

from ._hibf import (
    Config,
    HierarchicalInterleavedBloomFilter,
    HyperLogLog,
    InterleavedBloomFilter,
    Layout,
    compute_layout,
    library_version,
    read_layout_file,
    write_layout_file,
)

HIBF = HierarchicalInterleavedBloomFilter
IBF = InterleavedBloomFilter

try:
    __version__ = version("hibf")
except PackageNotFoundError:  # pragma: no cover
    __version__ = library_version

__all__ = [
    "HIBF",
    "IBF",
    "Config",
    "HierarchicalInterleavedBloomFilter",
    "HyperLogLog",
    "InterleavedBloomFilter",
    "Layout",
    "__version__",
    "compute_layout",
    "library_version",
    "read_layout_file",
    "write_layout_file",
]
