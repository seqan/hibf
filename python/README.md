<!--
SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

# HIBF Python bindings

Python bindings for the [Hierarchical Interleaved Bloom Filter (HIBF)](https://github.com/seqan/hibf), built with
[nanobind](https://github.com/wjakob/nanobind).

## Installation

Requirements: GCC >= 14 or Clang >= 20, CMake >= 3.20, Python >= 3.10.

```bash
pip install .                  # from the repository root
pip install ".[test]"          # additionally installs pytest
```

The build fetches the C++ dependencies (cereal, simde) via CPM. Set `CPM_SOURCE_CACHE` to reuse downloads across builds.

By default, the library is compiled with `-march=native`. The result is fast, but it may not run on other machines. For
portable wheels, disable this:

```bash
pip wheel . -C cmake.define.HIBF_NATIVE_BUILD=OFF
```

## Quick start

```python
import numpy as np
import hibf

# Each user bin is a set of unsigned 64-bit integers, e.g., hashed k-mers or minimisers.
user_bins = [
    np.arange(1, 11, dtype=np.uint64),
    np.array([1, 2, 3, 4, 5], dtype=np.uint64),
    np.array([3, 9, 11], dtype=np.uint64),
]

config = hibf.Config(user_bins, maximum_fpr=0.05, threads=4)
index = hibf.HIBF(config)

# Which user bins contain at least 2 of the query values?
index.membership_for([3, 9, 12, 14], threshold=2)  # array([2, 0], dtype=uint64); the order is unspecified
```

## Providing input

`Config.input` provides the values of each user bin. It can be:

* **A sequence**, where `input[i]` holds the values of user bin `i`. If `number_of_user_bins` is not given, it defaults
  to `len(input)`.
* **A callable**, where `input(i)` returns the values of user bin `i`. Use this when the data does not fit into memory,
  e.g., to read a file per user bin. `number_of_user_bins` is required.

The values of a user bin can be any iterable of non-negative integers. A contiguous NumPy `uint64` array is used without
copying. Other inputs, such as lists, generators, sets, or arrays with a different dtype, are converted first.

The library requests the values of each user bin **more than once** (once for the layout, once for building the
index). With `threads > 1`, requests come from several threads at once. The input must therefore return the same values
every time it is called. A callable that returns a new generator on each call works; a list of generators does not.

User bins must not be empty. Exceptions raised by the input, including `KeyboardInterrupt`, stop the construction and
are raised by the constructor.

## Querying

```python
# Single queries. The result is an unsorted array of user bin ids.
hits = index.membership_for(query, threshold)

# Agents avoid reallocating buffers. Use one agent per thread.
agent = index.membership_agent()
for query in queries:
    hits = agent.membership_for(query, threshold)

# Counts per user bin. `dtype` is uint16 (default), uint32, or uint64.
counts = index.counting_agent(dtype=np.uint32).bulk_count(query, threshold=1)

# Many queries in parallel. `threshold` is a single value or one value per query.
results = index.batch_membership_for(queries, [int(0.8 * len(q)) for q in queries], threads=8)
```

`batch_membership_for` releases the GIL and runs in parallel on its own threads. This is the fastest way to answer many
queries.

## Thread safety

Querying from several threads is safe. Give each thread its own agent; `batch_membership_for` manages its threads
itself.

Methods that modify an object are not synchronised with other methods on the same object:
`InterleavedBloomFilter.emplace`, `clear`, `increase_bin_number_to` and `try_increase_bin_number_to`, and
`HyperLogLog.add`, `merge` and `reset`. Do not call them while another thread uses that object. Some operations, such
as `batch_membership_for` and `save`, release the GIL, so another thread can run while they work. If such a thread
increases the number of bins of the IBF they are reading, the interpreter can crash.

## Serialisation

```python
index.save("index.hibf")
index = hibf.HIBF.load("index.hibf")

import pickle
index = pickle.loads(pickle.dumps(index))
```

`InterleavedBloomFilter` and `HyperLogLog` support the same methods. `Config` can be pickled, too. It is stored
together with its `input`.

Two independently built HIBFs may number their internal IBFs differently. They answer queries identically, but `==`
can still report them as different. `==` does return `True` for a saved and reloaded index.

## Layouts

The HIBF constructor computes a layout, which assigns user bins to IBFs, and then builds the index. You can also run the
two steps separately. For example, you can store a layout or reuse a layout computed by
[chopper](https://github.com/seqan/chopper):

```python
layout = hibf.compute_layout(config)
hibf.write_layout_file("layout.txt", config, layout)

config, layout = hibf.read_layout_file("layout.txt")
config.input = user_bins  # The input is not stored in layout files.
index = hibf.HIBF(config, layout)
```

## Interleaved Bloom Filter

```python
ibf = hibf.IBF(bin_count=8, bin_size=8192, hash_function_count=2)
ibf.emplace(np.array([126, 712], dtype=np.uint64), bin=0)
ibf.emplace(712, bin=3)

ibf.containment_agent().bulk_contains(712)  # boolean array with one entry per bin
ibf.membership_for([126, 712], threshold=2)  # array([0], dtype=uint64)

ibf = hibf.IBF(config)  # one bin per user bin, sized for config.maximum_fpr
```

## HyperLogLog

```python
sketch = hibf.HyperLogLog(num_bits=12)
sketch.add(hashed_values)  # values should be uniformly distributed, e.g., hashes
sketch.estimate()
```

## Development

```bash
pip install --no-build-isolation -e ".[test]"  # or: pip install ".[test]"
pytest
```

Type stubs (`_hibf.pyi`) are generated during the build. The C++ sources are in [`src/`](src), and the tests are in
[`tests/`](tests).
