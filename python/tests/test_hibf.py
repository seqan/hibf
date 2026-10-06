# SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

import copy
import gc
import pickle
import weakref

import numpy as np
import pytest

import hibf

README_DATA = [list(range(1, 11)), [1, 2, 3, 4, 5], [3, 9, 11]]


def random_user_bins(number_of_user_bins=100, seed=0):
    """User bins of varying sizes with disjoint values, so that each query has exactly one true hit."""
    rng = np.random.default_rng(seed)
    sizes = rng.integers(50, 2000, size=number_of_user_bins)
    values = rng.permutation(np.uint64(2**62) - np.arange(sizes.sum(), dtype=np.uint64))
    return np.split(values, np.cumsum(sizes)[:-1])


def assert_same_results(index, other, user_bins):
    """Independent builds may number their IBFs differently, so `==` is not suitable. Compare query results."""
    queries = [values[:100] for values in user_bins]
    expected = index.batch_membership_for(queries, 80)
    actual = other.batch_membership_for(queries, 80)
    assert [sorted(r) for r in actual] == [sorted(r) for r in expected]


@pytest.fixture(scope="module")
def user_bins():
    return random_user_bins()


@pytest.fixture(scope="module")
def index(user_bins):
    return hibf.HIBF(hibf.Config(user_bins, threads=2))


# --------------------------------------------------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------------------------------------------------


def test_config_defaults():
    config = hibf.Config()
    assert config.input is None
    assert config.number_of_user_bins == 0
    assert config.number_of_hash_functions == 2
    assert config.maximum_fpr == pytest.approx(0.05)
    assert config.relaxed_fpr == pytest.approx(0.3)
    assert config.threads == 1
    assert config.sketch_bits == 12
    assert config.tmax == 0
    assert not config.disable_rearrangement


def test_config_infers_number_of_user_bins():
    assert hibf.Config(README_DATA).number_of_user_bins == 3
    assert hibf.Config(README_DATA, number_of_user_bins=2).number_of_user_bins == 2
    assert hibf.Config(lambda i: [i]).number_of_user_bins == 0

    config = hibf.Config()
    config.input = README_DATA
    assert config.number_of_user_bins == 3


def test_config_number_of_user_bins_follows_input():
    data = random_user_bins(5)
    config = hibf.Config(data[:3])
    config.input = data
    assert config.number_of_user_bins == 5
    assert hibf.HIBF(config).number_of_user_bins == 5
    assert hibf.Config.from_string(config.to_string()).number_of_user_bins == 5

    restored = pickle.loads(pickle.dumps(config))
    restored.input = data[:4]
    assert restored.number_of_user_bins == 4

    config.number_of_user_bins = 2  # An explicit value takes precedence.
    config.input = data[:4]
    assert config.number_of_user_bins == 2

    config.number_of_user_bins = 0  # Back to inference.
    assert config.number_of_user_bins == 4


def test_config_rejects_invalid_input():
    with pytest.raises(TypeError):
        hibf.Config(42)


def test_config_validate():
    config = hibf.Config(README_DATA)
    config.validate_and_set_defaults()
    assert config.tmax == 64

    with pytest.raises(ValueError, match="number_of_hash_functions"):
        hibf.Config(README_DATA, number_of_hash_functions=6).validate_and_set_defaults()

    with pytest.raises(ValueError, match="input"):
        hibf.Config(number_of_user_bins=3).validate_and_set_defaults()


def test_config_equality():
    assert hibf.Config(README_DATA) == hibf.Config(lambda i: README_DATA[i], number_of_user_bins=3)
    assert hibf.Config(number_of_hash_functions=2) != hibf.Config(number_of_hash_functions=5)
    assert hibf.Config(track_occupancy=False) != hibf.Config(track_occupancy=True)


def test_config_string_roundtrip():
    config = hibf.Config(README_DATA, maximum_fpr=0.01, tmax=128)
    parsed = hibf.Config.from_string(config.to_string())
    assert parsed == config
    assert parsed.input is None


def test_config_pickle_and_copy():
    config = hibf.Config(README_DATA, maximum_fpr=0.01)
    for other in (pickle.loads(pickle.dumps(config)), copy.copy(config), copy.deepcopy(config)):
        assert other == config
        assert other.input == README_DATA
    assert copy.copy(config).input is config.input
    assert copy.deepcopy(config).input is not config.input


def test_config_reference_cycle_is_collected():
    """A bound method as input creates a cycle: object -> config -> input -> object."""

    class Indexer:
        def __init__(self):
            self.data = [np.arange(1, 6, dtype=np.uint64)]
            self.config = hibf.Config(self.read_bin, number_of_user_bins=1)

        def read_bin(self, user_bin_id):
            return self.data[user_bin_id]

    indexer = Indexer()
    hibf.HIBF(indexer.config)
    ref = weakref.ref(indexer)
    del indexer
    gc.collect()
    assert ref() is None


def test_config_repr():
    assert repr(hibf.Config(maximum_fpr=0.01)).startswith("Config(input=None, number_of_user_bins=0")


# --------------------------------------------------------------------------------------------------------------------
# HIBF construction
# --------------------------------------------------------------------------------------------------------------------


def test_readme_example():
    index = hibf.HIBF(hibf.Config(README_DATA))
    assert index.number_of_user_bins == 3
    assert sorted(index.membership_for([3, 9, 12, 14], 2)) == [0, 2]

    agent = index.membership_agent()
    assert sorted(agent.membership_for(np.array([1, 2, 3, 4, 5], dtype=np.uint64), 5)) == [0, 1]


@pytest.mark.parametrize(
    "make_input",
    [
        pytest.param(lambda data: data, id="list of arrays"),
        pytest.param(lambda data: [list(map(int, ub)) for ub in data], id="list of lists"),
        pytest.param(lambda data: tuple(ub.astype(np.int64) for ub in data), id="tuple of int64 arrays"),
        pytest.param(lambda data: (lambda i: data[i]), id="callable"),
        pytest.param(lambda data: (lambda i: (int(v) for v in data[i])), id="callable returning generator"),
        pytest.param(lambda data: (lambda i: set(map(int, data[i]))), id="callable returning set"),
    ],
)
def test_input_types(make_input):
    data = random_user_bins(20)
    config = hibf.Config(make_input(data), number_of_user_bins=len(data))
    index = hibf.HIBF(config)
    for user_bin_id, values in enumerate(data):
        assert user_bin_id in index.membership_for(values, len(values))


def test_input_2d_array():
    data = np.arange(1, 1001, dtype=np.uint64).reshape(10, 100)
    index = hibf.HIBF(hibf.Config(data))
    assert index.number_of_user_bins == 10
    assert 3 in index.membership_for(data[3], 100)


def test_no_false_negatives(index, user_bins):
    assert index.number_of_ibfs > 1
    for user_bin_id, values in enumerate(user_bins):
        assert user_bin_id in index.membership_for(values, len(values))


def test_threads_give_same_results(user_bins, index):
    assert_same_results(index, hibf.HIBF(hibf.Config(user_bins, threads=1)), user_bins)
    assert_same_results(index, hibf.HIBF(hibf.Config(user_bins, threads=4)), user_bins)


def test_missing_input():
    with pytest.raises(ValueError, match="input"):
        hibf.HIBF(hibf.Config(number_of_user_bins=3))


def test_empty_user_bin():
    with pytest.raises(ValueError, match="User bin 1 is empty"):
        hibf.HIBF(hibf.Config([[1, 2, 3], [], [4, 5]]))


@pytest.mark.parametrize("threads", [1, 4])
def test_callback_exception_is_propagated(threads):
    data = random_user_bins(50)

    def input_fn(user_bin_id):
        if user_bin_id == 17:
            raise KeyError("broken user bin")
        return data[user_bin_id]

    with pytest.raises(KeyError, match="broken user bin"):
        hibf.HIBF(hibf.Config(input_fn, number_of_user_bins=len(data), threads=threads))

    with pytest.raises(KeyError, match="broken user bin"):
        hibf.IBF(hibf.Config(input_fn, number_of_user_bins=len(data), threads=threads))


def test_callback_exception_during_build_is_propagated():
    """The layout stage succeeds; the error occurs when filling the index."""
    data = random_user_bins(20)
    calls = [0] * len(data)

    def input_fn(user_bin_id):
        calls[user_bin_id] += 1
        if user_bin_id == 5 and calls[user_bin_id] > 1:
            raise RuntimeError("second call fails")
        return data[user_bin_id]

    with pytest.raises(RuntimeError, match="second call fails"):
        hibf.HIBF(hibf.Config(input_fn, number_of_user_bins=len(data)))


@pytest.mark.parametrize("build", [hibf.HIBF, hibf.IBF, hibf.compute_layout])
def test_input_reassigned_during_construction(build):
    """The input drops the last reference to itself while the library still requests values."""
    data = random_user_bins(30)
    config = hibf.Config(number_of_user_bins=len(data))

    class Source:
        def __call__(self, user_bin_id):
            config.input = data
            return data[user_bin_id]

    config.input = Source()
    result = build(config)
    if not isinstance(result, hibf.Layout):
        for user_bin_id, values in enumerate(data):
            assert user_bin_id in result.membership_for(values, len(values))


def test_user_bin_empty_during_build():
    """The input returns different values when requested again."""
    data = random_user_bins(20)
    calls = [0] * len(data)

    def input_fn(user_bin_id):
        calls[user_bin_id] += 1
        return data[user_bin_id] if calls[user_bin_id] == 1 else []

    with pytest.raises(ValueError, match="is empty"):
        hibf.HIBF(hibf.Config(input_fn, number_of_user_bins=len(data)))


def test_ibf_allows_empty_user_bins():
    ibf = hibf.IBF(hibf.Config([[1, 2, 3], [], [4, 5]]))
    assert list(ibf.membership_for([4, 5], 2)) == [2]


def test_invalid_values():
    with pytest.raises(OverflowError):
        hibf.HIBF(hibf.Config([[1, 2], [-1]]))
    with pytest.raises(TypeError):
        hibf.HIBF(hibf.Config([[1, 2], [[1, 2], [3, 4]]]))
    with pytest.raises(ValueError, match="one-dimensional"):
        hibf.HIBF(hibf.Config([[1, 2], np.ones((2, 2), dtype=np.uint64)]))
    with pytest.raises(TypeError, match="dtype float64"):
        hibf.HIBF(hibf.Config([[1, 2], np.array([1.5, 2.5])]))


def test_value_conversion():
    ibf = hibf.IBF(bin_count=1, bin_size=1 << 16)

    for values in (np.array([1.5]), [1.0], 2.5, np.float64(2.0), np.array([True]), b"\x01", "1", [[1, 2]]):
        with pytest.raises(TypeError):
            ibf.emplace(values, 0)
        with pytest.raises(TypeError):
            ibf.membership_for(values, 1)
    for values in ([-1], -1, [2**64]):
        with pytest.raises(OverflowError):
            ibf.emplace(values, 0)
    with pytest.raises(ValueError, match="one-dimensional"):
        ibf.emplace(np.ones((2, 2), dtype=np.uint64), 0)

    ibf.emplace(np.array([-1], dtype=np.int64), 0)  # reinterpreted as 2**64 - 1
    ibf.emplace([2**64 - 2], 0)
    ibf.emplace(np.arange(10, 15, dtype=">u8"), 0)  # non-native byte order
    ibf.emplace(np.uint8(7), 0)
    ibf.emplace((v for v in [20, 21]), 0)
    assert list(ibf.membership_for([2**64 - 1, 2**64 - 2, 10, 14, 7, 21], 6)) == [0]


def test_timings(index):
    assert set(index.timings) == {
        "layout_compute_sketches",
        "layout_union_estimation",
        "layout_rearrangement",
        "layout_dp_algorithm",
        "index_allocation",
        "user_bin_io",
        "merge_kmers",
        "fill_ibf",
    }
    assert index.timings["layout_compute_sketches"] > 0


# --------------------------------------------------------------------------------------------------------------------
# HIBF queries
# --------------------------------------------------------------------------------------------------------------------


def test_membership_result_type(index, user_bins):
    result = index.membership_for(user_bins[0], 1)
    assert isinstance(result, np.ndarray)
    assert result.dtype == np.uint64


def test_membership_query_input_types(index, user_bins):
    values = user_bins[7]
    expected = sorted(index.membership_for(values, len(values)))
    assert 7 in expected
    for query in (list(map(int, values)), tuple(map(int, values)), values.astype(np.int64), values[::-1].copy()):
        assert sorted(index.membership_for(query, len(values))) == expected
    assert sorted(index.membership_for(values[::2], len(values[::2]))) == sorted(
        index.membership_for(values[::2].copy(), len(values[::2]))
    )


@pytest.mark.parametrize("dtype", ["uint16", "uint32", np.uint64])
def test_counting_agent(index, user_bins, dtype):
    agent = index.counting_agent(dtype)
    counts = agent.bulk_count(user_bins[3])
    assert counts.dtype == np.dtype(dtype)
    assert counts.shape == (index.number_of_user_bins,)
    assert counts[3] >= len(user_bins[3])


def test_counting_agent_invalid(index):
    with pytest.raises(ValueError, match="Unsupported dtype"):
        index.counting_agent("int16")
    with pytest.raises(ValueError, match="threshold"):
        index.counting_agent().bulk_count([1, 2, 3], 0)


def test_query_size_limits(index):
    values = np.arange(1, 65_537, dtype=np.uint64)  # 65536 values
    ibf = hibf.IBF(bin_count=1, bin_size=1 << 20)
    ibf.emplace(values, 0)

    for filter in (index, ibf):
        for query in (filter.membership_for, filter.membership_agent().membership_for):
            with pytest.raises(ValueError, match="at most 65535 values, got 65536"):
                query(values, 1)
            query(values[:-1], 1)
        with pytest.raises(ValueError, match="at most 65535 values"):
            filter.batch_membership_for([values[:10], values], 1)
        with pytest.raises(ValueError, match="at most 65535 values"):
            filter.counting_agent().bulk_count(values)
        filter.counting_agent("uint32").bulk_count(values)

    assert list(ibf.membership_for(values[:-1], len(values) - 1)) == [0]
    assert ibf.counting_agent("uint32").bulk_count(values)[0] == len(values)


def test_split_bin_sum_exceeds_uint16():
    """A query value may be a false positive in the other technical bins of a split bin; the sum exceeds 65535."""
    values = np.arange(1, 70_001, dtype=np.uint64)
    index = hibf.HIBF(hibf.Config([values], maximum_fpr=0.3, relaxed_fpr=0.3))
    query = values[:60_000]
    assert index.counting_agent("uint32").bulk_count(query)[0] > 65_535
    assert list(index.membership_for(query, 60_000)) == [0]
    assert index.counting_agent().bulk_count(query)[0] == 65_535


def test_batch_membership(index, user_bins):
    queries = user_bins[:10]
    thresholds = [len(q) for q in queries]
    expected = [sorted(index.membership_for(q, t)) for q, t in zip(queries, thresholds)]

    for threads in (1, 3):
        results = index.batch_membership_for(queries, thresholds, threads=threads)
        assert [sorted(r) for r in results] == expected

    assert len(index.batch_membership_for(queries, 1)) == 10
    assert index.batch_membership_for([], 1) == []

    with pytest.raises(ValueError, match="number of thresholds"):
        index.batch_membership_for(queries, [1, 2])
    with pytest.raises(ValueError, match="threads"):
        index.batch_membership_for(queries, 1, threads=0)


def test_unhashable(index):
    """The classes compare by value, so the identity-based default hash would break sets and dicts."""
    for obj in (hibf.Config(), hibf.Layout(), hibf.IBF(bin_count=4, bin_size=64), index):
        with pytest.raises(TypeError, match="unhashable"):
            hash(obj)


def test_agent_keeps_index_alive(user_bins):
    index = hibf.HIBF(hibf.Config(user_bins[:10]))
    agent = index.membership_agent()
    del index
    gc.collect()
    assert 4 in agent.membership_for(user_bins[4], len(user_bins[4]))


# --------------------------------------------------------------------------------------------------------------------
# Serialisation
# --------------------------------------------------------------------------------------------------------------------


def test_hibf_pickle_and_copy(index, user_bins):
    for other in (pickle.loads(pickle.dumps(index)), copy.copy(index), copy.deepcopy(index)):
        assert other == index
        assert other is not index
        assert 5 in other.membership_for(user_bins[5], len(user_bins[5]))


def test_hibf_save_load(index, tmp_path):
    path = tmp_path / "index.hibf"
    index.save(path)
    assert hibf.HIBF.load(path) == index
    assert hibf.HIBF.load(str(path)) == index


def test_load_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        hibf.HIBF.load(tmp_path / "does_not_exist")


def test_setstate_invalid_bytes():
    with pytest.raises(Exception):
        pickle.loads(pickle.dumps(hibf.HIBF(hibf.Config(README_DATA)))[:-20])


# --------------------------------------------------------------------------------------------------------------------
# Layout
# --------------------------------------------------------------------------------------------------------------------


def test_layout(user_bins, index, tmp_path):
    config = hibf.Config(user_bins)
    layout = hibf.compute_layout(config)
    assert layout.number_of_user_bins == len(user_bins)
    assert layout.number_of_max_bins == index.number_of_ibfs - 1
    assert hibf.Layout.from_string(layout.to_string()) == layout

    assert_same_results(index, hibf.HIBF(config, layout), user_bins)

    path = tmp_path / "layout.txt"
    hibf.write_layout_file(path, config, layout)
    read_config, read_layout = hibf.read_layout_file(path)
    assert read_layout == layout
    assert read_config == config
    assert read_config.input is None

    read_config.input = user_bins
    assert_same_results(index, hibf.HIBF(read_config, read_layout), user_bins)


def test_layout_must_match_config(user_bins):
    config = hibf.Config(user_bins[:3])
    with pytest.raises(ValueError, match="contains 0 user bins, but the config has 3"):
        hibf.HIBF(config, hibf.Layout())
    with pytest.raises(ValueError, match="contains 100 user bins, but the config has 3"):
        hibf.HIBF(config, hibf.compute_layout(hibf.Config(user_bins)))

    # Each line after the header is `user_bin_id<TAB>technical_bin_indices<TAB>number_of_technical_bins`.
    header, lines = hibf.compute_layout(config).to_string().split("#USER_BIN_IDX")
    column_names, *records = lines.rstrip("\n").split("\n")
    for replacement in ("7", "0"):  # out of range, duplicate
        changed = [replacement + record[1:] if record.startswith("2\t") else record for record in records]
        layout = hibf.Layout.from_string(header + "#USER_BIN_IDX" + "\n".join([column_names, *changed]) + "\n")
        with pytest.raises(ValueError, match=f"Found {replacement}"):
            hibf.HIBF(config, layout)


# --------------------------------------------------------------------------------------------------------------------
# IBF
# --------------------------------------------------------------------------------------------------------------------


def test_ibf_basic():
    ibf = hibf.IBF(bin_count=12, bin_size=8192, hash_function_count=3)
    assert (ibf.bin_count, ibf.bin_size, ibf.hash_function_count) == (12, 8192, 3)
    assert ibf.bit_size == 64 * 8192

    ibf.emplace([126, 712], 0)
    ibf.emplace(712, 3)
    ibf.emplace(np.array([1, 2, 3], dtype=np.uint64), 11)

    contained = ibf.containment_agent().bulk_contains(712)
    assert contained.dtype == np.bool_
    assert contained.shape == (12,)
    assert contained[0] and contained[3]

    assert list(ibf.membership_for([126, 712], 2)) == [0]
    counts = ibf.counting_agent().bulk_count([126, 712, 1])
    assert counts[0] == 2 and counts[3] == 1 and counts[11] == 1

    ibf.clear(0)
    assert not ibf.containment_agent().bulk_contains(126)[0]
    ibf.clear([3, 11])
    assert not ibf.containment_agent().bulk_contains(712).any()


def test_ibf_invalid():
    with pytest.raises(Exception, match="bins must be > 0"):
        hibf.IBF(bin_count=0, bin_size=10)
    ibf = hibf.IBF(bin_count=4, bin_size=64)
    with pytest.raises(IndexError):
        ibf.emplace(1, 4)
    with pytest.raises(IndexError):
        ibf.clear([0, 5])


def test_ibf_increase_bins():
    ibf = hibf.IBF(bin_count=4, bin_size=64)
    assert ibf.try_increase_bin_number_to(64)
    assert ibf.bin_count == 64
    assert not ibf.try_increase_bin_number_to(65)
    ibf.increase_bin_number_to(65)
    assert ibf.bin_count == 65


def test_ibf_agents_after_increasing_bins():
    ibf = hibf.IBF(bin_count=4, bin_size=1024)
    ibf.emplace([1, 2, 3], 0)
    containment, counting, membership = ibf.containment_agent(), ibf.counting_agent(), ibf.membership_agent()
    filled_bins = [0]

    # 64 bins fit into the allocated technical bins; 100 000 bins require a reallocation.
    for increase, new_bin_count in ((ibf.try_increase_bin_number_to, 64), (ibf.increase_bin_number_to, 100_000)):
        increase(new_bin_count)
        last_bin = new_bin_count - 1
        ibf.emplace([1, 2, 3], last_bin)
        filled_bins.append(last_bin)

        contained = containment.bulk_contains(2)
        assert contained.shape == (new_bin_count,)
        assert contained[0] and contained[last_bin]

        counts = counting.bulk_count([1, 2, 3])
        assert counts.shape == (new_bin_count,)
        assert counts[0] == 3 and counts[last_bin] == 3

        assert list(membership.membership_for([1, 2, 3], 3)) == filled_bins


def test_ibf_occupancy():
    ibf = hibf.IBF(bin_count=4, bin_size=1024, track_occupancy=True)
    assert ibf.track_occupancy
    ibf.emplace([1, 2, 3], 1)
    ibf.emplace([1, 2, 3], 2)
    assert ibf.occupancy[:4] == [0, 3, 3, 0]

    ibf.clear(1)
    ibf.clear([2])
    assert ibf.occupancy[:4] == [0, 0, 0, 0]


def test_ibf_from_config(user_bins):
    data = user_bins[:30]
    ibf = hibf.IBF(hibf.Config(data, maximum_fpr=0.01))
    assert ibf.bin_count == 30
    for user_bin_id, values in enumerate(data):
        assert user_bin_id in ibf.membership_for(values, len(values))

    results = ibf.batch_membership_for(data, [len(v) for v in data], threads=2)
    assert all(i in r for i, r in enumerate(results))


def test_ibf_serialisation(tmp_path):
    ibf = hibf.IBF(bin_count=8, bin_size=1024)
    ibf.emplace([1, 2, 3], 5)
    assert pickle.loads(pickle.dumps(ibf)) == ibf
    assert copy.deepcopy(ibf) == ibf
    ibf.save(tmp_path / "ibf")
    assert hibf.IBF.load(tmp_path / "ibf") == ibf


# --------------------------------------------------------------------------------------------------------------------
# HyperLogLog
# --------------------------------------------------------------------------------------------------------------------


def hashed(values):
    """A cheap 64-bit mixing function; HyperLogLog expects hashed values."""
    values = np.asarray(values, dtype=np.uint64)
    with np.errstate(over="ignore"):
        values = (values ^ (values >> np.uint64(33))) * np.uint64(0xFF51AFD7ED558CCD)
        values = (values ^ (values >> np.uint64(33))) * np.uint64(0xC4CEB9FE1A85EC53)
    return values ^ (values >> np.uint64(33))


def test_hyperloglog():
    sketch = hibf.HyperLogLog(12)
    assert sketch.data_size == 4096
    sketch.add(hashed(np.arange(10_000)))
    assert sketch.estimate() == pytest.approx(10_000, rel=0.05)

    other = hibf.HyperLogLog(12)
    other.add(hashed(np.arange(5_000, 20_000)))
    assert sketch.merge_and_estimate(other) == pytest.approx(20_000, rel=0.05)

    restored = pickle.loads(pickle.dumps(sketch))
    assert restored.estimate() == sketch.estimate()

    sketch.reset()
    assert sketch.estimate() == 0

    for num_bits in (4, 33, 40):  # 40 would allocate 1 TiB before the library checks the value
        with pytest.raises(ValueError, match="num_bits must be in"):
            hibf.HyperLogLog(num_bits)


def test_hyperloglog_merge_requires_same_size():
    small, big = hibf.HyperLogLog(5), hibf.HyperLogLog(16)
    big.add(hashed(np.arange(1000)))
    for sketch, other in ((big, small), (small, big)):
        for merge in (sketch.merge, sketch.merge_and_estimate):
            with pytest.raises(ValueError, match="same num_bits"):
                merge(other)
    assert big.estimate() == pytest.approx(1000, rel=0.05)
