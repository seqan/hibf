// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Python bindings for the HIBF library via nanobind.
 */

#include <atomic>      // for atomic
#include <cerrno>      // for errno
#include <concepts>    // for same_as
#include <cstddef>     // for size_t
#include <cstdint>     // for uint64_t, uint16_t, uint32_t
#include <exception>   // for exception_ptr, current_exception, rethrow_exception
#include <filesystem>  // for path
#include <fstream>     // for ifstream, ofstream
#include <memory>      // for addressof, make_unique, unique_ptr
#include <mutex>       // for mutex, lock_guard
#include <span>        // for span
#include <sstream>     // for ostringstream, istringstream
#include <stdexcept>   // for invalid_argument
#include <streambuf>   // for streambuf
#include <string>      // for string, to_string
#include <string_view> // for string_view
#include <tuple>       // for tuple
#include <utility>     // for move
#include <variant>     // for variant, get_if
#include <vector>      // for vector

#include <cereal/archives/binary.hpp> // for BinaryInputArchive, BinaryOutputArchive
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <hibf/config.hpp>                                // for config, insert_iterator
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter
#include <hibf/interleaved_bloom_filter.hpp>              // for interleaved_bloom_filter, bin_index, ...
#include <hibf/layout/compute_layout.hpp>                 // for compute_layout
#include <hibf/layout/layout.hpp>                         // for layout
#include <hibf/misc/iota_vector.hpp>                      // for iota_vector
#include <hibf/misc/timer.hpp>                            // for concurrent_timer
#include <hibf/sketch/compute_sketches.hpp>               // for compute_sketches
#include <hibf/sketch/estimate_kmer_counts.hpp>           // for estimate_kmer_counts
#include <hibf/sketch/hyperloglog.hpp>                    // for hyperloglog
#include <hibf/version.hpp>                               // for hibf_version_cstring

namespace nb = nanobind;
using namespace nb::literals;

namespace
{

using ibf_t = seqan::hibf::interleaved_bloom_filter;
using hibf_t = seqan::hibf::hierarchical_interleaved_bloom_filter;
using layout_t = seqan::hibf::layout::layout;
using hyperloglog_t = seqan::hibf::sketch::hyperloglog;

//!\brief A read-only, contiguous, one-dimensional array of uint64_t values.
using u64_array = nb::ndarray<uint64_t const, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

template <typename value_t>
using numpy_array = nb::ndarray<nb::numpy, value_t, nb::ndim<1>>;

bool accept_any(PyObject *) noexcept
{
    return true;
}

//!\brief Accepts any object. Arguments of this type are converted via as_u64_array; the name is used in the stubs.
class array_like : public nb::object
{
    NB_OBJECT_DEFAULT(array_like, object, "numpy.typing.ArrayLike", accept_any)
};

//!\brief Accepts any object. Arguments of this type are converted via numpy.dtype; the name is used in the stubs.
class dtype_like : public nb::object
{
    NB_OBJECT_DEFAULT(dtype_like, object, "numpy.typing.DTypeLike", accept_any)
};

// ---------------------------------------------------------------------------------------------------------------------
// Conversion helpers
// ---------------------------------------------------------------------------------------------------------------------

nb::module_ numpy()
{
    return nb::module_::import_("numpy");
}

//!\brief Wraps the values in an array that owns them.
u64_array to_u64_array(std::vector<uint64_t> values)
{
    auto owned = std::make_unique<std::vector<uint64_t>>(std::move(values));
    uint64_t const * const data = owned->data();
    size_t const size = owned->size();

    nb::capsule owner(owned.get(),
                      [](void * ptr) noexcept
                      {
                          delete static_cast<std::vector<uint64_t> *>(ptr);
                      });
    owned.release();

    return u64_array(data, {size}, owner);
}

//!\brief Converts a Python integer (or any object implementing `__index__`) to uint64_t. Raises for other values.
uint64_t index_to_u64(nb::handle obj)
{
    nb::object const index = nb::steal(PyNumber_Index(obj.ptr()));
    if (!index.is_valid())
        throw nb::python_error{};

    unsigned long long const value = PyLong_AsUnsignedLongLong(index.ptr());
    if (value == static_cast<unsigned long long>(-1) && PyErr_Occurred())
        throw nb::python_error{};

    return value;
}

/*!\brief Converts an arbitrary Python object into a one-dimensional uint64 array.
 * \details
 * A C-contiguous uint64 array (anything supporting the buffer protocol or DLPack) is used without copying.
 * Other arrays, buffers and NumPy scalars must have an integer dtype and are converted by NumPy. Signed values are
 * reinterpreted as two's complement, e.g., -1 becomes 2^64 - 1.
 * Python integers and iterables of them (e.g., lists, generators, sets) are converted value by value. Each value must
 * be an integer in [0, 2^64).
 * Floats are rejected in both cases instead of being truncated. Requires the GIL.
 */
u64_array as_u64_array(nb::handle obj)
{
    u64_array result;
    if (nb::try_cast(obj, result, /*convert=*/false))
        return result;

    if (PyObject_CheckBuffer(obj.ptr()) || nb::hasattr(obj, "__array__") || nb::hasattr(obj, "__array_interface__")
        || nb::hasattr(obj, "__dlpack__"))
    {
        nb::module_ np = numpy();
        nb::object const array = np.attr("asarray")(obj);
        nb::object const dtype = array.attr("dtype");
        std::string const kind = nb::cast<std::string>(dtype.attr("kind"));

        if (kind != "i" && kind != "u")
            throw nb::type_error(
                ("Expected an array of integers, got an array of dtype " + nb::cast<std::string>(nb::str(dtype)) + ".")
                    .c_str());

        if (nb::cast<size_t>(array.attr("ndim")) > 1u)
            throw std::invalid_argument{"Expected a one-dimensional array of integers."};

        return nb::cast<u64_array>(np.attr("ascontiguousarray")(array, "dtype"_a = "uint64"));
    }

    if (PyIndex_Check(obj.ptr()))
        return to_u64_array({index_to_u64(obj)});

    if (!nb::hasattr(obj, "__iter__") && !nb::hasattr(obj, "__getitem__"))
        throw nb::type_error(("Expected an integer, an array of integers, or an iterable of integers, got "
                              + nb::cast<std::string>(obj.type().attr("__name__")) + ".")
                                 .c_str());

    std::vector<uint64_t> values;
    if (nb::hasattr(obj, "__len__"))
        values.reserve(nb::len(obj));
    for (nb::handle item : obj)
        values.push_back(index_to_u64(item));

    return to_u64_array(std::move(values));
}

std::span<uint64_t const> as_span(u64_array const & array)
{
    return {array.data(), array.size()};
}

//!\brief Copies a sized range into a newly allocated NumPy array.
template <typename value_t, typename range_t>
numpy_array<value_t> to_numpy(range_t const & range)
{
    size_t const size = range.size();
    auto * data = new value_t[size];
    for (size_t i = 0; i < size; ++i)
        data[i] = static_cast<value_t>(range[i]);

    nb::capsule owner(data,
                      [](void * ptr) noexcept
                      {
                          delete[] static_cast<value_t *>(ptr);
                      });

    return numpy_array<value_t>(data, {size}, owner);
}

enum class count_type
{
    uint16,
    uint32,
    uint64
};

//!\brief Maps anything accepted by `numpy.dtype` to one of the supported counter types.
count_type parse_count_dtype(nb::handle dtype)
{
    std::string const name = nb::cast<std::string>(numpy().attr("dtype")(dtype).attr("name"));

    if (name == "uint16")
        return count_type::uint16;
    if (name == "uint32")
        return count_type::uint32;
    if (name == "uint64")
        return count_type::uint64;

    throw std::invalid_argument{"Unsupported dtype '" + name + "'. Supported: uint16, uint32, uint64."};
}

//!\brief Resolves a threshold that is either a single value or one value per query.
std::vector<uint16_t> resolve_thresholds(std::variant<uint16_t, std::vector<uint16_t>> const & threshold,
                                         size_t const number_of_queries)
{
    if (uint16_t const * value = std::get_if<uint16_t>(&threshold))
        return std::vector<uint16_t>(number_of_queries, *value);

    auto const & values = std::get<std::vector<uint16_t>>(threshold);
    if (values.size() != number_of_queries)
        throw std::invalid_argument{"The number of thresholds (" + std::to_string(values.size())
                                    + ") must match the number of queries (" + std::to_string(number_of_queries)
                                    + ")."};
    return values;
}

// ---------------------------------------------------------------------------------------------------------------------
// Serialisation helpers (cereal binary archives)
// ---------------------------------------------------------------------------------------------------------------------

//!\brief A read-only stream buffer over existing memory. Avoids copying the pickled bytes.
struct memory_buffer : public std::streambuf
{
    memory_buffer(char const * data, size_t const size)
    {
        char * begin = const_cast<char *>(data);
        setg(begin, begin, begin + size);
    }
};

[[noreturn]] void raise_os_error(std::filesystem::path const & path)
{
    PyErr_SetFromErrnoWithFilename(PyExc_OSError, path.string().c_str());
    throw nb::python_error{};
}

template <typename object_t>
nb::bytes to_bytes(object_t const & object)
{
    std::string buffer;
    {
        nb::gil_scoped_release nogil{};
        std::ostringstream stream{std::ios::binary};
        {
            cereal::BinaryOutputArchive archive{stream};
            archive(object);
        }
        buffer = std::move(stream).str();
    }
    return nb::bytes(buffer.data(), buffer.size());
}

template <typename object_t>
object_t from_bytes(nb::bytes const & data)
{
    object_t object{};
    char const * ptr = data.c_str();
    size_t const size = data.size();
    {
        nb::gil_scoped_release nogil{};
        memory_buffer buffer{ptr, size};
        std::istream stream{&buffer};
        cereal::BinaryInputArchive archive{stream};
        archive(object);
    }
    return object;
}

template <typename object_t>
void save_to_file(object_t const & object, std::filesystem::path const & path)
{
    std::ofstream stream{path, std::ios::binary};
    if (!stream.good())
        raise_os_error(path);

    nb::gil_scoped_release nogil{};
    cereal::BinaryOutputArchive archive{stream};
    archive(object);
}

template <typename object_t>
object_t load_from_file(std::filesystem::path const & path)
{
    std::ifstream stream{path, std::ios::binary};
    if (!stream.good())
        raise_os_error(path);

    object_t object{};
    {
        nb::gil_scoped_release nogil{};
        cereal::BinaryInputArchive archive{stream};
        archive(object);
    }
    return object;
}

//!\brief Adds pickle, copy, save and load support to a class that is serialisable with cereal.
template <typename object_t, typename... extra_t>
void add_serialisation(nb::class_<object_t, extra_t...> & cls)
{
    cls.def("__getstate__", &to_bytes<object_t>)
        .def("__setstate__",
             [](object_t & self, nb::bytes const & state)
             {
                 // Deserialise first, so that `self` is only initialised if loading succeeds.
                 new (&self) object_t{from_bytes<object_t>(state)};
             })
        .def("__copy__",
             [](object_t const & self)
             {
                 return object_t{self};
             })
        .def(
            "__deepcopy__",
            [](object_t const & self, nb::handle)
            {
                return object_t{self};
            },
            "memo"_a)
        .def("save",
             &save_to_file<object_t>,
             "path"_a,
             "Serialises the object into a binary file. The file can be read with ``load()``.")
        .def_static("load",
                    &load_from_file<object_t>,
                    "path"_a,
                    "Loads an object from a binary file that was written with ``save()``.");
}

// ---------------------------------------------------------------------------------------------------------------------
// Config and user input
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief The Python view of seqan::hibf::config.
 * \details
 * seqan::hibf::config::input_fn cannot be exposed directly, because a std::function owning a Python object may be
 * copied or destroyed by library threads that do not hold the GIL. Instead, the Python object is stored here and
 * wrapped by an input_adapter for the duration of a construction.
 *
 * A number_of_user_bins of 0 means that the number of user bins is inferred from the input whenever it is needed.
 * The inferred value is only stored by validate_and_set_defaults(), like other defaults.
 */
struct py_config : public seqan::hibf::config
{
    nb::object input{nb::none()};

    //!\brief Returns number_of_user_bins or, if unset and the input is a sized sequence, the length of the input.
    size_t inferred_number_of_user_bins() const
    {
        if (number_of_user_bins == 0u && !input.is_none() && !PyCallable_Check(input.ptr())
            && PyObject_HasAttrString(input.ptr(), "__len__"))
            return nb::len(input);
        return number_of_user_bins;
    }

    void infer_number_of_user_bins()
    {
        number_of_user_bins = inferred_number_of_user_bins();
    }

    //!\brief Returns a copy of the library config with the inferred number of user bins.
    seqan::hibf::config resolved() const
    {
        seqan::hibf::config config{static_cast<seqan::hibf::config const &>(*this)};
        config.number_of_user_bins = inferred_number_of_user_bins();
        return config;
    }
};

void validate_input(nb::handle input)
{
    if (input.is_none() || PyCallable_Check(input.ptr()))
        return;

    if (!PyObject_HasAttrString(input.ptr(), "__getitem__"))
        throw nb::type_error("Config.input must be a callable `f(user_bin_id) -> values` or an indexable sequence "
                             "`input[user_bin_id] -> values`.");
}

/*!\brief Provides the user's Python input to the library.
 * \details
 * The library calls seqan::hibf::config::input_fn from OpenMP worker threads and expects it not to throw: an
 * exception escaping an OpenMP region terminates the process. Hence, the adapter acquires the GIL and records the
 * first error. Afterwards, the input is no longer requested. Callers must check for a recorded error via
 * rethrow_if_failed() after the library returns.
 *
 * The HIBF construction requires non-empty user bins; inserting an empty user bin divides by zero. Unless
 * `allow_empty` is set, an empty user bin is recorded as an error and a placeholder value is inserted instead. This
 * also applies to the calls after an error. The resulting index is discarded.
 *
 * The adapter holds a reference to the Python object: the input may be reassigned (e.g., by the input itself) while
 * the library is running. The adapter must be destroyed while holding the GIL.
 */
class input_adapter
{
public:
    input_adapter(nb::handle source, bool const allow_empty) :
        source{nb::borrow(source)},
        is_callable{PyCallable_Check(source.ptr()) == 1},
        allow_empty{allow_empty}
    {}

    input_adapter(input_adapter const &) = delete;
    input_adapter & operator=(input_adapter const &) = delete;

    void operator()(size_t const user_bin_id, seqan::hibf::insert_iterator it) noexcept
    {
        size_t const inserted = failed.load(std::memory_order_relaxed) ? 0u : insert(user_bin_id, it);

        if (inserted == 0u && !allow_empty)
        {
            record(std::make_exception_ptr(
                std::invalid_argument{"User bin " + std::to_string(user_bin_id)
                                      + " is empty. Each user bin must contain at least one value."}));
            it = 0u;
        }
    }

    void rethrow_if_failed() const
    {
        if (error)
            std::rethrow_exception(error);
    }

private:
    //!\brief Inserts the values of a user bin and returns their number. Returns 0 on error.
    size_t insert(size_t const user_bin_id, seqan::hibf::insert_iterator & it) noexcept
    {
        try
        {
            nb::gil_scoped_acquire gil{};

            // Only has an effect on the main thread. Allows interrupting long constructions via Ctrl+C.
            if (PyErr_CheckSignals() != 0)
                throw nb::python_error{};

            nb::object data = is_callable ? source(user_bin_id) : source[nb::int_(user_bin_id)];
            u64_array const values = as_u64_array(data);
            size_t const size = values.size();

            nb::gil_scoped_release nogil{};
            for (uint64_t const value : as_span(values))
                it = value;

            return size;
        }
        catch (...)
        {
            record(std::current_exception());
            return 0u;
        }
    }

    //!\brief Records an error. Only the first error is kept.
    void record(std::exception_ptr exception) noexcept
    {
        std::lock_guard lock{mutex};
        if (!error)
            error = std::move(exception);
        failed.store(true, std::memory_order_relaxed);
    }

    nb::object source;
    bool is_callable{};
    bool allow_empty{};
    std::atomic<bool> failed{false};
    std::mutex mutex;
    std::exception_ptr error;
};

//!\brief Returns a copy of the config whose input_fn forwards to the adapter.
seqan::hibf::config make_cpp_config(py_config const & py_cfg, input_adapter & adapter)
{
    seqan::hibf::config config = py_cfg.resolved();

    config.input_fn = [&adapter](size_t const user_bin_id, seqan::hibf::insert_iterator && it)
    {
        adapter(user_bin_id, std::move(it));
    };

    return config;
}

nb::handle require_input(py_config const & config)
{
    if (config.input.is_none())
        throw std::invalid_argument{"[HIBF CONFIG ERROR] You did not set the required Config.input."};
    return config.input;
}

std::string config_to_string(seqan::hibf::config const & config)
{
    std::ostringstream stream;
    config.write_to(stream);
    return std::move(stream).str();
}

py_config config_from_string(std::string const & text)
{
    py_config config{};
    std::istringstream stream{text};
    config.read_from(stream);
    return config;
}

// ---------------------------------------------------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------------------------------------------------

struct layout_timers
{
    seqan::hibf::concurrent_timer compute_sketches{};
    seqan::hibf::concurrent_timer union_estimation{};
    seqan::hibf::concurrent_timer rearrangement{};
    seqan::hibf::concurrent_timer dp_algorithm{};
};

/*!\brief Computes the layout for the given config. Mirrors the first half of the HIBF constructor.
 * \details
 * Running the stages individually allows raising errors in the input, e.g., empty user bins, before the index is
 * built. Must be called without the GIL. `config` must have been validated.
 */
layout_t compute_layout_impl(seqan::hibf::config const & config, input_adapter & adapter, layout_timers & timers)
{
    std::vector<hyperloglog_t> sketches{};
    std::vector<size_t> kmer_counts{};

    timers.compute_sketches.start();
    seqan::hibf::sketch::compute_sketches(config, sketches);
    adapter.rethrow_if_failed();
    seqan::hibf::sketch::estimate_kmer_counts(sketches, kmer_counts);
    timers.compute_sketches.stop();

    timers.dp_algorithm.start();
    layout_t layout = seqan::hibf::layout::compute_layout(config,
                                                          kmer_counts,
                                                          sketches,
                                                          seqan::hibf::iota_vector(config.number_of_user_bins),
                                                          timers.union_estimation,
                                                          timers.rearrangement);
    timers.dp_algorithm.stop();

    return layout;
}

layout_t compute_layout_for(py_config const & py_cfg)
{
    input_adapter adapter{require_input(py_cfg), /*allow_empty=*/false};
    seqan::hibf::config config = make_cpp_config(py_cfg, adapter);
    config.validate_and_set_defaults();
    layout_timers timers{};

    nb::gil_scoped_release nogil{};
    return compute_layout_impl(config, adapter, timers);
}

/*!\brief Checks that the layout's user bins are exactly the config's user bins.
 * \details
 * The library trusts the layout. An empty layout crashes the construction, and user bin ids that are out of range are
 * reported by queries and written out of bounds by counting agents.
 */
void check_layout_matches(layout_t const & layout, size_t const number_of_user_bins)
{
    if (layout.user_bins.size() != number_of_user_bins)
        throw std::invalid_argument{"The layout contains " + std::to_string(layout.user_bins.size())
                                    + " user bins, but the config has " + std::to_string(number_of_user_bins) + "."};

    std::vector<bool> seen(number_of_user_bins);
    for (auto const & user_bin : layout.user_bins)
    {
        if (user_bin.idx >= number_of_user_bins || seen[user_bin.idx])
            throw std::invalid_argument{"The layout's user bin ids must be 0, ..., "
                                        + std::to_string(number_of_user_bins - 1) + ", each exactly once. Found "
                                        + std::to_string(user_bin.idx) + "."};
        seen[user_bin.idx] = true;
    }
}

hibf_t build_hibf(py_config const & py_cfg, layout_t const * layout)
{
    input_adapter adapter{require_input(py_cfg), /*allow_empty=*/false};
    seqan::hibf::config config = make_cpp_config(py_cfg, adapter);
    config.validate_and_set_defaults();

    if (layout != nullptr)
        check_layout_matches(*layout, config.number_of_user_bins);

    nb::gil_scoped_release nogil{};
    try
    {
        if (layout != nullptr)
        {
            hibf_t hibf{config, *layout};
            adapter.rethrow_if_failed();
            return hibf;
        }

        layout_timers timers{};
        layout_t const computed_layout = compute_layout_impl(config, adapter, timers);
        hibf_t hibf{config, computed_layout};
        adapter.rethrow_if_failed();

        hibf.layout_compute_sketches_timer = std::move(timers.compute_sketches);
        hibf.layout_union_estimation_timer = std::move(timers.union_estimation);
        hibf.layout_rearrangement_timer = std::move(timers.rearrangement);
        hibf.layout_dp_algorithm_timer = std::move(timers.dp_algorithm);
        return hibf;
    }
    catch (...)
    {
        // A failing callback is the root cause of any subsequent error.
        adapter.rethrow_if_failed();
        throw;
    }
}

ibf_t build_ibf(py_config const & py_cfg, size_t const max_bin_elements)
{
    input_adapter adapter{require_input(py_cfg), /*allow_empty=*/true};
    seqan::hibf::config config = make_cpp_config(py_cfg, adapter);

    nb::gil_scoped_release nogil{};
    try
    {
        ibf_t ibf{config, max_bin_elements};
        adapter.rethrow_if_failed();
        return ibf;
    }
    catch (...)
    {
        adapter.rethrow_if_failed();
        throw;
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Queries
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief An agent as exposed to Python.
 * \details
 * IBF agents size their buffers by the number of bins. Using such an agent after the number of bins changed, e.g., via
 * interleaved_bloom_filter::increase_bin_number_to(), accesses memory out of bounds. Hence, the agent is rebuilt if the
 * number of bins changed since it was created. An HIBF cannot change; its agents are wrapped for uniformity.
 */
template <typename filter_t, typename agent_t>
class py_agent
{
public:
    py_agent() = delete;
    py_agent(py_agent const &) = default;
    py_agent & operator=(py_agent const &) = default;
    py_agent(py_agent &&) = default;
    py_agent & operator=(py_agent &&) = default;
    ~py_agent() = default;

    explicit py_agent(filter_t const & filter) :
        filter{std::addressof(filter)},
        bin_count{number_of_bins(filter)},
        agent{filter}
    {}

    //!\brief Returns the agent. Rebuilds it if the number of bins of the filter changed.
    agent_t & get()
    {
        if (size_t const current = number_of_bins(*filter); current != bin_count)
        {
            agent = agent_t{*filter};
            bin_count = current;
        }
        return agent;
    }

private:
    static size_t number_of_bins(filter_t const & filter)
    {
        if constexpr (std::same_as<filter_t, ibf_t>)
            return filter.bin_count();
        else
            return 0u;
    }

    filter_t const * filter{nullptr};
    size_t bin_count{};
    agent_t agent;
};

using ibf_containment_agent = py_agent<ibf_t, ibf_t::containment_agent_type>;
using ibf_membership_agent = py_agent<ibf_t, ibf_t::membership_agent_type>;
using hibf_membership_agent = py_agent<hibf_t, hibf_t::membership_agent_type>;

template <typename filter_t, typename value_t>
using counting_agent = py_agent<filter_t, typename filter_t::template counting_agent_type<value_t>>;

template <typename agent_t>
numpy_array<uint64_t> membership_for(agent_t & agent, array_like const & values, uint16_t const threshold)
{
    u64_array const array = as_u64_array(values);
    return to_numpy<uint64_t>(agent.membership_for(as_span(array), threshold));
}

//!\brief Queries many value sets in parallel. Each thread uses its own membership agent.
template <typename filter_t>
nb::typed<nb::list, numpy_array<uint64_t>>
batch_membership_for(filter_t const & filter,
                     nb::typed<nb::sequence, array_like> const & queries,
                     std::variant<uint16_t, std::vector<uint16_t>> const & threshold,
                     size_t const threads)
{
    if (threads == 0u)
        throw std::invalid_argument{"threads must be > 0."};

    size_t const number_of_queries = nb::len(queries);
    std::vector<uint16_t> const thresholds = resolve_thresholds(threshold, number_of_queries);

    std::vector<u64_array> arrays;
    arrays.reserve(number_of_queries);
    for (nb::handle query : queries)
        arrays.push_back(as_u64_array(query));

    std::vector<std::vector<uint64_t>> results(number_of_queries);
    {
        nb::gil_scoped_release nogil{};
#pragma omp parallel num_threads(threads)
        {
            auto agent = filter.membership_agent();
#pragma omp for schedule(dynamic)
            for (size_t i = 0; i < number_of_queries; ++i)
            {
                auto const & result = agent.membership_for(as_span(arrays[i]), thresholds[i]);
                results[i].assign(result.begin(), result.end());
            }
        }
    }

    nb::typed<nb::list, numpy_array<uint64_t>> output;
    for (auto const & result : results)
        output.append(to_numpy<uint64_t>(result));
    return output;
}

//!\brief Registers `name` as a counting agent for a filter. Counting agents are templated on the counter type.
template <typename agent_t, typename value_t, typename... bulk_count_args_t>
void bind_counting_agent(nb::handle scope, char const * name, char const * doc)
{
    auto cls = nb::class_<agent_t>(scope, name, doc);

    if constexpr (sizeof...(bulk_count_args_t) == 0u)
    {
        cls.def(
            "bulk_count",
            [](agent_t & agent, array_like const & values)
            {
                u64_array const array = as_u64_array(values);
                return to_numpy<value_t>(agent.get().bulk_count(as_span(array)));
            },
            "values"_a,
            "Returns, for each bin, how many of the given values are contained in it.");
    }
    else
    {
        cls.def(
            "bulk_count",
            [](agent_t & agent, array_like const & values, size_t const threshold)
            {
                if (threshold == 0u)
                    throw std::invalid_argument{"threshold must be > 0."};

                u64_array const array = as_u64_array(values);
                return to_numpy<value_t>(agent.get().bulk_count(as_span(array), threshold));
            },
            "values"_a,
            "threshold"_a = 1u,
            "Returns, for each user bin, how many of the given values are contained in it. "
            "Counts below ``threshold`` may be reported as 0, because subtrees of the hierarchy are skipped if "
            "their count is below ``threshold``.");
    }
}

//!\brief Creates a counting agent with the requested counter type. The agent keeps the filter alive.
template <typename filter_t>
nb::object make_counting_agent(filter_t const & filter, nb::handle dtype)
{
    switch (parse_count_dtype(dtype))
    {
    case count_type::uint16:
        return nb::cast(counting_agent<filter_t, uint16_t>{filter});
    case count_type::uint32:
        return nb::cast(counting_agent<filter_t, uint32_t>{filter});
    default:
        return nb::cast(counting_agent<filter_t, uint64_t>{filter});
    }
}

void check_bin(ibf_t const & ibf, size_t const bin)
{
    if (bin >= ibf.bin_count())
        throw nb::index_error(("Bin index " + std::to_string(bin) + " is out of range for an IBF with "
                               + std::to_string(ibf.bin_count()) + " bins.")
                                  .c_str());
}

// ---------------------------------------------------------------------------------------------------------------------
// Bindings
// ---------------------------------------------------------------------------------------------------------------------

void bind_config(nb::module_ & m)
{
    nb::class_<py_config>(m,
                          "Config",
                          "The configuration used to build an HIBF or IBF.\n\n"
                          "``input`` provides the values of each user bin. It is either\n\n"
                          "* a sequence with ``input[user_bin_id]`` returning the values of a user bin, or\n"
                          "* a callable ``input(user_bin_id)`` returning the values of a user bin.\n\n"
                          "The values of a user bin are an iterable of unsigned 64-bit integers, ideally a NumPy "
                          "``uint64`` array. The input is requested multiple times per user bin and possibly from "
                          "multiple threads (see ``threads``), so it must return the same values on each request.\n\n"
                          "If ``input`` is a sequence and ``number_of_user_bins`` is not set, ``len(input)`` is "
                          "used.")
        .def(
            "__init__",
            [](py_config * self,
               nb::object input,
               size_t number_of_user_bins,
               size_t number_of_hash_functions,
               double maximum_fpr,
               double relaxed_fpr,
               size_t threads,
               uint8_t sketch_bits,
               size_t tmax,
               double empty_bin_fraction,
               bool track_occupancy,
               double alpha,
               double max_rearrangement_ratio,
               bool disable_estimate_union,
               bool disable_rearrangement)
            {
                validate_input(input);
                py_config * config = new (self) py_config{};
                config->input = std::move(input);
                config->number_of_user_bins = number_of_user_bins;
                config->number_of_hash_functions = number_of_hash_functions;
                config->maximum_fpr = maximum_fpr;
                config->relaxed_fpr = relaxed_fpr;
                config->threads = threads;
                config->sketch_bits = sketch_bits;
                config->tmax = tmax;
                config->empty_bin_fraction = empty_bin_fraction;
                config->track_occupancy = track_occupancy;
                config->alpha = alpha;
                config->max_rearrangement_ratio = max_rearrangement_ratio;
                config->disable_estimate_union = disable_estimate_union;
                config->disable_rearrangement = disable_rearrangement;
            },
            "input"_a = nb::none(),
            nb::kw_only(),
            "number_of_user_bins"_a = 0u,
            "number_of_hash_functions"_a = 2u,
            "maximum_fpr"_a = 0.05,
            "relaxed_fpr"_a = 0.3,
            "threads"_a = 1u,
            "sketch_bits"_a = 12u,
            "tmax"_a = 0u,
            "empty_bin_fraction"_a = 0.0,
            "track_occupancy"_a = false,
            "alpha"_a = 1.2,
            "max_rearrangement_ratio"_a = 0.5,
            "disable_estimate_union"_a = false,
            "disable_rearrangement"_a = false)
        .def_prop_rw(
            "input",
            [](py_config const & self)
            {
                return self.input;
            },
            [](py_config & self, nb::object input)
            {
                validate_input(input);
                self.input = std::move(input);
            },
            "A sequence or callable providing the values of each user bin.")
        .def_prop_rw(
            "number_of_user_bins",
            [](py_config const & self)
            {
                return self.inferred_number_of_user_bins();
            },
            [](py_config & self, size_t const number_of_user_bins)
            {
                self.number_of_user_bins = number_of_user_bins;
            },
            "The number of user bins. If set to 0 (the default) and ``input`` is a sequence, ``len(input)`` is "
            "used.")
        .def_rw("number_of_hash_functions",
                &py_config::number_of_hash_functions,
                "The number of hash functions for the underlying Bloom filters. Must be in [1, 5].")
        .def_rw("maximum_fpr", &py_config::maximum_fpr, "The desired maximum false positive rate of the IBFs.")
        .def_rw("relaxed_fpr",
                &py_config::relaxed_fpr,
                "The desired relaxed false positive rate of IBFs that are not the last level.")
        .def_rw("threads", &py_config::threads, "The number of threads to use for construction.")
        .def_rw("sketch_bits", &py_config::sketch_bits, "The number of bits for the HyperLogLog sketches.")
        .def_rw("tmax", &py_config::tmax, "The maximum number of technical bins on each level. 0: sqrt(user bins).")
        .def_rw("empty_bin_fraction",
                &py_config::empty_bin_fraction,
                "The fraction of empty bins in the top-level IBF (for later insertions).")
        .def_rw("track_occupancy", &py_config::track_occupancy, "Whether to track the occupancy of each bin.")
        .def_rw("alpha", &py_config::alpha, "The weight of merged bins in the layout algorithm.")
        .def_rw("max_rearrangement_ratio",
                &py_config::max_rearrangement_ratio,
                "The maximum ratio of user bins to rearrange. Only used if rearrangement is enabled.")
        .def_rw("disable_estimate_union",
                &py_config::disable_estimate_union,
                "Whether to disable union estimation in the layout. Also disables rearrangement.")
        .def_rw("disable_rearrangement",
                &py_config::disable_rearrangement,
                "Whether to disable rearrangement of user bins in the layout.")
        .def(
            "validate_and_set_defaults",
            [](py_config & self)
            {
                // input_fn is checked by the library but never called here.
                self.infer_number_of_user_bins();
                require_input(self);
                self.input_fn = [](size_t const, seqan::hibf::insert_iterator &&) {};
                try
                {
                    self.validate_and_set_defaults();
                }
                catch (...)
                {
                    self.input_fn = {};
                    throw;
                }
                self.input_fn = {};
            },
            "Checks the config for errors and sets defaults, e.g., for ``tmax``.")
        .def(
            "to_string",
            [](py_config const & self)
            {
                return config_to_string(self.resolved());
            },
            "Returns the config in the textual format used by layout files. ``input`` is not included.")
        .def_static("from_string", &config_from_string, "text"_a, "Parses a config written by ``to_string()``.")
        .def(
            "__eq__",
            [](py_config const & self, py_config const & other)
            {
                return self.resolved() == other.resolved();
            },
            nb::is_operator(),
            "Two configs are equal if all options, except ``input``, are equal.")
        .def("__getstate__",
             [](py_config const & self)
             {
                 // Without inferring number_of_user_bins: the restored config infers it from its input.
                 return std::make_pair(config_to_string(self), self.input);
             })
        .def("__setstate__",
             [](py_config & self, std::pair<std::string, nb::object> const & state)
             {
                 py_config config = config_from_string(state.first);
                 config.input = state.second;
                 new (&self) py_config{std::move(config)};
             })
        .def("__copy__",
             [](py_config const & self)
             {
                 return py_config{self};
             })
        .def(
            "__deepcopy__",
            [](py_config const & self, nb::handle memo)
            {
                py_config config{self};
                config.input = nb::module_::import_("copy").attr("deepcopy")(self.input, memo);
                return config;
            },
            "memo"_a)
        .def("__repr__",
             [](py_config const & self)
             {
                 std::string const input_repr = nb::cast<std::string>(nb::repr(self.input));
                 std::ostringstream stream;
                 stream << "Config(input=" << (input_repr.size() > 60 ? input_repr.substr(0, 57) + "..." : input_repr)
                        << ", number_of_user_bins=" << self.inferred_number_of_user_bins()
                        << ", number_of_hash_functions=" << self.number_of_hash_functions
                        << ", maximum_fpr=" << self.maximum_fpr << ", relaxed_fpr=" << self.relaxed_fpr
                        << ", threads=" << self.threads << ", sketch_bits=" << static_cast<int>(self.sketch_bits)
                        << ", tmax=" << self.tmax << ", empty_bin_fraction=" << self.empty_bin_fraction
                        << ", track_occupancy=" << (self.track_occupancy ? "True" : "False") << ", alpha=" << self.alpha
                        << ", max_rearrangement_ratio=" << self.max_rearrangement_ratio
                        << ", disable_estimate_union=" << (self.disable_estimate_union ? "True" : "False")
                        << ", disable_rearrangement=" << (self.disable_rearrangement ? "True" : "False") << ")";
                 return stream.str();
             });
}

void bind_layout(nb::module_ & m)
{
    nb::class_<layout_t>(m,
                         "Layout",
                         "The hierarchical layout of an HIBF, i.e., which user bins are stored in which IBFs.\n\n"
                         "Layouts can be computed with ``compute_layout()`` or read from layout files, e.g., those "
                         "written by chopper.")
        .def(nb::init<>())
        .def_ro("top_level_max_bin_id", &layout_t::top_level_max_bin_id)
        .def_prop_ro(
            "number_of_user_bins",
            [](layout_t const & self)
            {
                return self.user_bins.size();
            },
            "The number of user bins in the layout.")
        .def_prop_ro(
            "number_of_max_bins",
            [](layout_t const & self)
            {
                return self.max_bins.size();
            },
            "The number of merged bins in the layout, i.e., the number of lower-level IBFs.")
        .def(
            "to_string",
            [](layout_t const & self)
            {
                std::ostringstream stream;
                self.write_to(stream);
                return std::move(stream).str();
            },
            "Returns the layout in the textual format used by layout files.")
        .def_static(
            "from_string",
            [](std::string const & text)
            {
                layout_t layout{};
                std::istringstream stream{text};
                layout.read_from(stream);
                return layout;
            },
            "text"_a,
            "Parses a layout written by ``to_string()``.")
        .def(nb::self == nb::self)
        .def("__repr__",
             [](layout_t const & self)
             {
                 return "<Layout number_of_user_bins=" + std::to_string(self.user_bins.size())
                      + " number_of_max_bins=" + std::to_string(self.max_bins.size()) + ">";
             });

    m.def("compute_layout",
          &compute_layout_for,
          "config"_a,
          "Computes the HIBF layout for a config.\n\n"
          "``HierarchicalInterleavedBloomFilter(config)`` does this internally. Computing the layout separately "
          "allows storing it via ``write_layout_file()`` or passing it to "
          "``HierarchicalInterleavedBloomFilter(config, layout)``.");

    m.def(
        "read_layout_file",
        [](std::filesystem::path const & path)
        {
            std::ifstream stream{path};
            if (!stream.good())
                raise_os_error(path);

            py_config config{};
            config.read_from(stream);
            layout_t layout{};
            layout.read_from(stream);
            return std::make_pair(std::move(config), std::move(layout));
        },
        "path"_a,
        "Reads a layout file, e.g., one written by ``write_layout_file()`` or chopper.\n\n"
        "Returns a tuple ``(config, layout)``. ``config.input`` must be set before building an HIBF.");

    m.def(
        "write_layout_file",
        [](std::filesystem::path const & path, py_config const & config, layout_t const & layout)
        {
            std::ofstream stream{path};
            if (!stream.good())
                raise_os_error(path);

            config.resolved().write_to(stream);
            layout.write_to(stream);
        },
        "path"_a,
        "config"_a,
        "layout"_a,
        "Writes a config and layout to a layout file.");
}

void bind_ibf(nb::module_ & m)
{
    auto cls = nb::class_<ibf_t>(m,
                                 "InterleavedBloomFilter",
                                 "The Interleaved Bloom Filter (IBF): a set of equally sized Bloom filters (bins) "
                                 "that can be queried simultaneously.\n\n"
                                 "Methods that modify the IBF are not synchronised. Do not call them while another "
                                 "thread uses the IBF.");

    cls.def(
           "__init__",
           [](ibf_t * self,
              size_t bin_count,
              size_t bin_size,
              size_t hash_function_count,
              double empty_bin_fraction,
              bool track_occupancy)
           {
               new (self) ibf_t{seqan::hibf::bin_count{bin_count},
                                seqan::hibf::bin_size{bin_size},
                                seqan::hibf::hash_function_count{hash_function_count},
                                seqan::hibf::empty_bin_fraction{empty_bin_fraction},
                                seqan::hibf::track_occupancy{track_occupancy}};
           },
           "bin_count"_a,
           "bin_size"_a,
           "hash_function_count"_a = 2u,
           "empty_bin_fraction"_a = 0.0,
           "track_occupancy"_a = false,
           "Creates an empty IBF with ``bin_count`` bins of ``bin_size`` bits each.")
        .def(
            "__init__",
            [](ibf_t * self, py_config const & config, size_t max_bin_elements)
            {
                new (self) ibf_t{build_ibf(config, max_bin_elements)};
            },
            "config"_a,
            "max_bin_elements"_a = 0u,
            "Builds an IBF with one bin per user bin from a config.\n\n"
            "The bin size is chosen such that the bin with the most elements achieves ``config.maximum_fpr``. "
            "The number of elements in the biggest bin is determined from the input, unless "
            "``max_bin_elements`` is given.")
        .def(
            "emplace",
            [](ibf_t & self, array_like const & values, size_t bin)
            {
                check_bin(self, bin);
                u64_array const array = as_u64_array(values);
                nb::gil_scoped_release nogil{};
                for (uint64_t const value : as_span(array))
                    self.emplace(value, seqan::hibf::bin_index{bin});
            },
            "values"_a,
            "bin"_a,
            "Inserts a value or an array of values into a bin.")
        .def(
            "clear",
            [](ibf_t & self, std::variant<size_t, std::vector<size_t>> const & bins)
            {
                std::vector<seqan::hibf::bin_index> indices;
                if (size_t const * bin = std::get_if<size_t>(&bins))
                    indices.push_back({*bin});
                else
                    for (size_t const bin : std::get<std::vector<size_t>>(bins))
                        indices.push_back({bin});

                for (auto const & index : indices)
                    check_bin(self, index.value);

                nb::gil_scoped_release nogil{};
                self.clear(indices);
            },
            "bins"_a,
            "Removes all values from a bin or a list of bins.")
        .def(
            "try_increase_bin_number_to",
            [](ibf_t & self, size_t new_bin_count)
            {
                return self.try_increase_bin_number_to(seqan::hibf::bin_count{new_bin_count});
            },
            "new_bin_count"_a,
            "Increases the number of bins without reallocation, if possible. Returns whether the bins were "
            "increased. Must not be called while another thread uses this IBF.")
        .def(
            "increase_bin_number_to",
            [](ibf_t & self, size_t new_bin_count)
            {
                self.increase_bin_number_to(seqan::hibf::bin_count{new_bin_count});
            },
            "new_bin_count"_a,
            "Increases the number of bins. Requires reallocation if the number of technical bins grows. "
            "Existing agents adapt to the new number of bins.\n\n"
            "Must not be called while another thread uses this IBF: the reallocation frees the memory that thread "
            "may be reading.")
        .def_prop_ro("hash_function_count", &ibf_t::hash_function_count, "The number of hash functions.")
        .def_prop_ro("bin_count", &ibf_t::bin_count, "The number of bins.")
        .def_prop_ro("bin_size", &ibf_t::bin_size, "The size of each bin in bits.")
        .def_prop_ro("bit_size", &ibf_t::bit_size, "The total size of the IBF in bits.")
        .def_ro("occupancy", &ibf_t::occupancy, "The number of values inserted into each technical bin.")
        .def_ro("track_occupancy", &ibf_t::track_occupancy, "Whether occupancy is tracked.")
        .def(
            "containment_agent",
            [](ibf_t const & self)
            {
                return ibf_containment_agent{self};
            },
            nb::keep_alive<0, 1>(),
            "Returns an agent for single-value containment queries.")
        .def(
            "counting_agent",
            [](ibf_t const & self, dtype_like const & dtype)
            {
                return make_counting_agent(self, dtype);
            },
            "dtype"_a = "uint16",
            nb::keep_alive<0, 1>(),
            "Returns an agent that counts the occurrences of values per bin. ``dtype`` is the counter type: "
            "uint16, uint32, or uint64.")
        .def(
            "membership_agent",
            [](ibf_t const & self)
            {
                return ibf_membership_agent{self};
            },
            nb::keep_alive<0, 1>(),
            "Returns an agent that determines the bins containing at least ``threshold`` of the given values.")
        .def(
            "membership_for",
            [](ibf_t const & self, array_like const & values, uint16_t threshold)
            {
                auto agent = self.membership_agent();
                return membership_for(agent, values, threshold);
            },
            "values"_a,
            "threshold"_a,
            "Returns the bins containing at least ``threshold`` of the given values.\n\n"
            "Convenience for ``membership_agent().membership_for(values, threshold)``.")
        .def("batch_membership_for",
             &batch_membership_for<ibf_t>,
             "queries"_a,
             "threshold"_a,
             "threads"_a = 1u,
             "Answers ``membership_for`` for many queries in parallel.\n\n"
             "``threshold`` is either a single value or one value per query. Returns a list with one array of bin "
             "indices per query.")
        .def(nb::self == nb::self)
        .def("__repr__",
             [](ibf_t const & self)
             {
                 return "<InterleavedBloomFilter bin_count=" + std::to_string(self.bin_count())
                      + " bin_size=" + std::to_string(self.bin_size())
                      + " hash_function_count=" + std::to_string(self.hash_function_count()) + ">";
             });

    add_serialisation(cls);

    nb::class_<ibf_containment_agent>(cls, "ContainmentAgent", "Answers single-value containment queries.")
        .def(
            "bulk_contains",
            [](ibf_containment_agent & agent, uint64_t value)
            {
                return to_numpy<bool>(agent.get().bulk_contains(value));
            },
            "value"_a,
            "Returns a boolean array indicating which bins (may) contain the value.");

    bind_counting_agent<counting_agent<ibf_t, uint16_t>, uint16_t>(cls, "CountingAgentUInt16", "Counting agent.");
    bind_counting_agent<counting_agent<ibf_t, uint32_t>, uint32_t>(cls, "CountingAgentUInt32", "Counting agent.");
    bind_counting_agent<counting_agent<ibf_t, uint64_t>, uint64_t>(cls, "CountingAgentUInt64", "Counting agent.");

    nb::class_<ibf_membership_agent>(cls, "MembershipAgent", "Answers membership queries.")
        .def(
            "membership_for",
            [](ibf_membership_agent & agent, array_like const & values, uint16_t threshold)
            {
                return membership_for(agent.get(), values, threshold);
            },
            "values"_a,
            "threshold"_a,
            "Returns the bins containing at least ``threshold`` of the given values.");
}

void bind_hibf(nb::module_ & m)
{
    auto cls = nb::class_<hibf_t>(m,
                                  "HierarchicalInterleavedBloomFilter",
                                  "The Hierarchical Interleaved Bloom Filter (HIBF): an index for approximate "
                                  "membership queries over many user bins of vastly different sizes.");

    cls.def(
           "__init__",
           [](hibf_t * self, py_config const & config)
           {
               new (self) hibf_t{build_hibf(config, nullptr)};
           },
           "config"_a,
           "Computes a layout and builds the HIBF from a config.")
        .def(
            "__init__",
            [](hibf_t * self, py_config const & config, layout_t const & layout)
            {
                new (self) hibf_t{build_hibf(config, &layout)};
            },
            "config"_a,
            "layout"_a,
            "Builds the HIBF from a config and a precomputed layout.")
        .def_ro("number_of_user_bins", &hibf_t::number_of_user_bins, "The number of user bins.")
        .def_prop_ro(
            "number_of_ibfs",
            [](hibf_t const & self)
            {
                return self.ibf_vector.size();
            },
            "The number of IBFs in the hierarchy.")
        .def_prop_ro(
            "bit_size",
            [](hibf_t const & self)
            {
                size_t bits{};
                for (auto const & ibf : self.ibf_vector)
                    bits += ibf.bit_size();
                return bits;
            },
            "The total size of all IBFs in bits.")
        .def_prop_ro(
            "timings",
            [](hibf_t const & self)
            {
                nb::dict timings;
                timings["layout_compute_sketches"] = self.layout_compute_sketches_timer.in_seconds();
                timings["layout_union_estimation"] = self.layout_union_estimation_timer.in_seconds();
                timings["layout_rearrangement"] = self.layout_rearrangement_timer.in_seconds();
                timings["layout_dp_algorithm"] = self.layout_dp_algorithm_timer.in_seconds();
                timings["index_allocation"] = self.index_allocation_timer.in_seconds();
                timings["user_bin_io"] = self.user_bin_io_timer.in_seconds();
                timings["merge_kmers"] = self.merge_kmers_timer.in_seconds();
                timings["fill_ibf"] = self.fill_ibf_timer.in_seconds();
                return timings;
            },
            "Time spent in each construction step in seconds. Not preserved by serialisation.")
        .def(
            "membership_agent",
            [](hibf_t const & self)
            {
                return hibf_membership_agent{self};
            },
            nb::keep_alive<0, 1>(),
            "Returns an agent that determines the user bins containing at least ``threshold`` of the given values.")
        .def(
            "counting_agent",
            [](hibf_t const & self, dtype_like const & dtype)
            {
                return make_counting_agent(self, dtype);
            },
            "dtype"_a = "uint16",
            nb::keep_alive<0, 1>(),
            "Returns an agent that counts the occurrences of values per user bin. ``dtype`` is the counter type: "
            "uint16, uint32, or uint64.")
        .def(
            "membership_for",
            [](hibf_t const & self, array_like const & values, uint16_t threshold)
            {
                auto agent = self.membership_agent();
                return membership_for(agent, values, threshold);
            },
            "values"_a,
            "threshold"_a,
            "Returns the user bins containing at least ``threshold`` of the given values. The result is unsorted.\n\n"
            "Convenience for ``membership_agent().membership_for(values, threshold)``.")
        .def("batch_membership_for",
             &batch_membership_for<hibf_t>,
             "queries"_a,
             "threshold"_a,
             "threads"_a = 1u,
             "Answers ``membership_for`` for many queries in parallel.\n\n"
             "``threshold`` is either a single value or one value per query. Returns a list with one array of user "
             "bin ids per query.")
        .def(nb::self == nb::self)
        .def("__repr__",
             [](hibf_t const & self)
             {
                 return "<HierarchicalInterleavedBloomFilter number_of_user_bins="
                      + std::to_string(self.number_of_user_bins)
                      + " number_of_ibfs=" + std::to_string(self.ibf_vector.size()) + ">";
             });

    add_serialisation(cls);

    nb::class_<hibf_membership_agent>(cls, "MembershipAgent", "Answers membership queries.")
        .def(
            "membership_for",
            [](hibf_membership_agent & agent, array_like const & values, uint16_t threshold)
            {
                return membership_for(agent.get(), values, threshold);
            },
            "values"_a,
            "threshold"_a,
            "Returns the user bins containing at least ``threshold`` of the given values. The result is unsorted.");

    bind_counting_agent<counting_agent<hibf_t, uint16_t>, uint16_t, size_t>(cls,
                                                                            "CountingAgentUInt16",
                                                                            "Counting agent.");
    bind_counting_agent<counting_agent<hibf_t, uint32_t>, uint32_t, size_t>(cls,
                                                                            "CountingAgentUInt32",
                                                                            "Counting agent.");
    bind_counting_agent<counting_agent<hibf_t, uint64_t>, uint64_t, size_t>(cls,
                                                                            "CountingAgentUInt64",
                                                                            "Counting agent.");
}

//!\brief The library only asserts that merged sketches have the same size; otherwise, it reads out of bounds.
void check_same_size(hyperloglog_t const & self, hyperloglog_t const & other)
{
    if (self.data_size() != other.data_size())
        throw std::invalid_argument{"Cannot merge a sketch with " + std::to_string(other.data_size())
                                    + " registers into a sketch with " + std::to_string(self.data_size())
                                    + " registers. Both sketches must have the same num_bits."};
}

void bind_hyperloglog(nb::module_ & m)
{
    auto cls = nb::class_<hyperloglog_t>(m,
                                         "HyperLogLog",
                                         "A HyperLogLog sketch for estimating the number of distinct values.\n\n"
                                         "Methods that modify the sketch are not synchronised. Do not call them "
                                         "while another thread uses the sketch.");

    cls.def(
           "__init__",
           [](hyperloglog_t * self, uint8_t const num_bits)
           {
               // The library allocates 2^num_bits bytes before checking num_bits.
               if (num_bits < 5u || num_bits > 32u)
                   throw std::invalid_argument{"num_bits must be in [5, 32], got " + std::to_string(num_bits) + "."};
               new (self) hyperloglog_t{num_bits};
           },
           "num_bits"_a = 5u,
           "Creates an empty sketch with 2^num_bits registers. ``num_bits`` must be in [5, 32].")
        .def(
            "add",
            [](hyperloglog_t & self, array_like const & values)
            {
                u64_array const array = as_u64_array(values);
                nb::gil_scoped_release nogil{};
                for (uint64_t const value : as_span(array))
                    self.add(value);
            },
            "values"_a,
            "Adds a value or an array of values. Values should be hashed.")
        .def("estimate", &hyperloglog_t::estimate, "Estimates the number of distinct values.")
        .def(
            "merge",
            [](hyperloglog_t & self, hyperloglog_t const & other)
            {
                check_same_size(self, other);
                self.merge(other);
            },
            "other"_a,
            "Merges another sketch with the same ``num_bits`` into this.")
        .def(
            "merge_and_estimate",
            [](hyperloglog_t & self, hyperloglog_t const & other)
            {
                check_same_size(self, other);
                return self.merge_and_estimate(other);
            },
            "other"_a,
            "Merges another sketch with the same ``num_bits`` into this and returns the new estimate.")
        .def("reset", &hyperloglog_t::reset, "Removes all values.")
        .def_prop_ro("data_size", &hyperloglog_t::data_size, "The number of registers.");

    add_serialisation(cls);
}

} // namespace

NB_MODULE(_hibf, m)
{
    m.doc() = "Python bindings for the Hierarchical Interleaved Bloom Filter (HIBF) library.";
    m.attr("library_version") = seqan::hibf::hibf_version_cstring;

    bind_config(m);
    bind_layout(m);
    bind_ibf(m);
    bind_hibf(m);
    bind_hyperloglog(m);
}
