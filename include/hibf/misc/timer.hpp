// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan::hibf::timer.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <algorithm> // for max
#include <atomic>    // for atomic, __atomic_base, memory_order_relaxed
#include <cassert>   // for assert
#include <chrono>    // for duration, time_point, operator-, steady_clock
#include <cinttypes> // for uint64_t
#include <concepts>  // for same_as
#include <utility>   // for move

#include <hibf/platform.hpp>

namespace seqan::hibf
{

class concurrent_timer;

/*!\brief A timer.
 * \ingroup hibf
 */
class serial_timer
{
private:
    using steady_clock_t = std::chrono::steady_clock;

    steady_clock_t::time_point start_point{std::chrono::time_point<steady_clock_t>::max()};
    steady_clock_t::time_point stop_point{};

    steady_clock_t::rep ticks{};
    steady_clock_t::rep max{};
    uint64_t count{};

    friend class concurrent_timer;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    serial_timer() = default;                                 //!< Defaulted.
    serial_timer(serial_timer const &) = default;             //!< Defaulted.
    serial_timer & operator=(serial_timer const &) = default; //!< Defaulted.
    serial_timer(serial_timer &&) = default;                  //!< Defaulted.
    serial_timer & operator=(serial_timer &&) = default;      //!< Defaulted.
    ~serial_timer() = default;                                //!< Defaulted.

    //!\}

    /*!\name Modification
     * \{
     */
    //!\brief Starts the timer.
    void start()
    {
        start_point = steady_clock_t::now();
    }

    /*!\brief Stops the timer.
     * \details
     * In Debug mode, an assertion checks that `start()` has been called before.
     */
    void stop()
    {
        stop_point = steady_clock_t::now();
        assert(stop_point >= start_point);
        steady_clock_t::rep duration = (stop_point - start_point).count();

        ticks += duration;
        max = std::max(max, duration);
        ++count;
    }

    //!\brief Adds another timer.
    template <typename timer_t>
        requires (std::same_as<timer_t, serial_timer> || std::same_as<timer_t, concurrent_timer>)
    void operator+=(timer_t const & other)
    {
        steady_clock_t::rep ticks_to_add{};
        if constexpr (std::same_as<timer_t, concurrent_timer>)
            ticks_to_add = other.ticks.load();
        else
            ticks_to_add = other.ticks;

        ticks += ticks_to_add;
        max = std::max(max, ticks_to_add);
        ++count;
    }
    //!\}

    /*!\name Access
     * \{
     */
    //!\brief Returns the measured time in seconds.
    double in_seconds() const
    {
        return std::chrono::duration<double>(steady_clock_t::duration{ticks}).count();
    }

    /*!\brief Returns the maximum measured time interval in seconds.
     * \details
     * A time interval may be:
     * * The elpased time between `start()` and `stop()`.
     * * The elapsed time added via `operator+=()`.
     */
    double max_in_seconds() const
    {
        return std::chrono::duration<double>(steady_clock_t::duration{max}).count();
    }

    /*!\brief Returns the average measured time interval in seconds.
     * \details
     * A time interval may be:
     * * The elpased time between `start()` and `stop()`.
     * * The elapsed time added via `operator+=()`.
     *
     * The count used for averaging is the number of calls to `stop()` and `operator+=()`.
     * \warning Calling this function when neither `stop()` or `operator+=()` have been used is undefined behaviour.
     */
    double avg_in_seconds() const
    {
        assert(count > 0u);
        return in_seconds() / count;
    }
    //!\}

    /*!\name Comparison
     * \{
     */
    //!\brief Two timer are always equal.
    constexpr bool operator==(serial_timer const &) const
    {
        return true;
    }

    //!\brief Two timer are always equal.
    constexpr bool operator==(concurrent_timer const &) const
    {
        return true;
    }
    //!\}
};

/*!\brief A timer with a thread-safe `operator+=()`.
 * \ingroup hibf
 */
class concurrent_timer
{
private:
    using steady_clock_t = std::chrono::steady_clock;

    alignas(64) std::atomic<steady_clock_t::rep> ticks{};

    steady_clock_t::time_point start_point{std::chrono::time_point<steady_clock_t>::max()};
    steady_clock_t::time_point stop_point{};

    alignas(64) std::atomic<steady_clock_t::rep> max{};
    alignas(64) std::atomic<uint64_t> count{};

    friend class serial_timer;

    void update_max(steady_clock_t::rep const value)
    {
        for (steady_clock_t::rep previous_value = max;
             previous_value < value && !max.compare_exchange_weak(previous_value, value, std::memory_order_relaxed);)
            ;
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Defaulted.
    concurrent_timer() = default;
    //!\brief Defaulted.
    concurrent_timer(concurrent_timer const & other) :
        ticks{other.ticks.load()},
        start_point{other.start_point},
        stop_point{other.stop_point},
        max{other.max.load()},
        count{other.count.load()}
    {}
    //!\brief Defaulted.
    concurrent_timer & operator=(concurrent_timer const & other)
    {
        ticks = other.ticks.load();
        start_point = other.start_point;
        stop_point = other.stop_point;
        max = other.max.load();
        count = other.count.load();
        return *this;
    }
    //!\brief Defaulted.
    concurrent_timer(concurrent_timer && other) noexcept :
        ticks{other.ticks.load()},
        start_point{other.start_point},
        stop_point{other.stop_point},
        max{other.max.load()},
        count{other.count.load()}
    {}
    //!\brief Defaulted.
    concurrent_timer & operator=(concurrent_timer && other) noexcept
    {
        ticks = other.ticks.load();
        start_point = other.start_point;
        stop_point = other.stop_point;
        max = other.max.load();
        count = other.count.load();
        return *this;
    }
    //!\brief Defaulted.
    ~concurrent_timer() = default;
    //!\}

    /*!\name Modification
     * \{
     */
    //!\copydoc seqan::hibf::serial_timer::start
    void start()
    {
        start_point = steady_clock_t::now();
    }

    /*!\copydoc seqan::hibf::serial_timer::stop
     * \attention This function is **not** thread-safe.
     */
    void stop()
    {
        stop_point = steady_clock_t::now();
        assert(stop_point >= start_point);
        steady_clock_t::rep duration = (stop_point - start_point).count();

        ticks += duration;
        update_max(duration);
        ++count;
    }

    /*!\copydoc seqan::hibf::serial_timer::operator+=
     * This function is thread-safe.
     */
    template <typename timer_t>
        requires (std::same_as<timer_t, serial_timer> || std::same_as<timer_t, concurrent_timer>)
    void operator+=(timer_t const & other)
    {
        ticks += other.ticks;
        update_max(other.ticks);
        ++count;
    }
    //!\}

    /*!\name Access
     * \{
     */
    //!\copydoc seqan::hibf::serial_timer::in_seconds
    double in_seconds() const
    {
        return std::chrono::duration<double>(steady_clock_t::duration{ticks.load()}).count();
    }

    //!\copydoc seqan::hibf::serial_timer::max_in_seconds
    double max_in_seconds() const
    {
        return std::chrono::duration<double>(steady_clock_t::duration{max.load()}).count();
    }

    //!\copydoc seqan::hibf::serial_timer::avg_in_seconds
    double avg_in_seconds() const
    {
        assert(count.load() > 0u);
        return in_seconds() / count.load();
    }

    /*!\name Comparison
     * \{
     */
    //!\copydoc seqan::hibf::serial_timer::operator==
    constexpr bool operator==(serial_timer const &) const
    {
        return true;
    }

    //!\copydoc seqan::hibf::serial_timer::operator==
    constexpr bool operator==(concurrent_timer const &) const
    {
        return true;
    }
    //!\}
};

} // namespace seqan::hibf
