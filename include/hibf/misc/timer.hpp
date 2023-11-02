// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan::hibf::timer.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <atomic>      // for atomic
#include <cassert>     // for assert
#include <chrono>      // for steady_clock, duration, operator-, time_point
#include <type_traits> // for conditional_t

#include <hibf/platform.hpp>

namespace seqan::hibf
{

/*!\brief Whether the timer is concurrent.
 * \ingroup hibf_build
 */
enum class concurrent
{
    no, //!< Not concurrent.
    yes //!< Concurrent.
};

/*!\brief Timer.
 * \ingroup hibf_build
 */
template <concurrent concurrency>
class timer
{
private:
    static constexpr bool is_concurrent{concurrency == concurrent::yes};

    using rep_t =
        std::conditional_t<is_concurrent, std::atomic<std::chrono::steady_clock::rep>, std::chrono::steady_clock::rep>;

    using count_t = std::conditional_t<is_concurrent, std::atomic<uint64_t>, uint64_t>;

    static constexpr int alignment = is_concurrent ? 64 : 8;

    template <concurrent concurrency_>
    friend class timer;

public:
    timer() = default;
    timer(timer &&) = default;
    timer & operator=(timer &&) = default;
    ~timer() = default;

    timer(timer const & other)
        requires (!is_concurrent)
    = default;
    timer & operator=(timer const & other)
        requires (!is_concurrent)
    = default;

    timer(timer const & other)
        requires is_concurrent
        :
        start_{other.start_},
        stop_{other.stop_},
        ticks{other.ticks.load()},
        max{other.max.load()},
        count{other.count.load()}
    {}
    timer & operator=(timer const & other)
        requires is_concurrent
    {
        start_ = other.start_;
        stop_ = other.stop_;
        ticks = other.ticks.load();
        max = other.max.load();
        count = other.count.load();
        return *this;
    }

    void start()
    {
        start_ = std::chrono::steady_clock::now();
    }

    void stop()
    {
        stop_ = std::chrono::steady_clock::now();
        assert(stop_ >= start_);
        std::chrono::steady_clock::rep duration = (stop_ - start_).count();
        ticks += duration;
        update_max(duration);
        ++count;
    }

    template <concurrent concurrency_>
    void operator+=(timer<concurrency_> const & other)
    {
        ticks += other.ticks;
        update_max(other.ticks);
        ++count;
    }

    double in_seconds() const
        requires is_concurrent
    {
        return std::chrono::duration<double>(std::chrono::steady_clock::duration{ticks.load()}).count();
    }

    double max_in_seconds() const
        requires is_concurrent
    {
        return std::chrono::duration<double>(std::chrono::steady_clock::duration{max.load()}).count();
    }

    double avg_in_seconds() const
        requires is_concurrent
    {
        assert(count.load() > 0u);
        return in_seconds() / count.load();
    }

    // GCOVR_EXCL_START
    double in_seconds() const
        requires (!is_concurrent)
    {
        return std::chrono::duration<double>(std::chrono::steady_clock::duration{ticks}).count();
    }

    double max_in_seconds() const
        requires (!is_concurrent)
    {
        return std::chrono::duration<double>(std::chrono::steady_clock::duration{max}).count();
    }

    double avg_in_seconds() const
        requires (!is_concurrent)
    {
        assert(count > 0u);
        return in_seconds() / count;
    }
    // GCOVR_EXCL_STOP

    // Timer are always equal.
    constexpr bool operator==(timer const &) const
    {
        return true;
    }

private:
    std::chrono::steady_clock::time_point start_{std::chrono::time_point<std::chrono::steady_clock>::max()};
    std::chrono::steady_clock::time_point stop_{};

    alignas(alignment) rep_t ticks{};
    alignas(alignment) rep_t max{};
    alignas(alignment) count_t count{};

    void update_max(std::chrono::steady_clock::rep const value)
        requires is_concurrent
    {
        for (std::chrono::steady_clock::rep previous_value = max;
             previous_value < value && !max.compare_exchange_weak(previous_value, value, std::memory_order_relaxed);)
            ;
    }

    // GCOVR_EXCL_START
    void update_max(std::chrono::steady_clock::rep const value)
        requires (!is_concurrent)
    {
        max = std::max(max, value);
    }
    // GCOVR_EXCL_STOP
};

/*!\brief Alias for timer<concurrent::no>
 * \ingroup hibf_build
 */
using serial_timer = timer<concurrent::no>;
/*!\brief Alias for timer<concurrent::yes>
 * \ingroup hibf_build
 */
using concurrent_timer = timer<concurrent::yes>;

} // namespace seqan::hibf
