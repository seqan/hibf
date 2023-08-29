#pragma once

#include <cstddef> // for size_t

#include <hibf/config.hpp>            // for config
#include <hibf/detail/data_store.hpp> // for data_store

namespace seqan::hibf
{

size_t execute(config const &, data_store &);

} // namespace seqan::hibf
