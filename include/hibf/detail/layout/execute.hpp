#pragma once

#include <cstddef> // for size_t

#include <hibf/detail/configuration.hpp> // for configuration
#include <hibf/detail/data_store.hpp>    // for data_store

namespace hibf
{

size_t execute(configuration &, data_store &);

} // namespace hibf
