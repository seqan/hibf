#pragma once

#include <filesystem> // for path
#include <string>     // for string

#include <hibf/platform.hpp>

#include <cereal/macros.hpp> // for CEREAL_LOAD_FUNCTION_NAME, CEREAL_SAVE_FUNCTION_NAME

namespace cereal
{

template <typename archive_t>
void CEREAL_SAVE_FUNCTION_NAME(archive_t & archive, std::filesystem::path const & path)
{
    std::string const str{path.string()};
    archive(str);
}

template <typename archive_t>
void CEREAL_LOAD_FUNCTION_NAME(archive_t & archive, std::filesystem::path & path)
{
    std::string str;
    archive(str);
    path.assign(str);
}

} // namespace cereal
