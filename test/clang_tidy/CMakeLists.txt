cmake_minimum_required (VERSION 3.10)

project (hibf_tidy NONE)

find_program (CLANG_TIDY "clang-tidy-17" REQUIRED)
set (CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY})

add_subdirectory (../unit/ unit)
