// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides platform and dependency checks.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

// IWYU pragma: always_keep
// IWYU pragma: begin_exports

#include <version> // for __cpp_lib_constexpr_vector

// IWYU pragma: end_exports

// ============================================================================
//  Documentation
// ============================================================================

// Doxygen related
// this macro is a NO-OP unless doxygen parses it, in which case it resolves to the argument
#ifndef HIBF_DOXYGEN_ONLY
#    define HIBF_DOXYGEN_ONLY(x)
#endif

// ============================================================================
//  Compiler support general
// ============================================================================

/*!\def HIBF_COMPILER_IS_GCC
 * \brief Whether the current compiler is GCC.
 * \private
 * \details
 * __GNUC__ is also used to indicate the support for GNU compiler extensions. To detect the presence of the GCC
 * compiler, one has to rule out other compilers.
 *
 * \sa https://sourceforge.net/p/predef/wiki/Compilers
 */
#if defined(__GNUC__) && !defined(__llvm__) && !defined(__INTEL_COMPILER) && !defined(__INTEL_LLVM_COMPILER)
#    define HIBF_COMPILER_IS_GCC 1
#else
#    define HIBF_COMPILER_IS_GCC 0
#endif

/*!\def HIBF_HAS_AVX512
 * \brief Whether AVX512F and AVX512BW are available.
 * \private
 */
#ifndef HIBF_HAS_AVX512
#    if defined(__AVX512F__) && defined(__AVX512BW__)
#        define HIBF_HAS_AVX512 1
#    else
#        define HIBF_HAS_AVX512 0
#    endif
#endif

// ============================================================================
//  Compiler support
// ============================================================================

#if HIBF_COMPILER_IS_GCC && (__GNUC__ < 14)
#    error "At least GCC 14 is needed."
#endif

#if defined(__INTEL_LLVM_COMPILER) && (__INTEL_LLVM_COMPILER < 20250000)
#    error "At least Intel OneAPI 2025 is needed."
#endif

#if defined(__clang__) && defined(__clang_major__) && (__clang_major__ < 20) && !defined(__INTEL_LLVM_COMPILER)
#    error "At least Clang 20 is needed."
#endif

#if defined(_GLIBCXX_USE_CXX11_ABI) && _GLIBCXX_USE_CXX11_ABI == 0
#    pragma message "We do not actively support compiler that have -D_GLIBCXX_USE_CXX11_ABI=0 set."
#endif // _GLIBCXX_USE_CXX11_ABI == 0

// ============================================================================
//  Standard library support
// ============================================================================

#if defined(_LIBCPP_VERSION) && (_LIBCPP_VERSION < 200000)
#    error "At least libc++ 20 is required."
#endif

#if defined(_GLIBCXX_RELEASE) && (_GLIBCXX_RELEASE < 14)
#    error "At least libstdc++ 14 is needed."
#endif

// ============================================================================
//  C++ standard and features
// ============================================================================

// C++ standard [required]
#ifdef __cplusplus
#    if (__cplusplus < 202302L)
#        error "C++23 is required, make sure that you have set -std=c++23."
#    endif
#else
#    error "This is not a C++ compiler."
#endif

// ============================================================================
//  Dependencies
// ============================================================================

// HIBF [required]
#if __has_include(<hibf/version.hpp>)
#    include <hibf/version.hpp>
#else
#    error "HIBF include directory not set correctly. Forgot to add -I ${INSTALLDIR}/include to your CXXFLAGS?"
#endif
