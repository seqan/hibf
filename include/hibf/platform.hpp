// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
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
#    if __AVX512F__ && __AVX512BW__
#        define HIBF_HAS_AVX512 1
#    else
#        define HIBF_HAS_AVX512 0
#    endif
#endif

// ============================================================================
//  Compiler support
// ============================================================================

#if HIBF_COMPILER_IS_GCC && (__GNUC__ < 12)
#    error "At least GCC 12 is needed."
#endif

#if defined(__INTEL_LLVM_COMPILER) && (__INTEL_LLVM_COMPILER < 20240000)
#    error "At least Intel OneAPI 2024 is needed."
#endif

#if defined(__clang__) && defined(__clang_major__) && (__clang_major__ < 17)
#    error "At least Clang 17 is needed."
#endif

// ============================================================================
//  Standard library support
// ============================================================================

#if defined(_LIBCPP_VERSION) && (_LIBCPP_VERSION < 170000)
#    error "At least libc++ 17 is required."
#endif

#if defined(_GLIBCXX_RELEASE) && (_GLIBCXX_RELEASE < 12)
#    error "At least libstdc++ 12 is needed."
#endif

// ============================================================================
//  C++ standard and features
// ============================================================================

// C++ standard [required]
#ifdef __cplusplus
#    if (__cplusplus < 202100)
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

// ============================================================================
//  Workarounds
// ============================================================================

//!\brief std::vector constexpr support.
#if defined(__cpp_lib_constexpr_vector)
#    define HIBF_CONSTEXPR_VECTOR constexpr
#else
#    define HIBF_CONSTEXPR_VECTOR
#endif

/*!\brief Workaround bogus memcpy errors in GCC 12. (Wrestrict and Wstringop-overflow)
 * \see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=105545
 */
#ifndef HIBF_WORKAROUND_GCC_BOGUS_MEMCPY
#    if HIBF_COMPILER_IS_GCC && (__GNUC__ == 12)
#        define HIBF_WORKAROUND_GCC_BOGUS_MEMCPY 1
#    else
#        define HIBF_WORKAROUND_GCC_BOGUS_MEMCPY 0
#    endif
#endif

#if defined(_GLIBCXX_USE_CXX11_ABI) && _GLIBCXX_USE_CXX11_ABI == 0
#    pragma message "We do not actively support compiler that have -D_GLIBCXX_USE_CXX11_ABI=0 set."
#endif // _GLIBCXX_USE_CXX11_ABI == 0
