// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/library-template/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides platform and dependency checks.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cinttypes>
#include <ciso646> // makes _LIBCPP_VERSION available
#include <cstddef> // makes __GLIBCXX__ available

// macro cruft
//!\cond
#define LIBRARY_TEMPLATE_STR_HELPER(x) #x
#define LIBRARY_TEMPLATE_STR(x) LIBRARY_TEMPLATE_STR_HELPER(x)
//!\endcond

// ============================================================================
//  Documentation
// ============================================================================

// Doxygen related
// this macro is a NO-OP unless doxygen parses it, in which case it resolves to the argument
#ifndef LIBRARY_TEMPLATE_DOXYGEN_ONLY
#    define LIBRARY_TEMPLATE_DOXYGEN_ONLY(x)
#endif

// ============================================================================
//  Compiler support general
// ============================================================================

/*!\def LIBRARY_TEMPLATE_COMPILER_IS_GCC
 * \brief Whether the current compiler is GCC.
 * \ingroup core
 * \details
 * __GNUC__ is also used to indicate the support for GNU compiler extensions. To detect the presence of the GCC
 * compiler, one has to rule out other compilers.
 *
 * \sa https://sourceforge.net/p/predef/wiki/Compilers
 */
#if defined(__GNUC__) && !defined(__llvm__) && !defined(__INTEL_COMPILER)
#    define LIBRARY_TEMPLATE_COMPILER_IS_GCC 1
#else
#    define LIBRARY_TEMPLATE_COMPILER_IS_GCC 0
#endif

#if LIBRARY_TEMPLATE_DOXYGEN_ONLY(1) 0
//!\brief This disables the warning you would get if your compiler is not known to work.
#    define LIBRARY_TEMPLATE_DISABLE_COMPILER_CHECK
#endif // LIBRARY_TEMPLATE_DOXYGEN_ONLY(1)0

// ============================================================================
//  Compiler support GCC
// ============================================================================

#if LIBRARY_TEMPLATE_COMPILER_IS_GCC
#    if (__GNUC__ < 10)
#        error "At least GCC 10 is needed."
#    endif // (__GNUC__ < 10)

#    if (__GNUC__ == 10 && __GNUC_MINOR__ <= 3)
#        pragma GCC warning "Be aware that GCC < 10.4 might have bugs that cause compile failure."
#    endif // (__GNUC__ == 10 && __GNUC_MINOR__ <= 3)

#    if (__GNUC__ == 11 && __GNUC_MINOR__ <= 2)
#        pragma GCC warning "Be aware that GCC < 11.3 might have bugs that cause compile failure."
#    endif // (__GNUC__ == 11 && __GNUC_MINOR__ <= 2)

#    if (__GNUC__ == 12 && __GNUC_MINOR__ <= 1)
#        pragma GCC warning "Be aware that GCC < 12.2 might have bugs that cause compile failure."
#    endif // (__GNUC__ == 12 && __GNUC_MINOR__ <= 1)

#    if LIBRARY_TEMPLATE_DOXYGEN_ONLY(1) 0
//!\brief This disables the warning you would get if your compiler is newer than the latest supported version.
#        define LIBRARY_TEMPLATE_DISABLE_NEWER_COMPILER_DIAGNOSTIC
#    endif // LIBRARY_TEMPLATE_DOXYGEN_ONLY(1)0

#    ifndef LIBRARY_TEMPLATE_DISABLE_NEWER_COMPILER_DIAGNOSTIC
#        if (__GNUC__ > 12)
#            pragma message                                                                                            \
                "Your compiler is newer than the latest supported compiler version (gcc-12). It might be that compiling fails. You can disable this warning by setting -DLIBRARY_TEMPLATE_DISABLE_NEWER_COMPILER_DIAGNOSTIC."
#        endif // (__GNUC__ > 12)
#    endif     // LIBRARY_TEMPLATE_DISABLE_NEWER_COMPILER_DIAGNOSTIC

// ============================================================================
//  Compiler support other
// ============================================================================

#elif !defined(LIBRARY_TEMPLATE_DISABLE_COMPILER_CHECK)
#    error                                                                                                             \
        "Your compiler is not supported. Currently, only GCC is known to work. You can disable this error by setting -DLIBRARY_TEMPLATE_DISABLE_COMPILER_CHECK."
#endif // LIBRARY_TEMPLATE_COMPILER_IS_GCC

// ============================================================================
//  C++ standard and features
// ============================================================================

#if __has_include(<version>)
#    include <version>
#endif

// C++ standard [required]
// Note: gcc10 -std=c++20 still defines __cplusplus=201709
#ifdef __cplusplus
#    if (__cplusplus < 201709)
#        error "C++20 is required, make sure that you have set -std=c++20."
#    endif
#else
#    error "This is not a C++ compiler."
#endif

// ============================================================================
//  Dependencies
// ============================================================================

// LIBRARY_TEMPLATE [required]
#if __has_include(<library_template/version.hpp>)
#    include <library_template/version.hpp>
#else
#    error                                                                                                             \
        "LIBRARY_TEMPLATE include directory not set correctly. Forgot to add -I ${INSTALLDIR}/include to your CXXFLAGS?"
#endif

// ============================================================================
//  Deprecation Messages
// ============================================================================

//!\brief _Pragma requires a string-literal and # makes it a string
#ifndef LIBRARY_TEMPLATE_PRAGMA
#    define LIBRARY_TEMPLATE_PRAGMA(non_string_literal) _Pragma(#non_string_literal)
#endif

//!\brief Deprecation message for deprecated header.
#ifndef LIBRARY_TEMPLATE_DEPRECATED_HEADER
#    ifndef LIBRARY_TEMPLATE_DISABLE_DEPRECATED_WARNINGS
#        define LIBRARY_TEMPLATE_DEPRECATED_HEADER(message) LIBRARY_TEMPLATE_PRAGMA(GCC warning message)
#    else
#        define LIBRARY_TEMPLATE_DEPRECATED_HEADER(message) /**/
#    endif
#endif

//!\brief Deprecation message for release.
#ifndef LIBRARY_TEMPLATE_REMOVE_DEPRECATED_100
#    ifndef LIBRARY_TEMPLATE_DEPRECATED_100
#        ifndef LIBRARY_TEMPLATE_DISABLE_DEPRECATED_WARNINGS
#            define LIBRARY_TEMPLATE_DEPRECATED_100                                                                    \
                [[deprecated("This will be removed in version 1.0.0; please see the documentation.")]]
#        else
#            define LIBRARY_TEMPLATE_DEPRECATED_100 /**/
#        endif
#    endif
#endif

// ============================================================================
//  Workarounds
// ============================================================================

/*!\brief Workaround bogus memcpy errors in GCC 12.1 and 12.2. (Wrestrict and Wstringop-overflow)
 * \see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=105545
 */
#ifndef LIBRARY_TEMPLATE_WORKAROUND_GCC_BOGUS_MEMCPY
#    if LIBRARY_TEMPLATE_COMPILER_IS_GCC && (__GNUC__ == 12 && __GNUC_MINOR__ < 3)
#        define LIBRARY_TEMPLATE_WORKAROUND_GCC_BOGUS_MEMCPY 1
#    else
#        define LIBRARY_TEMPLATE_WORKAROUND_GCC_BOGUS_MEMCPY 0
#    endif
#endif

/*!\brief This is needed to support CentOS 7 or RHEL 7; Newer CentOS's include a more modern default-gcc version making
 *        this macro obsolete.
 *
 * In GCC 5 there was a bigger ABI change and modern systems compile with dual ABI, but some enterprise systems (those
 * where gcc 4 is the standard compiler) don't support dual ABI. This has the effect that even community builds of gcc
 * are build with --disable-libstdcxx-dual-abi. Only building the compiler yourself would solve this problem.
 *
 * \see https://github.com/seqan/library-template/issues/2244
 * \see https://gcc.gnu.org/onlinedocs/libstdc++/manual/using_dual_abi.html
 */
#ifndef LIBRARY_TEMPLATE_WORKAROUND_GCC_NO_CXX11_ABI
#    if defined(_GLIBCXX_USE_CXX11_ABI) && _GLIBCXX_USE_CXX11_ABI == 0
#        define LIBRARY_TEMPLATE_WORKAROUND_GCC_NO_CXX11_ABI 1
#    else
#        define LIBRARY_TEMPLATE_WORKAROUND_GCC_NO_CXX11_ABI 0
#    endif
#endif

//!\brief Our char literals returning std::vector should be constexpr if constexpr std::vector is supported.
#if defined(__cpp_lib_constexpr_vector) && __cpp_lib_constexpr_vector >= 201907L
#    define LIBRARY_TEMPLATE_WORKAROUND_LITERAL constexpr
#else
#    define LIBRARY_TEMPLATE_WORKAROUND_LITERAL inline
#endif

#if LIBRARY_TEMPLATE_DOXYGEN_ONLY(1) 0
//!\brief This disables the warning you would get if -D_GLIBCXX_USE_CXX11_ABI=0 is set.
#    define LIBRARY_TEMPLATE_DISABLE_LEGACY_STD_DIAGNOSTIC
#endif // LIBRARY_TEMPLATE_DOXYGEN_ONLY(1)0

#if defined(_GLIBCXX_USE_CXX11_ABI) && _GLIBCXX_USE_CXX11_ABI == 0
#    ifndef LIBRARY_TEMPLATE_DISABLE_LEGACY_STD_DIAGNOSTIC
#        pragma message                                                                                                \
            "We do not actively support compiler that have -D_GLIBCXX_USE_CXX11_ABI=0 set, and it might be that compiling fails. It is known that all compilers of CentOS 7 / RHEL 7 set this flag by default (and that it cannot be overridden!). Note that these versions of the OSes are community-supported. You can disable this warning by setting -DLIBRARY_TEMPLATE_DISABLE_LEGACY_STD_DIAGNOSTIC."
#    endif // LIBRARY_TEMPLATE_DISABLE_LEGACY_STD_DIAGNOSTIC
#endif     // _GLIBCXX_USE_CXX11_ABI == 0

// ============================================================================
//  Backmatter
// ============================================================================

// macro cruft undefine
#undef LIBRARY_TEMPLATE_STR
#undef LIBRARY_TEMPLATE_STR_HELPER
