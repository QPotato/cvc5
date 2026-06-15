###############################################################################
# This file is part of the cvc5 project.
#
# Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
# in the top-level source directory and their institutional affiliations.
# All rights reserved.  See the file COPYING in the top-level source
# directory for licensing information.
# #############################################################################
#
# Find Singular
# Singular_FOUND - system has Singular lib
# Singular_INCLUDE_DIRS - the Singular include directories
# Singular_LIBRARIES - Libraries needed to use Singular
##

include(deps-helper)

# Singular is only supported as a system dependency via pkg-config; there is no
# auto-download support for it.
find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(Singular Singular)
endif()

if(NOT Singular_FOUND)
  message(FATAL_ERROR
    "Could not find Singular via pkg-config. Please install the system "
    "package (Debian/Ubuntu: 'sudo apt-get install libsingular4-dev'; "
    "Fedora: 'Singular-devel'). Auto-download is not supported for Singular.")
endif()

# The pkg-config '--cflags' carry the mandatory defines -DSING_NDEBUG and
# -DOM_NDEBUG. Split the flags into compile definitions (without the leading
# -D) and include directories.
set(Singular_COMPILE_DEFINITIONS "")
set(Singular_INCLUDE_DIR "")
foreach(flag ${Singular_CFLAGS} ${Singular_CFLAGS_OTHER})
  if(flag MATCHES "^-D(.*)$")
    list(APPEND Singular_COMPILE_DEFINITIONS "${CMAKE_MATCH_1}")
  elseif(flag MATCHES "^-I(.*)$")
    list(APPEND Singular_INCLUDE_DIR "${CMAKE_MATCH_1}")
  endif()
endforeach()
list(APPEND Singular_INCLUDE_DIR ${Singular_INCLUDE_DIRS})
list(REMOVE_DUPLICATES Singular_INCLUDE_DIR)

add_library(Singular INTERFACE IMPORTED GLOBAL)
set_target_properties(Singular PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${Singular_INCLUDE_DIRS}"
)
set_target_properties(Singular PROPERTIES
  INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${Singular_INCLUDE_DIRS}"
)
set_target_properties(Singular PROPERTIES
  INTERFACE_LINK_LIBRARIES "${Singular_LIBRARIES}"
)
# pkg-config reports the library names (e.g. singular-Singular) plus a -L search
# directory separately; propagate the search directory so the linker can find
# them (needed when Singular lives outside the default linker search path).
set_target_properties(Singular PROPERTIES
  INTERFACE_LINK_DIRECTORIES "${Singular_LIBRARY_DIRS}"
)
set_target_properties(Singular PROPERTIES
  INTERFACE_COMPILE_DEFINITIONS "SING_NDEBUG;OM_NDEBUG"
)

mark_as_advanced(Singular_FOUND)
mark_as_advanced(Singular_INCLUDE_DIR)
mark_as_advanced(Singular_INCLUDE_DIRS)
mark_as_advanced(Singular_LIBRARIES)
mark_as_advanced(Singular_COMPILE_DEFINITIONS)

message(STATUS "Found Singular ${Singular_VERSION}: ${Singular_LIBRARIES}")
