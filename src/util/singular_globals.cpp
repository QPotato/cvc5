/******************************************************************************
 * This file is part of the cvc5 project.
 *
 * Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
 * in the top-level source directory and their institutional affiliations.
 * All rights reserved.  See the file COPYING in the top-level source
 * directory for licensing information.
 * ****************************************************************************
 *
 * Singleton libSingular global initializer
 */

#include "cvc5_public.h"

#ifdef CVC5_USE_SINGULAR

#include "util/singular_globals.h"

// external includes
#include <Singular/libsingular.h>
#include <reporter/reporter.h>
#include <resources/feFopen.h>

namespace cvc5::internal {

namespace {
/** Swallow Singular output so that it never reaches cvc5's streams. */
void swallow(const char*) {}
}  // namespace

void initSingular()
{
  static bool s_initialized = false;
  if (!s_initialized)
  {
    // Redirect all of Singular's output to a swallow function *before* siInit,
    // so the library never writes to cvc5's stdout/stderr (siInit itself emits
    // a standard.lib warning that we want to suppress).
    WerrorS_callback = swallow;
    WarnS_callback = swallow;
    PrintS_callback = swallow;
    // The argument is unused for kernel-only work.
    siInit((char*)"cvc5");
    s_initialized = true;
  }
}

}  // namespace cvc5::internal

#endif /* CVC5_USE_SINGULAR */
