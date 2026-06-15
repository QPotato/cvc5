/******************************************************************************
 * This file is part of the cvc5 project.
 *
 * Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
 * in the top-level source directory and their institutional affiliations.
 * All rights reserved.  See the file COPYING in the top-level source
 * directory for licensing information.
 * ****************************************************************************
 *
 * Singleton libSingular global initializer.
 *
 * libSingular must be initialized (via siInit) before most Singular kernel
 * operations are performed.
 *
 * It must be initialized **exactly once** per process. We guard this with a
 * static boolean so that repeated calls are no-ops.
 *
 * Before initializing, we redirect Singular's output callbacks to a function
 * that swallows the output, so that the library never writes to cvc5's stdout
 * or stderr.
 */

#include "cvc5_public.h"

#ifdef CVC5_USE_SINGULAR

#ifndef CVC5__UTIL__SINGULAR_GLOBALS_H
#define CVC5__UTIL__SINGULAR_GLOBALS_H

namespace cvc5::internal {

/**
 * Initializes libSingular if it has not been initialized already.
 *
 * On the first call, this redirects Singular's print/warn/error callbacks to a
 * swallow function and then calls siInit. Subsequent calls are no-ops.
 */
void initSingular();

}  // namespace cvc5::internal

#endif /* CVC5__UTIL__SINGULAR_GLOBALS_H */

#endif /* CVC5_USE_SINGULAR */
