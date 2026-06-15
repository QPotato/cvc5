/******************************************************************************
 * This file is part of the cvc5 project.
 *
 * Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
 * in the top-level source directory and their institutional affiliations.
 * All rights reserved.  See the file COPYING in the top-level source
 * directory for licensing information.
 * ****************************************************************************
 *
 * Black box testing of groebner basis core computation.
 *
 * The single test case that used to live here (DistinctRootsPoly) did not
 * actually test root finding: it exercised the CoCoA-based UNSAT-core
 * dependency tracer (theory/ff/core.h's `ff::Tracer`), which hooked into
 * CoCoA's Groebner-basis reduction callbacks.
 *
 * The Singular wrapper exposes no equivalent of that callback-based tracer
 * (libSingular does not provide the reduction hooks CoCoA did), so this
 * functionality has been removed. The test case is therefore dropped rather
 * than ported. See the report accompanying the CoCoA->Singular migration.
 */

#include "cvc5_private.h"

#ifdef CVC5_USE_SINGULAR

// No test cases: the only test here exercised the removed CoCoA Tracer/core.

#endif  // CVC5_USE_SINGULAR
