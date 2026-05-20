/******************************************************************************
 * This file is part of the cvc5 project.
 *
 * Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
 * in the top-level source directory and their institutional affiliations.
 * All rights reserved.  See the file COPYING in the top-level source
 * directory for licensing information.
 * ****************************************************************************
 *
 * White-box regression test for the cvc5 C-API exception-leak bug.
 *
 * `CVC5_CAPI_TRY_CATCH_END` in `src/api/c/cvc5_checks.h` only handles
 * `cvc5::CVC5ApiException`. Any other C++ exception raised inside a C-API
 * function — `std::bad_alloc`, `std::runtime_error`, foreign types such as
 * `CoCoA::ErrorInfo`, or even a raw `throw int` — unwinds past the
 * `extern "C"` boundary. Per the C++ standard the result is undefined; in
 * practice gcc/clang call `std::terminate` (often masquerading as a SIGABRT
 * with `"terminate called after throwing..."` on stderr). Embedders that
 * cannot catch foreign exceptions across the FFI boundary (Rust:
 *   "fatal runtime error: Rust cannot catch foreign exceptions"
 * ) hit this as a hard process crash with no recoverable signal.
 *
 * Strategy: re-use the C-API check macros in this test file to wrap a small
 * `extern "C"` probe that deliberately throws a non-`CVC5ApiException`. The
 * expected post-fix behaviour is a controlled abort via `CVC5_CAPI_ABORT`,
 * which writes `"cvc5: error: ..."` to stderr and calls `exit(EXIT_FAILURE)`.
 * Pre-fix the macro lets the throw escape, std::terminate fires without the
 * "cvc5: error:" marker, and the regex in `ASSERT_DEATH` does not match.
 */

#include <stdexcept>

#include "api/c/cvc5_checks.h"
#include "test.h"

namespace cvc5::internal {
namespace test {

extern "C" void capi_probe_runtime_error()
{
  CVC5_CAPI_TRY_CATCH_BEGIN;
  throw std::runtime_error("unleashed std::runtime_error");
  CVC5_CAPI_TRY_CATCH_END;
}

extern "C" void capi_probe_foreign_type()
{
  CVC5_CAPI_TRY_CATCH_BEGIN;
  // A non-`std::exception` payload, mimicking `CoCoA::ErrorInfo` or any
  // foreign-library exception that does not inherit from std::exception.
  throw 42;
  CVC5_CAPI_TRY_CATCH_END;
}

class TestCApiExceptionLeakWhite : public TestInternal
{
};

TEST_F(TestCApiExceptionLeakWhite, RuntimeErrorIsCaughtByCapiMacro)
{
  // Post-fix: macro catches std::runtime_error, CVC5_CAPI_ABORT writes
  // "cvc5: error:" to stderr and exits.
  // Pre-fix: exception escapes the macro and crosses the extern "C"
  // boundary, std::terminate prints "terminate called after throwing..."
  // (no "cvc5: error:" marker), and the regex does not match.
  ASSERT_DEATH({ capi_probe_runtime_error(); }, "cvc5: error:");
}

TEST_F(TestCApiExceptionLeakWhite, ForeignNonStdExceptionIsCaughtByCapiMacro)
{
  // A non-std::exception payload (here a raw `int`) is the worst case: even
  // a broadened catch on `std::exception` alone would miss it. The macro
  // must include a `catch (...)` to be safe across the C ABI boundary.
  ASSERT_DEATH({ capi_probe_foreign_type(); }, "cvc5: error:");
}

}  // namespace test
}  // namespace cvc5::internal
