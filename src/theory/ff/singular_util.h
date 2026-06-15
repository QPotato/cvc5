/******************************************************************************
 * This file is part of the cvc5 project.
 *
 * Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
 * in the top-level source directory and their institutional affiliations.
 * All rights reserved.  See the file COPYING in the top-level source
 * directory for licensing information.
 * ****************************************************************************
 *
 * libSingular utilities.
 *
 * A thin, value-semantics C++ wrapper over libSingular's polynomial kernel,
 * specialized to prime fields GF(p). Provides the primitives the finite-field
 * theory needs: polynomial arithmetic, Groebner bases, normal forms, ideal
 * dimension, minimal polynomials, univariate root finding, gcd, and
 * evaluation.
 *
 * This mirrors the public surface that the CoCoA-based code previously used
 * (see cocoa_util.h, uni_roots.h, multi_roots.h).
 *
 * NB: this header deliberately does NOT include <Singular/libsingular.h>.
 * Singular defines a global `typedef ... poly;` which clashes with libpoly's
 * global `namespace poly` (pulled in transitively by expr/node.h). To let the
 * FF solver files use both, the Singular handle types are kept opaque here
 * (forward-declared struct pointers) and libsingular.h is included only in
 * singular_util.cpp.
 */

#ifdef CVC5_USE_SINGULAR

#include "cvc5_private.h"

#ifndef CVC5__THEORY__FF__SINGULAR_UTIL_H
#define CVC5__THEORY__FF__SINGULAR_UTIL_H

// external includes
#include <gmp.h>

// std includes
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

// internal includes
#include "util/finite_field_value.h"
#include "util/integer.h"

// Forward declarations of Singular's opaque kernel structs (defined in
// <Singular/libsingular.h>, included only in the .cpp). These are the targets
// of Singular's `ring`/`coeffs`/`poly`/`number`/`ideal` typedefs; we use the
// struct-pointer spelling here so the typedef names never enter a translation
// unit that also pulls in libpoly's `namespace poly`.
struct ip_sring;
struct n_Procs_s;
struct spolyrec;
struct snumber;
struct sip_sideal;

namespace cvc5::internal {
namespace theory {
namespace ff {

/**
 * Owns a Singular polynomial ring over GF(p) and its coefficient field.
 *
 * Non-copyable. Calls rDelete on destruction. All Poly/Ideal objects living in
 * this ring must be destroyed *before* the SingularRing is destroyed.
 */
class SingularRing
{
 public:
  /**
   * Construct a polynomial ring over GF(prime) with the given variable names,
   * using a degree-lexicographic (dp) ordering.
   *
   * If prime <= FACTORY_MAX_PRIME a true small prime field (n_Zp) is used;
   * otherwise a custom GMP-backed large prime field (see singular_util.cpp) is
   * used -- Singular has no built-in large prime field, and its integers-mod-n
   * type (n_Zn) is a ring engine whose Groebner bases are incomplete here.
   */
  SingularRing(const Integer& prime, const std::vector<std::string>& varNames);
  ~SingularRing();

  SingularRing(const SingularRing&) = delete;
  SingularRing& operator=(const SingularRing&) = delete;

  /** The underlying Singular ring. */
  ip_sring* raw() const;
  /** The underlying Singular coefficient field. */
  n_Procs_s* cf() const;
  /** The number of indeterminates. */
  size_t nVars() const;
  /** The field characteristic / size. */
  const Integer& prime() const;
  /** True iff this uses the small-field (n_Zp) path. */
  bool isSmall() const;

  /** Convert a Singular number to its canonical representative in [0,prime). */
  Integer numberToInteger(snumber* n) const;
  /** Convert an Integer to a (newly-allocated) Singular number in [0,prime). */
  snumber* integerToNumber(const Integer& i) const;

 private:
  Integer d_prime;
  size_t d_nVars;
  bool d_isSmall;
  n_Procs_s* d_cf;
  ip_sring* d_ring;
  /** Backing storage for the n_Zn base; only allocated/owned for big primes. */
  mpz_t d_base;
};

/**
 * Value-semantics RAII wrapper over a Singular `poly` together with its ring.
 *
 * A default-constructed Poly is "null" (no ring); only assignment is valid on
 * it. Copying deep-copies the underlying polynomial; the destructor frees it.
 *
 * The referenced SingularRing must outlive every Poly in it.
 */
class Poly
{
 public:
  Poly();
  Poly(const Poly&);
  Poly(Poly&&) noexcept;
  Poly& operator=(const Poly&);
  Poly& operator=(Poly&&) noexcept;
  ~Poly();

  /** Adopt an existing raw poly (consumes it; util-internal). */
  static Poly adopt(const SingularRing& r, spolyrec* p);

  /** The zero polynomial. */
  static Poly zero(const SingularRing& r);
  /** The one (unit) polynomial. */
  static Poly one(const SingularRing& r);
  /** A constant polynomial equal to the field element c. */
  static Poly constant(const SingularRing& r, const Integer& c);
  /** The indeterminate with the given 0-based index. */
  static Poly indet(const SingularRing& r, size_t varIdx);

  /** Is this the zero polynomial? */
  bool isZero() const;
  /** Is this the constant 1? */
  bool isOne() const;
  /**
   * Is this constant? True iff no monomial has a positive exponent (so a
   * nonzero constant, OR the zero polynomial, both return true).
   */
  bool isConstant() const;
  /** Total degree. Returns 0 for a nonzero constant and -1 for zero. */
  long deg() const;
  /** Number of terms. */
  size_t numTerms() const;
  /** String representation (Singular's, which uses a symmetric coeff repr). */
  std::string str() const;
  /**
   * The 0-based index of the unique variable that occurs (with positive
   * degree) across all terms, or -1 if zero/constant or more than one variable
   * occurs.
   */
  int univariateIndetIndex() const;
  /**
   * Coefficient of the all-exponents-zero term, as a canonical [0,prime)
   * integer; 0 if there is no such term.
   */
  Integer constantCoeff() const;
  /** Leading coefficient, as a canonical [0,prime) integer (0 for zero). */
  Integer leadingCoeff() const;
  /** This divided by its leading coefficient (so the result is monic). */
  Poly monic() const;
  /** This raised to the power e. */
  Poly power(unsigned long e) const;

  Poly operator+(const Poly&) const;
  Poly operator-(const Poly&) const;
  Poly operator*(const Poly&) const;
  Poly operator-() const;
  Poly& operator+=(const Poly&);
  Poly& operator*=(const Poly&);

  /**
   * Exact polynomial division (the divisor is assumed to divide this exactly).
   * For small fields uses singclap_pdivide; for big fields uses a hand-rolled
   * univariate schoolbook division.
   */
  Poly exactDiv(const Poly& divisor) const;

  bool operator==(const Poly&) const;

  /** The ring this polynomial lives in. */
  const SingularRing& sring() const;
  /** The underlying raw poly (util-internal use). */
  spolyrec* raw() const;

 private:
  void clear();
  const SingularRing* d_ring;
  spolyrec* d_poly;
};

/** A field element is represented as a constant Poly. */
using Scalar = Poly;
/** A list of polynomials. */
using Polys = std::vector<Poly>;
/** An input (point/vector) to a polynomial. */
using Point = std::vector<Scalar>;
/** A partial input (point/vector with optional entries) to a polynomial. */
using PartialPoint = std::vector<std::optional<Scalar>>;

/**
 * An ideal with a lazily-computed, cached reduced Groebner basis.
 *
 * The referenced SingularRing must outlive the Ideal.
 *
 * NB: Groebner-basis computation (kStd) is uninterruptible. Timeouts are the
 * caller's responsibility: check the ResourceManager before calling gbasis()
 * (and between findZero iterations) and throw FfTimeoutException there.
 */
class Ideal
{
 public:
  Ideal(const SingularRing& r, std::vector<Poly> generators);

  /** The ring this ideal lives in. */
  const SingularRing& sring() const;
  /** The generators this ideal was constructed from. */
  const std::vector<Poly>& gens() const;

  /** Compute (and cache) the reduced Groebner basis. */
  const std::vector<Poly>& gbasis();

  /** Is this the whole ring (GB is a single nonzero constant)? */
  bool isWholeRing();
  /** Is the variety zero-dimensional? */
  bool isZeroDim();
  /** Normal form of p modulo the (cached) Groebner basis. */
  Poly normalForm(const Poly& p);

  /**
   * Minimal polynomial of the indeterminate varIdx (0-based) in the quotient
   * ring (the ideal must be zero-dimensional).
   *
   * Returned as a monic univariate polynomial in that indeterminate.
   *
   * Computed by the normal-form power-dependency method: reduce 1, x, x^2,
   * ... mod the GB and find the first power whose normal form is linearly
   * dependent (over GF(p)) on the earlier ones; that dependency (normalized so
   * the highest power has coefficient 1) is the minimal polynomial. The degree
   * is bounded by the vector-space dimension of the quotient (scMult0Int).
   */
  Poly minimalPolynomial(size_t varIdx);

  /** The Krull dimension of the variety (0 means zero-dimensional). */
  long dimension();
  /** The vector-space dimension of the quotient ring (scMult0Int). */
  long vectorSpaceDim();

 private:
  /** Build a fresh Singular ideal from the cached GB (caller frees it). */
  sip_sideal* makeGbIdeal();

  const SingularRing* d_ring;
  std::vector<Poly> d_gens;
  std::vector<Poly> d_gb;
  bool d_hasGb;
};

/**
 * Given a univariate f over GF(p), return a list of its roots in the field, as
 * constant Polys, sorted by string representation.
 *
 * For small fields this uses Singular's factorization; for big fields it uses
 * the Berlekamp-Rabin algorithm.
 */
std::vector<Poly> roots(const Poly& f);

/**
 * Compute a monic polynomial q of minimal degree with the same root set as f
 * (a product of distinct linear factors). Used by the Rabin root finder.
 */
Poly distinctRootsPoly(const Poly& f);

/**
 * Univariate gcd over GF(p), returned monic. For small fields uses
 * singclap_gcd; for big fields uses a hand-rolled Euclidean algorithm.
 */
Poly polyGcd(const Poly& a, const Poly& b);

/**
 * Partial evaluation: substitute the given (optional) field values for the
 * indeterminates. Returns an empty optional if some occurring variable has no
 * value.
 */
std::optional<Scalar> singularEval(const Poly& p, const PartialPoint& values);

/** Total evaluation: substitute the given field values for all variables. */
Scalar singularEval(const Poly& p, const Point& values);

/** Convert a constant Poly (a field element) to a FiniteFieldValue. */
FiniteFieldValue toFfVal(const Scalar& elem, const FfSize& size);

}  // namespace ff
}  // namespace theory
}  // namespace cvc5::internal

#endif /* CVC5__THEORY__FF__SINGULAR_UTIL_H */

#endif /* CVC5_USE_SINGULAR */
