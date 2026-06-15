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
 */

#ifdef CVC5_USE_SINGULAR

#include "theory/ff/singular_util.h"

// external includes
#include <Singular/libsingular.h>
#include <coeffs/numbers.h>
#include <coeffs/rmodulon.h>
#include <coeffs/si_gmp.h>
#include <kernel/combinatorics/stairc.h>
#include <polys/clapsing.h>

// The largest prime for which Singular's factory small-field path (n_Zp) is
// usable. It is defined in <coeffs/modulop.h>, but that header transitively
// pulls in NTL headers that may not be on the include path; we provide it
// directly to avoid that fragile dependency.
#ifndef FACTORY_MAX_PRIME
#define FACTORY_MAX_PRIME 536870909
#endif

// std includes
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <unordered_map>

// internal includes
#include "base/check.h"
#include "base/output.h"
#include "util/random.h"
#include "util/singular_globals.h"

// ===========================================================================
// Custom large-prime FIELD coefficient domain for Singular.
//
// Singular's built-in coefficient types cannot represent GF(p) for a large
// prime p: n_Zp caps at p < 2^31, and the integers-mod-n type (n_Zn) is a
// *ring* engine that does not perform the coefficient divisions a field
// Groebner basis needs (it returns incomplete bases -- e.g. NF(1)=1 for a
// system where 1 is in the ideal -- which would make the FF solver UNSOUND).
//
// Singular's GB kernel (kStd/kNF/scDimInt) is coefficient-domain-generic: it
// calls all arithmetic through the coeffs vtable. So we register a custom
// GMP-backed field domain that implements true mod-p arithmetic (with a real
// modular inverse) and reports itself as a field; the SAME kernel then runs
// the field algorithm correctly over GF(p) for an arbitrary prime.
//
// CRITICAL subtlety (the reason a naive heap-mpz coeffs FAILS): Singular's
// Buchberger engine (initenterpairs) tests pGetComp(p) = p->exp[pCompIndex]
// to decide whether two polynomials share a module component before creating
// their critical pair. For a plain ideal pCompIndex == -1, so the access
// aliases the polynomial's COEFFICIENT word and the engine compares the raw
// bits of the two leading coefficients. n_Zp/n_Q work because equal small
// coefficients are bit-identical immediates; a fresh heap mpz_ptr gives two
// different pointers for two equal coefficients, so the (bogus) test fails,
// the pair is silently dropped, and the basis is incomplete. The fix is to
// INTERN residues: one canonical mpz_ptr per residue value, so equal field
// elements are bit-identical. (This is an internal-layout assumption of
// Singular 4.3.2; the litmus checks in the FF unit tests guard against a
// future change silently reintroducing incomplete bases.)
// ===========================================================================
namespace {

/** Per-coeffs data: the prime and the residue-interning table. */
struct BigField
{
  mpz_t prime;
  /** decimal residue in [0,prime) -> the canonical mpz_ptr for that value. */
  std::unordered_map<std::string, mpz_ptr> intern;
};
/** Parameter passed to nInitChar to construct a BigField coeffs. */
struct BigFieldInfo
{
  mpz_ptr prime;
};

inline BigField* BF(const coeffs r) { return static_cast<BigField*>(r->data); }
inline mpz_ptr BF_P(const coeffs r) { return BF(r)->prime; }

/**
 * Given a freshly computed residue z (already reduced into [0,prime)), return
 * the unique interned mpz_ptr for that value, freeing z if a canonical copy
 * already exists. This interning is what makes equal field elements
 * bit-identical (see the header comment).
 */
number bf_canon(mpz_ptr z, const coeffs r)
{
  BigField* bf = BF(r);
  char* s = mpz_get_str(nullptr, 10, z);
  std::string key(s);
  void (*freefn)(void*, size_t);
  mp_get_memory_functions(nullptr, nullptr, &freefn);
  freefn(s, key.size() + 1);

  auto it = bf->intern.find(key);
  if (it != bf->intern.end())
  {
    mpz_clear(z);
    omFree((void*)z);
    return (number)it->second;
  }
  bf->intern.emplace(std::move(key), z);
  return (number)z;
}
inline mpz_ptr bf_tmp()
{
  mpz_ptr z = (mpz_ptr)omAlloc(sizeof(mpz_t));
  mpz_init(z);
  return z;
}

number bfInit(long i, const coeffs r)
{
  mpz_ptr z = bf_tmp();
  mpz_set_si(z, i);
  mpz_mod(z, z, BF_P(r));
  return bf_canon(z, r);
}
number bfInitMPZ(mpz_t m, const coeffs r)
{
  mpz_ptr z = bf_tmp();
  mpz_set(z, m);
  mpz_mod(z, z, BF_P(r));
  return bf_canon(z, r);
}
void bfDelete(number* a, const coeffs) { *a = nullptr; }  // interned: no free
number bfCopy(number a, const coeffs) { return a; }       // interned: share

BOOLEAN bfIsZero(number a, const coeffs) { return mpz_sgn((mpz_ptr)a) == 0; }
BOOLEAN bfIsOne(number a, const coeffs) { return mpz_cmp_si((mpz_ptr)a, 1) == 0; }
BOOLEAN bfIsMOne(number a, const coeffs r)
{
  mpz_t t;
  mpz_init(t);
  mpz_sub_ui(t, BF_P(r), 1);
  BOOLEAN res = (mpz_cmp((mpz_ptr)a, t) == 0);
  mpz_clear(t);
  return res;
}
BOOLEAN bfEqual(number a, number b, const coeffs)
{
  return mpz_cmp((mpz_ptr)a, (mpz_ptr)b) == 0;
}
number bfAdd(number a, number b, const coeffs r)
{
  mpz_ptr z = bf_tmp();
  mpz_add(z, (mpz_ptr)a, (mpz_ptr)b);
  mpz_mod(z, z, BF_P(r));
  return bf_canon(z, r);
}
number bfSub(number a, number b, const coeffs r)
{
  mpz_ptr z = bf_tmp();
  mpz_sub(z, (mpz_ptr)a, (mpz_ptr)b);
  mpz_mod(z, z, BF_P(r));
  return bf_canon(z, r);
}
number bfMult(number a, number b, const coeffs r)
{
  mpz_ptr z = bf_tmp();
  mpz_mul(z, (mpz_ptr)a, (mpz_ptr)b);
  mpz_mod(z, z, BF_P(r));
  return bf_canon(z, r);
}
number bfInpNeg(number a, const coeffs r)
{
  mpz_ptr z = bf_tmp();
  if (mpz_sgn((mpz_ptr)a) != 0)
    mpz_sub(z, BF_P(r), (mpz_ptr)a);
  else
    mpz_set_ui(z, 0);
  return bf_canon(z, r);
}
number bfInvers(number a, const coeffs r)
{
  mpz_ptr z = bf_tmp();
  if (mpz_sgn((mpz_ptr)a) == 0)
  {
    WerrorS(nDivBy0);
    mpz_set_ui(z, 0);
    return bf_canon(z, r);
  }
  mpz_invert(z, (mpz_ptr)a, BF_P(r));
  return bf_canon(z, r);
}
number bfDiv(number a, number b, const coeffs r)
{
  mpz_ptr z = bf_tmp();
  if (mpz_sgn((mpz_ptr)b) == 0)
  {
    WerrorS(nDivBy0);
    mpz_set_ui(z, 0);
    return bf_canon(z, r);
  }
  mpz_t inv;
  mpz_init(inv);
  mpz_invert(inv, (mpz_ptr)b, BF_P(r));
  mpz_mul(z, (mpz_ptr)a, inv);
  mpz_mod(z, z, BF_P(r));
  mpz_clear(inv);
  return bf_canon(z, r);
}
void bfInpMult(number& a, number b, const coeffs r) { a = bfMult(a, b, r); }
void bfInpAdd(number& a, number b, const coeffs r) { a = bfAdd(a, b, r); }

/** Symmetric representative in (-prime/2, prime/2], like n_Zp's npInt. */
long bfInt(number& a, const coeffs r)
{
  mpz_t t, half;
  mpz_init_set(t, (mpz_ptr)a);
  mpz_init(half);
  mpz_fdiv_q_2exp(half, BF_P(r), 1);
  if (mpz_cmp(t, half) > 0) mpz_sub(t, t, BF_P(r));
  long res = mpz_get_si(t);
  mpz_clear(t);
  mpz_clear(half);
  return res;
}
void bfMPZ(mpz_t result, number& a, const coeffs)
{
  mpz_set(result, (mpz_ptr)a);  // representative in [0,prime)
}
/** Field gcd: gcd of nonzeros is 1 (0 only if both are zero). */
number bfGcd(number a, number b, const coeffs r)
{
  mpz_ptr z = bf_tmp();
  if (mpz_sgn((mpz_ptr)a) == 0 && mpz_sgn((mpz_ptr)b) == 0)
    mpz_set_ui(z, 0);
  else
    mpz_set_ui(z, 1);
  return bf_canon(z, r);
}
number bfLcm(number a, number b, const coeffs r) { return bfGcd(a, b, r); }
int bfSize(number a, const coeffs)
{
  if (mpz_sgn((mpz_ptr)a) == 0) return 0;
  return (int)mpz_size((mpz_ptr)a) + 1;
}
void bfPower(number a, int i, number* result, const coeffs r)
{
  mpz_ptr z = bf_tmp();
  if (i < 0)
  {
    mpz_t inv;
    mpz_init(inv);
    mpz_invert(inv, (mpz_ptr)a, BF_P(r));
    mpz_powm_ui(z, inv, (unsigned long)(-(long)i), BF_P(r));
    mpz_clear(inv);
  }
  else
  {
    mpz_powm_ui(z, (mpz_ptr)a, (unsigned long)i, BF_P(r));
  }
  *result = bf_canon(z, r);
}
BOOLEAN bfGreaterZero(number a, const coeffs r)
{
  mpz_t half;
  mpz_init(half);
  mpz_fdiv_q_2exp(half, BF_P(r), 1);
  BOOLEAN res = (mpz_cmp((mpz_ptr)a, half) <= 0);
  mpz_clear(half);
  return res;
}
BOOLEAN bfGreater(number a, number b, const coeffs)
{
  return mpz_cmp((mpz_ptr)a, (mpz_ptr)b) > 0;
}
void bfWriteLong(number a, const coeffs r)
{
  mpz_t t, half;
  mpz_init_set(t, (mpz_ptr)a);
  mpz_init(half);
  mpz_fdiv_q_2exp(half, BF_P(r), 1);
  if (mpz_cmp(t, half) > 0) mpz_sub(t, t, BF_P(r));
  char* s = mpz_get_str(nullptr, 10, t);
  StringAppendS(s);
  void (*freefn)(void*, size_t);
  mp_get_memory_functions(nullptr, nullptr, &freefn);
  freefn(s, std::strlen(s) + 1);
  mpz_clear(t);
  mpz_clear(half);
}
void bfNormalize(number&, const coeffs) {}

const char* bfRead(const char* s, number* a, const coeffs r)
{
  mpz_t z;
  mpz_init(z);
  const char* p = s;
  bool any = false;
  while (*p >= '0' && *p <= '9')
  {
    mpz_mul_ui(z, z, 10);
    mpz_add_ui(z, z, *p - '0');
    ++p;
    any = true;
  }
  if (!any) mpz_set_ui(z, 1);
  mpz_ptr t = bf_tmp();
  mpz_set(t, z);
  mpz_mod(t, t, BF_P(r));
  mpz_clear(z);
  *a = bf_canon(t, r);
  return p;
}

BOOLEAN bfCoeffIsEqual(const coeffs r, n_coeffType n, void* parameter)
{
  if (n != r->type) return FALSE;
  BigFieldInfo* info = (BigFieldInfo*)parameter;
  return mpz_cmp(BF_P(r), info->prime) == 0;
}
char* bfCoeffName(const coeffs r)
{
  static char buf[512];
  char* s = mpz_get_str(nullptr, 10, BF_P(r));
  snprintf(buf, sizeof(buf), "ZZ/%s", s);
  void (*freefn)(void*, size_t);
  mp_get_memory_functions(nullptr, nullptr, &freefn);
  freefn(s, std::strlen(s) + 1);
  return buf;
}
void bfCoeffWrite(const coeffs r, BOOLEAN) { PrintS(bfCoeffName(r)); }
void bfKillChar(coeffs r)
{
  BigField* bf = BF(r);
  if (bf != nullptr)
  {
    for (auto& kv : bf->intern)
    {
      mpz_clear(kv.second);
      omFree((void*)kv.second);
    }
    mpz_clear(bf->prime);
    delete bf;
    r->data = nullptr;
  }
}

number bfMapIdent(number from, const coeffs, const coeffs) { return from; }
number bfMapFromZ(number from, const coeffs, const coeffs dst)
{
  mpz_ptr z = bf_tmp();
  mpz_set(z, (mpz_ptr)from);
  mpz_mod(z, z, BF_P(dst));
  return bf_canon(z, dst);
}
number bfMapFromQ(number from, const coeffs src, const coeffs dst)
{
  mpz_t num, den;
  mpz_init(num);
  mpz_init(den);
  number fn = n_GetNumerator(from, src), fd = n_GetDenom(from, src);
  n_MPZ(num, fn, src);
  n_MPZ(den, fd, src);
  n_Delete(&fn, src);
  n_Delete(&fd, src);
  mpz_mod(num, num, BF_P(dst));
  if (mpz_cmp_ui(den, 1) != 0)
  {
    mpz_t dinv;
    mpz_init(dinv);
    mpz_invert(dinv, den, BF_P(dst));
    mpz_mul(num, num, dinv);
    mpz_clear(dinv);
  }
  mpz_ptr z = bf_tmp();
  mpz_set(z, num);
  mpz_mod(z, z, BF_P(dst));
  mpz_clear(num);
  mpz_clear(den);
  return bf_canon(z, dst);
}
nMapFunc bfSetMap(const coeffs src, const coeffs dst)
{
  if (src == dst) return bfMapIdent;
  if (getCoeffType(src) == n_Z) return bfMapFromZ;
  if (getCoeffType(src) == n_Q) return bfMapFromQ;
  if (getCoeffType(src) == dst->type && mpz_cmp(BF_P(src), BF_P(dst)) == 0)
    return bfMapIdent;
  return nullptr;
}

/** The InitChar proc: fills the vtable and flags from a BigFieldInfo*. */
BOOLEAN bfInitChar(coeffs r, void* p)
{
  BigFieldInfo* info = (BigFieldInfo*)p;
  BigField* bf = new BigField();
  mpz_init_set(bf->prime, info->prime);
  r->data = (void*)bf;

  r->is_field = TRUE;  // routes kStd to the FIELD engine
  r->is_domain = TRUE;
  r->rep = n_rep_gmp;  // number == mpz_ptr
  r->ch = 0;           // 254-bit p doesn't fit int ch; harmless on the field
                       // path (is_field=TRUE keeps us off the ring code path).
  r->has_simple_Alloc = TRUE;    // interned & immutable: Copy shares, Delete no-op
  r->has_simple_Inverse = TRUE;  // inverse is cheap (mpz_invert)

  r->cfInit = bfInit;
  r->cfInitMPZ = bfInitMPZ;
  r->cfInt = bfInt;
  r->cfMPZ = bfMPZ;
  r->cfDelete = bfDelete;
  r->cfCopy = bfCopy;
  r->cfAdd = bfAdd;
  r->cfInpAdd = bfInpAdd;
  r->cfSub = bfSub;
  r->cfMult = bfMult;
  r->cfInpMult = bfInpMult;
  r->cfDiv = bfDiv;
  r->cfExactDiv = bfDiv;
  r->cfInpNeg = bfInpNeg;
  r->cfInvers = bfInvers;
  r->cfIsZero = bfIsZero;
  r->cfIsOne = bfIsOne;
  r->cfIsMOne = bfIsMOne;
  r->cfEqual = bfEqual;
  r->cfGreaterZero = bfGreaterZero;
  r->cfGreater = bfGreater;
  r->cfSize = bfSize;
  r->cfPower = bfPower;
  r->cfGcd = bfGcd;
  r->cfLcm = bfLcm;
  r->cfNormalize = bfNormalize;
  r->cfWriteLong = bfWriteLong;
  r->cfRead = bfRead;
  r->cfSetMap = bfSetMap;
  r->cfCoeffName = bfCoeffName;
  r->cfCoeffWrite = bfCoeffWrite;
  r->nCoeffIsEqual = bfCoeffIsEqual;
  r->cfKillChar = bfKillChar;

  return FALSE;
}

/** The (process-global) coefficient-type id for our big-prime field domain. */
n_coeffType bigFieldType()
{
  // C++11 guarantees thread-safe one-time initialization of this static.
  static n_coeffType s_type = nRegister(n_unknown, bfInitChar);
  return s_type;
}

}  // namespace

namespace cvc5::internal {
namespace theory {
namespace ff {

/* -------------------------------------------------------------------------- */
/* SingularRing                                                               */
/* -------------------------------------------------------------------------- */

SingularRing::SingularRing(const Integer& prime,
                           const std::vector<std::string>& varNames)
    : d_prime(prime), d_nVars(varNames.size())
{
  initSingular();
  d_isSmall = (prime <= Integer(FACTORY_MAX_PRIME));
  if (d_isSmall)
  {
    d_cf = nInitChar(n_Zp, (void*)(long)prime.getLong());
  }
  else
  {
    // A genuine prime FIELD over a large prime, via our custom GMP-backed
    // coeffs (see the BigField anonymous namespace above). Singular's built-in
    // big-modulus type (n_Zn) is a ring engine whose Groebner bases are
    // incomplete here, which would make the solver unsound.
    mpz_init_set_str(d_base, prime.toString().c_str(), 10);
    BigFieldInfo info;
    info.prime = d_base;
    d_cf = nInitChar(bigFieldType(), (void*)&info);
  }
  char** names = (char**)omAlloc(varNames.size() * sizeof(char*));
  for (size_t i = 0; i < varNames.size(); ++i)
  {
    names[i] = omStrDup(varNames[i].c_str());
  }
  // rDefault takes ownership of d_cf and names.
  d_ring = rDefault(d_cf, (int)varNames.size(), names, ringorder_dp);
  rChangeCurrRing(d_ring);
}

SingularRing::~SingularRing()
{
  // rDelete frees the coefficient field and the variable names too.
  rDelete(d_ring);
  if (!d_isSmall)
  {
    mpz_clear(d_base);
  }
}

ring SingularRing::raw() const { return d_ring; }
coeffs SingularRing::cf() const { return d_cf; }
size_t SingularRing::nVars() const { return d_nVars; }
const Integer& SingularRing::prime() const { return d_prime; }
bool SingularRing::isSmall() const { return d_isSmall; }

Integer SingularRing::numberToInteger(number n) const
{
  if (d_isSmall)
  {
    // n_Int gives a symmetric representative in (-p/2, p/2].
    long v = n_Int(n, d_cf);
    Integer r(v);
    if (r.sgn() < 0)
    {
      r += d_prime;
    }
    return r;
  }
  // n_MPZ gives a representative in [0, p) directly. We go through a decimal
  // string so as not to depend on cvc5's Integer being GMP-backed (it may be
  // CLN-backed), matching the impl-agnostic approach of the CoCoA code.
  mpz_t z;
  mpz_init(z);
  n_MPZ(z, n, d_cf);
  std::vector<char> buf(mpz_sizeinbase(z, 10) + 2);
  mpz_get_str(buf.data(), 10, z);
  mpz_clear(z);
  return Integer(std::string(buf.data()), 10);
}

number SingularRing::integerToNumber(const Integer& i) const
{
  // canonicalize into [0, prime)
  Integer v = i.floorDivideRemainder(d_prime);
  if (d_isSmall)
  {
    return n_Init(v.getLong(), d_cf);
  }
  mpz_t z;
  mpz_init_set_str(z, v.toString().c_str(), 10);
  number n = n_InitMPZ(z, d_cf);
  mpz_clear(z);
  return n;
}

/* -------------------------------------------------------------------------- */
/* Poly                                                                       */
/* -------------------------------------------------------------------------- */

Poly::Poly() : d_ring(nullptr), d_poly(nullptr) {}

Poly::Poly(const Poly& o) : d_ring(o.d_ring), d_poly(nullptr)
{
  if (d_ring != nullptr)
  {
    d_poly = p_Copy(o.d_poly, d_ring->raw());
  }
}

Poly::Poly(Poly&& o) noexcept : d_ring(o.d_ring), d_poly(o.d_poly)
{
  o.d_ring = nullptr;
  o.d_poly = nullptr;
}

Poly& Poly::operator=(const Poly& o)
{
  if (this != &o)
  {
    clear();
    d_ring = o.d_ring;
    if (d_ring != nullptr)
    {
      d_poly = p_Copy(o.d_poly, d_ring->raw());
    }
  }
  return *this;
}

Poly& Poly::operator=(Poly&& o) noexcept
{
  if (this != &o)
  {
    clear();
    d_ring = o.d_ring;
    d_poly = o.d_poly;
    o.d_ring = nullptr;
    o.d_poly = nullptr;
  }
  return *this;
}

Poly::~Poly() { clear(); }

void Poly::clear()
{
  if (d_ring != nullptr && d_poly != nullptr)
  {
    p_Delete(&d_poly, d_ring->raw());
  }
  d_poly = nullptr;
}

Poly Poly::adopt(const SingularRing& r, poly p)
{
  Poly out;
  out.d_ring = &r;
  out.d_poly = p;
  return out;
}

Poly Poly::zero(const SingularRing& r) { return adopt(r, nullptr); }

Poly Poly::one(const SingularRing& r) { return adopt(r, p_ISet(1, r.raw())); }

Poly Poly::constant(const SingularRing& r, const Integer& c)
{
  number n = r.integerToNumber(c);
  // p_NSet consumes n.
  return adopt(r, p_NSet(n, r.raw()));
}

Poly Poly::indet(const SingularRing& r, size_t varIdx)
{
  poly m = p_ISet(1, r.raw());
  p_SetExp(m, (int)varIdx + 1, 1, r.raw());  // 1-based
  p_Setm(m, r.raw());
  return adopt(r, m);
}

bool Poly::isZero() const { return d_poly == nullptr; }

bool Poly::isOne() const
{
  return d_poly != nullptr && p_IsOne(d_poly, d_ring->raw());
}

bool Poly::isConstant() const
{
  if (d_poly == nullptr)
  {
    return true;
  }
  ring R = d_ring->raw();
  int nv = rVar(R);
  for (poly q = d_poly; q != nullptr; q = pNext(q))
  {
    for (int v = 1; v <= nv; ++v)
    {
      if (p_GetExp(q, v, R) > 0)
      {
        return false;
      }
    }
  }
  return true;
}

long Poly::deg() const
{
  if (d_poly == nullptr)
  {
    return -1;  // <= 0 for the zero polynomial
  }
  return p_Deg(d_poly, d_ring->raw());
}

size_t Poly::numTerms() const { return (size_t)pLength(d_poly); }

std::string Poly::str() const
{
  return std::string(p_String(d_poly, d_ring->raw()));
}

int Poly::univariateIndetIndex() const
{
  if (d_poly == nullptr)
  {
    return -1;
  }
  ring R = d_ring->raw();
  int nv = rVar(R);
  int found = -1;
  for (poly q = d_poly; q != nullptr; q = pNext(q))
  {
    for (int v = 1; v <= nv; ++v)
    {
      if (p_GetExp(q, v, R) > 0)
      {
        if (found == -1)
        {
          found = v;
        }
        else if (found != v)
        {
          return -1;
        }
      }
    }
  }
  return found == -1 ? -1 : found - 1;  // 0-based
}

Integer Poly::constantCoeff() const
{
  ring R = d_ring->raw();
  int nv = rVar(R);
  for (poly q = d_poly; q != nullptr; q = pNext(q))
  {
    bool isConst = true;
    for (int v = 1; v <= nv; ++v)
    {
      if (p_GetExp(q, v, R) != 0)
      {
        isConst = false;
        break;
      }
    }
    if (isConst)
    {
      return d_ring->numberToInteger(pGetCoeff(q));
    }
  }
  return Integer(0);
}

Integer Poly::leadingCoeff() const
{
  if (d_poly == nullptr)
  {
    return Integer(0);
  }
  return d_ring->numberToInteger(pGetCoeff(d_poly));
}

Poly Poly::monic() const
{
  ring R = d_ring->raw();
  coeffs cf = d_ring->cf();
  poly p = p_Copy(d_poly, R);
  number lc = n_Copy(pGetCoeff(p), cf);
  number inv = n_Invers(lc, cf);
  p = p_Mult_nn(p, inv, R);
  n_Delete(&lc, cf);
  n_Delete(&inv, cf);
  return adopt(*d_ring, p);
}

Poly Poly::power(unsigned long e) const
{
  ring R = d_ring->raw();
  return adopt(*d_ring, p_Power(p_Copy(d_poly, R), (int)e, R));
}

Poly Poly::operator+(const Poly& o) const
{
  ring R = d_ring->raw();
  poly r = p_Add_q(p_Copy(d_poly, R), p_Copy(o.d_poly, R), R);
  return adopt(*d_ring, r);
}

Poly Poly::operator-(const Poly& o) const
{
  ring R = d_ring->raw();
  poly negO = p_Neg(p_Copy(o.d_poly, R), R);
  poly r = p_Add_q(p_Copy(d_poly, R), negO, R);
  return adopt(*d_ring, r);
}

Poly Poly::operator*(const Poly& o) const
{
  ring R = d_ring->raw();
  poly r = pp_Mult_qq(d_poly, o.d_poly, R);
  return adopt(*d_ring, r);
}

Poly Poly::operator-() const
{
  ring R = d_ring->raw();
  return adopt(*d_ring, p_Neg(p_Copy(d_poly, R), R));
}

Poly& Poly::operator+=(const Poly& o)
{
  ring R = d_ring->raw();
  d_poly = p_Add_q(d_poly, p_Copy(o.d_poly, R), R);
  return *this;
}

Poly& Poly::operator*=(const Poly& o)
{
  *this = (*this) * o;
  return *this;
}

namespace {

/**
 * Hand-rolled univariate schoolbook polynomial division over the (prime)
 * field, used on the big-field path where factory routines fail. Computes the
 * quotient if `wantQuotient`, else the remainder, of num / den. Both are
 * assumed univariate in the same indeterminate (or den constant).
 */
Poly univDivMod(const Poly& num, const Poly& den, bool wantQuotient)
{
  const SingularRing& r = num.sring();
  ring R = r.raw();
  coeffs cf = r.cf();
  int varIdx = den.univariateIndetIndex();  // 0-based; -1 if constant
  int v1 = varIdx + 1;                        // 1-based; 0 if constant

  long denDeg = den.raw() != nullptr ? p_Deg(den.raw(), R) : 0;
  number denLc = n_Copy(pGetCoeff(den.raw()), cf);
  number denLcInv = n_Invers(denLc, cf);

  // A nonzero constant divides everything in a field: the remainder is 0, and
  // the quotient is num scaled by the inverse of that constant. (Without this
  // special case the loop below would never reduce the degree of a non-constant
  // dividend, and would spin forever.)
  if (denDeg == 0)
  {
    n_Delete(&denLc, cf);
    if (wantQuotient)
    {
      poly q = p_Mult_nn(p_Copy(num.raw(), R), denLcInv, R);
      n_Delete(&denLcInv, cf);
      return Poly::adopt(r, q);
    }
    n_Delete(&denLcInv, cf);
    return Poly::zero(r);
  }

  poly rem = p_Copy(num.raw(), R);
  poly quot = nullptr;

  while (rem != nullptr && p_Deg(rem, R) >= denDeg)
  {
    long remDeg = p_Deg(rem, R);
    number remLc = pGetCoeff(rem);
    number coef = n_Mult(remLc, denLcInv, cf);
    poly term = p_NSet(n_Copy(coef, cf), R);
    if (v1 > 0 && remDeg - denDeg > 0)
    {
      p_SetExp(term, v1, (int)(remDeg - denDeg), R);
      p_Setm(term, R);
    }
    n_Delete(&coef, cf);
    if (wantQuotient)
    {
      quot = p_Add_q(quot, p_Copy(term, R), R);
    }
    poly sub = pp_Mult_qq(term, den.raw(), R);
    rem = p_Add_q(rem, p_Neg(sub, R), R);
    p_Delete(&term, R);
  }
  n_Delete(&denLc, cf);
  n_Delete(&denLcInv, cf);
  if (wantQuotient)
  {
    p_Delete(&rem, R);
    return Poly::adopt(r, quot);
  }
  return Poly::adopt(r, rem);
}

/** Univariate remainder: num mod den. */
Poly univRem(const Poly& num, const Poly& den)
{
  return univDivMod(num, den, /*wantQuotient=*/false);
}

}  // namespace

Poly Poly::exactDiv(const Poly& divisor) const
{
  ring R = d_ring->raw();
  if (d_ring->isSmall())
  {
    poly q = singclap_pdivide(p_Copy(d_poly, R), p_Copy(divisor.d_poly, R), R);
    return adopt(*d_ring, q);
  }
  // hand-rolled univariate schoolbook exact division over n_Zn
  return univDivMod(*this, divisor, /*wantQuotient=*/true);
}

bool Poly::operator==(const Poly& o) const
{
  if (d_poly == nullptr && o.d_poly == nullptr)
  {
    return true;
  }
  if (d_ring == nullptr || o.d_ring == nullptr)
  {
    return false;
  }
  return p_EqualPolys(d_poly, o.d_poly, d_ring->raw());
}

const SingularRing& Poly::sring() const { return *d_ring; }
poly Poly::raw() const { return d_poly; }

/* -------------------------------------------------------------------------- */
/* polyGcd, redMod, powerMod                                                  */
/* -------------------------------------------------------------------------- */

Poly polyGcd(const Poly& a, const Poly& b)
{
  const SingularRing& r = a.sring();
  ring R = r.raw();
  if (r.isSmall())
  {
    poly g = singclap_gcd(p_Copy(a.raw(), R), p_Copy(b.raw(), R), R);
    Poly result = Poly::adopt(r, g);
    if (result.isZero())
    {
      return result;
    }
    return result.monic();
  }
  // Euclidean algorithm over the field.
  Poly x = a;
  Poly y = b;
  while (!y.isZero())
  {
    Poly t = univRem(x, y);
    x = y;
    y = t;
  }
  if (x.isZero())
  {
    return x;
  }
  return x.monic();
}

namespace {

/** Reduce b modulo m (treating {m} as the modulus), via univariate remainder. */
Poly redMod(const Poly& b, const Poly& m) { return univRem(b, m); }

/** Compute b^e mod m, by repeated squaring with reductions by m each step. */
Poly powerMod(const Poly& b, const Integer& e, const Poly& m)
{
  const SingularRing& r = b.sring();
  Poly acc = Poly::one(r);
  Poly bPower = b;
  Integer ee = e;
  Integer two(2);
  while (ee.sgn() != 0)
  {
    if (ee.isBitSet(0))
    {
      acc = acc * bPower;
      acc = redMod(acc, m);
    }
    bPower = bPower * bPower;
    bPower = redMod(bPower, m);
    ee = ee.floorDivideQuotient(two);
  }
  return acc;
}

/** A random Integer in [0, hi], built from cvc5's Random. */
Integer randomBigInt(const Integer& hi)
{
  Integer range = hi + Integer(1);
  // Build a random non-negative Integer with enough bits, then reduce mod
  // range. The slight modulo bias is irrelevant for the Rabin sampler.
  size_t bits = range.length() + 64;
  Integer acc(0);
  Integer shift = Integer(2).pow(64);
  for (size_t produced = 0; produced < bits; produced += 64)
  {
    uint64_t word = Random::getRandom()();
    acc = acc * shift + Integer(word);
  }
  return acc.floorDivideRemainder(range);
}

/** Sort polys by their string representation (stable), mirroring CoCoA. */
std::vector<Poly> sortHack(const std::vector<Poly>& values)
{
  std::vector<std::string> strs;
  std::unordered_map<std::string, size_t> origIndices;
  for (const auto& v : values)
  {
    std::string s = v.str();
    origIndices.emplace(s, strs.size());
    strs.push_back(s);
  }
  std::sort(strs.begin(), strs.end());
  std::vector<Poly> output;
  for (const auto& s : strs)
  {
    output.push_back(values[origIndices[s]]);
  }
  return output;
}

}  // namespace

/* -------------------------------------------------------------------------- */
/* roots, distinctRootsPoly                                                   */
/* -------------------------------------------------------------------------- */

Poly distinctRootsPoly(const Poly& f)
{
  const SingularRing& r = f.sring();
  int idx = f.univariateIndetIndex();
  Assert(idx >= 0);
  Poly x = Poly::indet(r, (size_t)idx);
  // Prime field, so q == characteristic and LogCardinality == 1.
  Integer q = r.prime();
  Poly fieldPoly = powerMod(x, q, f) - x;
  return polyGcd(f, fieldPoly);
}

std::vector<Poly> roots(const Poly& f)
{
  const SingularRing& r = f.sring();
  ring R = r.raw();
  int idx = f.univariateIndetIndex();
  Assert(idx >= 0);
  Poly x = Poly::indet(r, (size_t)idx);
  Integer q = r.prime();
  std::vector<Poly> output;

  if (r.isSmall())
  {
    // Use Singular's factorization (only valid for small fields).
    intvec* v = nullptr;
    ideal facs = singclap_factorize(p_Copy(f.raw(), R), &v, 0, R);
    for (int i = 0; i < IDELEMS(facs); ++i)
    {
      Poly factor = Poly::adopt(r, p_Copy(facs->m[i], R));
      if (factor.deg() == 1)
      {
        // monic x - a => root a = -constantCoeff
        Poly m = factor.monic();
        Integer cc = m.constantCoeff();
        Integer root = (-cc).floorDivideRemainder(r.prime());
        output.push_back(Poly::constant(r, root));
      }
    }
    id_Delete(&facs, R);
    if (v != nullptr)
    {
      delete v;
    }
  }
  else
  {
    // Berlekamp-Rabin root finding (big prime field).
    Assert(q.isBitSet(0)) << "Rabin requires an odd field size";
    Integer s = q.floorDivideQuotient(Integer(2));
    std::vector<Poly> toFactor;
    toFactor.push_back(distinctRootsPoly(f));

    while (!toFactor.empty())
    {
      Poly p = toFactor.back();
      toFactor.pop_back();
      Trace("ff::roots") << "toFactor " << p.str() << std::endl;
      if (p.deg() == 0)
      {
        // constant: no factors
      }
      else if (p.constantCoeff().sgn() == 0)
      {
        // It has a zero root.
        output.push_back(Poly::constant(r, Integer(0)));
        toFactor.push_back(p.exactDiv(x));
      }
      else if (p.deg() == 1)
      {
        // It is linear; root = -constantCoeff (after making monic).
        Poly m = p.monic();
        Integer cc = m.constantCoeff();
        Integer root = (-cc).floorDivideRemainder(r.prime());
        output.push_back(Poly::constant(r, root));
      }
      else
      {
        // Super-linear, no zero-root: split with random delta.
        while (true)
        {
          Integer deltaInt = randomBigInt(q - Integer(1));
          Poly delta = Poly::constant(r, deltaInt);
          // split(X) = (X - delta)^(q/2) mod p - 1
          Poly base = x - delta;
          Poly split = powerMod(base, s, p) - Poly::one(r);
          Poly h = polyGcd(p, split);
          if (h.deg() > 0 && h.deg() < p.deg())
          {
            toFactor.push_back(h);
            toFactor.push_back(p.exactDiv(h));
            break;
          }
          // else: guess a new delta
        }
      }
    }
  }
  return sortHack(output);
}

/* -------------------------------------------------------------------------- */
/* Ideal                                                                      */
/* -------------------------------------------------------------------------- */

Ideal::Ideal(const SingularRing& r, std::vector<Poly> generators)
    : d_ring(&r), d_gens(std::move(generators)), d_hasGb(false)
{
}

const SingularRing& Ideal::sring() const { return *d_ring; }
const std::vector<Poly>& Ideal::gens() const { return d_gens; }

const std::vector<Poly>& Ideal::gbasis()
{
  if (d_hasGb)
  {
    return d_gb;
  }
  ring R = d_ring->raw();
  int n = (int)d_gens.size();
  ideal I = idInit(n < 1 ? 1 : n, 1);
  for (int i = 0; i < n; ++i)
  {
    I->m[i] = p_Copy(d_gens[i].raw(), R);
  }
  si_opt_1 |= Sy_bit(OPT_REDSB);
  ideal G = kStd(I, NULL, testHomog, NULL);
  idSkipZeroes(G);
  d_gb.clear();
  for (int i = 0; i < IDELEMS(G); ++i)
  {
    d_gb.push_back(Poly::adopt(*d_ring, p_Copy(G->m[i], R)));
  }
  id_Delete(&I, R);
  id_Delete(&G, R);
  d_hasGb = true;
  return d_gb;
}

ideal Ideal::makeGbIdeal()
{
  ring R = d_ring->raw();
  int n = (int)d_gb.size();
  ideal G = idInit(n < 1 ? 1 : n, 1);
  for (int i = 0; i < n; ++i)
  {
    G->m[i] = p_Copy(d_gb[i].raw(), R);
  }
  return G;
}

bool Ideal::isWholeRing()
{
  const auto& g = gbasis();
  return g.size() == 1 && !g[0].isZero() && g[0].deg() <= 0;
}

bool Ideal::isZeroDim() { return dimension() == 0; }

Poly Ideal::normalForm(const Poly& p)
{
  gbasis();
  ring R = d_ring->raw();
  ideal G = makeGbIdeal();
  poly nf = kNF(G, NULL, p_Copy(p.raw(), R));
  Poly result = Poly::adopt(*d_ring, nf);
  id_Delete(&G, R);
  return result;
}

long Ideal::dimension()
{
  gbasis();
  ring R = d_ring->raw();
  ideal G = makeGbIdeal();
  long d = d_ring->isSmall() ? scDimInt(G, NULL) : scDimIntRing(G, NULL);
  id_Delete(&G, R);
  return d;
}

long Ideal::vectorSpaceDim()
{
  gbasis();
  ring R = d_ring->raw();
  ideal G = makeGbIdeal();
  long d = scMult0Int(G, NULL);
  id_Delete(&G, R);
  return d;
}

namespace {

/**
 * Build a coefficient map for a polynomial, keyed by a string encoding of the
 * monomial's exponent vector. Coefficients are canonical [0,prime) integers.
 */
std::map<std::string, Integer> polyToCoeffMap(const Poly& p)
{
  std::map<std::string, Integer> out;
  const SingularRing& r = p.sring();
  ring R = r.raw();
  int nv = rVar(R);
  for (poly q = p.raw(); q != nullptr; q = pNext(q))
  {
    std::string key;
    for (int v = 1; v <= nv; ++v)
    {
      key += std::to_string(p_GetExp(q, v, R));
      key += ",";
    }
    out[key] = r.numberToInteger(pGetCoeff(q));
  }
  return out;
}

/** row -= factor * other (mod p), over a sparse coefficient map. */
void subtractScaled(std::map<std::string, Integer>& row,
                    const std::map<std::string, Integer>& other,
                    const Integer& factor,
                    const Integer& p)
{
  for (const auto& [mon, c] : other)
  {
    Integer cur(0);
    auto it = row.find(mon);
    if (it != row.end())
    {
      cur = it->second;
    }
    // cur - factor*c (mod p)
    Integer prod = factor.modMultiply(c, p);
    Integer neg = (p - prod).floorDivideRemainder(p);
    Integer nv = cur.modAdd(neg, p);
    if (nv.sgn() == 0)
    {
      if (it != row.end())
      {
        row.erase(it);
      }
    }
    else
    {
      row[mon] = nv;
    }
  }
}

}  // namespace

Poly Ideal::minimalPolynomial(size_t varIdx)
{
  gbasis();
  const SingularRing& r = *d_ring;
  const Integer& p = r.prime();
  long vdim = vectorSpaceDim();  // degree bound for the minimal polynomial
  Poly x = Poly::indet(r, varIdx);

  // We collect NF(x^0), NF(x^1), ... and look for the first linear dependence
  // over GF(p). Each NF is expressed as a sparse coefficient vector over the
  // monomials it touches; we run incremental Gaussian elimination over the
  // field.
  //
  // rowsRed[i] is the reduced (pivot-normalized) coefficient vector of the
  // i-th independent NF; pivotMon[i] is its pivot monomial; combos[i] records
  // the combination of the original NF(x^j) that produced rowsRed[i]. When a
  // new NF reduces to zero, the combination gives the dependency, hence the
  // minimal polynomial.
  std::vector<std::map<std::string, Integer>> rowsRed;
  std::vector<std::vector<Integer>> combos;
  std::vector<std::string> pivotMon;

  Poly cur = Poly::one(r);  // x^0
  for (long k = 0; k <= vdim; ++k)
  {
    Poly nf = normalForm(cur);
    std::map<std::string, Integer> row = polyToCoeffMap(nf);
    std::vector<Integer> combo(k + 1, Integer(0));
    combo[k] = Integer(1);  // this row == 1 * NF(x^k)

    // Reduce row against existing pivots.
    for (size_t i = 0; i < rowsRed.size(); ++i)
    {
      auto it = row.find(pivotMon[i]);
      if (it != row.end() && it->second.sgn() != 0)
      {
        Integer factor = it->second;  // pivot of rowsRed[i] normalized to 1
        subtractScaled(row, rowsRed[i], factor, p);
        for (size_t j = 0; j < combos[i].size(); ++j)
        {
          Integer prod = factor.modMultiply(combos[i][j], p);
          Integer neg = (p - prod).floorDivideRemainder(p);
          combo[j] = combo[j].modAdd(neg, p);
        }
      }
    }

    // Did row reduce to zero? Then we have a dependency.
    std::string newPivot;
    Integer pivotVal(0);
    for (const auto& [mon, c] : row)
    {
      if (c.sgn() != 0)
      {
        newPivot = mon;
        pivotVal = c;
        break;
      }
    }
    if (newPivot.empty())
    {
      // Dependency: sum_j combo[j] * x^j == 0 (mod ideal).
      Poly minPoly = Poly::zero(r);
      for (long j = 0; j <= k; ++j)
      {
        if (combo[j].sgn() != 0)
        {
          Poly term = Poly::constant(r, combo[j]) * x.power((unsigned long)j);
          minPoly += term;
        }
      }
      return minPoly.monic();
    }

    // Normalize the new pivot row so the pivot coefficient is 1.
    Integer inv = pivotVal.modInverse(p);
    for (auto& [mon, c] : row)
    {
      c = c.modMultiply(inv, p);
    }
    for (auto& c : combo)
    {
      c = c.modMultiply(inv, p);
    }
    rowsRed.push_back(std::move(row));
    combos.push_back(std::move(combo));
    pivotMon.push_back(newPivot);

    cur = cur * x;  // x^(k+1)
  }
  Unreachable() << "minimalPolynomial: no dependency found within vdim";
}

/* -------------------------------------------------------------------------- */
/* singularEval, toFfVal                                                      */
/* -------------------------------------------------------------------------- */

std::optional<Scalar> singularEval(const Poly& p, const PartialPoint& values)
{
  const SingularRing& r = p.sring();
  ring R = r.raw();
  int nv = rVar(R);
  Poly out = Poly::zero(r);
  for (poly q = p.raw(); q != nullptr; q = pNext(q))
  {
    Poly term = Poly::constant(r, r.numberToInteger(pGetCoeff(q)));
    for (int v = 1; v <= nv; ++v)
    {
      int e = p_GetExp(q, v, R);
      if (e != 0)
      {
        if (!values[v - 1].has_value())
        {
          return {};
        }
        term *= values[v - 1].value().power((unsigned long)e);
      }
    }
    out += term;
  }
  return {out};
}

Scalar singularEval(const Poly& p, const Point& values)
{
  const SingularRing& r = p.sring();
  ring R = r.raw();
  int nv = rVar(R);
  Poly out = Poly::zero(r);
  for (poly q = p.raw(); q != nullptr; q = pNext(q))
  {
    Poly term = Poly::constant(r, r.numberToInteger(pGetCoeff(q)));
    for (int v = 1; v <= nv; ++v)
    {
      int e = p_GetExp(q, v, R);
      if (e != 0)
      {
        term *= values[v - 1].power((unsigned long)e);
      }
    }
    out += term;
  }
  return out;
}

FiniteFieldValue toFfVal(const Scalar& elem, const FfSize& size)
{
  return {elem.constantCoeff(), size};
}

}  // namespace ff
}  // namespace theory
}  // namespace cvc5::internal

#endif /* CVC5_USE_SINGULAR */
