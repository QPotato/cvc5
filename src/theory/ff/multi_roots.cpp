/******************************************************************************
 * This file is part of the cvc5 project.
 *
 * Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
 * in the top-level source directory and their institutional affiliations.
 * All rights reserved.  See the file COPYING in the top-level source
 * directory for licensing information.
 * ****************************************************************************
 *
 * Multivariate root finding. Implements "FindZero" from [OKTB23].
 *
 * [OKTB23]: https://doi.org/10.1007/978-3-031-37703-7_8
 */

#ifdef CVC5_USE_SINGULAR

#include "theory/ff/multi_roots.h"

// std includes
#include <algorithm>
#include <memory>
#include <unordered_map>

// internal includes
#include "theory/ff/singular_util.h"
#include "theory/ff/util.h"
#include "util/resource_manager.h"

namespace cvc5::internal {
namespace theory {
namespace ff {

AssignmentEnumerator::~AssignmentEnumerator() = default;

ListEnumerator::ListEnumerator(std::vector<Poly>&& options)
    : d_remainingOptions(std::move(options))
{
  std::reverse(d_remainingOptions.begin(), d_remainingOptions.end());
}

ListEnumerator::~ListEnumerator() {};

std::optional<Poly> ListEnumerator::next()
{
  if (d_remainingOptions.empty())
  {
    return {};
  }
  else
  {
    Poly v = d_remainingOptions.back();
    d_remainingOptions.pop_back();
    return v;
  }
}

std::string ListEnumerator::name() { return "list"; }

std::unique_ptr<ListEnumerator> factorEnumerator(Poly univariatePoly)
{
  int varIdx = univariatePoly.univariateIndetIndex();
  Assert(varIdx >= 0);
  Trace("ff::model::factor")
      << "roots for: " << univariatePoly.str() << std::endl;
  std::vector<Poly> theRoots = roots(univariatePoly);
  std::vector<Poly> linears{};
  Poly var = Poly::indet(univariatePoly.sring(), varIdx);
  for (const auto& r : theRoots)
  {
    linears.push_back(var - r);
  }
  return std::make_unique<ListEnumerator>(std::move(linears));
}

RoundRobinEnumerator::RoundRobinEnumerator(const std::vector<Poly>& vars,
                                           const SingularRing& ring)
    : d_vars(vars),
      d_ring(ring),
      d_idx(0),
      // LogCardinality is 1, so the field has prime() elements
      d_maxIdx(ring.prime() * Integer(static_cast<unsigned long>(vars.size())))
{
}

RoundRobinEnumerator::~RoundRobinEnumerator() {}

std::optional<Poly> RoundRobinEnumerator::next()
{
  std::optional<Poly> ret{};
  if (d_idx != d_maxIdx)
  {
    Integer numVars(static_cast<unsigned long>(d_vars.size()));
    size_t whichVar = d_idx.floorDivideRemainder(numVars).getUnsignedLong();
    Integer whichVal = d_idx.floorDivideQuotient(numVars);
    Poly val = Poly::constant(d_ring, whichVal);
    ret = d_vars[whichVar] - val;
    d_idx += Integer(1);
  }
  return ret;
}

std::string RoundRobinEnumerator::name() { return "round-robin"; }

bool isUnsat(Ideal& ideal)
{
  const auto& g = ideal.gbasis();
  return g.size() == 1 && !g[0].isZero() && g[0].deg() <= 0;
}

std::pair<size_t, Poly> extractAssignment(const Poly& elem)
{
  Assert(elem.deg() == 1);
  Assert(elem.numTerms() <= 2);
  const Poly m = elem.monic();
  int varNumber = elem.univariateIndetIndex();
  Assert(varNumber >= 0);
  const Integer& prime = elem.sring().prime();
  // value is -ConstantCoeff(m), normalized into [0, prime)
  Integer c = m.constantCoeff();
  Integer val = (prime - c).floorDivideRemainder(prime);
  return {static_cast<size_t>(varNumber), Poly::constant(elem.sring(), val)};
}

std::unordered_set<std::string> assignedVars(Ideal& ideal)
{
  std::unordered_set<std::string> ret{};
  for (const auto& g : ideal.gbasis())
  {
    if (g.deg() == 1)
    {
      int varNumber = g.univariateIndetIndex();
      if (varNumber >= 0)
      {
        ret.insert(Poly::indet(ideal.sring(), varNumber).str());
      }
    }
  }
  return ret;
}

bool allVarsAssigned(Ideal& ideal)
{
  return assignedVars(ideal).size() == ideal.sring().nVars();
}

std::unique_ptr<AssignmentEnumerator> applyRule(Ideal& ideal,
                                                FfStatistics* stats)
{
  const SingularRing& ring = ideal.sring();
  Assert(!isUnsat(ideal));
  // first, we look for super-linear univariate polynomials.
  for (const auto& p : ideal.gbasis())
  {
    int varNumber = p.univariateIndetIndex();
    if (varNumber >= 0 && p.deg() > 1)
    {
      return factorEnumerator(p);
    }
  }
  // now, we check the dimension
  if (ideal.isZeroDim())
  {
    if (stats) ++stats->d_idealMinPoly;
    // If zero-dimensional, we compute a minimal polynomial in some unset
    // variable.
    std::unordered_set<std::string> alreadySet = assignedVars(ideal);
    for (size_t i = 0, n = ring.nVars(); i < n; ++i)
    {
      if (!alreadySet.count(Poly::indet(ring, i).str()))
      {
        Poly minPoly = ideal.minimalPolynomial(i);
        return factorEnumerator(minPoly);
      }
    }
    Unreachable()
        << "There should be no unset variables in zero-dimensional ideal";
  }
  else
  {
    if (stats) ++stats->d_idealPosDim;
    // If positive dimensional, we make a list of unset variables and
    // round-robin guess.
    //
    // TODO(aozdemir): better model construction (cvc5-wishues/issues/138)
    std::unordered_set<std::string> alreadySet = assignedVars(ideal);
    std::vector<Poly> toGuess{};
    for (size_t i = 0, n = ring.nVars(); i < n; ++i)
    {
      Poly var = Poly::indet(ring, i);
      if (!alreadySet.count(var.str()))
      {
        toGuess.push_back(var);
      }
    }
    return std::make_unique<RoundRobinEnumerator>(toGuess, ring);
  }
}

std::vector<Poly> findZero(Ideal& initialIdeal,
                           const Env& env,
                           FfStatistics* stats)
{
  const SingularRing& ring = initialIdeal.sring();
  // We maintain two stacks:
  // * one of ideals
  // * one of branchers
  //
  // If brancher B has the same index as ideal I, then B represents possible
  // expansions of ideal I (equivalently, restrictions of I's variety).
  //
  // NB: FindZero of [OKTB23] also takes a partial map M as input. GB(I)
  // implicitly represents M: GB(I) contains a univariate linear polynomial
  // Xi - k, if and only iff M[Xi] = k.
  //
  // NB: FindZero of [OKTB23] is recursive. That recursion is flattened here
  // using the two stacks. The stack of ideals represents the input to
  // recursive FindZero: GB(I). The stack of branchers represents the
  // continuation context (which iteration of the for loop to return to).

  // goal: find a zero for any ideal in the stack.
  std::vector<Ideal> ideals{initialIdeal};
  if (TraceIsOn("ff::model::branch"))
  {
    Trace("ff::model::branch") << "init polys: " << std::endl;
    for (const auto& p : initialIdeal.gens())
    {
      Trace("ff::model::branch") << " * " << p.str() << std::endl;
    }
  }

  std::vector<std::unique_ptr<AssignmentEnumerator>> branchers{};
  // while some ideal might have a zero.
  while (!ideals.empty())
  {
    // check for timeout
    if (env.getResourceManager()->outOfTime())
    {
      throw FfTimeoutException("findZero");
    }
    // choose one ideal
    Ideal& ideal = ideals.back();
    // make sure we have a GBasis:
    ideal.gbasis();
    // If the ideal is UNSAT, drop it.
    if (isUnsat(ideal))
    {
      ideals.pop_back();
    }
    // If the ideal has a linear polynomial in each variable, we've found a
    // variety element (a model).
    else if (allVarsAssigned(ideal))
    {
      std::unordered_map<size_t, Poly> varNumToValue{};
      const auto& gens = ideal.gbasis();
      size_t numIndets = ring.nVars();
      Assert(gens.size() == numIndets);
      for (const auto& g : gens)
      {
        varNumToValue.insert(extractAssignment(g));
      }
      std::vector<Poly> values{};
      for (size_t i = 0; i < numIndets; ++i)
      {
        values.push_back(varNumToValue.at(i));
      }
      return values;
    }
    // If there are more ideals than branchers, branch
    else if (ideals.size() > branchers.size())
    {
      Assert(ideals.size() == branchers.size() + 1);
      branchers.push_back(applyRule(ideal, stats));
      Trace("ff::model::branch")
          << "brancher: " << branchers.back()->name() << std::endl;
      if (TraceIsOn("ff::model::branch"))
      {
        Trace("ff::model::branch") << "ideal polys: " << std::endl;
        for (const auto& p : ideal.gens())
        {
          Trace("ff::model::branch") << " * " << p.str() << std::endl;
        }
      }
    }
    // Otherwise, this ideal should have a brancher; get the next branch
    else
    {
      Assert(ideals.size() == branchers.size());
      std::optional<Poly> choicePoly = branchers.back()->next();
      // construct a new ideal from the branch
      if (choicePoly.has_value())
      {
        Trace("ff::model::branch")
            << "level: " << branchers.size()
            << ", brancher: " << branchers.back()->name()
            << ", branch: " << choicePoly.value().str() << std::endl;
        std::vector<Poly> newGens = ideal.gbasis();
        newGens.push_back(choicePoly.value());
        ideals.push_back(Ideal(ideal.sring(), newGens));
      }
      // or drop this ideal & brancher if we're out of branches.
      else
      {
        branchers.pop_back();
        ideals.pop_back();
      }
    }
  }
  // Could not find any solution; return empty.
  return {};
}

}  // namespace ff
}  // namespace theory
}  // namespace cvc5::internal

#endif /* CVC5_USE_SINGULAR */
