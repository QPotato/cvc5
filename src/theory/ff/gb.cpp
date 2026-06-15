/******************************************************************************
 * This file is part of the cvc5 project.
 *
 * Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
 * in the top-level source directory and their institutional affiliations.
 * All rights reserved.  See the file COPYING in the top-level source
 * directory for licensing information.
 * ****************************************************************************
 *
 */

#ifdef CVC5_USE_SINGULAR

#include "theory/ff/gb.h"

// internal includes
#include "options/ff_options.h"
#include "smt/env.h"
#include "theory/ff/multi_roots.h"
#include "theory/ff/singular_encoder.h"
#include "theory/ff/singular_util.h"
#include "util/resource_manager.h"

namespace cvc5::internal {
namespace theory {
namespace ff {

FfResult gb(const std::vector<Node>& facts,
            const FfSize& size,
            const Env& env,
            FfStatistics* stats)
{
  SingularEncoder enc(env.getNodeManager(), size);
  // collect leaves
  for (const Node& node : facts)
  {
    enc.addFact(node);
  }
  enc.endScan();
  // assert facts
  for (const Node& node : facts)
  {
    enc.addFact(node);
  }

  // compute a GB
  Polys generators;
  generators.insert(generators.end(), enc.polys().begin(), enc.polys().end());
  generators.insert(
      generators.end(), enc.bitsumPolys().begin(), enc.bitsumPolys().end());
  if (env.getOptions().ff.ffFieldPolys)
  {
    // x^p - x; field polys are only used on small fields, where the exponent
    // fits in an unsigned long.
    Assert(enc.ring().isSmall());
    unsigned long p = enc.ring().prime().getUnsignedLong();
    for (size_t i = 0, n = enc.ring().nVars(); i < n; ++i)
    {
      Poly var = Poly::indet(enc.ring(), i);
      generators.push_back(var.power(p) - var);
    }
  }
  Ideal ideal(enc.ring(), generators);
  if (stats) ++stats->d_numGbRuns;
  if (env.getResourceManager()->outOfTime())
  {
    throw FfTimeoutException("GBasis");
  }
  {
    CodeTimer timer(stats ? &stats->d_timeGbRuns : nullptr);
    ideal.gbasis();
  }

  // if it is trivial, create a conflict
  if (ideal.isWholeRing())
  {
    Trace("ff::gb") << "Trivial GB" << std::endl;
    if (stats) ++stats->d_numTrivialUnsat;
    // coarse conflict core
    return facts;
  }
  else
  {
    Trace("ff::gb") << "Non-trivial GB" << std::endl;

    // common root (vec of base ring elements)
    Point root;
    {
      CodeTimer timer(stats ? &stats->d_modelConstructionTime : nullptr);
      root = findZero(ideal, env, stats);
    }

    if (root.empty())
    {
      // set trivial conflict
      return facts;
    }
    else
    {
      FfModel model;
      for (const auto& [idx, node] : enc.nodeIndets())
      {
        if (isFfLeaf(node))
        {
          model.emplace(node, toFfVal(root[idx], size));
        }
      }
      return model;
    }
  }
}

}  // namespace ff
}  // namespace theory
}  // namespace cvc5::internal

#endif
