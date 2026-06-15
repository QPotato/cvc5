/******************************************************************************
 * This file is part of the cvc5 project.
 *
 * Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
 * in the top-level source directory and their institutional affiliations.
 * All rights reserved.  See the file COPYING in the top-level source
 * directory for licensing information.
 * ****************************************************************************
 *
 * encoding Nodes as singular ring elements.
 */

#include "cvc5_private.h"

#ifdef CVC5_USE_SINGULAR

#ifndef CVC5__THEORY__FF__SINGULAR_ENCODER_H
#define CVC5__THEORY__FF__SINGULAR_ENCODER_H

// std includes
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>

// internal includes
#include "expr/node.h"
#include "theory/ff/singular_util.h"
#include "theory/ff/util.h"

namespace cvc5::internal {
namespace theory {
namespace ff {

/**
 *  Create a Singular identifier, sanitizing varName.
 *  If index is given, fold it into the name as a "_<index>" suffix.
 *
 *  Singular identifiers must start with a letter and contain only letters,
 *  numbers, and underscores.
 */
std::string singularSym(const std::string& varName,
                        std::optional<size_t> index = {});

/**
 * Class for encoding a Node as a Poly.
 *
 * Requires two passes over the nodes. On the first pass it collects variables,
 * !=s, and bitsums. On the second, it encodes. The first stage is called
 * "Stage::Scan", the second is "Stage::Encode".
 *
 * Two stages are necessary because when creating a Singular polynomial ring,
 * one must declare all the variables up-front. So, before we create any
 * polynomials (to encode terms), we must know all the (Singular) variables.
 * Singular variables are used to encode cvc5 variables, bitsums, and witnesses
 * of disequality (a != b is encoded as (a - b)w = 1, where w is the witness).
 */
class SingularEncoder : public FieldObj
{
 public:
  /** Create a new encoder, for this field. */
  SingularEncoder(NodeManager* nm, const FfSize& size);
  /** Add a fact (one must call this twice per fact, once per stage). */
  void addFact(const Node& fact);
  /** Start Stage::Encode. */
  void endScan();
  /**
   * Get the polys who's common zero we are finding (excluding bitsums).
   * Available in Stage::Encode.
   */
  const Polys& polys() const { return d_polys; }
  /**
   * Get the bitsum polys.
   * These have form: x - b0 - 2*b1 - 4b2 ... - 2^n*b_n.
   * Available in Stage::Encode.
   */
  const Polys& bitsumPolys() const { return d_bitsumPolys; }
  /**
   * Get the poly for this term
   * Available in Stage::Encode.
   */
  const Poly& getTermEncoding(const Node& t) const { return d_cache.at(t); }
  /**
   * Get the bitsum terms (for the bitsumPolys).
   * Available in Stage::Encode.
   */
  std::vector<Node> bitsums() const;
  /**
   * The poly ring we've encoded into.
   * Available in Stage::Encode.
   */
  const SingularRing& ring() const { return *d_ring; }
  /**
   * A list of (indeterminant num, Node) pairs. Useful for extracting a model.
   * Available in Stage::Encode.
   */
  std::vector<std::pair<size_t, Node>> nodeIndets() const;
  /**
   * Does some fact that imply this poly?
   */
  bool polyHasFact(const Poly& poly) const;
  /**
   * Get the fact that implies this poly.
   */
  const Node& polyFact(const Poly& poly) const;

 private:
  /**
   * Get a fresh symbol that starts with varName.
   * If index is given, subscript the symbol by it.
   */
  std::string freshSym(const std::string& varName,
                       std::optional<size_t> index = {});
  /** a bitsum or a var */
  const Node& symNode(const std::string& s) const;
  /** have we assigned this symbol to some Node? */
  bool hasNode(const std::string& s) const;
  /** get the poly for this symbol */
  const Poly& symPoly(const std::string& s) const;
  /** encode this term as a poly (caching) */
  void encodeTerm(const Node& t);
  /** encode this fact as a poly that must be zero (caching) */
  void encodeFact(const Node& f);

  /** Which pass we're in. */
  enum class Stage
  {
    /** collecting: variable, !=, bitsum */
    Scan,
    /** encoding terms */
    Encode,
  };

  // configuration

  /**
   * The polynomial ring (populated at the end of Stage::Scan).
   *
   * Declared as the FIRST data member so that it is destroyed LAST: every Poly
   * stored in the maps/vectors below borrows this ring, and a Poly's destructor
   * frees its polynomial against the ring (p_Delete). The ring must therefore
   * outlive all of them.
   */
  std::unique_ptr<SingularRing> d_ring{};

  /** the stage that we're in; initially scanning */
  Stage d_stage{Stage::Scan};

  // populated during Stage::Scan

  /** all nodes scanned */
  std::unordered_set<Node> d_scanned{};
  /** all variables seen */
  std::unordered_set<std::string> d_vars{};
  /** map: bitsum term to its symbol */
  std::unordered_map<Node, std::string> d_bitsumSyms{};
  /** map: variable term to its symbol */
  std::unordered_map<Node, std::string> d_varSyms{};
  /** map: term (a != b) to the symbol for the inverse of (a - b) */
  std::unordered_map<Node, std::string> d_diseqSyms{};
  /** all symbols */
  std::vector<std::string> d_syms{};
  /** map: symbol name to polynomial */
  std::unordered_map<std::string, Poly> d_symPolys{};
  /** map: symbol name to term */
  std::unordered_map<std::string, Node> d_symNodes{};

  // populated during Stage::Encode

  /** encoding cache */
  std::unordered_map<Node, Poly> d_cache{};
  /** polynomials that must be zero (except bitsums) */
  Polys d_polys{};
  /** bitsum polynomials that must be zero */
  Polys d_bitsumPolys{};
  /** polys to the facts that imply them */
  std::unordered_map<std::string, Node> d_polyFacts{};
};

}  // namespace ff
}  // namespace theory
}  // namespace cvc5::internal

#endif /* CVC5__THEORY__FF__SINGULAR_ENCODER_H */

#endif /* CVC5_USE_SINGULAR */
