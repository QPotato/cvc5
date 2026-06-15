/******************************************************************************
 * This file is part of the cvc5 project.
 *
 * Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
 * in the top-level source directory and their institutional affiliations.
 * All rights reserved.  See the file COPYING in the top-level source
 * directory for licensing information.
 * ****************************************************************************
 *
 * Black box testing of ff univariate root finding.
 */

#include "cvc5_private.h"

#ifdef CVC5_USE_SINGULAR

#include <memory>
#include <utility>
#include <vector>

#include "test_smt.h"
#include "theory/ff/singular_util.h"
#include "util/integer.h"
#include "util/singular_globals.h"

namespace cvc5::internal {

using namespace context;
using namespace theory;

namespace test {

class TestTheoryFfRootsBlack : public TestSmt
{
  void SetUp() override
  {
    TestSmt::SetUp();
    initSingular();
  }
};

#define TINY_MODULUS "7"
#define BIG_MODULUS                                                            \
  "57896044618658097711785492504343953926634992332820282019728792003956564819" \
  "949"

TEST_F(TestTheoryFfRootsBlack, DistinctRootsPoly)
{
  {
    ff::SingularRing ring(Integer(TINY_MODULUS), {"x", "y", "z"});
    ff::Poly x = ff::Poly::indet(ring, 0);
    ff::Poly one = ff::Poly::one(ring);
    {
      ff::Poly f = x;
      ff::Poly ex = x;
      EXPECT_EQ(ff::distinctRootsPoly(f), ex);
    }

    {
      ff::Poly f = x * x * x;
      ff::Poly ex = x;
      EXPECT_EQ(ff::distinctRootsPoly(f), ex);
    }

    {
      ff::Poly f = x * x * x * (x - one);
      ff::Poly ex = x * (x - one);
      EXPECT_EQ(ff::distinctRootsPoly(f), ex);
    }

    {
      ff::Poly f = x * x * x * (x - one) * (x * x + one);
      ff::Poly ex = x * (x - one);
      EXPECT_EQ(ff::distinctRootsPoly(f), ex);
    }
  }

  {
    // 2^255-19 (prime)
    ff::SingularRing ring(Integer(BIG_MODULUS), {"x", "y", "z"});
    ff::Poly x = ff::Poly::indet(ring, 0);
    ff::Poly one = ff::Poly::one(ring);
    ff::Poly two = ff::Poly::constant(ring, Integer(2));
    {
      ff::Poly f = x;
      ff::Poly ex = x;
      EXPECT_EQ(ff::distinctRootsPoly(f), ex);
    }

    {
      ff::Poly f = x * x * x;
      ff::Poly ex = x;
      EXPECT_EQ(ff::distinctRootsPoly(f), ex);
    }

    {
      ff::Poly f = x * x * x * (x - one);
      ff::Poly ex = x * (x - one);
      EXPECT_EQ(ff::distinctRootsPoly(f), ex);
    }

    {
      ff::Poly f = x * x * x * (x - one) * (x * x + two);
      ff::Poly ex = x * (x - one);
      EXPECT_EQ(ff::distinctRootsPoly(f), ex);
    }
  }
}

TEST_F(TestTheoryFfRootsBlack, RootsZero)
{
  {
    ff::SingularRing ring(Integer(TINY_MODULUS), {"x", "y", "z"});
    ff::Poly x = ff::Poly::indet(ring, 0);
    ff::Poly one = ff::Poly::one(ring);
    ff::Poly z = ff::Poly::zero(ring);
    {
      ff::Poly f = x;
      std::vector<ff::Poly> roots = {ff::Poly::constant(ring, Integer(0))};
      EXPECT_EQ(ff::roots(f), roots);
    }

    {
      ff::Poly f = x * x * x;
      std::vector<ff::Poly> roots = {ff::Poly::constant(ring, Integer(0))};
      EXPECT_EQ(ff::roots(f), roots);
    }

    {
      ff::Poly f = x * (x * x + one);
      std::vector<ff::Poly> roots = {ff::Poly::constant(ring, Integer(0))};
      EXPECT_EQ(ff::roots(f), roots);
    }
  }

  {
    ff::SingularRing ring(Integer(BIG_MODULUS), {"x", "y", "z"});
    ff::Poly x = ff::Poly::indet(ring, 0);
    ff::Poly two = ff::Poly::constant(ring, Integer(2));
    {
      ff::Poly f = x;
      std::vector<ff::Poly> roots = {ff::Poly::constant(ring, Integer(0))};
      EXPECT_EQ(ff::roots(f), roots);
    }

    {
      ff::Poly f = x * x * x;
      std::vector<ff::Poly> roots = {ff::Poly::constant(ring, Integer(0))};
      EXPECT_EQ(ff::roots(f), roots);
    }

    {
      ff::Poly f = x * (x * x + two);
      std::vector<ff::Poly> roots = {ff::Poly::constant(ring, Integer(0))};
      EXPECT_EQ(ff::roots(f), roots);
    }
  }
}

TEST_F(TestTheoryFfRootsBlack, RootsFull)
{
  {
    ff::SingularRing ring(Integer(TINY_MODULUS), {"x", "y", "z"});
    ff::Poly x = ff::Poly::indet(ring, 0);
    ff::Poly one = ff::Poly::one(ring);

    {
      ff::Poly f = x * (x - one);
      std::vector<ff::Poly> roots = {ff::Poly::constant(ring, Integer(0)),
                                     ff::Poly::constant(ring, Integer(1))};
      EXPECT_EQ(ff::roots(f), roots);
    }

    {
      ff::Poly f = x * (x - one) * x * (x - one);
      std::vector<ff::Poly> roots = {ff::Poly::constant(ring, Integer(0)),
                                     ff::Poly::constant(ring, Integer(1))};
      EXPECT_EQ(ff::roots(f), roots);
    }

    {
      ff::Poly f =
          x * (x - one) * x * (x - one) * (x * x + one) * (x * x + one);
      std::vector<ff::Poly> roots = {ff::Poly::constant(ring, Integer(0)),
                                     ff::Poly::constant(ring, Integer(1))};
      EXPECT_EQ(ff::roots(f), roots);
    }

    {
      ff::Poly f = (x * x + one) * (x * x + one);
      std::vector<ff::Poly> roots = {};
      EXPECT_EQ(ff::roots(f), roots);
    }

    {
      ff::Poly f = x * x - x + one;
      // roots over GF(7) are 5 and 3; sorted by Singular's (symmetric) string
      // representation that is ["-2"(=5), "3"], i.e. {5, 3}.
      std::vector<ff::Poly> roots = {ff::Poly::constant(ring, Integer(5)),
                                     ff::Poly::constant(ring, Integer(3))};
      EXPECT_EQ(ff::roots(f), roots);
    }
  }

  {
    ff::SingularRing ring(Integer(BIG_MODULUS), {"x", "y", "z"});
    ff::Poly x = ff::Poly::indet(ring, 0);
    ff::Poly one = ff::Poly::one(ring);
    ff::Poly two = ff::Poly::constant(ring, Integer(2));

    {
      ff::Poly f = x * (x - one);
      std::vector<ff::Poly> roots = {ff::Poly::constant(ring, Integer(0)),
                                     ff::Poly::constant(ring, Integer(1))};
      EXPECT_EQ(ff::roots(f), roots);
    }

    {
      ff::Poly f = x * (x - one) * x * (x - one);
      std::vector<ff::Poly> roots = {ff::Poly::constant(ring, Integer(0)),
                                     ff::Poly::constant(ring, Integer(1))};
      EXPECT_EQ(ff::roots(f), roots);
    }

    {
      ff::Poly f =
          x * (x - one) * x * (x - one) * (x * x + two) * (x * x + two);
      std::vector<ff::Poly> roots = {ff::Poly::constant(ring, Integer(0)),
                                     ff::Poly::constant(ring, Integer(1))};
      EXPECT_EQ(ff::roots(f), roots);
    }

    {
      ff::Poly f = (x * x + two) * (x * x + two);
      std::vector<ff::Poly> roots = {};
      EXPECT_EQ(ff::roots(f), roots);
    }

    {
      ff::Poly f = x * x - x + one;
      // The two roots of x^2-x+1 mod (2^255-19), as canonical [0,p) integers.
      // roots() sorts by Singular's symmetric string representation; the larger
      // canonical value (> p/2) prints with a leading '-', so it sorts first.
      std::vector<ff::Poly> roots = {
          ff::Poly::constant(
              ring,
              Integer("3251576818157896011469325613977277291600281449988881"
                      "3854556049534830466505397")),
          ff::Poly::constant(
              ring,
              Integer("2538027643707913759709223636457118101063217783293146"
                      "8165172742469126098314553")),
      };
      EXPECT_EQ(ff::roots(f), roots);
    }
  }
}

}  // namespace test
}  // namespace cvc5::internal
#endif  // CVC5_USE_SINGULAR
