/******************************************************************************
 * This file is part of the cvc5 project.
 *
 * Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
 * in the top-level source directory and their institutional affiliations.
 * All rights reserved.  See the file COPYING in the top-level source
 * directory for licensing information.
 * ****************************************************************************
 *
 * Black box testing of ff multivariate roots.
 */

#include "cvc5_private.h"

#ifdef CVC5_USE_SINGULAR

#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "test_env.h"
#include "theory/ff/multi_roots.h"
#include "theory/ff/singular_util.h"
#include "util/integer.h"
#include "util/singular_globals.h"

namespace cvc5::internal {

using namespace context;
using namespace theory;

namespace test {

class TestTheoryFfModelBlack : public TestEnv
{
  void SetUp() override
  {
    TestEnv::SetUp();
    initSingular();
  }
};

TEST_F(TestTheoryFfModelBlack, UnivariateEnumerator)
{
  ff::SingularRing ring(Integer(7), {"a"});
  {
    ff::Poly var = ff::Poly::indet(ring, 0);
    std::vector<ff::Poly> assignments = {var - ff::Poly::zero(ring),
                                         var - ff::Poly::one(ring)};
    std::unique_ptr<ff::AssignmentEnumerator> a =
        std::make_unique<ff::ListEnumerator>(std::move(assignments));
    std::optional<ff::Poly> first = {var - ff::Poly::zero(ring)};
    std::optional<ff::Poly> second = {var - ff::Poly::one(ring)};
    std::optional<ff::Poly> third = {};
    EXPECT_EQ(a->next(), first);
    EXPECT_EQ(a->next(), second);
    EXPECT_EQ(a->next(), third);
    EXPECT_EQ(a->next(), third);
    EXPECT_EQ(a->next(), third);
  }
}

TEST_F(TestTheoryFfModelBlack, RoundRobinEnumerator)
{
  ff::SingularRing ring(Integer(3), {"a", "b", "c"});
  {
    std::vector<ff::Poly> vars{ff::Poly::indet(ring, 0),
                               ff::Poly::indet(ring, 1),
                               ff::Poly::indet(ring, 2)};
    std::unique_ptr<ff::AssignmentEnumerator> a =
        std::make_unique<ff::RoundRobinEnumerator>(vars, ring);
    std::optional<ff::Poly> next{};
    ff::Poly zero = ff::Poly::constant(ring, Integer(0));
    ff::Poly one = ff::Poly::constant(ring, Integer(1));
    ff::Poly two = ff::Poly::constant(ring, Integer(2));
    next = vars[0] - zero;
    EXPECT_EQ(a->next(), next);
    next = vars[1] - zero;
    EXPECT_EQ(a->next(), next);
    next = vars[2] - zero;
    EXPECT_EQ(a->next(), next);
    next = vars[0] - one;
    EXPECT_EQ(a->next(), next);
    next = vars[1] - one;
    EXPECT_EQ(a->next(), next);
    next = vars[2] - one;
    EXPECT_EQ(a->next(), next);
    next = vars[0] - two;
    EXPECT_EQ(a->next(), next);
    next = vars[1] - two;
    EXPECT_EQ(a->next(), next);
    next = vars[2] - two;
    EXPECT_EQ(a->next(), next);
    next = {};
    EXPECT_EQ(a->next(), next);
    next = {};
    EXPECT_EQ(a->next(), next);
  }
}

TEST_F(TestTheoryFfModelBlack, IsUnsat)
{
  ff::SingularRing ring(Integer(3), {"a", "b", "c"});
  {
    ff::Poly a = ff::Poly::indet(ring, 0);
    ff::Poly b = ff::Poly::indet(ring, 1);
    ff::Poly c = ff::Poly::indet(ring, 2);
    ff::Poly one = ff::Poly::one(ring);
    ff::Poly two = ff::Poly::constant(ring, Integer(2));

    ff::Ideal i1(ring, {a * (a - one)});
    EXPECT_EQ(ff::isUnsat(i1), false);
    ff::Ideal i2(ring, {a});
    EXPECT_EQ(ff::isUnsat(i2), false);
    ff::Ideal i3(ring, {a, b - one});
    EXPECT_EQ(ff::isUnsat(i3), false);
    ff::Ideal i4(ring, {a, b - one, c});
    EXPECT_EQ(ff::isUnsat(i4), false);
    ff::Ideal i5(ring, {a, a - one});
    EXPECT_EQ(ff::isUnsat(i5), true);
    // false b/c no field poly
    ff::Ideal i6(ring, {a * (a - one) * (a - two) - one});
    EXPECT_EQ(ff::isUnsat(i6), false);
    ff::Ideal i7(ring,
                 {a * (a - one) * (a - two) - one, a * a * a - a});
    EXPECT_EQ(ff::isUnsat(i7), true);
    ff::Ideal i8(ring, {(a - b) * c - one, a - b});
    EXPECT_EQ(ff::isUnsat(i8), true);
  }
}

TEST_F(TestTheoryFfModelBlack, ExtractAssignment)
{
  ff::SingularRing ring(Integer(3), {"a", "b", "c"});
  {
    ff::Poly a = ff::Poly::indet(ring, 0);
    ff::Poly c = ff::Poly::indet(ring, 2);
    ff::Poly one = ff::Poly::one(ring);
    std::pair<size_t, ff::Poly> r;
    r = {0, ff::Poly::one(ring)};
    EXPECT_EQ(ff::extractAssignment(a - one), r);
    r = {0, ff::Poly::zero(ring)};
    EXPECT_EQ(ff::extractAssignment(a), r);
    r = {0, ff::Poly::one(ring)};
    EXPECT_EQ(ff::extractAssignment((-a) + one), r);
    r = {2, ff::Poly::zero(ring)};
    EXPECT_EQ(ff::extractAssignment(c), r);
  }
}

TEST_F(TestTheoryFfModelBlack, CommonRoot)
{
  ff::SingularRing ring(Integer(3), {"a", "b"});
  ff::Poly a = ff::Poly::indet(ring, 0);
  ff::Poly b = ff::Poly::indet(ring, 1);
  ff::Poly one = ff::Poly::one(ring);
  ff::Poly z = ff::Poly::zero(ring);

  {
    ff::Ideal ideal(ring, {a * a - a, b * b - b, a - b, a});
    std::vector<ff::Poly> values = {z, z};
    EXPECT_EQ(ff::findZero(ideal, *d_env), values);
  }

  {
    ff::Ideal ideal(ring, {a * a - a, b * b - b, a + b - one, a});
    std::vector<ff::Poly> values = {z, z + one};
    EXPECT_EQ(ff::findZero(ideal, *d_env), values);
  }

  {
    ff::Ideal ideal(ring, {a, a - one});
    std::vector<ff::Poly> values = {};
    EXPECT_EQ(ff::findZero(ideal, *d_env), values);
  }

  {
    ff::Ideal ideal(ring,
                    {a * (a - one) * (a - ff::Poly::constant(ring, Integer(2)))
                     - one});
    std::vector<ff::Poly> values = {};
    EXPECT_EQ(ff::findZero(ideal, *d_env), values);
  }

  {
    ff::Ideal ideal(ring, {a * b - one});
    std::vector<ff::Poly> values = ff::findZero(ideal, *d_env);
    EXPECT_EQ(values[0] * values[1], z + one);
  }

  {
    ff::Ideal ideal(ring, {a * b - one, b});
    std::vector<ff::Poly> values = ff::findZero(ideal, *d_env);
    EXPECT_EQ(values.size(), 0);
  }

  {
    ff::Ideal ideal(ring, {a * b - one, b - ff::Poly::constant(ring, Integer(2))});
    std::vector<ff::Poly> values = {z + ff::Poly::constant(ring, Integer(2)),
                                    z + ff::Poly::constant(ring, Integer(2))};
    EXPECT_EQ(ff::findZero(ideal, *d_env), values);
  }
}

TEST_F(TestTheoryFfModelBlack, CommonRootBig)
{
  ff::SingularRing ring(Integer(17), {"a", "b", "c", "d"});
  ff::Poly a = ff::Poly::indet(ring, 0);
  ff::Poly b = ff::Poly::indet(ring, 1);
  ff::Poly c = ff::Poly::indet(ring, 2);
  ff::Poly d = ff::Poly::indet(ring, 3);
  ff::Poly one = ff::Poly::one(ring);
  ff::Poly z = ff::Poly::zero(ring);

  ff::Ideal ideal(ring, {a * a - a, b * b - b, a - b, a, c * d - one});
  std::vector<ff::Poly> values = ff::findZero(ideal, *d_env);
  EXPECT_EQ(values[0], z);
  EXPECT_EQ(values[1], z);
  EXPECT_EQ(values[2] * values[3], z + one);
}

TEST_F(TestTheoryFfModelBlack, CommonRootCosntraints)
{
  ff::SingularRing ring(Integer(17), {"a", "b", "c"});
  ff::Poly a = ff::Poly::indet(ring, 0);
  ff::Poly b = ff::Poly::indet(ring, 1);
  ff::Poly c = ff::Poly::indet(ring, 2);
  ff::Poly one = ff::Poly::one(ring);
  ff::Poly z = ff::Poly::zero(ring);
  // b is a perfect square
  // c is its inverse
  ff::Ideal ideal(ring, {a * a - b, b * c - one});
  std::vector<ff::Poly> values = ff::findZero(ideal, *d_env);
  EXPECT_EQ(values[0] * values[0], values[1]);
  EXPECT_EQ(values[1] * values[2], z + one);
}

}  // namespace test
}  // namespace cvc5::internal
#endif  // CVC5_USE_SINGULAR
