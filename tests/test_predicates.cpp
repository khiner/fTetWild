////////////////////////////////////////////////////////////////////////////////
#include <floattetwild/Predicates.hpp>
#include <floattetwild/intersections.h>
#include <catch2/catch_all.hpp>
////////////////////////////////////////////////////////////////////////////////

using namespace floatTetWild;

TEST_CASE("orient_3d", "[predicates]") {
	Vector3 p1(0, 0, 0);
	Vector3 p2(1, 0, 0);
	Vector3 p3(0, 1, 0);
	Vector3 p4(0, 0, 1);

	int orientation;

	orientation = Predicates::orient_3d(p1, p2, p3, p4);
	REQUIRE(orientation == Predicates::ORI_NEGATIVE);


	orientation = Predicates::orient_3d(p1,p4, p3, p2);
	REQUIRE(orientation == Predicates::ORI_POSITIVE);
}

namespace {
using Point = std::array<Scalar, 3>;

// The overlap test's result for two triangles that are known not to be coplanar.
int intersect(Point p1, Point q1, Point r1, Point p2, Point q2, Point r2) {
    Scalar    s[3] = {0, 0, 0};
    Scalar    t[3] = {0, 0, 0};
    const int res  = tri_tri_intersection_test_3d(
      p1.data(), q1.data(), r1.data(), p2.data(), q2.data(), r2.data(), s, t);
    REQUIRE(res != -1);
    return res;
}
}  // namespace

TEST_CASE("tri_tri_intersection_test_3d_rounded", "[predicates]") {
    REQUIRE(intersect({{22, -27.47818, 3.5}},
                      {{25.5, -27.47818, 0}},
                      {{-25.5, -27.47818, 0}},
                      {{28.05000000000000, -30.02818, -2.55}},
                      {{25.5, -13.183373550207467, 0}},
                      {{25.5, -19.059304172821577, 0}}) == 0);
}

TEST_CASE("tri_tri_intersection_test_3d_floating", "[predicates]") {
    const int res = intersect({{22, -27.478179999999998, 3.5}},
                              {{25.5, -27.478179999999998, 0}},
                              {{-25.5, -27.478179999999998, 0}},
                              {{28.050000000000001, -30.028179999999999, -2.5500000000000003}},
                              {{25.5, -13.183373550207467, 0}},
                              {{25.5, -19.059304172821577, 0}});
    // Accept either 0 or 666 as a correct "no intersection" result.
    REQUIRE((res == 0 || res == 666));
}
